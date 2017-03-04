/*

$Id: ungrout_main.c,v 1.2 2012/12/05 10:50:56 jim Exp $ 

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <fitsio.h>
#include <math.h>

#include <tools.h>

/* Command line arguments and defaults */

enum {SRCPATH_ARG,
      VERB_ARG,
      NOVERB_ARG};

static struct option myoptions [] = {
    {"srcpath",required_argument,NULL,SRCPATH_ARG},
    {"verbose",no_argument,NULL,VERB_ARG},
    {"noverbose",no_argument,NULL,NOVERB_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: ungrout incatalogue outcatalogue\n"
    "[--srcpath=%s] [--(no)verbose (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define SRCPATH_DEF "."
#define VERB_DEF 0

/* This is specifically for VIRCAM files. Only correct the first 7 apertures */

#define NPAWS 6
#define NEXTN 16
#define NAPCOR 7

/* imcore catalogue type 6 column definitions we need */

#define COL_X           3
#define COL_Y           5
#define COL_PEAKHEIGHT 18
#define COL_PKHTERR    19
#define COL_APFLUX1    20
#define COL_APFLUX1ERR 21

typedef struct {
    struct wcsprm *wcs;
    int           nxout;
    int           nyout;
    float         apcors[NAPCOR];
    float         magzpt;
    float         seeing;
    float         ellipt;
} tables;

static int gettabinfo(char *cat, tables *tab, char *errmsg); 
static int getinputinfo(char *infile, char *srcpath, int verbose,
			char pawcats[][BUFSIZ], char pawconfs[][BUFSIZ], 
			char *errmsg);
static void tidytabs(tables maintab, tables *tabs);

int main (int argc, char *argv[]) {
    char *infile,*outfile,errmsg[BUFSIZ],pawcats[NPAWS][BUFSIZ],sval[16];
    char pawconfs[NPAWS][BUFSIZ],catextn[BUFSIZ],msg[BUFSIZ],colname[16];
    int c,option_index,nhdui,retval,m,i,j,status,hdutype,anynul;
    int npts,k,ix,iy,ind,k2,res_col;
    long nrows,naxes[2];
    short int *confmap;
    fitsfile *iptr=NULL,*optr=NULL;
    tables maintab,pawtabs[NPAWS*NEXTN];
    double *dspace,*x,*y,*x2,*y2,*coords;
    double *xin,*yin,*xout,*yout;
    float *fspace,*cnumb,*fluxes,*fluxerrs,*fluxcor,*pk,*pkerr,delta_mag;
    float confwt,expon,scale,xmax,ymax,work[NEXTN];

    /* Set up defaults for command line switches */

    char srcpath[BUFSIZ];
    (void)strcpy(srcpath,SRCPATH_DEF);
    int verbose = VERB_DEF;

    /* First get the command line arguments */

    if (argc < 3) {
	fprintf(stderr,USAGE,srcpath,YESNO(verbose));
	exit(1);
    } else {
	infile = argv[1];
	outfile = argv[2];
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case SRCPATH_ARG:
	    (void)strcpy(srcpath,optarg);
	    break;
	case VERB_ARG:
	    verbose = 1;
	    break;
	case NOVERB_ARG:
	    verbose = 0;
	    break;
	default:
	    fprintf(stderr,USAGE,SRCPATH_DEF,YESNO(VERB_DEF));
	    exit(1);
	}
    }

    /* If doing verbose output set stdout to have no buffering */

    if (verbose)
	setvbuf(stdout,NULL,_IONBF,0);

    /* Do some sanity checks. First make sure we can open the input file */

    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(iptr,&nhdui,&status);
    if (status != 0) {
	fprintf(stderr,"UNGROUT: Can't open file %s\n",infile);
	closefits(iptr);
	exit(1);
    } else if (nhdui != 2) {
	fprintf(stderr,"UNGROUT: Input file %s should have only 1 extension\n",
		infile);
	closefits(iptr);
	exit(1);
    }

    /* See if this has already been grouted. */

    (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
    (void)fits_read_key(iptr,TLOGICAL,"GROUTED",&k,NULL,&status);
    if (status != 0) {
	status = 0;
	if (verbose)
	    fprintf(stdout,"UNGROUT: File %s hasn't been grouted\n",infile);
	closefits(iptr);
	exit(0);
    }

    /* Initialise wcs structures */

    maintab.wcs = NULL;
    for (i = 0; i < NPAWS*NEXTN; i++) 
	pawtabs[i].wcs = NULL;

    /* Read WCS, apcors etc from input table */

    (void)sprintf(catextn,"%s[1]",infile);
    retval = gettabinfo(catextn,&maintab,errmsg);
    if (retval != CIR_OK) {
        fprintf(stderr,"UNGROUT: Error reading input table %s -- %s\n",catextn,
		errmsg);
	exit(1);
    }

    /* Get PROV info and find the input files */

    retval = getinputinfo(infile,srcpath,verbose,pawcats,pawconfs,errmsg);
    if (retval != CIR_OK) {
	fprintf(stderr,"UNGROUT: Error finding source files %s\n",errmsg);
	exit(1);
    }

    /* Read object information from input table. First get some memory */

    (void)fits_get_num_rows(iptr,&nrows,&status);
    dspace = cir_malloc(4*nrows*sizeof(double));
    x = dspace;
    y = dspace + nrows;
    x2 = dspace + 2*nrows;
    y2 = dspace + 3*nrows;
    fspace = cir_calloc((3+3*NAPCOR)*nrows,sizeof(float));
    cnumb = fspace;
    pk = fspace + nrows;
    pkerr = fspace + 2*nrows;
    fluxcor = fspace + 3*nrows;
    fluxes = fspace + (3+NAPCOR)*nrows;
    fluxerrs = fspace + (3+2*NAPCOR)*nrows;

    /* Read the x,y coordinates, fluxes and flux errors. NB: flux + errs
       are stored by aperture rather than by object so you have to skip 
       over nrows entries to get the next value for a particular object. 
       Not particularly good visually, but helps with the cfitsio reads 
       and writes */
    
    (void)fits_read_col(iptr,TDOUBLE,COL_X,1,1,nrows,NULL,x,&anynul,&status);
    (void)fits_read_col(iptr,TDOUBLE,COL_Y,1,1,nrows,NULL,y,&anynul,&status);
    (void)fits_read_col(iptr,TFLOAT,COL_PEAKHEIGHT,1,1,nrows,NULL,pk,
			&anynul,&status);
    (void)fits_read_col(iptr,TFLOAT,COL_PKHTERR,1,1,nrows,NULL,pkerr,
			&anynul,&status);
    for (i = 0; i < NAPCOR; i++) {
	(void)fits_read_col(iptr,TFLOAT,COL_APFLUX1+2*i,1,1,nrows,NULL,
			    fluxes+i*nrows,&anynul,&status);
	(void)fits_read_col(iptr,TFLOAT,COL_APFLUX1ERR+2*i,1,1,nrows,NULL,
			    fluxerrs+i*nrows,&anynul,&status);
    }
    if (status != 0) {
	(void)fits_get_errstatus(status,msg);
	fprintf(stderr,"UNGROUT: Error reading input catalogue %s -- %s",
		infile,msg);
	closefits(iptr);
	freespace(fspace);
	freespace(dspace);
	return(CIR_FATAL);
    }

    /* Read WCS, apcors etc from paw tables */

    m = 0;
    for (i = 0; i < NPAWS; i++) {
	for (j = 0; j < NEXTN; j++) {
	    (void)sprintf(catextn,"%s[%d]",pawcats[i],j+1);
	    retval = gettabinfo(catextn,&(pawtabs[m]),errmsg);
	    if (retval != CIR_OK) {
		freespace(fspace);
		freespace(dspace);
		closefits(iptr);
		fprintf(stderr,"UNGROUT: Error reading paw table %s -- %s\n",
			catextn,errmsg);
		exit(1);
	    }
	    m++;
	}
    }

    /* Get rid of the output file if it already exists and make a new one */

    if (access(outfile,F_OK) == 0)
	remove(outfile);
    (void)fits_create_file(&optr,outfile,&status);
    (void)fits_copy_file(iptr,optr,1,1,1,&status);
    closefits(iptr);
    closefits(optr);
    if (status != 0) {
	fprintf(stderr,"UNGROUT: Error creating output table %s\n",outfile);
	freespace(fspace);
	freespace(dspace);
	exit(1);
    }

    /* Friendly message */

    if (verbose) 
	fprintf(stdout,"File: %s --> %s\n",infile,outfile);

    /* Loop for each stack confidence map image extension */

    m = -1;
    for (i = 0; i < NPAWS; i++) {
	(void)fits_open_file(&iptr,pawconfs[i],READONLY,&status);
	if (verbose) 
	    fprintf(stdout,"Doing %s",pawcats[i]);
	for (j = 0; j < NEXTN; j++) {
	    (void)fits_movabs_hdu(iptr,j+2,&hdutype,&status);
	    if (verbose)
		fprintf(stdout,".");
	    m++;
	    npts = pawtabs[m].nxout*pawtabs[m].nyout;
	    confmap = cir_malloc(npts*sizeof(short int));
	    (void)fits_read_img(iptr,TSHORT,1,npts,NULL,confmap,&anynul,
				&status);
	    delta_mag = maintab.magzpt - pawtabs[m].magzpt;

	    /* Shift the coordinates to the paw wcs */

	    cir_xytoxy_list(maintab.wcs,pawtabs[m].wcs,nrows,x,y,NULL,NULL,
			    x2,y2);

	    /* Loop through each object now and see if it's on this image
	       If it is, then work out the contribution of the flux correction
	       that you get from the current image */

	    xmax = (float)pawtabs[m].nxout - 0.5;
	    ymax = (float)pawtabs[m].nyout;
	    for (k = 0; k < nrows; k++) {
		if (x2[k] >= 0.5 && x2[k] < xmax && 
		    y2[k] >= 0.5 && y2[k] < ymax) {
		    ix = (int)(x2[k] + 0.5) - 1;
		    iy = (int)(y2[k] + 0.5) - 1;
		    ind = iy*pawtabs[m].nxout + ix;
		    confwt = (float)max(1.0,confmap[ind]);
		    cnumb[k] += confwt;
		    for (k2 = 0; k2 < NAPCOR; k2++) {
			expon = 0.4*(pawtabs[m].apcors[k2] - 
				     maintab.apcors[k2] + delta_mag);
			fluxcor[k2*nrows+k] += confwt*pow(10.0,expon);
		    }
		}
	    }
	    freespace(confmap);
	}
	closefits(iptr);
	if (verbose)
	    fprintf(stdout,"\n");
    }

    /* Now normalise the flux corrections and apply them */

    for (k = 0; k < nrows; k++) {
	for (k2 = 0; k2 < NAPCOR; k2++) {
	    if (cnumb[k] != 0.0) 
		scale = fluxcor[k2*nrows+k]/cnumb[k];
	    else
		scale = 1.0;
	    fluxes[k2*nrows+k] /= scale;
	    fluxerrs[k2*nrows+k] /= scale;
	    if (k2 == 0) {
		pk[k] /= scale;
	        pkerr[k] /= scale;
	    }
	}
    }

    /* Re-open the output file now and write the updated fluxes */

    (void)fits_open_file(&iptr,outfile,READWRITE,&status);
    (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
    k = 1;
    (void)fits_delete_key(iptr,"GROUTED",&status);
    (void)fits_delete_key(iptr,"MJD_DAY",&status);
    (void)fits_write_col(iptr,TFLOAT,COL_PEAKHEIGHT,1,1,nrows,pk,&status);
    (void)fits_write_col(iptr,TFLOAT,COL_PKHTERR,1,1,nrows,pkerr,&status);
    for (i = 0; i < NAPCOR; i++) {
	(void)fits_write_col(iptr,TFLOAT,COL_APFLUX1+2*i,1,1,nrows,
			    fluxes+i*nrows,&status);
	(void)fits_write_col(iptr,TFLOAT,COL_APFLUX1ERR+2*i,1,1,nrows,
			    fluxerrs+i*nrows,&status);
    }

    /* Delete the MJD column if it exists */

    (void)fits_get_colnum(iptr,CASEINSEN,"MJDoff",&res_col,&status);
    if (status != 0) {
	status = 0;
    } else {
	status = 0;   
	if (verbose) 
	    fprintf(stdout,"MJD column removed %d\n",res_col);
	(void)sprintf(colname,"TTYPE%d",res_col);
	(void)sprintf(sval,"Blank%d",res_col);
	(void)fits_update_key(iptr,TSTRING,colname,sval,NULL,&status);
	(void)sprintf(colname,"TUNIT%d",res_col);
	(void)fits_update_key(iptr,TSTRING,colname,"",NULL,&status);
	(void)memset(pk,0,nrows*sizeof(float));
	(void)fits_write_col(iptr,TFLOAT,res_col,1,1,nrows,pk,&status);
    }

    /* Delete some cards */

    (void)fits_delete_key(iptr,"PAWMAGZP",&status);
    status = 0;
    (void)fits_delete_key(iptr,"PAWELLPT",&status);
    status = 0;
    (void)fits_delete_key(iptr,"PAWSEENG",&status);
    status = 0;

    /* Free some memory */

    freespace(fspace);
    freespace(dspace);
    tidytabs(maintab,pawtabs);

    /* Stamp the primary */

    (void)fits_movabs_hdu(iptr,1,&hdutype,&status);
    retval = casu_stamp(iptr,"ungrout");
    closefits(iptr);
    if (verbose) 
	fprintf(stdout,"Done\n");

    /* Get out of here */

    exit(0);
}

static int gettabinfo(char *cat, tables *tab, char *errmsg) {
    int status,hdutype,i,ival;
    float val;
    double dval,cd11,cd12,cd21,cd22;
    fitsfile *iptr;
    struct wcsprm *wcs;
    char key[FLEN_KEYWORD];

    /* Open the file */

    status = 0;
    (void)fits_open_file(&iptr,cat,READONLY,&status);
    if (status != 0) {
	(void)sprintf(errmsg,"Unable to open table: %s\n",cat);
	return(CIR_FATAL);
    }
    (void)fits_get_hdu_type(iptr,&hdutype,&status);
    if (hdutype != BINARY_TBL) {
	(void)sprintf(errmsg,"FITS file %s is not a binary table\n",cat);
	closefits(iptr);
	return(CIR_FATAL);
    }

    /* Now get the WCS */

    if (cir_wcsopen(cat,&wcs,errmsg) != 0) {
	closefits(iptr);
	return(CIR_FATAL);
    }
    tab->wcs = wcs;

    /* Get the aperture corrections */

    for (i = 1; i <= NAPCOR; i++) {
	(void)snprintf(key,FLEN_KEYWORD,"APCOR%d",i);
	(void)fits_read_key(iptr,TFLOAT,key,&val,NULL,&status);
	if (status != 0) {
	    (void)sprintf(errmsg,"Keyword %s missing in %s\n",key,cat);
	    closefits(iptr);
	    return(CIR_FATAL);
	}
	tab->apcors[i-1] = val;
    }

    /* And the magnitude zero point */

    (void)fits_read_key(iptr,TFLOAT,"MAGZPT",&val,NULL,&status);
    if (status != 0) {
	(void)sprintf(errmsg,"Keyword MAGZPT missing in %s\n",cat);
	closefits(iptr);
	return(CIR_FATAL);
    }
    tab->magzpt = val;

    /* And the seeing */

    (void)fits_read_key(iptr,TFLOAT,"SEEING",&val,NULL,&status);
    if (status != 0) {
	(void)sprintf(errmsg,"Keyword SEEING missing in %s\n",cat);
	closefits(iptr);
	return(CIR_FATAL);
    }
    tab->seeing = val;

    /* Convert the seeing to arcseconds */

    (void)fits_read_key(iptr,TDOUBLE,"TC3_3",&cd11,NULL,&status);
    (void)fits_read_key(iptr,TDOUBLE,"TC3_5",&cd12,NULL,&status);
    (void)fits_read_key(iptr,TDOUBLE,"TC5_3",&cd21,NULL,&status);
    (void)fits_read_key(iptr,TDOUBLE,"TC5_5",&cd22,NULL,&status);
    if (status != 0) {
	(void)sprintf(errmsg,"CD matrix keyword missing in %s\n",cat);
	closefits(iptr);
	return(CIR_FATAL);
    }
    tab->seeing *= sqrt(cd11*cd22 - cd21*cd12)*3600.0;

    /* And the ellipticity */

    (void)fits_read_key(iptr,TFLOAT,"ELLIPTIC",&val,NULL,&status);
    if (status != 0) {
	(void)sprintf(errmsg,"Keyword ELLIPTIC missing in %s\n",cat);
	closefits(iptr);
	return(CIR_FATAL);
    }
    tab->ellipt = val;

    /* Finally the image size */

    (void)fits_read_key(iptr,TINT,"NXOUT",&ival,NULL,&status);
    if (status != 0) {
	(void)sprintf(errmsg,"Keyword NXOUT missing in %s\n",cat);
	closefits(iptr);
	return(CIR_FATAL);
    }
    tab->nxout = ival;
    (void)fits_read_key(iptr,TINT,"NYOUT",&ival,NULL,&status);
    if (status != 0) {
	(void)sprintf(errmsg,"Keyword NYOUT missing in %s\n",cat);
	closefits(iptr);
	return(CIR_FATAL);
    }
    tab->nyout = ival;
    
    /* Get out of here */

    closefits(iptr);
    return(CIR_OK);
}

static int getinputinfo(char *infile, char *srcpath, int verbose, 
			char pawcats[][BUFSIZ], char pawconfs[][BUFSIZ], 
			char *errmsg) {
    int status,i,ist,n,j,nerr,status1,status2,nxout,nyout,hdutype;
    long naxis[2];
    char key[FLEN_KEYWORD],prov[FLEN_VALUE],root[BUFSIZ],*extn_string;
    char msg[BUFSIZ],msg2[BUFSIZ];
    fitsfile *iptr,*confptr,*catptr;

    /* Read the provenance records in the input catalogue */

    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    nerr = 0;
    msg[0] = '\0';
    for (i = 0; i < NPAWS; i++) {
	(void)snprintf(key,FLEN_KEYWORD,"PROV%04d",i+1);
	(void)fits_read_key(iptr,TSTRING,key,prov,NULL,&status);
	if (status != 0) {
	    (void)sprintf(errmsg,"Unable to find keyword %s in %s\n",
			  key,infile);
	    closefits(iptr);
	    return(CIR_FATAL);
	}
	
	/* If this has been generated using nebulised files, then the
	   string 'tmp_' will appear. We want to ignore that bit... */

	ist = (strncmp(prov,"tmp_",4) ? 0 : 4);
	
	/* Mark the location of the file extension name */

	extn_string = strstr(prov,".fit");

	/* Now form the root name */

	n = extn_string - prov - ist;
	strncpy(root,prov+ist,n);
	root[n] = '\0';

	/* Ok, here are the input confidence map and catalogue names */

	(void)sprintf(pawcats[i],"%s/%s_cat.fits",srcpath,root);
	(void)sprintf(pawconfs[i],"%s/%s_conf%s",srcpath,root,extn_string);

	/* If this is verbose, then send a friendly message */

	if (verbose) {
	    if (i == 0) 
		fprintf(stdout,"Catalogue %s catalogue and conf maps: \n",
			infile);
	    fprintf(stdout,"    %s     %s\n",pawcats[i],pawconfs[i]);
	}

	/* Double check to make sure you can open both. Then check the
	   image dimensions */

	status1 = 0;
	(void)fits_open_file(&catptr,pawcats[i],READONLY,&status1);
	if (status1 != 0) {
	    (void)sprintf(msg2,"Can't open catalogue %s\n",pawcats[i]);
	    nerr++;
	    strcat(msg,msg2);
	}
	status2 = 0;
	(void)fits_open_file(&confptr,pawconfs[i],READONLY,&status2);
	if (status2 != 0) {
	    (void)sprintf(msg2,"Can't open conf map %s\n",pawconfs[i]);
	    nerr++;
	    strcat(msg,msg2);
	}
	if (status1 == 0 && status2 == 0) {
	    for (j = 1; j <= NEXTN; j++) {
		status1 = 0;
		(void)fits_movabs_hdu(catptr,j+1,&hdutype,&status1);
		if (status1 != 0) {
		    (void)sprintf(msg2,"Can't move to %s[%d]\n",pawcats[i],j);
		    nerr++;
		    strcat(msg,msg2);
		}
		status2 = 0;
		(void)fits_movabs_hdu(confptr,j+1,&hdutype,&status2);
		if (status2 != 0) {
		    (void)sprintf(msg2,"Can't move to %s[%d]\n",pawconfs[i],j);
		    nerr++;
		    strcat(msg,msg2);
		}
		if (status1 == 0 && status2 == 0) {
		    (void)fits_read_key(catptr,TINT,"NXOUT",&nxout,NULL,
					&status1);
		    (void)fits_read_key(catptr,TINT,"NYOUT",&nyout,NULL,
					&status1);
		    if (status1 != 0) {
			(void)sprintf(msg2,
				      "Missing NXOUT,NYOUT info in %s[%d]\n",
				      pawcats[i],j);
			nerr++;
			strcat(msg,msg2);
		    }
		    (void)fits_get_img_size(confptr,2,naxis,&status2);
		    if (status2 != 0) {
			(void)sprintf(msg2,"Couldn't get image size %s[%d]\n",
				      pawconfs[i],j);
			nerr++;
			strcat(msg,msg2);
		    }
		    if (status1 == 0 && status2 == 0 && (nxout != naxis[0] ||
							 nyout != naxis[1])) {
			(void)sprintf(msg2,"%s[%d] and %s[%d] don't match\n",
				      pawcats[i],j,pawconfs[i],j);
			nerr++;
			strcat(msg,msg2);
		    }
		}
	    }
	}
	closefits(confptr);
	closefits(catptr);
    }

    /* If there were any problems, then report them */

    closefits(iptr);
    if (nerr != 0) {
	sprintf(errmsg,"Problems with provenance: %s\n",msg);
	return(CIR_FATAL);
    } else {
	return(CIR_OK);
    }
}

static void tidytabs(tables maintab, tables *tabs) {
    int i;

    freewcs(maintab.wcs);
    for (i = 0; i < NEXTN*NPAWS; i++) 
	freewcs(tabs[i].wcs);
}

    
/* 

$Log: ungrout_main.c,v $
Revision 1.2  2012/12/05 10:50:56  jim
Added traps so that we just exit if the grouting or ungrouting is
unnecessary

Revision 1.1  2012-01-19 13:08:14  jim
New entry


*/
