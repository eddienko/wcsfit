/*

$Id: grout_main.c,v 1.6 2012/12/05 10:50:56 jim Exp $

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
      NOVERB_ARG,
      MJDMAP_ARG};

static struct option myoptions [] = {
    {"srcpath",required_argument,NULL,SRCPATH_ARG},
    {"verbose",no_argument,NULL,VERB_ARG},
    {"noverbose",no_argument,NULL,NOVERB_ARG},
    {"mjdmap",required_argument,NULL,MJDMAP_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: grout incatalogue outcatalogue\n"
    "[--srcpath=%s] [--(no)verbose (%s)] [--mjdmap=%s]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define SRCPATH_DEF "."
#define VERB_DEF 0
#define MJDMAP_DEF ""

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
    float         mjdobs;
} tables;

static int gettabinfo(char *cat, tables *tab, char *errmsg); 
static int getmjd(char *pawcat, int extn, char *srcpath, double *mjd, 
		  char *errmsg);
static int getinputinfo(char *infile, char *srcpath, int verbose,
			char pawcats[][BUFSIZ], char pawconfs[][BUFSIZ], 
			char *errmsg);
static void tidytabs(tables maintab, tables *tabs);

int main (int argc, char *argv[]) {
    char *infile,*outfile,errmsg[BUFSIZ],pawcats[NPAWS][BUFSIZ],sval[16];
    char pawconfs[NPAWS][BUFSIZ],catextn[BUFSIZ],msg[BUFSIZ],colname[16];
    int c,option_index,nhdui,retval,m,i,j,status,hdutype,anynul;
    int npts,k,ix,iy,ind,k2,res_col,mjd_day;
    long nrows,naxes[2];
    short int *confmap;
    fitsfile *iptr=NULL,*optr=NULL;
    tables maintab,pawtabs[NPAWS*NEXTN];
    double *dspace,*x,*y,*x2,*y2,mjd,*coords;
    double *xin,*yin,*xout,*yout;
    float *fspace,*cnumb,*fluxes,*fluxerrs,*fluxcor,*pk,*pkerr,delta_mag;
    float confwt,expon,scale,*mjds,*mjdarray,*mjdw,xmax,ymax,work[NEXTN];

    /* Set up defaults for command line switches */

    char srcpath[BUFSIZ],mjdmap[BUFSIZ];
    (void)strcpy(srcpath,SRCPATH_DEF);
    int verbose = VERB_DEF;
    (void)strcpy(mjdmap,MJDMAP_DEF);

    /* First get the command line arguments */

    if (argc < 3) {
	fprintf(stderr,USAGE,srcpath,YESNO(verbose),MJDMAP_DEF);
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
	case MJDMAP_ARG:
	    (void)strcpy(mjdmap,optarg);
	    break;
	default:
	    fprintf(stderr,USAGE,SRCPATH_DEF,YESNO(VERB_DEF),MJDMAP_DEF);
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
	fprintf(stderr,"GROUT: Can't open file %s\n",infile);
	closefits(iptr);
	exit(1);
    } else if (nhdui != 2) {
	fprintf(stderr,"GROUT: Input file %s should have only 1 extension\n",
		infile);
	closefits(iptr);
	exit(1);
    }

    /* See if this has already been grouted. */

    (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
    (void)fits_read_key(iptr,TLOGICAL,"GROUTED",&k,NULL,&status);
    if (status != 0) {
	status = 0;
    } else if (k == 1) {
	if (verbose)
	    fprintf(stdout,"GROUT: File %s is already grouted\n",infile);
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
        fprintf(stderr,"GROUT: Error reading input table %s -- %s\n",catextn,
		errmsg);
	exit(1);
    }

    /* Get PROV info and find the input files */

    retval = getinputinfo(infile,srcpath,verbose,pawcats,pawconfs,errmsg);
    if (retval != CIR_OK) {
	fprintf(stderr,"GROUT: Error finding source files %s\n",errmsg);
	exit(1);
    }

    /* Read object information from input table. First get some memory */

    (void)fits_get_num_rows(iptr,&nrows,&status);
    dspace = cir_malloc(4*nrows*sizeof(double));
    x = dspace;
    y = dspace + nrows;
    x2 = dspace + 2*nrows;
    y2 = dspace + 3*nrows;
    fspace = cir_calloc((4+3*NAPCOR)*nrows,sizeof(float));
    cnumb = fspace;
    pk = fspace + nrows;
    pkerr = fspace + 2*nrows;
    mjds = fspace + 3*nrows;
    fluxcor = fspace + 4*nrows;
    fluxes = fspace + (4+NAPCOR)*nrows;
    fluxerrs = fspace + (4+2*NAPCOR)*nrows;

    /* If you are doing an mjd map, then get some space for the output 
       data array */

    mjdarray = NULL;
    if (mjdmap[0] != '\0') {
	mjdarray = cir_calloc(2*maintab.nxout*maintab.nyout,sizeof(float));
	mjdw = mjdarray + maintab.nxout*maintab.nyout;
    }

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
	fprintf(stderr,"GROUT: Error reading input catalogue %s -- %s",
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
		fprintf(stderr,"GROUT: Error reading paw table %s -- %s\n",
			catextn,errmsg);
		exit(1);
	    }
	    m++;
	}
    }

    /* Read MJD info. We do this for each extension in case there were
       extensions rejected during the original stacking */

    m = 0;
    for (i = 0; i < NPAWS; i++) {
	for (j = 0; j < NEXTN; j++) {
	    retval = getmjd(pawcats[i],j+1,srcpath,&mjd,errmsg);
	    if (retval != CIR_OK) {
		freespace(fspace);
		freespace(dspace);
		closefits(iptr);
		fprintf(stderr,"GROUT: Error getting MJD info from %s -- %s\n",
			pawcats[i],errmsg);
		exit(1);
	    }
	    if (i == 0 && j == 0) 
		mjd_day = (int)mjd;
	    pawtabs[m++].mjdobs = (float)(1440.0*(mjd - (double)mjd_day));
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
	fprintf(stderr,"GROUT: Error creating output table %s\n",outfile);
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
		    mjds[k] += confwt*pawtabs[m].mjdobs;
		    cnumb[k] += confwt;
		    for (k2 = 0; k2 < NAPCOR; k2++) {
			expon = 0.4*(pawtabs[m].apcors[k2] - 
				     maintab.apcors[k2] + delta_mag);
			fluxcor[k2*nrows+k] += confwt*pow(10.0,expon);
		    }
		}
	    }

	    /* If you are doing the mjd map, then get some space for 
	       the xy arrays and then do the calculation for the current
	       confidence map */

	    if (mjdmap[0] != '\0') {
		npts = pawtabs[m].nxout*pawtabs[m].nyout;
		coords = cir_malloc(4*npts*sizeof(double));
		xin = coords;
		yin = coords + npts;
		xout = coords + 2*npts;
		yout = coords + 3*npts;
		ind = -1;
		for (iy = 1; iy <= pawtabs[m].nyout; iy++) {
		    for (ix = 1; ix <= pawtabs[m].nxout; ix++) {
			ind++;
			xin[ind] = (double)ix;
			yin[ind] = (double)iy;
		    }
		}

		/* Now shift these coordinates to the tile reference frame */

		cir_xytoxy_list(pawtabs[m].wcs,maintab.wcs,npts,xin,yin,NULL,
				NULL,xout,yout);

		/* Do the calculation now */
		
		for (k = 0; k < npts; k++) {
		    if (xout[k] >= 0.5 && xout[k] <= maintab.nxout - 0.5 &&
			yout[k] >= 0.5 && yout[k] <= maintab.nyout - 0.5) {
			ix = (int)(xin[k] + 0.5) - 1;
			iy = (int)(yin[k] + 0.5) - 1;
			ind = iy*pawtabs[m].nxout + ix;
		        confwt = (float)max(1.0,confmap[ind]);
			ix = (int)(xout[k] + 0.5) - 1;
			iy = (int)(yout[k] + 0.5) - 1;
			ind = iy*maintab.nxout + ix;
			mjdarray[ind] += confwt*pawtabs[m].mjdobs;
			mjdw[ind] += confwt;
		    }
		}
		freespace(coords);
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
	    fluxes[k2*nrows+k] *= scale;
	    fluxerrs[k2*nrows+k] *= scale;
	    if (k2 == 0) {
		pk[k] *= scale;
	        pkerr[k] *= scale;
	    }
	}
	if (cnumb[k] != 0.0)
	    mjds[k] /= cnumb[k];
    }

    /* Re-open the output file now and write the updated fluxes */

    (void)fits_open_file(&iptr,outfile,READWRITE,&status);
    (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
    k = 1;
    (void)fits_update_key(iptr,TLOGICAL,"GROUTED",&k,"Table has been grouted",
			  &status);
    (void)fits_update_key(iptr,TINT,"MJD_DAY",&mjd_day,"MJD day",&status);
    (void)fits_write_col(iptr,TFLOAT,COL_PEAKHEIGHT,1,1,nrows,pk,&status);
    (void)fits_write_col(iptr,TFLOAT,COL_PKHTERR,1,1,nrows,pkerr,&status);
    for (i = 0; i < NAPCOR; i++) {
	(void)fits_write_col(iptr,TFLOAT,COL_APFLUX1+2*i,1,1,nrows,
			    fluxes+i*nrows,&status);
	(void)fits_write_col(iptr,TFLOAT,COL_APFLUX1ERR+2*i,1,1,nrows,
			    fluxerrs+i*nrows,&status);
    }

    /* Find a blank column and write out the mjd info */

    (void)fits_get_colnum(iptr,CASEINSEN,"Blank#",&res_col,&status);
    if (status != 0 && status != COL_NOT_UNIQUE) {
	fprintf(stderr,"GROUT: Can't find a blank column in %s for MJD info\n",
		outfile);
	status = 0;
    } else {
	status = 0;   
	if (verbose) 
	    fprintf(stdout,"MJD information written to column %d\n",res_col);
	(void)sprintf(colname,"TTYPE%d",res_col);
	(void)fits_update_key(iptr,TSTRING,colname,"MJDoff",NULL,&status);
	(void)sprintf(colname,"TUNIT%d",res_col);
	(void)fits_update_key(iptr,TSTRING,colname,"Minutes",NULL,&status);
	(void)fits_write_col(iptr,TFLOAT,res_col,1,1,nrows,mjds,&status);
    }

    /* Add a card with the median magzpt for each of the pawprints */

    m = 0;
    (void)fits_movabs_hdu(iptr,1,&hdutype,&status);
    msg[0] = '\0';
    for (i = 0; i < NPAWS; i++) {
	for (j = 0; j < NEXTN; j++) 
	    work[j] = pawtabs[m++].magzpt;
	(void)cir_med(work,NULL,NEXTN,&delta_mag,errmsg);
	(void)sprintf(sval,"%0.3f ",delta_mag);
	(void)strcat(msg,sval);
    }
    (void)fits_update_key(iptr,TSTRING,"PAWMAGZP",msg,"Median paw MAGZPT",
			  &status);

    /* Add a card with the median ellipticity for each of the pawprints */

    m = 0;
    msg[0] = '\0';
    for (i = 0; i < NPAWS; i++) {
	for (j = 0; j < NEXTN; j++) 
	    work[j] = pawtabs[m++].ellipt;
	(void)cir_med(work,NULL,NEXTN,&delta_mag,errmsg);
	(void)sprintf(sval,"%0.3f ",delta_mag);
	(void)strcat(msg,sval);
    }
    (void)fits_update_key(iptr,TSTRING,"PAWELLPT",msg,
			  "Median paw ELLIPTIC",&status);

    /* Add a card with the median seeing for each of the pawprints */

    m = 0;
    msg[0] = '\0';
    for (i = 0; i < NPAWS; i++) {
	for (j = 0; j < NEXTN; j++) 
	    work[j] = pawtabs[m++].seeing;
	(void)cir_med(work,NULL,NEXTN,&delta_mag,errmsg);
	(void)sprintf(sval,"%0.3f ",delta_mag);
	(void)strcat(msg,sval);
    }
    (void)fits_update_key(iptr,TSTRING,"PAWSEENG",msg,
			  "[arcsec] Median paw SEEING",&status);

    /* If writing out the mjd map, the do that now */

    if (mjdmap[0] != '\0') {
	naxes[0] = maintab.nxout;
	naxes[1] = maintab.nyout;
	npts = naxes[0]*naxes[1];
	for (k = 0; k < npts; k++) 
	    if (mjdw[k] != 0)
		mjdarray[k] /= mjdw[k];
	    else
		mjdarray[k] = 0.0;
	retval = cir_open_output(mjdmap,pawcats[0],&optr,NEWDATASZ,FLOAT_IMG,
				 2,naxes,errmsg);
	if (retval != CIR_OK) {
	    fprintf(stderr,"GROUT: Unable to create output MJD file %s -- %s\n",
		    mjdmap,errmsg);
	} else {
	    (void)fits_write_img(optr,TFLOAT,1,npts,mjdarray,&status);
	    closefits(optr);
	}
	freespace(mjdarray);
    }

    /* Free some memory */

    freespace(fspace);
    freespace(dspace);
    tidytabs(maintab,pawtabs);

    /* Stamp the primary */

    (void)fits_movabs_hdu(iptr,1,&hdutype,&status);
    retval = casu_stamp(iptr,"grout");
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

static int getmjd(char *pawcat, int extn, char *srcpath, double *mjd, 
		  char *errmsg) {
    int status,hdutype,ip;
    double mjdsum;
    fitsfile *iptr,*pptr;
    char key[FLEN_KEYWORD],prov[FLEN_VALUE],fname[BUFSIZ];

    /* Read the prov information from the correct extension of the catalogue */

    status = 0;
    (void)fits_open_file(&iptr,pawcat,READONLY,&status);
    (void)fits_movabs_hdu(iptr,extn+1,&hdutype,&status);
    mjdsum = 0.0;
    ip = 1;
    while (1) {
	(void)snprintf(key,FLEN_KEYWORD,"PROV%04d",ip);
	(void)fits_read_key(iptr,TSTRING,key,prov,NULL,&status);
	if (status != 0) {
	    status = 0;
	    break;
	}
	(void)snprintf(fname,BUFSIZ,"%s/%s",srcpath,prov);
	(void)fits_open_file(&pptr,fname,READONLY,&status);
	if (status != 0) {
	    (void)sprintf(errmsg,"Unable to get MJD info from %s\n",prov);
	    closefits(pptr);
	    closefits(iptr);
	    return(CIR_FATAL);
	}
	(void)fits_read_key(pptr,TDOUBLE,"MJD-OBS",mjd,NULL,&status);
	if (status != 0) {
	    (void)sprintf(errmsg,"Unable to get MJD info from %s\n",prov);
	    closefits(pptr);
	    closefits(iptr);
	    return(CIR_FATAL);
	}
	closefits(pptr);
	mjdsum += *mjd;
	ip++;
    }
    closefits(iptr);
    ip--;

    /* Was there any prov info? */

    if (ip == 0) {
	(void)sprintf(errmsg,"No provenance info in %s[%d]\n",pawcat,extn);
	return(CIR_FATAL);
    }

    /* Normalise the sum and get out of here */

    *mjd = mjdsum/(double)ip;
    return(CIR_OK);
}
	
static void tidytabs(tables maintab, tables *tabs) {
    int i;

    freewcs(maintab.wcs);
    for (i = 0; i < NEXTN*NPAWS; i++) 
	freewcs(tabs[i].wcs);
}

    
/* 

$Log: grout_main.c,v $
Revision 1.6  2012/12/05 10:50:56  jim
Added traps so that we just exit if the grouting or ungrouting is
unnecessary

Revision 1.5  2011-06-09 12:10:57  jim
Adds pawprint QC information into the header now

Revision 1.4  2011-02-11 16:09:07  jim
Fixed memory bug

Revision 1.3  2011-01-31 15:03:44  jim
now adds an MJD column and writes an MJD map if requested

Revision 1.2  2011-01-12 13:12:05  jim
Fixed bug where objects were falling off the edges of stacks

Revision 1.1  2010-11-22 10:43:47  jim
new entry


*/
