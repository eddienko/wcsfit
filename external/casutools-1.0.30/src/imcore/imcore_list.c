/*

$Id: imcore_list.c,v 1.4 2014/08/08 09:21:38 jim Exp $

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fitsio.h>
#include <casu_wcs.h>

#include "errcodes.h"
#include "ap.h"
#include "util.h"
#include "imcore.h"
#include "floatmath.h"


static float *indata = NULL;
static short int *confdata = NULL;
static float *confsqrt = NULL;
static unsigned char *mflag = NULL;
static ap_t ap;
static objstruct *ob = NULL;

static int nobjects;

static int imcore_rdlist(char *, char *, char *, char *);
static float covariance(ap_t *, float);
static void tidy();

extern int imcore_list(char *infile, char *conf, char *listfile, char *trans, 
		       float threshold, float rcore, int nbsize, char *outfile,
		       char *ellfile, int verb, int cattyp, char *errmsg) {

    int i,retval,nbit,j;
    float fconst,nullval,skymed,skysig,thresh,offset,junk,covcor;
    long npix,nx,ny;
    char errstr[ERRSTR_LEN];
    objstruct *oblist[IMNUM];

    /* Useful constants */

    fconst = 1.0/M_LN2;
    nullval = 0.0;
    tptr = NULL;
    cattype = cattyp;
    verbose = verb;
    dribble = 0;

    /* Open the list */

    errstr[0] = '\0';
    retval = imcore_rdlist(listfile,infile,trans,errstr);
    if (retval != ERRCODE_OK) {
	sprintf(errmsg,"imcore_rdlist: %s",errstr);
	tidy();
	return(retval);
    }

    /* Open input image and associated confidence map */

    retval = imcore_rdbuf_mef(infile,TFLOAT,(void *)&indata,&nx,&ny,verbose,
			      errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"imcore_rdbuf_mef: %s",errstr);
	tidy();
	return(retval);
    }
    retval = imcore_rdbuf_conf(conf,&confdata,&confsqrt,nx,ny,verbose,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"imcore_rdbuf_conf: %s",errstr);
	tidy();
	return(retval);
    }

    /* Get mflag array for flagging saturated pixels */

    npix = nx*ny;
    mflag = calloc(npix,sizeof(*mflag));
    if(mflag == NULL) {
        sprintf(errmsg,"unable to allocate workspace (mflag) in imcore_conf");
	tidy();
	return(ERRCODE_MALLOC);
    }    

    /* Open the ap structure and define some stuff in it */

    ap.lsiz = nx;
    ap.csiz = ny;
    apinit(&ap);
    ap.multiply = 1;
    ap.ipnop = 0;
    ap.data = indata;
    ap.conf = confdata;
    ap.mflag = mflag;
    ap.rcore = rcore;
    ap.icrowd = 0;
    ap.fconst = fconst;
    ap.filtfwhm = 0.0;
    ap.nbsize = nbsize;

    /* Open the output catalogue FITS table */

    retval = tabinit(&ap,infile,outfile,errstr);
    if (retval != ERRCODE_OK) {
	sprintf(errmsg,"tabinit: %s",errstr);
	tidy();
	return(retval);
    }

    /* Set up ellipse drawing file for ds9.  It's not an error if you can't
       open one */

    ellfp = fopen(ellfile,"w");

    /* Set up the data flags */

    for (i = 0; i < npix; i++) 
	if (confdata[i] == 0)
	    mflag[i] = MF_ZEROCONF;
        else if (indata[i] < STUPID_VALUE)
	    mflag[i] = MF_STUPID_VALUE;
	else 
	    mflag[i] = MF_CLEANPIX;

    /* Compute the background variation and remove it from the data*/

    if (nbsize > 0) {
        if (verbose)
            printf("Computing background....\n");
        retval = imcore_background(&ap,nbsize,nullval,verbose,errstr);
        if (retval != ERRCODE_OK) {
            sprintf(errmsg,"imcore_background: %s",errstr);
            tidy();
            return(retval);
        }
    }

    /* Compute background statistics */

    if (verbose) 
        printf("Computing background statistics....\n");
    retval = imcore_backstats(&ap,nullval,0,&skymed,&skysig,&junk,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"imcore_backstats: %s",errstr);
	tidy();
	return(retval);
    }

    /* Take mean sky level out of data */

    for (i = 0; i < nx*ny; i++)
	indata[i] -= skymed;
    
    /* Work out average covariance if the input map has been interpolated
       in some way that we know of... */

    if (dribble) {
	covcor = covariance(&ap,skysig);
	skysig = MAX(1.0,MIN(2.0,sqrt(covcor)))*skysig;
    }

    /* Work out isophotal detection threshold levels */

    thresh = threshold*skysig;
    if (verbose) 
	printf("\nSky level = %8.2f\nNoise level = %8.2f\nThreshold = %8.2f\n",
	       skymed,skysig,thresh);
    
    /* Actual areal profile levels: T, 2T, 4T, 8T,...but written wrt T
       i.e. threshold as a power of 2 */

    offset = logf(thresh)*fconst;

    /* Define a few things more things in ap structure */

    ap.mulpix = 0;
    ap.areal_offset = offset;
    ap.thresh = thresh;
    ap.xintmin = 0;
    ap.sigma = skysig;
    ap.background = skymed;
    ap.saturation = 0;

    /* Right, now for the extraction loop.  Begin by defining a group of
       three rows of data and confidence */

    if (verbose) 
	printf("Extracting images....\n");
    nbit = 0;
    for (j = 0; j < nobjects; j++) {

	/* Add the current object into the list */

        oblist[nbit++] = ob + j;

	/* If this is the last object in the input list, or this image
	   is not part of a blend, then terminate */

	if (j == nobjects - 1 || ((ob+j)->areal7 != -1.0)) {
	    process_results_list(&ap,oblist,nbit,errstr);
	    nbit = 0;

        /* If this is part of a blend and the next object isn't part of 
	   the same blend, then terminate */

	} else if (((ob+j+1)->areal8 == 0.0) || ((ob+j+1)->areal7 != -1.0)) {
	    process_results_list(&ap,oblist,nbit,errstr);
	    nbit = 0;
	}
    }

    /* Post process */

    retval = do_seeing(&ap,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"do_seeing: %s",errstr);
	tidy();
	return(retval);
    }
    retval = tabclose(&ap,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"tabclose: %s",errstr);
	tidy();
	return(retval);
    }

    /* Tidy and exit */  

    tidy();
    return(ERRCODE_OK);
}

static int imcore_rdlist(char *listfile, char *image, char *trans, 
			 char *errstr) {
    fitsfile *fptr;
    int status,hdutype,icat,nr,nalloc,ialloc=256,i,ncol,ctype,tdone;
    double xx,yy,ra,dec,xt[3],yt[3];
    struct wcsprm *wcs;
    FILE *fd;
    char lineofdata[BUFSIZ];

    /* Right, first of all see what sort of list we have here. Try and 
       open it as a FITS file. If it is a FITS file, then see if it is
       a binary or ascii table */

    status = 0;
    nobjects = 0;
    (void)fits_open_file(&fptr,listfile,READONLY,&status);
    if (status == 0) {
	(void)fits_get_hdu_type(fptr,&hdutype,&status);

	/* If this is some sort of FITS table, then determine whether this
	   is a catalogue as generated by imcore */

	if (hdutype == ASCII_TBL || hdutype == BINARY_TBL) {
            (void)fits_read_key(fptr,TINT,"CATTYPE",&ctype,NULL,&status);
            if (status != 0) {
                status = 0;
                (void)fits_get_num_cols(fptr,&ncol,&status);
                ctype = -1;
                for (icat = 1; icat < 4; icat++) {
                    if (get_ncol(icat) == ncol) {
                        ctype = icat;
                        break;
                    }
                }
            }
	    closefits(fptr);
            if (ctype < 0) {
		(void)sprintf(errstr,"Unrecognised catalogue type");
		return(ERRCODE_FILE_IO);
	    }

	    /* If it's the right sort of catalogue, then read it */

	    if (readtab_list(listfile,&ob,&nr,ctype,errstr) != ERRCODE_OK) 
		return(ERRCODE_FILE_IO);

	    /* See if there is a transformation file defined. This is
	       used to transform the x,y coordinates in the table to 
	       the x,y coordinates in the image */

	    tdone = 1;
	    if ((fd = fopen(trans,"r")) != NULL) {
		if ((fscanf(fd,"%lg %lg %lg\n",&(xt[0]),&(xt[1]),&(xt[2])) != 3) ||
		    (fscanf(fd,"%lg %lg %lg\n",&(yt[0]),&(yt[1]),&(yt[2])) != 3)) {
		    tdone = 0;
		} else {
		    for (i = 0; i < nr; i++) {
			xx = (ob+i)->x_t;
			yy = (ob+i)->y_t;
			(ob+i)->x_i = (float)(xt[0]*xx + xt[1]*yy + xt[2]);
			(ob+i)->y_i = (float)(yt[0]*xx + yt[1]*yy + yt[2]);
		    }
		}
		fclose(fd);
	    } else {
		tdone = 0;
	    }
			
	    /* Get a WCS from the table header if it exists. If the WCS does 
	       exist, then convert the x,y values to RA,Dec. This gives a 
	       better value than using the single precision values from the 
	       table */

	    if (cir_wcsopen(listfile,&wcs,errstr) == 0) {
		for (i = 0; i < nr; i++) {
		    xx = (double)((ob+i)->x_t);
		    yy = (double)((ob+i)->y_t);
		    cir_xytoradec(wcs,xx,yy,&ra,&dec);
		    (ob+i)->ra = ra;
		    (ob+i)->dec = dec;
		}
		cir_wcsclose(wcs);
	    }

  	    /* If the transformation hasn't already been done by a trans file
	       then use the image WCS to transform the coordinates. If that
	       fails then just copy the table coordinates and hope for the
	       best */

	    if (tdone == 0) {
		if (cir_wcsopen(image,&wcs,errstr) == 0) {
		    for (i = 0; i < nr; i++) {
			ra = (ob+i)->ra;
			dec = (ob+i)->dec;
			cir_radectoxy(wcs,ra,dec,&xx,&yy);
			(ob+i)->x_i = (float)xx;
			(ob+i)->y_i = (float)yy;
		    }
		    cir_wcsclose(wcs);
		} else {
		    for (i = 0; i < nr; i++) {
			(ob+i)->x_i = (ob+i)->x_t;
			(ob+i)->y_i = (ob+i)->y_t;
		    }
		}
	    }
	} else {
	    closefits(fptr);
	    (void)sprintf(errstr,"List file is not a FITS table");
	    return(ERRCODE_FILE_IO);
	}

    /* Ok, this isn't a FITS table. Assume it's a text file and 
       try to read the RA and Dec values from it. Start by trying to
       open the file */

    } else {
	if ((fd = fopen(listfile,"r")) == NULL) {
	    (void)sprintf(errstr,"Unable to open list file text file %s",
			  listfile);
	    return(ERRCODE_FILE_IO);
	}

	/* Open the image WCS */

	if (cir_wcsopen(image,&wcs,errstr) != 0) {
	    (void)sprintf(errstr,"Unable to get image WCS: %s",image);
	    fclose(fd);
	    return(ERRCODE_FILE_IO);
	}

	/* Get some workspace */

	nalloc = ialloc;
	ob = malloc(nalloc*sizeof(objstruct));
	nr = 0;
	while (fscanf(fd,"%[^\n] ",lineofdata) != EOF) {
	    sscanf(lineofdata,"%lg %lg\n",&ra,&dec);
	    if (nr == nalloc) {
		nalloc += ialloc;
		ob = realloc(ob,nalloc*sizeof(objstruct));
	    }
	    (ob+nr)->ra = ra;
	    (ob+nr)->dec = dec;
	    (ob+nr)->x_t = 0.0;
	    (ob+nr)->y_t = 0.0;
	    (ob+nr)->class = 0.0;
	    (ob+nr)->stat = 0.0;
	    (ob+nr)->areal7 = 0.0;
	    (ob+nr)->areal8 = 0.0;
	    (ob+nr)->coreflux = 0.0;
	    cir_radectoxy(wcs,ra,dec,&xx,&yy);
	    (ob+nr)->x_i = (float)xx;
	    (ob+nr)->y_i = (float)yy;
	    nr++;
	}
	fclose(fd);
	ob = realloc(ob,nr*sizeof(objstruct));
	cir_wcsclose(wcs);
    }
    nobjects = nr;

    /* Get out of here */

    return(ERRCODE_OK);
}
	
static float covariance(ap_t *ap, float skysig) {
    float *map,sum,clip,*autoc,sumcor,av,norm;
    int nx,ny,ii,iarg,i,j,nx1,nx2,ny1,ny2,jarg,ncorr=5,nco;
    int jauto,iyauto,ixauto,jargoff,iargoff;
    unsigned char *mf;

    /* Set some convenience variables */
    
    map = ap->data;
    mf = ap->mflag;
    nx = ap->lsiz;
    ny = ap->csiz;
    nx1 = nx/4;
    nx2 = nx - nx1;
    ny1 = ny/4;
    ny2 = ny - ny1;
    clip = 3.0*skysig;
    nco = ncorr/2;

    /* Work out what the average is well away from the edges */
	
    sum = 0.0;
    ii = 0;
    for (j = ny1; j < ny2; j++) {
	jarg = j*nx;
	for (i = nx1; i < nx2; i++) {
	    iarg = jarg + i;
	    if (map[iarg] < clip && mf[iarg] == MF_CLEANPIX) {
		sum += map[iarg];
		ii++;
            }
	}
    }
    av = sum/(float)ii;

    /* Get workspace for the auto-correlation array */

    autoc = calloc(ncorr*ncorr,sizeof(*autoc));
    
    /* Do the correlation matrix */

    jauto = -1;
    for (iyauto = -nco; iyauto <= nco; iyauto++) {
	for (ixauto = -nco; ixauto <= nco; ixauto++) {
	    jauto++;
	    ii = 0;
	    sum = 0.0;
	    for (j = ny1; j < ny2; j++) {
		jarg = j*nx;
		jargoff = (j+iyauto)*nx;
		for (i = nx1; i < nx2; i++) {
		    iarg = jarg + i;
		    iargoff = jargoff + i + ixauto;
		    if (map[iarg] - av < clip && mf[iarg] == MF_CLEANPIX &&
			map[iargoff] - av < clip && 
			mf[iargoff] == MF_CLEANPIX) {
			ii++;
			sum += (map[iarg]-av)*(map[iargoff]-av);
		    }
		}
	    }
	    autoc[jauto] = sum/(float)ii;
	}
    }
    
    /* Normalise the array by the central pixel */

    norm = autoc[nco*ncorr+nco];
    sumcor = 0.0;
    jauto = -1;
    for (iyauto = -nco; iyauto <= nco; iyauto++) {
	for (ixauto = -nco; ixauto <= nco; ixauto++) {
	    jauto++;
	    autoc[jauto] /= norm;
	    if (abs(ixauto) == nco || abs(iyauto) == nco)
		continue;
	    sumcor += autoc[jauto];
	}
    }
    
    /* Get out of here */

    free(autoc);
    return(sumcor);
}   

static void tidy() {
    int status = 0;

    freespace(indata);
    freespace(confdata);
    freespace(confsqrt);
    freespace(mflag);
    closefits(tptr);
    closefile(ellfp);
    apclose(&ap);
    freespace(ob);
}

/* 

$Log: imcore_list.c,v $
Revision 1.4  2014/08/08 09:21:38  jim
Added nbsize <= 0 option

Revision 1.3  2014/07/30 08:12:45  jim
Fixed long standing bug that means that in list mode the output equatorial
coordinates are written to the FITS tables in degrees rather than in radians
as advertised

Revision 1.2  2010/09/06 08:58:32  jim
Added NBSIZE to ap structure

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.4  2010/02/11 21:55:33  jim
changed a few routine declarations

Revision 1.3  2008/04/15 19:04:00  jim
Modifed the way that data values are flagged

Revision 1.2  2007/06/14 06:50:26  jim
Fixed bug in WCS reading
CVI: ----------------------------------------------------------------------

Revision 1.1  2007/06/04 10:35:07  jim
Initial entry


*/
