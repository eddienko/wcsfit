/*

$Id: imcore_conf.c,v 1.4 2014/07/31 12:45:16 jim Exp $

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fitsio.h>
#include <tools.h>

#include "errcodes.h"
#include "ap.h"
#include "util.h"
#include "imcore.h"
#include "floatmath.h"

#define NW 5

static float *indata = NULL;
static short int *confdata = NULL;
static float *confsqrt = NULL;
static float *smoothed = NULL;
static float *smoothedc = NULL;
static unsigned char *mflag = NULL;
static fitsfile *iptr = NULL;
static ap_t ap;

static float weights[NW*NW];
static float weightc[NW*NW];
static long nx;
static long ny;

static void crweights(float);
static void convolve(int);
static float covariance(ap_t *, float);
static void tidy();

extern int imcore_conf(char *infile, char *conf, int ipix, float threshold,
		       int icrowd, float rcore, int nbsize, float filtfwhm,
		       char *outfile, char *ellfile, int verb, int cattyp, 
		       char *errmsg) {

    int i,retval,mulpix,j,status,nw2,hdunum;
    float fconst,nullval,skymed,skysig,thresh,xintmin,offset;
    float isat,isatbc,*current,junk,covcor;
    float *currentc;
    long npix;
    char errstr[ERRSTR_LEN],outextn[BUFSIZ];

    /* Useful constants */

    fconst = 1.0/M_LN2;
    nullval = 0.0;
    tptr = NULL;
    cattype = cattyp;
    verbose = verb;
    dribble = 0;

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
    if (mflag == NULL) {
        sprintf(errmsg,"unable to allocate workspace (mflag) in imcore_conf");
	tidy();
	return(ERRCODE_MALLOC);
    }    

    /* Open the ap structure and define some stuff in it */

    ap.lsiz = nx;
    ap.csiz = ny;
    apinit(&ap);
    ap.multiply = 1;
    ap.ipnop = ipix;
    ap.data = indata;
    ap.conf = confdata;
    ap.mflag = mflag;
    ap.rcore = rcore;
    ap.icrowd = icrowd;
    ap.fconst = fconst;
    ap.filtfwhm = filtfwhm;
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

    /* Compute a saturation level before background correction. */

    retval = imcore_backstats(&ap,nullval,1,&skymed,&skysig,&isatbc,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"imcore_backstats: %s",errstr);
	tidy();
	return(retval);
    }    

    /* Define saturated pixels */

    for (i = 0; i < npix; i++)
	if (mflag[i] == MF_CLEANPIX && indata[i] > isatbc)
	    mflag[i] = MF_SATURATED;

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

    /* Compute a saturation level after correction. */

    retval = imcore_backstats(&ap,nullval,1,&skymed,&skysig,&isat,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"imcore_backstats: %s",errstr);
	tidy();
	return(retval);
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

    /* Write the sky level and sigma to the header of the input file */

    status = 0;
    (void)fits_open_file(&iptr,infile,READWRITE,&status);
    (void)fits_update_key_fixflt(iptr,"SKYLEVEL",skymed,2,
                                 "Median sky brightness (counts/pixel)",
				 &status);
    (void)fits_update_key_fixflt(iptr,"SKYNOISE",skysig,2,
                                 "Pixel noise at sky level (counts)",&status);
    closefits(iptr);

    /* Take mean sky level out of data */

    for (i = 0; i < nx*ny; i++) {
	indata[i] -= skymed;
	if (indata[i]*confsqrt[i] > 3.0*skysig && 
	    mflag[i] != MF_SATURATED && mflag[i] != MF_STUPID_VALUE) 
	    mflag[i] = MF_3SIG;
    }

    /* Work out average covariance if the input map has been interpolated
       in some way that we know of... */

    if (dribble) {
	covcor = covariance(&ap,skysig);
	skysig = MAX(1.0,MIN(2.0,sqrt(covcor)))*skysig;
    }
    for (i = 0; i < nx*ny; i++)
	if (mflag[i] == MF_3SIG)
	    mflag[i] = MF_CLEANPIX;
    
    /* Work out isophotal detection threshold levels */

    thresh = threshold*skysig;
    if (verbose) 
	printf("\nSky level = %8.2f\nNoise level = %8.2f\nThreshold = %8.2f\n"
	       "Sat clips = %7.1f %7.1f\n",skymed,skysig,thresh,isat,isatbc);
    
    /* Minimum intensity for consideration */

    xintmin = 1.5*thresh*((float)ipix);

    /* Minimum size for considering multiple images */

    mulpix = MAX(8,2*ipix);

    /* Actual areal profile levels: T, 2T, 4T, 8T,...but written wrt T
       i.e. threshold as a power of 2 */

    offset = logf(thresh)*fconst;

    /* Get a bit of workspace for buffers */

    smoothed = malloc(nx*sizeof(*smoothed));
    smoothedc = malloc(nx*sizeof(*smoothedc));

    /* Define a few things more things in ap structure */

    ap.mulpix = mulpix;
    ap.areal_offset = offset;
    ap.thresh = thresh;
    ap.xintmin = xintmin;
    ap.sigma = skysig;
    ap.background = skymed;
    ap.saturation = isat;

    /* Set the weights */

    crweights(filtfwhm);
    nw2 = NW/2;

    /* Right, now for the extraction loop.  Begin by defining a group of
       three rows of data and confidence */

    if (verbose) 
	printf("Extracting images....\n");
    for (j = nw2; j < ny-nw2; j++) {
	current = indata + j*nx;
	currentc = confsqrt + j*nx;
	convolve(j);

        /* Do the detection now */

        apline(&ap,current,currentc,smoothed,smoothedc,j,NULL);

        /* Make sure we're not overruning the stacks */

        if (ap.ibstack > (ap.maxbl - ap.lsiz)) {
            if (verbose) 
		printf("APFU!: ibstack = %d\n",ap.ibstack);
	    apfu(&ap);
	}
	if (ap.ipstack > (ap.maxpa*3/4)) {
            if (verbose) 
		printf("APFU!: ipstack = %d\n",ap.ipstack);
            for (i = 0; i < ap.maxpa*3/8; i++)
		apfu(&ap);
	}

	/* See if there are any images to terminate */

	if (ap.ipstack > 1)
	    terminate(&ap,NULL);
    }

    /* Post process */

    retval = do_seeing(&ap,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"do_seeing: %s",errstr);
	tidy();
	return(retval);
    }
    (void)fits_get_hdu_num(tptr,&hdunum);
    (void)sprintf(outextn,"%s[%d]",outfile,hdunum-1);
    retval = tabclose(&ap,errstr);
    if (retval != ERRCODE_OK) {
        sprintf(errmsg,"tabclose: %s",errstr);
	tidy();
	return(retval);
    }
    cir_catcoord(infile,outextn,errstr);

    /* Tidy and exit */  

    tidy();
    return(ERRCODE_OK);
}

static void crweights(float filtfwhm) {
    int i,j,nw2,n;
    double gsigsq,di,dj;
    float renorm;

    /* Get the kernel size */

    nw2 = NW/2;
    
    /* Set the normalisation constants */

    gsigsq = 1.0/(2.0*pow(MAX(1.0,(double)filtfwhm)/2.35,2.0));
    renorm = 0.0;

    /* Now work out the weights */

    n = -1;
    for (i = -nw2; i <= nw2; i++) {
	di = (double)i;
	di *= gsigsq*di;
	for (j = -nw2; j <= nw2; j++) {
	    dj = (double)j;
	    dj *= gsigsq*dj;
	    n++;
	    weights[n] = (float)exp(-(di + dj));
	    renorm += weights[n];
	}
    }

    /* Now normalise the weights */

    n = -1;
    for (i = -nw2; i <= nw2; i++) {
	for (j = -nw2; j <= nw2; j++) {
	    n++;
	    weights[n] /= renorm;
	    /* weightc[n] = 0.01*weights[n]; */
	    weightc[n] = weights[n];
	}
    }
}

static void convolve(int ir) {
    int i,nw2,ix,jx,jy,n;
    float *idata,*cdata;

    /* Zero the summations */

    for (i = 0; i < nx; i++) {
	smoothed[i] = 0.0;
	smoothedc[i] = 0.0;
    }

    /* Now big is the smoothing kernel? */

    nw2 = NW/2;

    /* Now loop for each column */

    for (ix = nw2; ix < nx-nw2; ix++) {
	n = -1;
	for (jy = ir-nw2; jy <= ir+nw2; jy++) {
	    idata = indata + jy*nx;
	    cdata = confsqrt + jy*nx;
	    for (jx = ix-nw2; jx <= ix+nw2; jx++) {
		n++;
		smoothed[ix] += weights[n]*idata[jx];
		smoothedc[ix] += weightc[n]*idata[jx]*cdata[jx];
	    }
	}
    }
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
    freespace(smoothed);
    freespace(smoothedc);
    freespace(mflag);
    closefits(iptr);
    closefits(tptr);
    closefile(ellfp);
    apclose(&ap);
}

/* 

$Log: imcore_conf.c,v $
Revision 1.4  2014/07/31 12:45:16  jim
Modified so that nbsize <= 0 gives a constant background

Revision 1.3  2012/03/11 19:07:30  jim
Fixed a typo

Revision 1.2  2010-09-06 08:58:04  jim
NBSIZE added to ap structure

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.20  2010/02/11 21:54:46  jim
detection is done by weighting to sqrt of confidence

Revision 1.19  2009/12/17 11:36:07  jim
Makes correction to background noise for input images that have been
interpolated. This affects the detection threshold.

Revision 1.18  2009/01/22 13:50:39  jim
Add pixel flag to mark 3sigma points

Revision 1.17  2008/04/15 19:02:56  jim
Changed code to assign data flags

Revision 1.16  2007/12/18 15:23:17  jim
Fixed typo in comment

Revision 1.15  2006/07/31 13:21:09  jim
Modified imcore now allows for a smoothing kernel with variable FWHM

Revision 1.14  2006/07/24 11:41:25  jim
Fixed some problems with the background estimation

Revision 1.13  2006/07/06 12:25:27  jim
removed redundant declaration

Revision 1.12  2006/06/29 13:29:45  jim
Modifications to smoothing kernel

Revision 1.11  2006/06/26 15:13:51  jim
Background stats routines improved to deal with small noise estimates

Revision 1.10  2006/06/22 15:07:05  jim
Fixed to flag pixels that very low values

Revision 1.9  2006/06/05 11:21:03  jim
Fixed apline so that it takes confidence maps into account better

Revision 1.8  2005/09/23 08:32:03  jim
fixed swapped call to fits_close_file to closefits

Revision 1.7  2005/08/12 11:57:47  jim
Modified to write skylevel and skynoise estimates to input image header

Revision 1.6  2005/05/11 13:15:08  jim
Fixed bug in smoothing kernal where normalisation was being multiplied
instead of divided

Revision 1.5  2005/05/09 08:44:29  jim
Changed smoothing parameters depending on size of rcore

Revision 1.4  2005/03/06 19:41:20  jim
Removed unnecessary diagnostic things

Revision 1.3  2004/09/07 14:18:57  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.2  2004/04/05 11:25:42  jim
Small modifications and bug fixes

Revision 1.1  2004/04/02 10:54:58  jim
New version for rewrite of imcore


*/
