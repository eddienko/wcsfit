/*

$Id: imcore_opm.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>

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
static float *incopy = NULL;
static short int *ccopy = NULL;

static float weights[NW*NW];
static float weightc[NW*NW];
static long nx;
static long ny;

static void crweights(float);
static void convolve(int);
static void tidy();


extern int imcore_opm(char *infile, char *conf, int ipix, float threshold,
		      int nbsize, float filtfwhm, char *outfile, int verb, 
		      int niter, char *errmsg) {

    int i,retval,j,status,nw2,iter,nclip;
    float fconst,nullval,skymed,skysig,thresh,xintmin,offset;
    float isat,isatbc,*current,junk;
    float *currentc;
    long npix;
    char errstr[ERRSTR_LEN];

    /* Useful constants */

    fconst = 1.0/M_LN2;
    nullval = 0.0;
    tptr = NULL;
    cattype = 4;
    verbose = verb;

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

    /* Make a copy of each */

    npix = nx*ny;
    incopy = malloc(npix*sizeof(*incopy));
    ccopy = malloc(npix*sizeof(*ccopy));
    memmove(incopy,indata,npix*sizeof(*incopy));
    memmove(ccopy,confdata,npix*sizeof(*ccopy));

    /* Get mflag array for flagging saturated pixels */

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
    ap.ipnop = ipix;
    ap.data = indata;
    ap.conf = confdata;
    ap.mflag = mflag;
    ap.fconst = fconst;
    ap.filtfwhm = filtfwhm;

    /* Open the output catalogue FITS table */

    retval = tabinit(&ap,infile,outfile,errstr);
    if (retval != ERRCODE_OK) {
	sprintf(errmsg,"tabinit: %s",errstr);
	tidy();
	return(retval);
    }

    /* Set up the data flags */

    for (i = 0; i < npix ; i++) 
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

    /* Flag saturated pixels */

    for (i = 0; i < npix ; i++) 
	if (mflag[i] == MF_CLEANPIX && indata[i] > isatbc)
	    mflag[i] = MF_SATURATED;

    /* Get a bit of workspace for buffers */

    smoothed = malloc(nx*sizeof(*smoothed));
    smoothedc = malloc(nx*sizeof(*smoothedc));

    /* Set the weights */

    crweights(filtfwhm);
    nw2 = NW/2;

    /* Iteration loop */

    for (iter = 0; iter < niter; iter++) {

        /* Compute the background variation and remove it from the data*/

        if (verbose)
	    printf("Computing background....\n");
	retval = imcore_background(&ap,nbsize,nullval,verbose,errstr);
	if (retval != ERRCODE_OK) {
	    sprintf(errmsg,"imcore_background: %s",errstr);
	    tidy();
	    return(retval);
	}

	/* Compute a saturation level before background correction. */

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
				     "Pixel noise at sky level (counts)",
				     &status);
	closefits(iptr);

	/* Take mean sky level out of data */

	for (i = 0; i < nx*ny; i++)
	    indata[i] -= skymed;

	/* Work out isophotal detection threshold levels */

	thresh = threshold*skysig;
	if (verbose) 
	    printf("\nSky level = %8.2f\nNoise level = %8.2f\nThreshold = %8.2f\n"
		   "Sat clips = %7.1f %7.1f\n",skymed,skysig,thresh,isat,isatbc);

	/* Minimum intensity for consideration */

	xintmin = 1.5*thresh*((float)ipix);

	/* Actual areal profile levels: T, 2T, 4T, 8T,...but written wrt T
	   i.e. threshold as a power of 2 */

	offset = logf(thresh)*fconst;

	/* Define a few things more things in ap structure */

	ap.areal_offset = offset;
	ap.thresh = thresh;
	ap.xintmin = xintmin;
	ap.sigma = skysig;
	ap.background = skymed;
	ap.saturation = isat;

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
    
	/* Restore input data to its former glory. Update confidence map */
    
	memmove(indata,incopy,npix*sizeof(*indata));
	nclip = 0;
	for (i = 0; i < npix; i++) {
	    if ((ap.opm)[i]) {
		confdata[i] = 0;
		(ap.opm)[i] = 0;
		nclip++;
	    }
	}
	if (nclip == 0)
	    break;
    }
    for (i = 0; i < npix; i++)
	(ap.opm)[i] = (confdata[i] == 0);
    memmove(confdata,ccopy,npix*sizeof(*confdata));
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
/* 	    weightc[n] = 0.01*weights[n]; */
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


static void tidy() {
    int status = 0;

    freespace(indata);
    freespace(confdata);
    freespace(confsqrt);
    freespace(smoothed);
    freespace(smoothedc);
    freespace(mflag);
    freespace(incopy);
    freespace(ccopy);
    closefits(iptr);
    closefits(tptr);
    closefile(ellfp);
    apclose(&ap);
}

/* 

$Log: imcore_opm.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.5  2010/03/10 10:58:38  jim
Does detection based on sqrt(conf)

Revision 1.4  2010/02/11 21:55:33  jim
changed a few routine declarations

Revision 1.3  2008/04/25 12:22:47  jim
Modified to take out pixels with zero confidence from the calculation
of the number of object pixels

Revision 1.2  2008/04/15 19:02:56  jim
Changed code to assign data flags

Revision 1.1  2007/07/31 12:04:53  jim
new routine


*/
