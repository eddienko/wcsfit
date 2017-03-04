/*

$Id: cir_vstpickup.c,v 1.1 2012/12/08 07:31:00 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <tools.h>

typedef struct {
    float *data;
    float sky;
    float skynoise;
    float skynorm;
} dstrct;

static fitsfile *iptr = NULL;
static fitsfile *fptr = NULL;
static dstrct *fileinfo = NULL;
static dstrct *medianinfo = NULL;

static long naxis[2];

static void skyest(float *data, long npts, float thresh, float *skymed, 
		   float *skynoise);
static void average(void);
static void sort1(float *a, int n);
static void tidy();

#define DEFNEXTN 32
#define DEFNX_TRIMMED 2048
#define DEFNY_TRIMMED 4100

#define SCMIN 0.5
#define SCMAX 1.5
#define SCGAP 0.1
#define NSCALE ((int)((SCMAX - SCMIN)/SCGAP) + 1)

extern int cir_vstpickup(char *infile, char *flat, int checkonly, 
			 int *haspickup, char *errmsg) {
    int status,nhdu,hdutype,i,anynul,nmask,jrem,n,k,useflat,j,nfu,ind;
    long npts;
    dstrct *dd;
    char msg[BUFSIZ];
    float *data,*mdata,value,wmin,*work,*profile,*sprofile,*scores,statistic;
    float scalebest[DEFNEXTN],wmed,wmad,a,b,c,offset,scoremin,scale,*fdata;
    float skyref,modflat,gaincor;
    unsigned char *mask;

    /* Open the file and make sure we have all the extensions */

    status = 0;
    (void)fits_open_file(&iptr,infile,READWRITE,&status);
    if (status != 0) {
        (void)sprintf(errmsg,"VSTPICKUP: Can't open input file %s",infile);
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_num_hdus(iptr,&nhdu,&status);
    nhdu--;
    if (nhdu != DEFNEXTN) {
        (void)sprintf(errmsg,
		      "VSTPICKUP: File %s has %d extensions must have %d",
		      infile,nhdu,DEFNEXTN);
	tidy();
	return(CIR_FATAL);
    }

    /* Loop through each of the detectors and make sure the x,y dimensions
       are appropriate for the trimmed files */

    for (i = 1; i <= nhdu; i++) {
	(void)fits_movabs_hdu(iptr,i+1,&hdutype,&status);
	(void)fits_get_img_size(iptr,2,naxis,&status);
	if (naxis[0] != DEFNX_TRIMMED || naxis[1] != DEFNY_TRIMMED) {
	    (void)sprintf(errmsg,
			  "VSTPICKUP: Extension %s[%d] incorrect dimensions [%ld,%ld]",
			  infile,i,naxis[0],naxis[1]);
	    tidy();
	    return(CIR_FATAL);
	}
    }

    /* If there is a flat field file specified and we can open it, then
       do that now. Otherwise don't do flat field correction */

    status = 0;
    useflat = 1;
    (void)fits_open_file(&fptr,flat,READONLY,&status);
    if (status != 0) {
	useflat = 0;
	status = 0;
    } else {
	(void)fits_get_num_hdus(fptr,&nhdu,&status);
	nhdu--;
	if (nhdu != DEFNEXTN) {
	    (void)sprintf(errmsg,
			  "VSTPICKUP: Flat %s has %d extensions must have %d\n",
			  flat,nhdu,DEFNEXTN);
	    tidy();
	    return(CIR_FATAL);
	}
	for (i = 1; i <= nhdu; i++) {
	    (void)fits_movabs_hdu(fptr,i+1,&hdutype,&status);
	    (void)fits_get_img_size(fptr,2,naxis,&status);
	    if (naxis[0] != DEFNX_TRIMMED || naxis[1] != DEFNY_TRIMMED) {
		(void)sprintf(errmsg,
			      "VSTPICKUP: Flat extension %s[%d] incorrect dimensions [%ld,%ld]",
			      flat,i,naxis[0],naxis[1]);
		tidy();
		return(CIR_FATAL);
	    }
	}
    }

    /* Get some memory for the input file info */

    fileinfo = cir_malloc(DEFNEXTN*sizeof(*fileinfo));
    for (i = 0; i < DEFNEXTN; i++) 
	(fileinfo+i)->data = NULL;

    /* Now read the data arrays and work out the background */

    npts = naxis[0]*naxis[1];
    for (i = 1; i <= nhdu; i++) {
	dd = fileinfo + i - 1;
	(void)fits_movabs_hdu(iptr,i+1,&hdutype,&status);
	dd->data = cir_malloc(npts*sizeof(float));
	(void)fits_read_img(iptr,TFLOAT,1,npts,NULL,dd->data,&anynul,&status);
	skyest(dd->data,npts,3.0,&(dd->sky),&(dd->skynoise));
    }

    /* Create an output stacked image with the appropriate extensions flipped
       Get the sky background of the stack */

    medianinfo = cir_malloc(sizeof(*medianinfo));
    medianinfo->data = cir_malloc(npts*sizeof(float));
    average();
    skyest(medianinfo->data,npts,3.0,&(medianinfo->sky),
	   &(medianinfo->skynoise));

    /* Get some memory for the 1d profiles */

    work = cir_malloc(naxis[0]*sizeof(*work));
    profile = cir_malloc(naxis[1]*sizeof(*profile));
    sprofile = cir_malloc(naxis[1]*sizeof(*sprofile));

    /* Form the profile from the MAD of each row from the median background.
       Compare this to the average sky in the stack */

    for (i = 0; i < naxis[1]; i++) {
	data = medianinfo->data + i*naxis[0];
	for (j = 0; j < naxis[0]; j++)
	    work[j] = (float)fabs(data[j] - medianinfo->sky);
	(void)cir_med(work,NULL,(int)naxis[0],profile+i,msg);
	profile[i] *= (1.48/medianinfo->skynoise);
	profile[i] = max(0.001,log10(profile[i]));
	sprofile[i] = profile[i];
    }
    freespace(work);
    
    /* Filter the copy of the profile */

    cir_filt1d(sprofile,(int)naxis[1],75,25,0.0);

    /* Look for 3-sigma events and clear out some workspace */

    nfu = 0;
    for (i = 0; i < naxis[1]; i++) {
	statistic = (float)pow(10.0,(double)(profile[i] - sprofile[i])) + 1.0;
	if (statistic > 3.0)
	    nfu++;
    }
    freespace(profile);
    freespace(sprofile);
    
    /* If there is no sign of any problems or if we are only checking then
       get out of here now */

    *haspickup = nfu;
    if (nfu == 0 || checkonly) {
	tidy();
	return(CIR_OK);
    }

    /* Create a pixel mask where the stacked image has bright pixels */

    mdata = medianinfo->data;
    mask = cir_calloc(npts,sizeof(*mask));
    nmask = 0;
    for (i = 0; i < npts; i++) {
	value = (float)fabs(mdata[i] - medianinfo->sky);
	if (value > 3.0*medianinfo->skynoise) {
	    mask[i] = 1;
	    nmask++;
	}
    }

    /* Get some workspace for all the bright pixels */

    work = cir_malloc(nmask*sizeof(*work));
    scores = cir_malloc(NSCALE*sizeof(*scores));

    /* Loop for each of the input extensions */

    for (i = 1; i <= DEFNEXTN; i++) {
        dd = fileinfo + i - 1;
	data = dd->data;

	/* Search between 0.5 and 1.5 for the best scale factor */

	wmin = 1.0e6;
	jrem = -1;
	for (j = 0; j < NSCALE; j++) {
	    scale = SCMIN + SCGAP*(float)j;
	    n = 0;
	    for (k = 0; k < npts; k++) {
		if (! mask[k])
		    continue;
		ind = k;
		if ((i >= 9 && i <= 16) || (i >= 25 && i <= 32)) 
		    ind = npts - ind - 1;
		work[n] = data[ind] - scale*mdata[k];
		n++;
	    }
	    (void)cir_medmad(work,NULL,n,&wmed,&wmad,msg);
	    wmad *= 1.48;
	    if (wmad < wmin) {
		scalebest[i-1] = scale;
		wmin = wmad;
		jrem = j;
	    }
	    scores[j] = wmad;
	}

	/* Parabolic interpolation about the minimum */

	if (jrem >= 1 && jrem < NSCALE-1) {
	    a = scores[jrem];
	    b = 0.5*(scores[jrem+1] - scores[jrem-1]);
	    c = 0.5*(scores[jrem+1] + scores[jrem-1] - 2.0*a);
	    offset = max(-1.0,min(1.0,-0.5*b/c));
	    scoremin = a + b*offset + c*offset*offset;
	    if (scoremin > scores[jrem])
		offset = 0.0;
	} else {
	    offset = 0.0;
	    scoremin = scores[jrem];
	}
	scalebest[i-1] += offset*SCGAP;
    }
    freespace(work);
    freespace(mask);
    freespace(scores);

    /* Right, now correct each of the images. Use the flat if it's been
       requested */

    fdata = NULL;
    if (useflat) 
	fdata = cir_malloc(npts*sizeof(*fdata));
    for (i = 1; i <= DEFNEXTN; i++) {
	dd = fileinfo + i - 1;
	data = dd->data;
	mdata = medianinfo->data;
	skyref = medianinfo->sky;
	if (useflat) {
	    (void)fits_movabs_hdu(fptr,i+1,&hdutype,&status);
	    (void)fits_read_img(fptr,TFLOAT,1,npts,NULL,fdata,&anynul,&status);
	    (void)fits_read_key(fptr,TFLOAT,"GAINCOR",&gaincor,NULL,&status);
	    if (status != 0) {
		status = 0;
		gaincor = 1.0;
	    }
	    for (j = 0; j < npts; j++) {
		ind = j;
		if ((i >= 9 && i <= 16) || (i >= 25 && i <= 32)) 
		    ind = npts - ind - 1;
		modflat = max(0.5,min(2.0,fdata[ind]));
		data[ind] -= ((mdata[j] - skyref)*gaincor*scalebest[i-1]/modflat);
	    }
	} else {
	    for (j = 0; j < npts; j++) {
		ind = j;
		if ((i >= 9 && i <= 16) || (i >= 25 && i <= 32)) 
		    ind = npts - ind - 1;
		data[ind] -= ((mdata[j] - skyref)*scalebest[i-1]);
	    }
	}
    }
    freespace(fdata);

    /* Ok, now write the corrected data back to the input file */

    for (i = 1; i <= DEFNEXTN; i++) {
	data = (fileinfo + i - 1)->data;
	(void)fits_movabs_hdu(iptr,i+1,&hdutype,&status);
	(void)fits_write_img(iptr,TFLOAT,1,npts,data,&status);
	(void)fits_update_key(iptr,TFLOAT,"PICKCOR",scalebest+i-1,
			      "Scale factor for pickup correction",&status);
    }

    /* Right, get out of here */

    tidy();
    return(CIR_OK);
}
	

static void average(void) {
    int nf1,nf2,nrejmax,i,j,ind,nrej;
    float work[DEFNEXTN],work2[DEFNEXTN],avsky,avnoise,value,cliplev;
    char msg[BUFSIZ];
    dstrct *dd;
    long npts;

    /* Define where in the buffer the medians are calculated */

    nf1 = (DEFNEXTN/2) - 1;
    nf2 = nf1 + 1;
    nrejmax = DEFNEXTN/4;

    /* Work out median sky and sky noise from the input extensions */

    for (i = 0; i < DEFNEXTN; i++) {
	work[i] = (fileinfo+i)->sky;
	work2[i] = (fileinfo+i)->skynoise;
    }
    (void)cir_med(work,NULL,DEFNEXTN,&avsky,msg);
    (void)cir_med(work2,NULL,DEFNEXTN,&avnoise,msg);
    for (i = 0; i < DEFNEXTN; i++) 
	(fileinfo+i)->skynorm = avsky - (fileinfo+i)->sky;
    
    /* Do pixel by pixel median. Reverse the readout for half the
       array. Do a 3-sigma upper clip */

    npts = naxis[0]*naxis[1];
    for (i = 0; i < npts; i++) {
	for (j = 1; j <= DEFNEXTN; j++) {
	    dd = fileinfo + j - 1;
	    ind = i;
	    if ((j >= 9 && j <= 16) || (j >= 25 && j <= 32)) 
		ind = npts - ind - 1;
	    work[j-1] = dd->data[ind] + dd->skynorm;
    	}
	sort1(work,DEFNEXTN);
	value = 0.5*(work[nf1] + work[nf2]);
	nrej = 0;
	cliplev = value + 3.0*avnoise;
	while (nrej < nrejmax && work[DEFNEXTN-nrej-1] > cliplev) 
	    nrej++;
	if (nrej > 0) 
	    (void)cir_med(work,NULL,DEFNEXTN-nrej,&value,msg);
	medianinfo->data[i] = value;
    }
}

static void skyest(float *data, long npts, float thresh, float *skymed, 
		   float *skynoise) {
    unsigned char *bpm;
    char msg[BUFSIZ];

    /* Get a dummy bad pixel mask */

    bpm = cir_calloc(npts,sizeof(*bpm));

    /* Get the stats */

    (void)cir_qmedsig(data,bpm,npts,thresh,2,-1000.0,65535.0,skymed,skynoise,
		      msg);

    /* Clean up */

    freespace(bpm);
}   

static void sort1(float *a, int n) {
    int iii,ii,i,ifin,j;
    float b;

    iii = 4;
    while (iii < n)
        iii *= 2;
    iii = min(n,(3*iii)/4 - 1);

    while (iii > 1) {
        iii /= 2;
        ifin = n - iii;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + iii;
            if (a[i] > a[j]) {
                b = a[j];
                while (1) {
                    a[j] = a[i];
                    j = i;
                    i = i - iii;
                    if (i < 0 || a[i] <= b) 
                        break;
                }
                a[j] = b;
            }
        }
    }
}

static void tidy() {
    int i;

    closefits(fptr);
    closefits(iptr);
    if (medianinfo != NULL) {
	freespace(medianinfo->data);
	freespace(medianinfo);
    }
    if (fileinfo != NULL) {
	for (i = 0; i < DEFNEXTN; i++) 
	    freespace((fileinfo+i)->data);
	freespace(fileinfo);
    }
}

/*

$Log: cir_vstpickup.c,v $
Revision 1.1  2012/12/08 07:31:00  jim
New entry


*/
