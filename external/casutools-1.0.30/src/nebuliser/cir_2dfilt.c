/*

$Id: cir_2dfilt.c,v 1.1 2010/09/06 09:03:50 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <tools.h>
#include <nebuliser.h>
#include <misc.h>

static float *idata = NULL;
static unsigned char *bdata = NULL;
static short int *cdata = NULL;
static float *odata = NULL;
static fitsfile *iptr = NULL;
static fitsfile *cpmptr = NULL;
static fitsfile *optr = NULL;
static char cvsid[] = "$Id: cir_2dfilt.c,v 1.1 2010/09/06 09:03:50 jim Exp $";

static void cir_conf2bpm(short int *conf, unsigned char *bpm, int nelm);
static void tidy();

/*+
 *  Name:
 *      cir_2dfilt
 *
 *  Purpose:
 *      Filter the background from a 2d image.
 *
 *  Description:
 *      An iterated sliding median filter is done to an image to remove
 *      large scale variations in the backgroud.  
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infile = char * (Given)
 *          The name of the input FITS image. 
 *      confmap = char * (Given) If a confidence map is specified then it will
 *          also be used as a bad pixel mask. If you don't want to use a
 *          confidence map the the value 'noconf' should be given.
 *      nfiltmed = int (Given)
 *          The size of the smoothing box for the median filter window.
 *      nfiltlin = int (Given)
 *          The size of the smoothing box for the linear filter window.
 *      backmap = char * (Given)
 *          An optional FITS file name for an output background map
 *      niter = int (Given)
 *          The number of clipping iteration loops
 *      axis = int (Given)
 *          Can take a value of either 1 or 2, indicating which axis will
 *          be filtered first.
 *      twod = int (Given)
 *          If set, then a full 2d filter will be done
 *      takeout_sky = int (Given)
 *          If set, then the median sky is also removed
 *      inorm = int (Given)
 *          If unset, the output will be the original data minus the background
 *          filtered map. If set, then the output will be the original data
 *          divided by the background map
 *      signeg = float (Given)
 *          The lower clipping threshold in units of background sigma
 *      sigpos = float (Given)
 *          The upper clipping threshold in units of background sigma
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      None
 *              
 *  Dependencies:
 *      cfitsio, cir_weight.c, cir_mikesubs.c
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2010-2013 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_2dfilt(char *infile, char *confmap, int nfiltmed, int nfiltlin,
		      char *backmap, int niter, int axis, int twod, 
		      int takeout_sky, int inorm, float signeg, float sigpos, 
		      char *errmsg) {
    int status,usebpm,anynul,nx,ny,wantback;
    long naxis[2],npts;
    char msg[BUFSIZ];

    /* Open the input file */

    status = 0;
    (void)fits_open_file(&iptr,infile,READWRITE,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"2DFILT: Can't open input file %s -- %s",infile,
		      msg);
        tidy();
        return(CIR_FATAL);
    }
    (void)fits_get_img_size(iptr,2,naxis,&status);
    nx = naxis[0];
    ny = naxis[1];
    npts = nx*ny;

    /* If there is a confidence map, then open that now...*/

    if (strcmp(confmap,"noconf") == 0) {
	usebpm = 0;
    } else {
	(void)fits_open_file(&cpmptr,confmap,READONLY,&status);
        if (status != 0) {
            fits_get_errstatus(status,msg);
            (void)sprintf(errmsg,"2DFILT: Can't open CPM file %s -- %s",
			  confmap,msg);
            tidy();
            return(CIR_FATAL);
        }
	usebpm = 1;
    }

    /* Read the input image into memory. */

    idata = cir_malloc(npts*sizeof(*idata));
    (void)fits_read_img(iptr,TFLOAT,1,npts,NULL,idata,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"2DFILT: Can't read image -- %s",msg);
        tidy();
        return(CIR_FATAL);
    }

    /* If a confidence map is available, then read it and convert it to a 
       bad pixel mask. If not, then just create one */

    if (usebpm) {
	cdata = cir_malloc(npts*sizeof(*cdata));
	bdata = cir_malloc(npts*sizeof(*bdata));
	(void)fits_read_img(cpmptr,TSHORT,1,npts,NULL,cdata,&anynul,&status);
        if (status != 0) {
            fits_get_errstatus(status,msg);
            (void)sprintf(errmsg,"2DFILT: Can't read conf map -- %s",msg);
            tidy();
            return(CIR_FATAL);
        }
        cir_conf2bpm(cdata,bdata,npts);
        freespace(cdata);
    } else {
	bdata = cir_calloc(npts,sizeof(*bdata));
    }

    /* Open output file for the background map if you want one */

    if (strlen(backmap)) {
	if (cir_open_output(backmap,infile,&optr,NEWBITPIX,FLOAT_IMG,0,NULL,
			    msg) != CIR_OK) {
	    (void)sprintf(errmsg,"2DFILT: Can't open output file -- %s",msg);
	    tidy();
	    return(CIR_FATAL);
	}
	wantback = 1;
    } else {
	optr = NULL;
	wantback = 0;
    }

    /* Do the filtering */

    twodfilt(idata,bdata,nx,ny,nfiltmed,nfiltlin,niter,axis,twod,takeout_sky,
	     inorm,wantback,signeg,sigpos,&odata);

    /* Write the data out */

    (void)fits_write_img(iptr,TFLOAT,1,npts,idata,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"2DFILT: Can't write file -- %s",msg);
        tidy();
        return(CIR_FATAL);
    }
    if (wantback) {
	(void)fits_write_img(optr,TFLOAT,1,npts,odata,&status);
	if (status != 0) {
	    fits_get_errstatus(status,msg);
	    (void)sprintf(errmsg,"2DFILT: Can't write file -- %s",msg);
	    tidy();
	    return(CIR_FATAL);
	}
    }

    /* Close up files and free workspace */

    tidy();

    /* Exit */

    return(CIR_OK);
}

static void cir_conf2bpm(short int *conf, unsigned char *bpm, int nelm) {
    int i;

    /* If the confidence is zero it's a bad pixel (bpm == 1) otherwise it's
       a good pixel (bpm == 0) */

    for (i = 0; i < nelm; i++)
        bpm[i] = ((conf[i] == 0) ? 1 : 0);
}
             
static void tidy() {
    int status;

    freespace(idata);
    freespace(bdata);
    freespace(cdata);
    freespace(odata);
    closefits(iptr);
    closefits(cpmptr);
    closefits(optr);
}

/*

$Log: cir_2dfilt.c,v $
Revision 1.1  2010/09/06 09:03:50  jim
New entry

Revision 1.1  2010/04/27 11:23:41  jim
New routine


*/
