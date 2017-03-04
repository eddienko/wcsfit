/*

$Id: cir_imstack_cat.c,v 1.11 2013/02/07 16:05:22 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>

#include "imstack.h"
#include <tools.h>

#define M_PI  3.14159265358979323846
#define TWOPI 2*M_PI

typedef struct {
    double     crval[2];
    double     crpix[2];
    double     cd[4];
    double     secd;
    double     tand;
    double     pv21;
    double     pv23;
    double     pv25;
    int        nx;
    int        ny;
} mywcs;

typedef struct {
    char           *infile;
    char           *conf;
    int            nx;
    int            ny;
    fitsfile       *fptr;
    fitsfile       *cptr;
    unsigned char  *bpm;
    mywcs          *vwcs;
    float          xoff;
    float          yoff;
    double         *trans;
    float          sky;
    float          skydiff;
    float          noise;
    float          expscale;
    float          weight;
    float          magzptscale;
} dstrct;

typedef struct {
    char   *fname;
    int    nrows;
    float  *x;
    float  *y;
} catstrct;

typedef struct {
    int    nrows;
    double *x1;
    double *y1;
    double *x2;
    double *y2;
} matchstrct;

typedef struct {
    char      *outname;
    char      *outcname;
    fitsfile  *optr;
    fitsfile  *ocptr;
    int       nxo;
    int       nyo;
    mywcs     *outwcs;
    float     *outdata;
    short int *outcdata;
} outstrct;

static float lsig;
static float hsig;
static float sumweight;

static void output_files(char *out, char *outc, dstrct *fileptrs,
			 int nimages, outstrct *outstr);
static void backgrnd_ov(dstrct *fileptrs, int nimages);
static void seeing_wt(dstrct *fileptrs, int nimages);
static void skyest(float *data, short int *cdata, long npts, float thresh, 
		   float *skymed, float *skynoise);
static int stack_nn(dstrct *fileptrs, int nim, outstrct *outstr);
static int stack_lin(dstrct *fileptrs, int nim, outstrct *outstr);
static void do_averages(dstrct *fileptrs, int ncontrib, float *data, 
			float *wconf, float *conf, unsigned char *id, 
			float lclip, float hclip, float extra, 
			short int *lowcl, short int *highcl, float *outval, 
			float *outvalc);
static void diffxy(mywcs *wref, mywcs *wprog, float *dx, float *dy);
static void outloc(mywcs *win, double xin, double yin, mywcs *wout, 
		   double *tdata, double *xout, double *yout);
static mywcs *getwcsinfo(fitsfile *fptr);
static double *transinit(void);
static void platexy(matchstrct *matchedxy, int nconst, double **cc, 
		    int *status);
static int plate6(double *xpos, double *ypos, double *xi, double *eta,
		  unsigned char *flag, int npts, double *a, double *b, 
		  double *c, double *d, double *e, double *f);
static int plate4(double *xpos, double *ypos, double *xi, double *eta,
		  unsigned char *flag, int npts, double *a, double *b,
		  double *c, double *d, double *e, double *f);
static void freematch (matchstrct *m);
static void freecat(catstrct *c);
static void shiftcat(mywcs *pwcs, mywcs *rwcs, catstrct *c);
static int matchcat(catstrct *ref, catstrct *prog, int nx, int ny, 
		    matchstrct *match);
static catstrct readcat(char *catfile);
static void sorty(float *x, float *y, int n);
static void cir_prov(int method, fitsfile *outstr, char *fname, char **infiles, 
		     int nimages);
static void writewcs(fitsfile *optr, mywcs *outwcs);

/*+
 *  Name:
 *      cir_imstack_cat
 *
 *  Purpose:
 *      Stack list of images with WCS information and optional catalogues
 *
 *  Description:
 *      A list of WCS corrected images are stacked into an output image 
 *      using the on-board WCS information encoded in the FITS header. If
 *      object catalogues are available these can be employed to refine
 *      the WCS solution. Input images are single extensions of FITS files.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infiles = char ** (Given)
 *          A list of input FITS extensions to be stacked. These must have
 *          a WCS defined in the header.
 *      confs = char ** (Given)
 *          A list of input FITS confidence maps. There must either be one
 *          for each image in the infiles list, or just a single confidence
 *          map which will be used by all of the input images.
 *      cats = char ** (Given)
 *          An optional list of input object catalogues in FITS tables. These
 *          can come from any source, but must have the columns "X_coordinate" 
 *          and "Y_coordinate". If these are used there must be one for each
 *          input image
 *      nimages = int (Given)
 *          The number of images in the infiles list
 *      nconfs = int (Given)
 *          The number of images in the confs list
 *      ncats = int (Given)
 *          The number of tables in the cats list
 *      lthr = float (Given)
 *          The low rejection threshold in units of background sigma
 *      hthr = float (Given)
 *          The high rejection threshold in units of background sigma
 *      method = int (Given)
 *          The interpolation method to be used. Currently this is either
 *          0: nearest neighbour or 1: bi-linear interpolation
 *      nplate = int (Given)
 *          The number of plate constants to be fit if using the catalogues
 *          to refine the WCS solution. This can either be 4 (two offsets,
 *          scale and rigid rotation angle) or 6 (two offsets, two scales and
 *          two non-rigid rotation angles). 
 *      expkey = char * (Given)
 *          The FITS keyword for the exposure time in the input header
 *      seeing = int (Given)
 *          If set then the images are weighted by their mean seeing. This
 *          is read from the keyword "SEEING" in the image headers
 *      magzptref = float (Given)
 *          If greater than zero, the reference magnitude zeropoint to
 *          which to scale the input images.
 *      out = char * (Given)
 *          The name of the FITS file for the output stack. 
 *      outc = char * (Given)
 *          The name of the FITS file for the output confidence map.
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      Only TAN and ZPN projections are currently supported 
 *
 *  Dependencies:
 *      cfitsio, wcslib, tools
 *
 *  Authors:
 *      Mike Irwin (CASU, IoA)
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2009-2011 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_imstack_cat(char **infiles, char **confs, char **cats,
			   int nimages, int nconfs, int ncats, float lthr, 
			   float hthr, int method, int nplate, char *expkey, 
			   int seeing, float magzptref, char *out, char *outc, 
			   char *errmsg) {
    int i,status,xcol,ycol,npts,domagscale;
    long naxis[2],naxisc[2];
    float exptime,expref,magzpt;
    dstrct *dd;
    char errs[BUFSIZ],errstr[BUFSIZ];
    outstrct outstr;
    dstrct *fileptrs = NULL;
    matchstrct matched;
    fitsfile *fptr;
    catstrct catref,catprog;

    /* Is there any point in being here? */

    if (nimages == 0) {
	(void)sprintf(errmsg,"IMSTACK_CAT: No images to combine\n");
	return(CIR_FATAL);
    }

    /* Check value of nconfs. It has to be either one or the same as nimages */
    
    if (nconfs != 1 && nconfs != nimages) {
	sprintf(errmsg,
		"IMSTACK_CAT: Number of conf maps must be 1 or nimages\n");
	return(CIR_FATAL);
    }

    /* Check value of ncats. It has to be either zero or the same as nimages */
    
    if (ncats != 0 && ncats != nimages) {
	sprintf(errmsg,
		"IMSTACK_CAT: Number of catalogues must be 0 or nimages\n");
	return(CIR_FATAL);
    }

    /* Test to see whether we can even open all of these files. If we can't 
       then get out of here right now */

    domagscale = (magzptref > 0.0);
    errs[0] = '\0';
    for (i = 0; i < nimages; i++) {
	status = 0;

	/* Check the input image */

	(void)fits_open_file(&fptr,infiles[i],READONLY,&status);
	if (status != 0) {
	    (void)sprintf(errstr,"%s can't be opened\n",infiles[i]);
	    strcat(errs,errstr);
	    continue;
        }
	(void)fits_get_img_size(fptr,2,naxis,&status);
	if (status != 0) {
	    (void)sprintf(errstr,"%s isn't an image\n",infiles[i]);
	    strcat(errs,errstr);
	    status = 0;
  	    (void)fits_close_file(fptr,&status);
	    continue;
        }
	if (domagscale) {
	    (void)fits_read_key(fptr,TFLOAT,"MAGZPT",&magzpt,NULL,&status);
	    if (status != 0) {
		(void)sprintf(errstr,"%s has no MAGZPT keyword\n",infiles[i]);
		strcat(errs,errstr);
		status = 0;
		(void)fits_close_file(fptr,&status);
		continue;
	    }
	}
	(void)fits_close_file(fptr,&status);

	/* Check confidence map image sizes match */

        if (i == 0 || nconfs == nimages) {
	    (void)fits_open_file(&fptr,confs[i],READONLY,&status);
	    if (status != 0) {
		(void)sprintf(errstr,"%s can't be opened\n",confs[i]);
		strcat(errs,errstr);
		continue;
	    }
   	    (void)fits_get_img_size(fptr,2,naxisc,&status);
  	    if (status != 0) {
	        (void)sprintf(errstr,"%s isn't an image\n",infiles[i]);
	        strcat(errs,errstr);
		naxisc[0] = 0;
		naxisc[1] = 0;
		status = 0;
	        (void)fits_close_file(fptr,&status);
	        continue;
            }
	    (void)fits_close_file(fptr,&status);
	}
	if (naxis[0] != naxisc[0] && naxis[1] != naxisc[1]) {
	    (void)sprintf(errstr,"%s dims don't match with conf map",
			  infiles[i]);
	    continue;
	}
    }

    /* Look at catalogues. Make sure they have the right columns */

    if (ncats != 0) {
	for (i = 0; i < ncats; i++) {
	    status = 0;
	    (void)fits_open_file(&fptr,cats[i],READONLY,&status);
	    if (status != 0) {
		(void)sprintf(errstr,"%s can't be opened\n",cats[i]);
		strcat(errs,errstr);
		continue;
	    }
	    (void)fits_get_colnum(fptr,CASEINSEN,"X_coordinate",&xcol,&status);
	    (void)fits_get_colnum(fptr,CASEINSEN,"Y_coordinate",&ycol,&status);
	    if (status != 0) {
	        (void)sprintf(errstr,"%s isn't a catalogue\n",cats[i]);
	        strcat(errs,errstr);
		status = 0;
	        (void)fits_close_file(fptr,&status);
		continue;
	    }
	    (void)fits_close_file(fptr,&status);
	}
    }
    if (strlen(errs)) {
	(void)sprintf(errmsg,"IMSTACK_CAT: Errors in input files:%s\n",errs);
	return(CIR_FATAL);
    }

    /* Initialise somc global variables */

    lsig = lthr;
    hsig = hthr;

    /* Allocate file struct array and fill it in. Start by opening the
       the input FITS images. We've already test that they exist and will
       open, so there's no need to do that again */

    fileptrs = cir_malloc(nimages*sizeof(*fileptrs));
    expref = 1.0;
    for (i = 0; i < nimages; i++) {
	dd = fileptrs + i;
	dd->infile = infiles[i];
	status = 0;
	(void)fits_open_file(&fptr,infiles[i],READONLY,&status);
	dd->fptr = fptr;

	/* Get the data array */

	(void)fits_get_img_size(fptr,2,naxis,&status);
	dd->nx = (int)naxis[0];
	dd->ny = (int)naxis[1];
	npts = dd->nx*dd->ny;
	dd->bpm = cir_calloc(npts,sizeof(unsigned char));

	/* Get the exposure time */

	(void)fits_read_key(fptr,TFLOAT,expkey,&exptime,NULL,&status);
	if (status != 0) {
	    status = 0;
	    exptime = 1.0;
	}
	if (i == 0)
	    expref = exptime;
	dd->expscale = exptime/expref;

	/* Get the magnitude zeropoint and work out the scale factor */

	dd->magzptscale = 1.0;
	if (domagscale) {
	    (void)fits_read_key(fptr,TFLOAT,"MAGZPT",&magzpt,NULL,&status);
	    dd->magzptscale = (float)pow(10.0,0.4*(double)(magzptref - magzpt));
	}

	/* Read the WCS from the header */

	dd->vwcs = getwcsinfo(fptr);
    
	/* Get a rough shift between the current frame and the reference */
	
	if (i != 0) {
   	    diffxy(fileptrs->vwcs,dd->vwcs,&(dd->xoff),&(dd->yoff));
	} else {
	    dd->xoff = 0.0;
	    dd->yoff = 0.0;
	}

	/* Get the confidence map. If there is only one, then make a reference
	   to the information in the reference image structure */

	if ((i == 0 && nconfs == 1) || (nconfs > 1)) {
  	    (void)fits_open_file(&fptr,confs[i],READONLY,&status);
	    dd->cptr = fptr;
	    dd->conf = confs[i];
	} else {
	    dd->cptr = fileptrs->cptr;
	    dd->conf = fileptrs->conf;
	}

	/* If catalogues exist, then work out a transformation correction
	   based on matched object positions */
	
	if (ncats == 0) {
	    dd->trans = NULL;
        } else {
            if (i == 0) {
		dd->trans = transinit();
		catref = readcat(cats[i]);
	    } else {
		catprog = readcat(cats[i]);
		shiftcat(dd->vwcs,fileptrs->vwcs,&catprog);
		if (matchcat(&catref,&catprog,fileptrs->nx,fileptrs->ny,
			     &matched) != CIR_OK) {
		    dd->trans = transinit();
		} else {
		    status = CIR_OK;
  		    platexy(&matched,nplate,&(dd->trans),&status);
		    freematch(&matched);
		}
		freecat(&catprog);
		if (i == nimages-1)
		    freecat(&catref);
            }
	}
    }

    /* Get the background levels in the overlap regions */

    backgrnd_ov(fileptrs,nimages);

    /* Do seeing weighting if you want to */

    if (seeing) 
	seeing_wt(fileptrs,nimages);

    /* Create the output file */

    output_files(out,outc,fileptrs,nimages,&outstr);

    /* Ok time to do the stacking. Our first job is to do a nearest neighbour
       stack in order to find the objects that will be rejected. If we 
       only want the NN algorithm, then we're done. Otherwise we can use the
       rejection information to restack using any other algorithm */

    stack_nn(fileptrs,nimages,&outstr);

    /* Switch for final stacking method */

    switch (method) {
    case 0:
	break;
    case 1:
	stack_lin(fileptrs,nimages,&outstr);
	status = 0;
	(void)fits_update_key(outstr.optr,TSTRING,"DRIBBLE","bilinear",
			      "Interpolation method",&status);
	break;
    default:
	break;
    }

    /* Add provenance keywords */

    cir_prov(method,outstr.optr,outstr.outname,infiles,nimages);
    cir_prov(method,outstr.ocptr,outstr.outcname,confs,nconfs);
    
    /* Write the output file and tidy the output structure away */

    status = 0;
    (void)fits_write_img(outstr.optr,TFLOAT,1,outstr.nxo*outstr.nyo,
			 outstr.outdata,&status);
    (void)fits_update_key(outstr.optr,TSTRING,"CIR_CPM",outstr.outcname,NULL,
			  &status);
    (void)fits_update_key(outstr.optr,TINT,"NICOMB",&nimages,
			  "Number of images in this stack",&status);
    if (domagscale) 
	(void)fits_update_key(outstr.optr,TFLOAT,"MAGZPT",&magzptref,
			      NULL,&status);
    writewcs(outstr.optr,outstr.outwcs);
    (void)fits_close_file(outstr.optr,&status);
    (void)fits_write_img(outstr.ocptr,TSHORT,1,outstr.nxo*outstr.nyo,
			 outstr.outcdata,&status);
    (void)fits_update_key(outstr.ocptr,TINT,"NICOMB",&nimages,
			  "Number of images in this stack",&status);
    writewcs(outstr.ocptr,outstr.outwcs);
    (void)fits_close_file(outstr.ocptr,&status);
    freespace(outstr.outdata);
    freespace(outstr.outcdata);
    freespace(outstr.outwcs);

    /* Tidy up the rest */

    for (i = 0; i < nimages; i++) {
	status = 0;
	dd = fileptrs + i;
	closefits(dd->fptr);
	if ((i == 0) || (nconfs > 1)) 
	    closefits(dd->cptr);
	freespace(dd->bpm);
	freespace(dd->vwcs);
	freespace(dd->trans);
    }
    freespace(fileptrs);
    
    /* Get out of here */

    return(CIR_OK);
}

static void output_files(char *out, char *outc, dstrct *fileptrs,
			 int nimages, outstrct *outstr) {
    double lowerleft[2],upperleft[2],lowerright[2],upperright[2],xout,yout;
    double xmin,ymin;
    int i,ixo,iyo;
    long naxis[2];
    dstrct *dd;
    mywcs *outwcs;
    char errmsg[BUFSIZ];

    /* Set some preliminary info */

    outstr->outname = out;
    outstr->outcname = outc;

    /* Work out the output image size. Initialise the first element of 
       the results arrays since we're doing all this relative to the first 
       image */

    lowerleft[0] = 1.0;
    lowerleft[1] = 1.0;
    upperleft[0] = 1.0;
    upperleft[1] = (double)fileptrs->ny;
    lowerright[0] = (double)fileptrs->nx;
    lowerright[1] = 1.0;
    upperright[0] = (double)fileptrs->nx;
    upperright[1] = (double)fileptrs->ny;

    /* Loop for all images relative to the first one */

    for (i = 1; i < nimages; i++) {
	dd = fileptrs + i;
	outloc(dd->vwcs,(double)1.0,(double)1.0,fileptrs->vwcs,dd->trans,
	       &xout,&yout);
	lowerleft[0] = min(lowerleft[0],xout);
	lowerleft[1] = min(lowerleft[1],yout);
	outloc(dd->vwcs,(double)1.0,(double)dd->ny,fileptrs->vwcs,dd->trans,
	       &xout,&yout);
	upperleft[0] = min(upperleft[0],xout);
	upperleft[1] = max(upperleft[1],yout);
	outloc(dd->vwcs,(double)dd->nx,(double)1.0,fileptrs->vwcs,dd->trans,
	       &xout,&yout);
	lowerright[0] = max(lowerright[0],xout);
	lowerright[1] = min(lowerright[1],yout);
	outloc(dd->vwcs,(double)dd->nx,(double)dd->ny,fileptrs->vwcs,dd->trans,
	       &xout,&yout);
	upperright[0] = max(upperright[0],xout);
	upperright[1] = max(upperright[1],yout);
    }

    /* Ok, what are the limits? */

    ixo = cir_nint(max(lowerright[0]-lowerleft[0],upperright[0]-upperleft[0])) + 1;
    iyo = cir_nint(max(upperright[1]-lowerright[1],upperleft[1]-lowerleft[1])) + 1;
    xmin = min(lowerleft[0],upperleft[0]);
    ymin = min(lowerleft[1],lowerright[1]);
    outstr->nxo = ixo;
    outstr->nyo = iyo;
    naxis[0] = (long)ixo;
    naxis[1] = (long)iyo;
    
    /* Create the images */

    cir_open_output(out,fileptrs->infile,&(outstr->optr),NEWDATASZ,FLOAT_IMG,2,
		    naxis,errmsg);
    cir_open_output(outc,fileptrs->infile,&(outstr->ocptr),NEWDATASZ,SHORT_IMG,
		    2,naxis,errmsg);
    outstr->outdata = cir_malloc(ixo*iyo*sizeof(float));
    outstr->outcdata = cir_malloc(ixo*iyo*sizeof(short int));

    /* Update the reference point for the WCS */

    outwcs = getwcsinfo(fileptrs->fptr);
    outwcs->crpix[0] -= (xmin - 1.0);
    outwcs->crpix[1] -= (ymin - 1.0);
    outstr->outwcs = outwcs;
}

static void seeing_wt(dstrct *fileptrs, int nimages) {
    int i,ierr,status;
    float *seeing,sval;
    dstrct *dd;

    /* Small array for seeing estimates */

    seeing = cir_malloc(nimages*sizeof(*seeing));

    /* Read the seeing keyword. Get out of here if one is missing and just
       use the background weights that we already have */

    ierr = 0;
    for (i = 0; i < nimages; i++) {
	dd = fileptrs + i;
	status = 0;
	(void)fits_read_key(dd->fptr,TFLOAT,"SEEING",&sval,NULL,&status);
	if (status != 0) {
	    fprintf(stderr,"IMSTACK: File %s missing SEEING keyword\n",
		    dd->infile);
	    ierr = 1;
	    break;
	} else if (sval == 0.0) {
	    fprintf(stderr,"IMSTACK: File %s has SEEING = %g\n",dd->infile,
		    sval);
	    ierr = 1;
	    break;
	} else {
	    seeing[i] = sval;
	}
    }

    /* If there weren't errors, then modify the image weights */

    if (! ierr) {
	sumweight = fileptrs->weight;
	for (i = 1; i < nimages; i++) {
	    dd = fileptrs + i;
	    dd->weight *= seeing[0]/seeing[i];
	    sumweight += dd->weight;
	}
    } else {
	fprintf(stderr,"         No seeing weighting will be done\n");
    }
    freespace(seeing);
}

static void backgrnd_ov(dstrct *fileptrs, int nimages) {
    int i,dx,dy,n,jj,ind1,nx,ii,ind,status,anynul;
    short int *cdata1,*cdata2,*cdata;
    long npts1,npts2,npts,start[2],end[2];
    dstrct *dd1,*dd2;
    float *data1,*data2,*data,sky1,skynoise1,sky2,skynoise2,frac,scfac1,scfac2;

    /* Initialise a few things */

    dd1 = fileptrs;
    npts1 = dd1->nx*dd1->ny;
    data1 = cir_malloc(npts1*sizeof(*data1));
    cdata1 = cir_malloc(npts1*sizeof(*cdata1));
    status = 0;
    (void)fits_read_img(dd1->fptr,TFLOAT,1,npts1,NULL,data1,&anynul,&status);
    (void)fits_read_img(dd1->cptr,TSHORT,1,npts1,NULL,cdata1,&anynul,&status);

    /* Start by working out the full background of the first frame */

    skyest(data1,cdata1,npts1,3.0,&sky1,&skynoise1);
    scfac1 = (dd1->magzptscale)/(dd1->expscale);
    dd1->sky = sky1*scfac1;
    dd1->noise = skynoise1*scfac1;
    dd1->skydiff = 0.0;
    dd1->weight = 1.0;
    sumweight = 1.0;

    /* Now loop for all the others */    

    for (i = 1; i < nimages; i++) {
	dd2 = fileptrs + i;
	scfac1 = dd1->magzptscale/dd1->expscale;
	scfac2 = dd2->magzptscale/dd2->expscale;
	npts2 = dd2->nx*dd2->ny;
	data2 = cir_malloc(npts2*sizeof(*data1));
	cdata2 = cir_malloc(npts2*sizeof(*cdata1));
	status = 0;
	(void)fits_read_img(dd2->fptr,TFLOAT,1,npts2,NULL,data2,&anynul,
			    &status);
	(void)fits_read_img(dd2->cptr,TSHORT,1,npts2,NULL,cdata2,&anynul,
			    &status);

	/* Offset differences */

	dx = cir_nint(dd1->xoff - dd2->xoff);
	dy = cir_nint(dd1->yoff - dd2->yoff);

	/* Ok, here's the tricky bit. Work out the range of x and y for
	   each of the images. The shift moves an image to the left and
	   and down. Hence positive values of dx or dy moves the second
	   frame further to the left of down compared to the first frame.
	   Here are the coordinate ranges for the first image */

	start[0] = max(0,-dx);
	start[1] = max(0,-dy);
	end[0] = min(dd1->nx-1,(dd1->nx)-dx-1);
	end[1] = min(dd1->ny-1,(dd1->ny)-dy-1);
	nx = dd1->nx;
	
	/* How much does this cover the first image? */

	npts = (end[0]-start[0]+1)*(end[1]-start[1]+1);
	frac = (float)npts/(float)npts1;

	/* If the coverage is more than 50% then calculate the background
	   in the coverage region for the first image */

	if (frac >= 0.5) {
 	    data = cir_malloc(npts*sizeof(*data));
 	    cdata = cir_malloc(npts*sizeof(*cdata));
	    n = 0;
	    for (jj = start[1]; jj <= end[1]; jj++) {
		ind1 = jj*nx;
		for (ii = start[0]; ii <= end[0]; ii++) {
		    ind = ind1 + ii;
		    data[n] = data1[ind];
		    cdata[n++] = cdata1[ind];
		}
	    }
	    skyest(data,cdata,npts,3.0,&sky1,&skynoise1);
	    sky1 *= scfac1;
            freespace(data);
	    freespace(cdata);
	} else {
	    sky1 = dd1->sky;
	}

	/* And here are the coordinate ranges for the second image if 
	   the coverage is more than 50%. Get the subset you need. */

	if (frac > 0.5) {
  	    start[0] = max(1,dx);
 	    start[1] = max(1,dy);
	    end[0] = min(dd2->nx-1,(dd2->nx)+dx-1);
	    end[1] = min(dd2->ny-1,(dd2->ny)+dy-1);
	    nx = dd2->nx;
 	    npts = (end[0]-start[0]+1)*(end[1]-start[1]+1);
 	    data = cir_malloc(npts*sizeof(*data));
	    cdata = cir_malloc(npts*sizeof(*cdata));
	    n = 0;
	    for (jj = start[1]; jj <= end[1]; jj++) {
		ind1 = jj*nx;
		for (ii = start[0]; ii <= end[0]; ii++) {
		    ind = ind1 + ii;
		    data[n] = data2[ind];
		    cdata[n++] = cdata2[ind];
		}
	    }
	    skyest(data,cdata,npts,3.0,&sky2,&skynoise2);
	    freespace(data);
	    freespace(cdata);
	} else {		
	    npts = (dd2->nx)*(dd2->ny);
	    skyest(data2,cdata2,npts,3.0,&sky2,&skynoise2);
	}

	/* Store info away */

	dd2->sky = sky2*scfac2;
	dd2->skydiff = sky1 - dd2->sky;
	dd2->noise = skynoise2*scfac2;
	if (dd2->noise != 0.0 && dd2->sky != 0.0) 
  	    dd2->weight = (float)(pow(dd1->noise,2.0)/pow(dd2->noise,2.0)); 
	else 
	    dd2->weight = 0.0;
	sumweight += dd2->weight;
	freespace(data2);
	freespace(cdata2);
    }
    freespace(data1);
    freespace(cdata1);
}

static void skyest(float *data, short int *cdata, long npts, float thresh, 
		   float *skymed, float *skynoise) {
    unsigned char *bpm;
    int i;
    char msg[BUFSIZ];

    /* Get a dummy bad pixel mask */

    bpm = cir_calloc(npts,sizeof(*bpm));
    for (i = 0; i < npts; i++)
	bpm[i] = (cdata[i] == 0);

    /* Get the stats */

    (void)cir_qmedsig(data,bpm,npts,thresh,2,-1000.0,65535.0,skymed,skynoise,
		      msg);

    /* Clean up */

    freespace(bpm);
}

static int stack_nn(dstrct *fileptrs, int nim, outstrct *outstr) {
    float *workbuf,*workbufc,*workbufw,lclip,hclip,*data,wt,med,noise,outvalc;
    float value,*dbuf,*wbuf,*cbuf,outval,avlev,avvar,*confwork,scfac;
    float *outdata,scale_avvar;
    double oxf,oyf;
    int *iloc,i,npts,ix,iy,ind1,ind,ox,oy,nn,bufloc,ncontrib,*lbuf;
    int j,nx,ny,n,nz,nxo,nyo,status,anynul;
    short int clipped_low,clipped_high,*cdata,*nbuf,*outdatac;
    unsigned char *id,*clipmap,*cc,*cclast,*ibuf;
    dstrct *dd;
    char msg[BUFSIZ];
    mywcs *outwcs;

    /* Useful constant */

    scale_avvar = sqrt(M_PI/2.0)/9.0;

    /* Output data info */

    outdata = outstr->outdata;
    outdatac = outstr->outcdata;
    outwcs = outstr->outwcs;
    nxo = outstr->nxo;
    nyo = outstr->nyo;

    /* Get workspace for collecting info over the whole output map */

    nz = nim + 2;
    workbuf = cir_calloc(nz*nxo*nyo,sizeof(*workbuf));
    workbufc = cir_calloc(nz*nxo*nyo,sizeof(*workbufc));
    workbufw = cir_calloc(nz*nxo*nyo,sizeof(*workbufw));
    id = cir_calloc(nz*nxo*nyo,sizeof(*id));
    iloc = cir_calloc(nz*nxo*nyo,sizeof(*iloc));
    nbuf = cir_calloc(nxo*nyo,sizeof(*nbuf));
    confwork = cir_calloc(nxo*nyo,sizeof(*confwork));

    /* Workspace for marking output rows with clipped pixels */

    clipmap = cir_calloc(2*nxo,sizeof(*clipmap));

    /* Buffers for averaging an individual pixel */

    dbuf = cir_calloc(2*nz,sizeof(*dbuf));
    wbuf = cir_calloc(2*nz,sizeof(*wbuf));
    cbuf = cir_calloc(2*nz,sizeof(*cbuf));
    ibuf = cir_calloc(2*nz,sizeof(*ibuf));
    lbuf = cir_calloc(2*nz,sizeof(*lbuf));

    /* Rejection thresholds */

    lclip = fileptrs->sky - lsig*(fileptrs->noise);
    hclip = fileptrs->sky + hsig*(fileptrs->noise);

    /* Ok, loop for each of the input files and collect all the information
       for each pixel */

    for (i = 0; i < nim; i++) {
	dd = fileptrs + i;
	nx = dd->nx;
	ny = dd->ny;
	npts = nx*ny;
        data = cir_malloc(npts*sizeof(*data));
	cdata = cir_malloc(npts*sizeof(*cdata));
        status = 0;
        (void)fits_read_img(dd->fptr,TFLOAT,1,npts,NULL,data,&anynul,&status);
        (void)fits_read_img(dd->cptr,TSHORT,1,npts,NULL,cdata,&anynul,&status);
	wt = dd->weight;
	scfac = (dd->magzptscale)/(dd->expscale);

	/* Loop for each pixel in the data array and work out where it belongs
	   in the output array */

	n = 0;
	for (iy = 0; iy < ny; iy++) {
	    ind1 = iy*nx;
	    for (ix = 0; ix < nx; ix++) {
		ind = ind1 + ix;
		outloc(dd->vwcs,(double)ix,(double)iy,outwcs,dd->trans,
		       &oxf,&oyf);
		ox = cir_nint(oxf) - 1;
		oy = cir_nint(oyf) - 1;
		if (ox < 0 || ox >= nxo || oy < 0 || oy >= nyo || 
		    cdata[ind] == 0)
		    continue;
		nn = oy*nxo + ox;
		bufloc = oy*nxo*nz + ox*nz + nbuf[nn];
		value = data[ind]*scfac + dd->skydiff;
		workbuf[bufloc] = value;
		workbufw[bufloc] = wt*(float)cdata[ind];
		workbufc[bufloc] = (float)cdata[ind];
		id[bufloc] = i;
		nbuf[nn] += 1;
		iloc[bufloc] = ind;
	    }
	}
	freespace(data);
	freespace(cdata);
    }

    /* Now loop through the output arrays and form an initial estimate */

    for (iy = 0; iy < nyo; iy++) {
	ind1 = iy*nxo;
	cc = clipmap + (iy % 2)*nxo;
	cclast = clipmap + ((iy-1) % 2)*nxo;
	for (ix = 0; ix < nxo; ix++) {
	    ind = ind1 + ix;
	    bufloc = iy*nxo*nz + ix*nz;
	    ncontrib = nbuf[ind];
	    if (ncontrib == 0) {
		outdata[ind] = fileptrs->sky;
		outdatac[ind] = 0;
		continue;
	    } 
	    cc[ix] = 0;

	    /* Put some stuff in buffers for the averaging routine */
	    
	    for (i = 0; i < ncontrib; i++) {
		dbuf[i] = workbuf[bufloc+i];
		wbuf[i] = workbufw[bufloc+i];
		cbuf[i] = workbufc[bufloc+i];
		ibuf[i] = id[bufloc+i];
		lbuf[i] = iloc[bufloc+i];
	    }

	    /* Do the averages with the nominal clipping */

            do_averages(fileptrs,ncontrib,dbuf,wbuf,cbuf,ibuf,lclip,hclip,0.0,
			&clipped_low,&clipped_high,&outval,&outvalc);
	    
	    /* Store these away */

	    outdata[ind] = outval;
	    confwork[ind] = outvalc;
	    cc[ix] = (clipped_high >= 0);
	    if (clipped_low >= 0) 
		((fileptrs+ibuf[clipped_low])->bpm)[lbuf[clipped_low]] = 2;
	    if (clipped_high >= 0) 
		((fileptrs+ibuf[clipped_high])->bpm)[lbuf[clipped_high]] = 1;
	}

	/* Look at the non-edge pixels and see if the local variation can
	   show whether any clipping that has been done is justified. NB:
	   we are looking at the _previous_ row here */

	if (iy < 2)
	    continue;
	for (ix = 1; ix < nxo-1; ix++) {
	    if (! cclast[ix])
		continue;
	    ind = (iy-1)*nxo + ix;
	    bufloc = (iy-1)*nxo*nz + ix*nz;
	    ncontrib = nbuf[ind];
	    if (ncontrib == 0)
		continue;
	    avlev = 0.0;
	    for (i = -1; i <= 1; i++) 
		for (j = -1; j <= 1; j++)  
		    avlev += outdata[ind + i*nxo + j];
	    avlev /= 9.0;
	    avvar = 0.0;
	    for (i = -1; i <= 1; i++) 
		for (j = -1; j <= 1; j++) 
		    avvar += (float)fabs(avlev - outdata[ind + i*nxo + j]);
	    avvar *= scale_avvar;
	    if (avlev <= hclip || avvar <= (fileptrs->noise))
		continue;

	    /* Put some stuff in buffers for the averaging routine */
	    
	    for (i = 0; i < ncontrib; i++) {
		dbuf[i] = workbuf[bufloc+i];
		wbuf[i] = workbufw[bufloc+i];
		cbuf[i] = workbufc[bufloc+i];
		ibuf[i] = id[bufloc+i];
		lbuf[i] = iloc[bufloc+i];
		(fileptrs+ibuf[i])->bpm[lbuf[i]] = 0;
	    }

	    /* Redo the averages with the nominal clipping */

            do_averages(fileptrs,ncontrib,dbuf,wbuf,cbuf,ibuf,lclip,hclip,
			hsig*avvar,&clipped_low,&clipped_high,&outval,
			&outvalc);
	    
	    /* Store these away */

	    outdata[ind] = outval;
	    confwork[ind] = outvalc;
	    if (clipped_low >= 0) 
		((fileptrs+ibuf[clipped_low])->bpm)[lbuf[clipped_low]] = 2;
	    if (clipped_high >= 0) 
		((fileptrs+ibuf[clipped_high])->bpm)[lbuf[clipped_high]] = 1;
	}    
    }

    /* Now renormalise the confidence map */

    npts = nxo*nyo;
    cir_qmedsig(confwork,NULL,npts,3.0,2,25.0,65535.0,&med,&noise,msg);
    for (i = 0; i < npts; i++) {
	confwork[i] *= (100.0/med);
	outdatac[i] = max(0,min(1000,cir_nint(confwork[i])));
    }

    /* Free up some workspace */

    freespace(workbuf);
    freespace(workbufc);
    freespace(workbufw);
    freespace(confwork);
    freespace(id);
    freespace(iloc);
    freespace(nbuf);
    freespace(clipmap);
    freespace(dbuf);
    freespace(wbuf);
    freespace(cbuf);
    freespace(ibuf);
    freespace(lbuf);
    return(CIR_OK);
}

static int stack_lin(dstrct *fileptrs, int nim, outstrct *outstr) {
    double oxf,oyf;
    int npts,nptso,i,ix,iy,ind1,ind,ox,oy,m,jj,ii,nn,nx,ny,n,nxo,nyo;
    int status,anynul;
    short int *outdatac,*cdata;
    float *workbufw,*outdata,*data,wt,dx,dy,w[4],value,cval;
    float med,noise,scfac;
    dstrct *dd;
    char msg[BUFSIZ];
    mywcs *outwcs;

    /* Output file info */

    outdata = outstr->outdata;
    outdatac = outstr->outcdata;
    outwcs = outstr->outwcs;
    nxo = outstr->nxo;
    nyo = outstr->nyo;
    nptso = nxo*nyo;

    /* Get workspace for collecting info over the whole output map */

    workbufw = cir_calloc(nptso,sizeof(*workbufw));

    /* Initialise the output map. This is because the output data array 
       probably still has the results of the nearest neighbour stack */

    for (i = 0; i < nptso; i++) {
	outdata[i] = 0.0;
	outdatac[i] = 0;
    }

    /* Ok, loop for each of the input files and collect all the information
       for each pixel */

    for (i = 0; i < nim; i++) {
	dd = fileptrs + i;
	nx = dd->nx;
	ny = dd->ny;
	npts = nx*ny;
        data = cir_malloc(npts*sizeof(*data));
	cdata = cir_malloc(npts*sizeof(*cdata));
        status = 0;
        (void)fits_read_img(dd->fptr,TFLOAT,1,npts,NULL,data,&anynul,&status);
        (void)fits_read_img(dd->cptr,TSHORT,1,npts,NULL,cdata,&anynul,&status);
	wt = dd->weight;
	scfac = (dd->magzptscale)/(dd->expscale);

	/* Loop for each pixel in the data array and work out where it belongs
	   in the output array */

	n = 0;
	for (iy = 0; iy < ny; iy++) {
	    ind1 = iy*nx;
	    for (ix = 0; ix < nx; ix++) {
		ind = ind1 + ix;
		outloc(dd->vwcs,(double)ix,(double)iy,outwcs,dd->trans,
		       &oxf,&oyf);
 		if (oxf < 0.0 || oyf < 0.0)
 		    continue;
		ox = (int)oxf;
		oy = (int)oyf;
		if ((dd->bpm)[ind] != 0)
		    continue;
		if (ox < 0 || ox >= nxo || oy < 0 || oy >= nyo ||
		    cdata[ind] == 0) 
		    continue;
		dx = oxf - ox;
		dy = oyf - oy;
		w[0] = (1.0-dx)*(1.0-dy);
		w[1] = dx*(1.0-dy);
		w[2] = dy*(1.0-dx);
		w[3] = dx*dy;
		m = 0;
		value = data[ind]*scfac + dd->skydiff;
		cval = (float)cdata[ind];
		for (jj = 0; jj <= 1; jj++) {
		    for (ii = 0; ii <= 1; ii++) {
			nn = (oy+jj)*nxo + ox + ii;
			if (nn >= nptso || ox+ii >= nxo || oy+jj >= nyo)
			    continue;
			outdata[nn] += w[m]*value*wt*cval;
			workbufw[nn] += w[m]*wt*cval;
			m++;
		    }
		}
	    }
	}
	freespace(data);
	freespace(cdata);
    }

    /* Now do the final averaging */

    for (i = 0; i < nptso; i++) {
	if (workbufw[i] != 0.0)
   	    outdata[i] /= workbufw[i];
	else
	    outdata[i] = fileptrs->sky;
 	workbufw[i] /= sumweight;
   }
    cir_qmedsig(workbufw,NULL,nptso,3.0,2,0.0,65535.0,&med,&noise,msg);
    for (i = 0; i < nptso; i++) {
	workbufw[i] *= (100.0/med);
	outdatac[i] = max(0,min(1000,cir_nint(workbufw[i])));
    }
    freespace(workbufw);
    return(CIR_OK);
}

static void do_averages(dstrct *fileptrs, int ncontrib, float *data, 
			float *wconf, float *conf, unsigned char *id, 
			float lclip, float hclip, float extra, 
			short int *lowcl, short int *highcl, float *outval, 
			float *outvalc) {
    float minval,maxval,sum,wsum,csum,reflev,clipval,clipscale;
    int imax,imin,i;

    /* Do the summation and keep track of the minimum and
       maximum values */

    minval = 1.0e10;
    maxval = -1.0e10;
    sum = 0.0;
    wsum = 0.0;
    csum = 0.0;
    imin = -1;
    imax = -1;
    for (i = 0; i < ncontrib; i++) {
	sum += data[i]*wconf[i];
	wsum += wconf[i];
	csum += conf[i];
	if (data[i] > maxval) {
	    maxval = data[i];
	    imax = i;
	}
	if (data[i] < minval) {
	    minval = data[i];
	    imin = i;
	}
    }

    /* First crack at output value */

    if (wsum > 0.0) 
	*outval = sum/wsum;
    else 
	*outval = fileptrs->sky;

    /* Look at high clipping for of cosmic rays */

    *highcl = -1;
    if (maxval > hclip && wsum > 150.0 && csum > 150.0) {
	reflev = (sum - data[imax]*wconf[imax])/(wsum - wconf[imax]);
	clipscale = max(1.0,reflev)/max(1.0,(fileptrs+id[imax])->sky);
	clipval = reflev + clipscale*hsig*(fileptrs+id[imax])->noise;
	clipval += extra;

	/* Clip the maximum point */

	if (maxval > clipval) {
	    sum -= data[imax]*wconf[imax];
	    wsum -= wconf[imax];
	    *outval = reflev;
	    *highcl = imax;
	}
    }

    /* Clip low points (presumably serious cosmetic issues on detectors) */

    *lowcl = -1;
    if (minval < lclip  && wsum > 150.0 && csum > 150.0) {
	reflev = (sum - data[imin]*wconf[imin])/(wsum - wconf[imin]);
	clipscale = max(1.0,reflev)/max(1.0,(fileptrs+id[imin])->sky);
	clipval = reflev + clipscale*hsig*(fileptrs+id[imin])->noise;
	if (minval < clipval) {
	    sum -= data[imin]*wconf[imin];
	    wsum -= wconf[imin];
	    *outval = reflev;
	    *lowcl = imin;
	}
    }

    /* Output confidence (not normalised) */

    *outvalc = wsum/sumweight;
}

static void diffxy(mywcs *wref, mywcs *wprog, float *dx, float *dy) {
    double xref,yref,xprog,yprog;

    /* Find the middle of the reference frame */

    xref = 0.5*(double)wref->nx;
    yref = 0.5*(double)wref->ny;

    /* Work out the output position */

    outloc(wref,xref,yref,wprog,NULL,&xprog,&yprog);

    /* Now the offsets */

    *dx = (float)(xref - xprog);
    *dy = (float)(yref - yprog);
}

static void outloc(mywcs *win, double xin, double yin, mywcs *wout, 
		   double *tdata, double *xout, double *yout) {
    double xt,yt,xi,eta,r,rfac,aa,tandec,denom,rp;
    int i;

    /* Do the conversion. First to standard coordinates in the frame of
       the input image */

    xt = xin - win->crpix[0];
    yt = yin - win->crpix[1];
    xi = win->cd[0]*xt + win->cd[1]*yt;
    eta = win->cd[2]*xt + win->cd[3]*yt;
    if (fabs(xt) < 1.0e-6 && fabs(yt) < 1.0e-6) {
	rfac = 1.0;
    } else {
	rp = sqrt(xi*xi + eta*eta);
	r = rp;
	for (i = 0; i < 3; i++) {
	    rfac = win->pv21 + win->pv23*pow(r,2.0) + win->pv25*pow(r,4.0);
	    r = rp/rfac;
	}
	rfac = tan(r)/rp;
    }
    xi *= rfac;
    eta *= rfac;
    aa = atan(xi*win->secd/(1.0-eta*win->tand));
    if (xi != 0.0) 
	tandec = (eta+win->tand)*sin(aa)/(xi*win->secd);
    else
	tandec = (eta+win->tand)/(1.0 - eta*win->tand);

    /* Now from standard coordinates in the frame of the output image */

    aa += (win->crval[0] - wout->crval[0]);
    denom = wout->tand*tandec + cos(aa);
    xi = wout->secd*sin(aa)/denom;
    eta = (tandec - wout->tand*cos(aa))/denom;
    rp = sqrt(xi*xi + eta*eta);
    r = atan(rp);
    rfac = wout->pv21 + wout->pv23*pow(rp,2.0) + wout->pv25*pow(rp,4.0);
    rfac *= r/rp;
    xi *= rfac;
    eta *= rfac; 
    denom = wout->cd[0]*wout->cd[3] - wout->cd[1]*wout->cd[2];
    xt = (xi*wout->cd[3] - eta*wout->cd[1])/denom + wout->crpix[0];
    yt = (eta*wout->cd[0] - xi*wout->cd[2])/denom + wout->crpix[1];
    
    /* If a fine tune transformation is defined, then do that now */

    if (tdata != NULL) {
	*xout = xt*tdata[0] + yt*tdata[1] + tdata[2];
	*yout = xt*tdata[3] + yt*tdata[4] + tdata[5];
    } else {
	*xout = xt;
	*yout = yt;
    }
}

static mywcs *getwcsinfo(fitsfile *fptr) {
    mywcs *w;
    int i,status;

    /* Copy over the relevant info */

    w = cir_malloc(sizeof(mywcs));
    status = 0;
    (void)fits_read_key(fptr,TDOUBLE,"CRVAL1",&(w->crval[0]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CRVAL2",&(w->crval[1]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CRPIX1",&(w->crpix[0]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CRPIX2",&(w->crpix[1]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD1_1",&(w->cd[0]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD1_2",&(w->cd[1]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD2_1",&(w->cd[2]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD2_2",&(w->cd[3]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"PV2_1",&(w->pv21),NULL,&status);
    if (status != 0) {
	w->pv21 = 1.0;
	status = 0;
    }
    (void)fits_read_key(fptr,TDOUBLE,"PV2_3",&(w->pv23),NULL,&status);
    if (status != 0) {
	w->pv23 = -0.3333333;
	status = 0;
    }
    (void)fits_read_key(fptr,TDOUBLE,"PV2_5",&(w->pv25),NULL,&status);
    if (status != 0) {
	w->pv25 = 0.0;
	status = 0;
    }
    (void)fits_read_key(fptr,TINT,"NAXIS1",&(w->nx),NULL,&status);
    (void)fits_read_key(fptr,TINT,"NAXIS2",&(w->ny),NULL,&status);

    /* Change angles to radians */

    for (i = 0; i < 2; i++) 
	w->crval[i] *= DEGRAD;
    for (i = 0; i < 4; i++)
	w->cd[i] *= DEGRAD;

    /* Add a couple of convenience values */

    w->tand = tan(w->crval[1]);
    w->secd = 1.0/cos(w->crval[1]);

    /* Get out of here */

    return(w);
}

static double *transinit(void) {
    double *td;

    td = cir_malloc(6*sizeof(double));
    td[0] = 1.0;
    td[1] = 0.0;
    td[2] = 0.0;
    td[3] = 0.0;
    td[4] = 1.0;
    td[5] = 0.0;
    return(td);
}

static void platexy(matchstrct *matchedxy, int nconst, double **cc, 
		    int *status) {
    int nxy,nc2,i,niter,nrej,ngood;
    unsigned char *isbad,*wptr2,*iwork;
    double *xptr1,*xptr2,*yptr1,*yptr2,*wptr,a,b,c,d,e,f,xfit,yfit;
    double dx,dy,averr;
    char msg[BUFSIZ];

    /* How many objects are in the input table? */

    *cc = NULL;
    nxy = matchedxy->nrows;
    nc2 = nconst/2;
    if (nxy < nc2) {
	*cc = transinit();
	return;
    }

    /* Get some workspace now */

    wptr = cir_malloc(2*nxy*sizeof(*wptr));
    iwork = cir_calloc(3*nxy,sizeof(*iwork));
    isbad = iwork;
    wptr2 = iwork + nxy;
    
    /* Get the data from the table and put it all into double precision
       arrays */

    xptr1 = matchedxy->x1;
    yptr1 = matchedxy->y1;
    xptr2 = matchedxy->x2;
    yptr2 = matchedxy->y2;

    /* Right, now loop for maximum number of iterations or until
       convergence */

    niter = 0;
    while (niter >= 0) {

        /* Do a plate solution */

        switch (nconst) {
        case 6:
            *status = plate6(xptr2,yptr2,xptr1,yptr1,isbad,nxy,&a,&b,
			     &c,&e,&d,&f);
            break;
        case 4:
            *status = plate4(xptr2,yptr2,xptr1,yptr1,isbad,nxy,&a,&b,
			     &c,&e,&d,&f);
            break;
        default:
            *status = plate6(xptr2,yptr2,xptr1,yptr1,isbad,nxy,&a,&b,
			     &c,&e,&d,&f);
            break;
        }
        if (*status != CIR_OK) {
            freespace(wptr);
	    freespace(iwork);
	    return;
        }

        /* Now look at the residuals and see if any should be rejected */

        for (i = 0; i < nxy; i++) {
            xfit = xptr2[i]*a + yptr2[i]*b + c;
            yfit = xptr2[i]*d + yptr2[i]*e + f;
            dx = fabs(xfit - xptr1[i]);
            dy = fabs(yfit - yptr1[i]);
            wptr[i*2] = dx;
            wptr[i*2+1] = dy;
            wptr2[i*2] = isbad[i];
            wptr2[i*2+1] = isbad[i];
        }
        (void)cir_dmed(wptr,wptr2,2*nxy,&averr,msg);
        averr *= 1.48;
        if (niter == 3)
            break;
        nrej = 0;
        ngood = 0;
        for (i = 0; i < nxy; i++) {
            if (!isbad[i] && (wptr[i*2] > 3.0*averr || wptr[i*2+1] > 3.0*averr))                nrej++;
            if (!isbad[i])
                ngood++;
        }
        ngood -= nrej;
        if (nrej == 0 || ngood < nconst)
            break;
        for (i = 0; i < nxy; i++) {
            if (!isbad[i] && (wptr[i*2] > 3.0*averr || wptr[i*2+1] > 3.0*averr))                isbad[i] = 1;
        }
        niter++;
    }

    /* Set the output array */

    *cc = cir_malloc(6*sizeof(double));
    (*cc)[0] = a;
    (*cc)[1] = b;
    (*cc)[2] = c;
    (*cc)[3] = d;
    (*cc)[4] = e;
    (*cc)[5] = f;

    /* Tidy and exit */

    freespace(wptr);
    freespace(iwork);
    return;
}

static int plate6(double *xpos, double *ypos, double *xi, double *eta,
		  unsigned char *flag, int npts, double *a, double *b,
		  double *c, double *d, double *e, double *f) {
    double sx1sq,sy1sq,sx1y1,sx1x2,sy1x2;
    double sy1y2,sx1y2,xposmean,yposmean,ximean,etamean,xx1,yy1,xx2,yy2;
    int i,ngood,nbad;
    char msg[BUFSIZ];

    /* Is it worthwhile even being here? */

    (void)cir_sumbpm(flag,npts,&nbad,msg);
    ngood = npts - nbad;
    if (ngood < 2)
        return(CIR_FATAL);

    /* Initialise all the counters and summations */

    sx1sq = 0.0;
    sy1sq = 0.0;
    sx1y1 = 0.0;
    sx1x2 = 0.0;
    sy1x2 = 0.0;
    sy1y2 = 0.0;
    sx1y2 = 0.0;
    xposmean = 0.0;
    yposmean = 0.0;
    ximean = 0.0;
    etamean = 0.0;

    /* Find means in each coordinate system */

    (void)cir_dmean(xpos,flag,npts,&xposmean,msg);
    (void)cir_dmean(ypos,flag,npts,&yposmean,msg);
    (void)cir_dmean(xi,flag,npts,&ximean,msg);
    (void)cir_dmean(eta,flag,npts,&etamean,msg);

    /* Now accumulate the sums */

    for (i = 0; i < npts; i++) {
        if (!flag[i]) {
	    xx1 = xpos[i] - xposmean;
	    yy1 = ypos[i] - yposmean;
	    xx2 = xi[i] - ximean;
	    yy2 = eta[i] - etamean;
            sx1sq += xx1*xx1;
            sy1sq += yy1*yy1;
            sx1y1 += xx1*yy1;
            sx1x2 += xx1*xx2;
            sy1x2 += yy1*xx2;
            sy1y2 += yy1*yy2;
            sx1y2 += xx1*yy2;
        }
    }

    /* Solution for X */

    *a = (sx1y1*sy1x2 - sx1x2*sy1sq)/(sx1y1*sx1y1 - sx1sq*sy1sq);
    *b = (sx1x2*sx1y1 - sx1sq*sy1x2)/(sx1y1*sx1y1 - sx1sq*sy1sq);
    *c = -xposmean*(*a) - yposmean*(*b) + ximean;

    /* Now the solution for Y */

    *d = (sx1y1*sx1y2 - sy1y2*sx1sq)/(sx1y1*sx1y1 - sy1sq*sx1sq);
    *e = (sy1y2*sx1y1 - sy1sq*sx1y2)/(sx1y1*sx1y1 - sy1sq*sx1sq);
    *f = -xposmean*(*e) - yposmean*(*d) + etamean;

    /* Get outta here */

    return (CIR_OK);
}

static int plate4(double *xpos, double *ypos, double *xi, double *eta,
		  unsigned char *flag, int npts, double *a, double *b,
		  double *c, double *d, double *e, double *f) {
    double sx1sq,sy1sq,sx1x2,sy1x2,sy1y2,sx1y2,xposmean,yposmean;
    double ximean,etamean,xx1,yy1,xx2,yy2,det,num,denom,theta,mag;
    double stheta,ctheta;
    int i,ngood,nbad;
    char msg[BUFSIZ];

    /* Is it worthwhile even being here? */

    (void)cir_sumbpm(flag,npts,&nbad,msg);
    ngood = npts - nbad;
    if (ngood < 2)
        return(CIR_FATAL);

    /* Initialise all the counters and summations */

    sx1sq = 0.0;
    sy1sq = 0.0;
    sx1x2 = 0.0;
    sy1x2 = 0.0;
    sy1y2 = 0.0;
    sx1y2 = 0.0;
    xposmean = 0.0;
    yposmean = 0.0;
    ximean = 0.0;
    etamean = 0.0;

    /* Find means in each coordinate system */

    (void)cir_dmean(xpos,flag,npts,&xposmean,msg);
    (void)cir_dmean(ypos,flag,npts,&yposmean,msg);
    (void)cir_dmean(xi,flag,npts,&ximean,msg);
    (void)cir_dmean(eta,flag,npts,&etamean,msg);

    /* Now accumulate the sums */

    for (i = 0; i < npts; i++) {
        if (!flag[i]) {
	    xx1 = xpos[i] - xposmean;
	    yy1 = ypos[i] - yposmean;
	    xx2 = xi[i] - ximean;
	    yy2 = eta[i] - etamean;
            sx1sq += xx1*xx1;
            sy1sq += yy1*yy1;
            sx1x2 += xx1*xx2;
            sy1x2 += yy1*xx2;
            sy1y2 += yy1*yy2;
            sx1y2 += xx1*yy2;
        }
    }

    /* Compute the rotation angle */

    det = sx1x2*sy1y2 - sy1x2*sx1y2;
    if (det < 0.0) {
        num = sy1x2 + sx1y2;
        denom = -sx1x2 + sy1y2;
    } else {
        num = sy1x2 - sx1y2;
        denom = sx1x2 + sy1y2;
    }
    if (num == 0.0 && denom == 0.0) {
        theta = 0.0;
    } else {
        theta = atan2(num,denom);
        if (theta < 0.0)
            theta += TWOPI;
    }

    /* Compute magnification factor */

    ctheta = cos(theta);
    stheta = sin(theta);
    num = denom*ctheta  + num*stheta;
    denom = sx1sq + sy1sq;
    if (denom <= 0.0) {
        mag = 1.0;
    } else {
        mag = num/denom;
    }

    /* Compute coeffs */

    if (det < 0.0) {
        *a = -mag*ctheta;
        *e = mag*stheta;
    } else {
        *a = mag*ctheta;
        *e = -mag*stheta;
    }
    *b = mag*stheta;
    *d = mag*ctheta;
    *c = -xposmean*(*a) - yposmean*(*b) + ximean;
    *f = -xposmean*(*e) - yposmean*(*d) + etamean;

    /* Get outta here */

    return (CIR_OK);
}

static void freematch (matchstrct *m) {
    freespace(m->x1);
    freespace(m->x2);
    freespace(m->y1);
    freespace(m->y2);
}
    
static void freecat(catstrct *c) {
    freespace(c->x);
    freespace(c->y);
}

static void shiftcat(mywcs *pwcs, mywcs *rwcs, catstrct *c) {
    int i,n;
    float *x,*y;
    double xx,yy,xx2,yy2;

    /* Reference a few things for convenience */

    x = c->x;
    y = c->y;
    n = c->nrows;

    /* Loop for each object */

    for (i = 0; i < n; i++) {
	xx = (double)x[i];
	yy = (double)y[i];
	outloc(pwcs,xx,yy,rwcs,NULL,&xx2,&yy2);
	x[i] = (float)xx2;
	y[i] = (float)yy2;
    }
}

static int matchcat(catstrct *ref, catstrct *prog, int nx, int ny, 
		    matchstrct *match) {
    int nref,nprog,ngrid,ngrid2,ibest,ig,j,nmatch,k,jm;
    float *xref,*yref,*xprog,*yprog,srad=5.0,aveden,errlim,xoffbest,yoffbest;
    float xoff,yoff,x,y,*xoffs,*yoffs,xoffset,yoffset;
    char errmsg[BUFSIZ];

    /* Reference a few variables */

    xref = ref->x;
    yref = ref->y;
    nref = ref->nrows;
    xprog = prog->x;
    yprog = prog->y;
    nprog = prog->nrows;

    /* Sort the reference list by y */

    sorty(xref,yref,nref);

    /* Calculate error limit */

    aveden = (float)nref/(float)(nx*ny);
    errlim = 1.0/sqrt(4.0*M_PI*aveden);
    errlim = min(errlim,5.0);
    ngrid = (int)(srad/errlim);
    ngrid = (ngrid/2)*2 + 1;    
    ngrid = max(5,min(21,ngrid));
    ngrid2 = ngrid/2 + 1;

    /* Now search for the best solution */

    ibest = 0;
    xoffbest = 0.0;
    yoffbest = 0.0;
    for (ig = -ngrid2; ig <= ngrid2; ig++) {
	xoff = (float)ig*errlim*M_SQRT2;
	for (j = -ngrid2; j <= ngrid2; j++) {
	    yoff = (float)j*errlim*M_SQRT2;
	    nmatch = 0;
	    for (k = 0; k < nprog; k++) {                  
		x = xprog[k] + xoff;
		y = yprog[k] + yoff;
		jm = fndmatch(x,y,xref,yref,nref,errlim);
		if (jm > -1) {
		    nmatch++;
		}
	    }
	    if (nmatch > ibest) {
		ibest = nmatch;
		xoffbest = xoff;
		yoffbest = yoff;
	    }
	}
    }

    /* Get some workspace */

    xoffs = cir_malloc(nprog*sizeof(*xoffs));
    yoffs = cir_malloc(nprog*sizeof(*yoffs));

    /* Ok, go through once more and find a median difference */

    nmatch = 0;
    for (k = 0; k < nprog; k++) {
	x = xprog[k] + xoffbest;
	y = yprog[k] + yoffbest;
	jm = fndmatch(x,y,xref,yref,nref,errlim);
	if (jm > -1) {
	    xoffs[nmatch] = (xref[jm] - xprog[k]);
	    yoffs[nmatch] = (yref[jm] - yprog[k]);
	    nmatch++;
	}
    }
    if (nmatch > 0) {
        (void)cir_med(xoffs,NULL,nmatch,&xoffset,errmsg);
        (void)cir_med(yoffs,NULL,nmatch,&yoffset,errmsg);
    } else {
        xoffset = 0.0;
        yoffset = 0.0;
    }
    freespace(xoffs);
    freespace(yoffs);

    /* Now make a list of the matching objects using the best offset */

    match->x1 = cir_malloc(nprog*sizeof(double));
    match->y1 = cir_malloc(nprog*sizeof(double));
    match->x2 = cir_malloc(nprog*sizeof(double));
    match->y2 = cir_malloc(nprog*sizeof(double));
    nmatch = 0;
    for (k = 0; k < nprog; k++) {
	x = xprog[k] + xoffset;
	y = yprog[k] + yoffset;
	jm = fndmatch(x,y,xref,yref,nref,errlim);
	if (jm > -1) {
	    (match->x1)[nmatch] = (double)xref[jm];
	    (match->y1)[nmatch] = (double)yref[jm];
	    (match->x2)[nmatch] = (double)xprog[k];
	    (match->y2)[nmatch] = (double)yprog[k];
	    nmatch++;
	}
    }
    match->nrows = nmatch;
    return(CIR_OK);
}
    
static catstrct readcat(char *catfile) {
    int status,xcol,ycol,anynul;
    long nrows;
    fitsfile *fptr;
    catstrct c;

    /* Open the file and get the number of rows */

    status = 0;
    (void)fits_open_file(&fptr,catfile,READONLY,&status);
    (void)fits_get_num_rows(fptr,&nrows,&status);
    (void)fits_get_colnum(fptr,CASEINSEN,"X_coordinate",&xcol,&status);
    (void)fits_get_colnum(fptr,CASEINSEN,"Y_coordinate",&ycol,&status);

    /* Get the workspace you need */

    c.x = cir_malloc(nrows*sizeof(float));
    c.y = cir_malloc(nrows*sizeof(float));
    c.nrows = nrows;
    c.fname = catfile;

    /* Read the xy coordinates */

    (void)fits_read_col(fptr,TFLOAT,xcol,1,1,nrows,NULL,c.x,&anynul,&status);
    (void)fits_read_col(fptr,TFLOAT,ycol,1,1,nrows,NULL,c.y,&anynul,&status);
    (void)fits_close_file(fptr,&status);
    return(c);
}

static void cir_prov(int method, fitsfile *optr, char *fname, char **infiles, 
		     int nimages) {
    int ncard,status,i,nkey;
    char key[FLEN_KEYWORD],value[FLEN_VALUE],comment[FLEN_COMMENT];
    char card[FLEN_CARD];
    char *inclist[] = {"PROV*"};
    char *algorithm[] = {"Nearest neighbour","Bi-linear interpolation"};

    /* Delete any provenance keywords that are already there */

    status = 0;
    (void)fits_find_nextkey(optr,inclist,1,NULL,0,card,&status);
    while (status == 0) {
	(void)fits_get_keyname(card,key,&nkey,&status);
	(void)fits_delete_key(optr,key,&status);
        (void)fits_find_nextkey(optr,inclist,1,NULL,0,card,&status);
    }
    status = 0;

    /* Write the zeroth card... */

    ncard = 0;
    snprintf(key,FLEN_KEYWORD,"PROV%04d",ncard);
    snprintf(value,FLEN_VALUE,"%s: formed from %s of:",basename(fname),
	     algorithm[method]);
    snprintf(comment,FLEN_COMMENT,"Output file name and combination algorithm");
    (void)fits_update_key(optr,TSTRING,key,value,comment,&status);

    /* Now write a card for each image in the file */

    for (i = 0; i < nimages; i++) {
        ncard++;
        snprintf(key,FLEN_KEYWORD,"PROV%04d",ncard);
        snprintf(value,FLEN_VALUE,"%s",basename(infiles[i]));
        snprintf(comment,FLEN_COMMENT,"Card # %d",ncard);
        (void)fits_update_key(optr,TSTRING,key,value,comment,&status);
    }
}

static void writewcs(fitsfile *optr, mywcs *outwcs) {
    int status;
    double val;

    /* We don't need to update the projection information because we're
       just inheritiing it from the old image. So only update the CR and
       CD keywords */

    status = 0;
    (void)fits_update_key(optr,TDOUBLE,"CRPIX1",&(outwcs->crpix[0]),NULL,
			  &status);
    (void)fits_update_key(optr,TDOUBLE,"CRPIX2",&(outwcs->crpix[1]),NULL,
			  &status);
    val = (outwcs->crval[0])/DEGRAD;
    (void)fits_update_key(optr,TDOUBLE,"CRVAL1",&val,NULL,&status);
    val = (outwcs->crval[1])/DEGRAD;
    (void)fits_update_key(optr,TDOUBLE,"CRVAL2",&val,NULL,&status);
    val = (outwcs->cd[0])/DEGRAD;
    (void)fits_update_key(optr,TDOUBLE,"CD1_1",&val,NULL,&status);
    val = (outwcs->cd[1])/DEGRAD;
    (void)fits_update_key(optr,TDOUBLE,"CD1_2",&val,NULL,&status);
    val = (outwcs->cd[2])/DEGRAD;
    (void)fits_update_key(optr,TDOUBLE,"CD2_1",&val,NULL,&status);
    val = (outwcs->cd[3])/DEGRAD;
    (void)fits_update_key(optr,TDOUBLE,"CD2_2",&val,NULL,&status);
}

static void sorty(float *x, float *y, int n) {
    int iii,ii,i,ifin,j;
    float b1,b2;

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
            if (y[i] > y[j]) {
		b1 = x[j];
		b2 = y[j];
                while (1) {
		    x[j] = x[i];
		    y[j] = y[i];
                    j = i;
                    i = i - iii;
                    if (i < 0 || y[i] <= b2)
                        break;
                }
		x[j] = b1;
		y[j] = b2;
            }
        }
    }
}

/*

$Log: cir_imstack_cat.c,v $
Revision 1.11  2013/02/07 16:05:22  jim
Fixed bug in outloc caused by CRPIX values that are exact integers

Revision 1.10  2012-12-08 07:29:26  jim
Modified to allow for magnitude zeropoint scaling

Revision 1.9  2012-08-13 09:53:36  jim
Fixed routine outloc so that zpn is calculated correctly

Revision 1.8  2012-07-20 09:34:44  jim
fixed minor bug in sort routine

Revision 1.7  2012-01-04 10:18:09  jim
Fixed a small bug in the overlap region calculation

Revision 1.6  2011-09-16 11:56:35  jim
Fixed little scaling problem

Revision 1.5  2010-10-29 11:19:31  jim
Modified provenance routines to fix overruning FITS key values

Revision 1.4  2010-10-01 09:55:29  jim
Added NICOMB to output headers

Revision 1.3  2010-09-06 09:01:17  jim
Tidied up documentation

Revision 1.2  2010-08-06 08:27:11  jim
Fixed little bug where the weightsum wasn't re-calculated if doing seeing
weighting

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.10  2010/04/27 11:36:22  jim
Fixed missing declarations

Revision 1.9  2010/04/27 11:24:47  jim
Partially rewrote to cut down on memory requirements. Added seeing weighting,
but it's not been enabled yet.

Revision 1.8  2010/03/10 10:55:37  jim
Puts provenance keywords into headers of confidence maps

Revision 1.7  2010/02/22 10:33:01  jim
Fixed lots of little bugs

Revision 1.6  2010/02/11 21:58:38  jim
Fixed bug in rejection algorithm

Revision 1.5  2010/02/05 10:04:14  jim
Adds name of confidence map to header of stack

Revision 1.4  2009/12/23 18:55:08  jim
Now adds DRIBBLE header keyword

Revision 1.3  2009/09/29 09:22:21  jim
Added provenance

Revision 1.2  2009/09/14 12:55:32  jim
Fixed declaration of nx,ny in mywcs structure

Revision 1.1  2009/09/09 10:42:07  jim
First entry


*/
