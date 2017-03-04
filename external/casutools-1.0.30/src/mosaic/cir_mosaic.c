/*

$Id: cir_mosaic.c,v 1.10 2013/02/07 16:05:22 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>

#include "mosaic.h"
#include <tools.h>

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
    int        istan;
} mywcs;

typedef struct {
    char           *infile;
    char           *conf;
    float          global_sky;
    float          global_noise;
    double         tpa;
    double         tpd;
    float          exptime;
    int            nhdus;

    float          sky_extn;
    mywcs          *wcs_extn;
    int            nx_extn;
    int            ny_extn;
    float          *data_extn;
    short int      *conf_extn;
    int            cur_extn;
    float          magzpt;
    float          scalefac;
    float          skydiff;
} dstrct;

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
    float     *work;
} outstrct;

static int nfptrs = 0;
static dstrct *fileptrs = NULL;

static void cir_prov(int method, fitsfile *optr, char *fname, char **infiles, 
		     int nimages);
static void free_output(outstrct *o);
static void average(dstrct *dd, outstrct *outstr, int interp, int conflim);
static void normal(outstrct *o);
static void read_extn(dstrct *dd, fitsfile *iptr, fitsfile *cptr);
static void free_extn(dstrct *dd);
static void outloc(mywcs *win, double xin, double yin, mywcs *wout, 
		   double *xout, double *yout);
static void xytord(mywcs *win, double xin, double yin, double *ra, 
		   double *dec);
static void rdtoxy(mywcs *win, double ra, double dec, double *x, double *y);
static void output_files(char *out, char *outc, outstrct *outstr);
static void writewcs(fitsfile *fptr, mywcs *o);
static mywcs *getwcsinfo(fitsfile *fptr);
static int read_global(dstrct *dd, char *expkey, int skyflag, char *errmsg);
static float distort_corr(double x, double y, mywcs *w);
static void tidy();

/*+
 *  Name:
 *      cir_mosaic
 *
 *  Purpose:
 *      Mosaic a list of image extensions with WCS information into a 
 *      single output image
 *
 *  Description:
 *      A list of WCS corrected images extensions are mosaicked into a
 *      single output image using the on-board WCS information encoded in the 
 *      FITS header. Confidence maps (either one per image or a single one
 *      that will be used for all images) must be supplied. The user can
 *      modify the output sky in a number of ways (see below).
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infiles = char ** (Given)
 *          A list of input FITS extensions to be tiled. These must have
 *          a WCS defined in the header.
 *      confs = char ** (Given)
 *          A list of input FITS confidence maps. There must either be one
 *          for each image in the infiles list, or just a single confidence
 *          map which will be used by all of the input images.
 *      nimages = int (Given)
 *          The number of images in the infiles list
 *      nconfs = int (Given)
 *          The number of images in the confs list
 *      interp = int (Given)
 *          The interpolation method to be used. Currently this is either
 *          0: nearest neighbour or 1: bi-linear interpolation
 *      skyflag = int (Given)
 *          A flag for how you want to deal with the background values. 
 *          The following values are currently supported:
 *              0: Don't do anything. Just use the backgrounds as they are
 *                 in the input files
 *              1: Set the sky background in the final image to a value given
 *                 in the 'skywish_in' argument
 *              2: Set the sky background to the mean of the backgrounds for
 *                 the input images
 *      skywish_in = float (Given)
 *          If you want to define the final sky background level (skyflag==1)
 *          then this is the value you want it to have.
 *      expkey = char * (Given)
 *          The FITS keyword for the exposure time in the input header
 *      conflim = int (Given)
 *          A confidence limit so that any pixels with a confidence value 
 *          below this one don't get used.
 *      out = char * (Given)
 *          The name of the FITS file for the output tile.
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
 *      Copyright (C) 2010-2013 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_mosaic(char **infiles, char **confs, int nimages, int nconfs,
		      int interp, int skyflag, float skywish_in, char *expkey, 
		      int conflim, char *out, char *outc, int verbose,
		      char *errmsg) {
    int i,j,status,hdutype,ind1,ind,ii,jj,nhdu;
    char errs[BUFSIZ],errstr[BUFSIZ];
    fitsfile *fptr,*cptr;
    long naxis[2],naxisc[2];
    float expref,skywish,zpt1,exp1,sky1,scalefac,ds,*data;
    double x,y,distort;
    dstrct *dd;
    outstrct outstr;

    /* Is there any point in being here? */

    if (nimages == 0) {
	(void)sprintf(errmsg,"MOSAIC: No images to combine\n");
	return(CIR_FATAL);
    }

    /* Check value of nconfs. It has to be either one or the same as nimages */
    
    if (nconfs != 1 && nconfs != nimages) {
	sprintf(errmsg,
		"MOSAIC: Number of conf maps must be 1 or nimages\n");
	return(CIR_FATAL);
    }

    /* Test to see whether we can even open all of these files. If we can't 
       then get out of here right now */

    errs[0] = '\0';
    for (i = 0; i < nimages; i++) {
	status = 0;

	/* Check the input image */

	(void)fits_open_file(&fptr,infiles[i],READONLY,&status);
	if (status != 0) {
	    (void)sprintf(errstr,"%s can't be opened\n",infiles[i]);
	    strcat(errs,errstr);
	    closefits(fptr);
	    continue;
        }
	(void)fits_get_num_hdus(fptr,&nhdu,&status);
	nhdu--;

	/* Check the confidence map(s) */

	if (nconfs == 1) 
	    (void)fits_open_file(&cptr,confs[0],READONLY,&status);
	else
	    (void)fits_open_file(&cptr,confs[i],READONLY,&status);
	if (status != 0) {
	    (void)sprintf(errstr,"%s can't be opened\n",confs[i]);
	    strcat(errs,errstr);
	    closefits(cptr);
	    continue;
	}

	/* Loop for the image extensions */

	for (j = 1; j <= nhdu; j++) {
	    
	    /* The input image extension */

	    (void)fits_movabs_hdu(fptr,j+1,&hdutype,&status);
	    (void)fits_get_img_size(fptr,2,naxis,&status);
	    if (status != 0) {
		(void)sprintf(errstr,"%s[%d] isn't an image\n",infiles[i],j);
		naxis[0] = 0;
		naxis[1] = 0;
		strcat(errs,errstr);
		status = 0;
		continue;
	    }

	    /* The confidence map */

	    (void)fits_movabs_hdu(cptr,j+1,&hdutype,&status);
	    (void)fits_get_img_size(cptr,2,naxisc,&status);
	    if (status != 0) {
		(void)sprintf(errstr,"%s[%d] isn't an image\n",infiles[i],j);
		naxisc[0] = 0;
		naxisc[1] = 0;
		strcat(errs,errstr);
		status = 0;
		continue;
	    }
    
	    /* Do the axes match? */
    
	    if (naxis[0] != naxisc[0] && naxis[1] != naxisc[1]) {
		if (nconfs == 1) 
		    (void)sprintf(errstr,"%s[%d], %s[%d] dims don't match",
				  infiles[i],j,confs[0],j);
		else
		    (void)sprintf(errstr,"%s[%d], %s[%d] dims don't match",
				  infiles[i],j,confs[i],j);
		continue;
	    }
	}
	closefits(fptr);
	closefits(cptr);
    }
    if (strlen(errs)) {
	(void)sprintf(errmsg,"MOSAIC: Errors in input files:%s\n",errs);
	return(CIR_FATAL);
    }


    /* Allocate file struct array and fill it in. Start by opening the
       the input FITS images. We've already tested that they exist and will
       open, so there's no need to do that again */

    fileptrs = cir_malloc(nimages*sizeof(*fileptrs));
    nfptrs = nimages;
    expref = 1.0;
    for (i = 0; i < nimages; i++) {
	dd = fileptrs + i;
	dd->infile = infiles[i];
	if (verbose)
	    fprintf(stdout,"Reading %s...\n",dd->infile);
	if (nconfs == 1) 
	    dd->conf = confs[0];
	else
	    dd->conf = confs[i];
	if (read_global(dd,expkey,skyflag,errmsg) != CIR_OK) {
	    tidy();
	    return(CIR_FATAL);
	}
    }

    /* If you are using the global sky values as a means of defining the
       output sky level, then do that now (skyflag == 2). If skyflag == 1
       then the output sky level should be a value that you've given. If
       skywish == 0 then don't change the sky levels at all */

    if (skyflag == 2) {
	skywish = 0.0;
	for (i = 0; i < nimages; i++)
	    skywish += (fileptrs+i)->global_sky;
	skywish /= (float)nimages;
    } else if (skyflag == 1) {
	skywish = skywish_in;
    } else {
	skywish = 0.0;
    } 
    if (skyflag != 0) {
	for (i = 0; i < nimages; i++)
	    (fileptrs+i)->global_sky = skywish;
    }

    /* Set up the output files */

    if (verbose)
	fprintf(stdout,"Setting up output files...\n");
    output_files(out,outc,&outstr);

    /* Right, now loop for each file and extension and add them into 
       the output buffer. Start by reading the current image extension
       into memory */

    zpt1 = 0.0;
    exp1 = 1.0;
    sky1 = 0.0;
    for (i = 0; i < nimages; i++) {
	dd = fileptrs + i;
	status = 0;
	(void)fits_open_file(&fptr,dd->infile,READONLY,&status);
	(void)fits_open_file(&cptr,dd->conf,READONLY,&status);
	for (j = 1; j <= dd->nhdus; j++) {
	    if (verbose) 
		fprintf(stdout,"Adding in image %s[%d]...\n",dd->infile,j);
	    (void)fits_movabs_hdu(fptr,j+1,&hdutype,&status);
	    (void)fits_movabs_hdu(cptr,j+1,&hdutype,&status);
	    dd->cur_extn = j;
	    read_extn(dd,fptr,cptr);

	    /* Keep a record of the exposure time, the photometric
	       zero point and the sky level for the first image */

	    if (i == 0 && j == 1) {
		zpt1 = dd->magzpt;
		exp1 = dd->exptime;
		sky1 = dd->global_sky;
	    }
	    scalefac = (exp1/dd->exptime)*pow(10.0,0.4*(zpt1-(dd->magzpt)));
	    ds = sky1 - scalefac*(dd->global_sky);
	    dd->scalefac = scalefac;
	    dd->skydiff = ds;

	    /* Modify the data array so that they are all on the same
	       system as the first image. Also take out the astrometric
	       distortion */
	    
	    data = dd->data_extn;
	    if (skyflag != 0) {
		for (jj = 0; jj < dd->ny_extn; jj++) {
		    ind1 = jj*dd->nx_extn;
		    y = (double)(jj + 1);
		    for (ii = 0; ii < dd->nx_extn; ii++) {
			ind = ind1 + ii;
			x = (double)(ii + 1);
			distort = distort_corr(x,y,dd->wcs_extn);
			data[ind] = (data[ind] - dd->sky_extn)*distort + 
			    skywish;
		    }
		}
	    }

	    /* Now average it in */

	    average(dd,&outstr,interp,conflim);

	    /* Free up the workspace for this extension */

	    free_extn(dd);
	}

	/* Close up the FITS files */

	closefits(fptr);
	closefits(cptr);
    }

    /* Normalise the output maps */
    
    normal(&outstr);

    /* Write out the output FITS files */

    status = 0;
    (void)fits_write_img(outstr.optr,TFLOAT,1,outstr.nxo*outstr.nyo,
			 outstr.outdata,&status);
    (void)fits_write_img(outstr.ocptr,TSHORT,1,outstr.nxo*outstr.nyo,
			 outstr.outcdata,&status);
    writewcs(outstr.optr,outstr.outwcs);
    writewcs(outstr.ocptr,outstr.outwcs);
    if (interp) 
	(void)fits_update_key(outstr.optr,TSTRING,"DRIBBLE","bilinear",
			      "Interpolation method",&status);
    cir_prov(interp,outstr.optr,outstr.outname,infiles,nimages);
    cir_prov(interp,outstr.ocptr,outstr.outcname,confs,nconfs);
    (void)fits_update_key(outstr.optr,TSTRING,"CIR_CPM",outc,NULL,&status);

    /* Free up output workspace */

    free_output(&outstr);

    /* Free up file descriptors */

    tidy();
    return(CIR_OK);
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
	    

static void free_output(outstrct *o) {
    closefits(o->optr);
    closefits(o->ocptr);
    freespace(o->outwcs);
    freespace(o->outdata);
    freespace(o->outcdata);
    freespace(o->work);
}

static void average(dstrct *dd, outstrct *outstr, int interp, int conflim) {
    int i,j,npo,ind1,ind,outind,j1,j2,i1,i2,jj,ii,ind1o,indo,ixo,iyo;
    float *data,*odata,*work,wi,wj,scf,fcd,ds;
    double x_in,y_in,x_out,y_out;
    short int *cdata;

    /* Set up some convenience variables */

    data = dd->data_extn;
    cdata = dd->conf_extn;
    odata = outstr->outdata;
    scf = dd->scalefac;
    ds = dd->skydiff;

    /* Get some workspace for the sum of the confidence */

    npo = outstr->nxo*outstr->nyo;
    work = outstr->work;

    /* Have separate loops for nearest neighbour and interpolation to
       speed things up */

    if (interp == 0) {
        for (j = 0; j < dd->ny_extn; j++) {
	    ind1 = j*dd->nx_extn;
  	    y_in = (double)(j+1);
 	    for (i = 0; i < dd->nx_extn; i++) {
	        ind = ind1 + i;
		if (cdata[ind] <= conflim)
		    continue;
	        x_in = (double)(i+1);
		outloc(dd->wcs_extn,x_in,y_in,outstr->outwcs,&x_out,&y_out);
		ixo = cir_nint(x_out);
		iyo = cir_nint(y_out);
		if (ixo < 1 || ixo > outstr->nxo || iyo < 1 || 
		    iyo > outstr->nyo)
		    continue;
		outind = (iyo-1)*outstr->nxo + ixo - 1;
		fcd = (float)cdata[ind];
		odata[outind] += fcd*(data[ind]*scf + ds);
		work[outind] += fcd;
	    }
	}
    } else {
        for (j = 0; j < dd->ny_extn; j++) {
	    ind1 = j*dd->nx_extn;
  	    y_in = (double)(j+1);
 	    for (i = 0; i < dd->nx_extn; i++) {
	        ind = ind1 + i;
		if (cdata[ind] <= conflim)
		    continue;
	        x_in = (double)(i+1);
		outloc(dd->wcs_extn,x_in,y_in,outstr->outwcs,&x_out,&y_out);
		if (x_out < 1.0 || y_out < 1.0 || x_out >= outstr->nxo || y_out >= outstr->nyo) {
/* 		    fprintf(stderr,"%d %d %g %g %g %g\n",i,j,x_in,y_in,x_out,y_out); */
		    
/* 		    xytord(dd->wcs_extn,x_in,y_in,&x_out,&y_out); */
/* 		    fprintf(stderr,"%g %g***\n",x_out,y_out); */
/* 		    exit(1); */
		    continue;
		}
		j1 = (int)y_out;
		j2 = j1 + 1;
		i1 = (int)x_out;
		i2 = i1 + 1;
		fcd = (float)cdata[ind];
		for (jj = j1; jj <= j2; jj++) {
		    wj = 1.0 - fabs(y_out - (double)jj);
		    ind1o = (jj - 1)*outstr->nxo;
		    for (ii = i1; ii <= i2; ii++) {
			indo = ind1o + ii - 1;
		        wi = 1.0 - fabs(x_out - (double)ii);
			odata [indo] += wi*wj*fcd*(scf*data[ind] + ds);
			work[indo] += wi*wj*fcd;
		    }
		}
	    }
	}
    }
}

static void normal(outstrct *o) {
    int npo,i;
    float *odata,*work,renorm,val,junk;
    short int *ocdata;
    char errmsg[BUFSIZ];

    /* Get some convenience variables */

    npo = o->nxo*o->nyo;
    odata = o->outdata;
    work = o->work;
    ocdata = o->outcdata;

    /* Normalise the output image */

    for (i = 0; i < npo; i++)  
	if (work[i] != 0.0) 
	    odata[i] /= work[i];

    /* Work out a mean background that you can use to fill in the blank spots */

    (void)cir_qmedsig(odata,NULL,npo,3.0,3,-1000.0,65535.0,&val,&junk,errmsg);
    for (i = 0; i < npo; i++)
	if (work[i] == 0.0)
	    odata[i] = val;

    /* Now normalise the confidence map */

    (void)cir_qmedsig(work,NULL,npo,3.0,3,-1000.0,65535.0,&val,&junk,errmsg);
    renorm = 100.0/val;
    for (i = 0; i < npo; i++)
	ocdata[i] = max(0,min(1000,cir_nint(work[i]*renorm)));

}			
		 
static void read_extn(dstrct *dd, fitsfile *iptr, fitsfile *cptr) {
    int status,npts,anynul,i;
    float percorr,magzpt,junk;
    long naxis[2];
    char errmsg[BUFSIZ];

    /* Get the image size */

    status = 0;
    (void)fits_get_img_size(iptr,2,naxis,&status);
    dd->nx_extn = naxis[0];
    dd->ny_extn = naxis[1];

    /* Get some information from the header */

    (void)fits_read_key(iptr,TFLOAT,"PERCORR",&percorr,NULL,&status);
    if (status != 0) {
	status = 0;
	percorr = 1.0;
    } else {
	percorr = pow(10.0,(double)0.4*percorr);
    }
    (void)fits_read_key(iptr,TFLOAT,"MAGZPT",&magzpt,NULL,&status);
    if (status != 0) {
	status = 0;
	magzpt = 0.0;
    }
    dd->magzpt = magzpt;
    
    /* Read the image data and the confidence map */

    npts = naxis[0]*naxis[1];
    dd->data_extn = cir_malloc(npts*sizeof(float));
    dd->conf_extn = cir_malloc(npts*sizeof(short int));
    (void)fits_read_img(iptr,TFLOAT,1,npts,NULL,dd->data_extn,&anynul,
			&status);
    (void)fits_read_img(cptr,TSHORT,1,npts,NULL,dd->conf_extn,&anynul,
			&status);

    /* Modify the data array to take the pedestal offset into account */

    if (percorr != 1.0) {
	for (i = 0; i < npts; i++)
	    (dd->data_extn)[i] *= percorr;
    }
    
    /* Get the on-board WCS */

    dd->wcs_extn = getwcsinfo(iptr);

    /* Now work out the median sky background for this extension */

    (void)cir_qmedsig(dd->data_extn,NULL,npts,3.0,3,-1000.0,65535.0,
		      &(dd->sky_extn),&junk,errmsg);
}

static void free_extn(dstrct *dd) {

    freespace(dd->data_extn);
    freespace(dd->conf_extn);
    freespace(dd->wcs_extn);
}
    
static void outloc(mywcs *win, double xin, double yin, mywcs *wout, 
		   double *xout, double *yout) {
    double xt,yt,xi,eta,r,rfac,aa,tandec,denom,rp;
    int i;

    /* Do the conversion. First to standard coordinates in the frame of
       the input image */

    xt = xin - win->crpix[0];
    yt = yin - win->crpix[1];
    xi = win->cd[0]*xt + win->cd[1]*yt;
    eta = win->cd[2]*xt + win->cd[3]*yt;
    if (! win->istan) {
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
    } else {
	rfac = 1.0;
    }
    xi *= rfac;
    eta *= rfac;
    aa = atan(xi*win->secd/(1.0-eta*win->tand));
    if (xi != 0.0) 
	tandec = (eta+win->tand)*sin(aa)/(xi*win->secd);
    else
	tandec = (eta+win->tand)/(1.0 - eta*win->tand);

    /* Now form standard coordinates in the frame of the output image */

    aa += (win->crval[0] - wout->crval[0]);
    denom = wout->tand*tandec + cos(aa);
    xi = wout->secd*sin(aa)/denom;
    eta = (tandec - wout->tand*cos(aa))/denom;
    if (! wout->istan) {
	rp = sqrt(xi*xi + eta*eta);
	r = atan(rp);
        rfac = wout->pv21 + wout->pv23*pow(rp,2.0) + wout->pv25*pow(rp,4.0);
	rfac *= r/rp;
    } else {
	rfac = 1.0;
    }
    xi *= rfac;
    eta *= rfac; 
    denom = wout->cd[0]*wout->cd[3] - wout->cd[1]*wout->cd[2];
    *xout = (xi*wout->cd[3] - eta*wout->cd[1])/denom + wout->crpix[0];
    *yout = (eta*wout->cd[0] - xi*wout->cd[2])/denom + wout->crpix[1];
}

static void xytord(mywcs *win, double xin, double yin, double *ra, 
		   double *dec) {
    double xt,yt,xi,eta,r,rfac,aa;

    xt = xin - win->crpix[0];
    yt = yin - win->crpix[1];
    xi = win->cd[0]*xt + win->cd[1]*yt;
    eta = win->cd[2]*xt + win->cd[3]*yt;
    r = sqrt(xi*xi + eta*eta);
    if (! win->istan) {
	rfac = win->pv21 + win->pv23*pow(r,2.0) + win->pv25*pow(r,4.0);
	r /= rfac;
	rfac = win->pv21 + win->pv23*pow(r,2.0) + win->pv25*pow(r,4.0);
	xi /= rfac;
	eta /= rfac;
    } else {
	if (r == 0.0) 
	    rfac = 1.0;
	else
	    rfac = tan(r)/r;
	xi *= rfac;
	eta *= rfac;
    }
    aa = atan(xi*win->secd/(1.0-eta*win->tand));
    *ra = aa + win->crval[0];
    if (xi != 0.0) {
	*dec = atan((eta + win->tand)*sin(aa)/(xi*win->secd));
    } else {
	*dec = atan((eta + win->tand)/(1.0 - eta*win->tand));
    }
    if (*ra > 2.0*M_PI)
	*ra -= 2.0*M_PI;
    else if (*ra < 0.0)
	*ra += 2.0*M_PI;
    *ra /= DEGRAD;
    *dec /= DEGRAD;
}

static void rdtoxy(mywcs *win, double ra, double dec, double *x, double *y) {
    double xi,eta,rra,ddec,denom,rp,rfac;

    rra = ra*DEGRAD - win->crval[0];
    ddec = dec*DEGRAD;
    denom = win->tand*tan(ddec) + cos(rra);
    xi = win->secd*sin(rra)/denom;
    eta = (tan(ddec) - win->tand*cos(rra))/denom;
    rp = sqrt(xi*xi + eta*eta);
    if (win->istan) {
	if (rp == 0) 
	    rfac = 1.0;
	else 
	    rfac = atan(rp)/rp;
    } else {
	rfac = win->pv21 + win->pv23*pow(rp,2.0) + win->pv25*pow(rp,4.0);
    }
    xi *= rfac;
    eta *= rfac;
    denom = win->cd[0]*win->cd[3] - win->cd[1]*win->cd[2];
    *x = (xi*win->cd[3] - eta*win->cd[1])/denom + win->crpix[0];
    *y = (eta*win->cd[0] - xi*win->cd[2])/denom + win->crpix[1];
}

static void output_files(char *out, char *outc, outstrct *outstr) {
    double lowerleft[2],upperleft[2],lowerright[2],upperright[2],xout,yout;
    double xmin,ymin,ra,dec;
    fitsfile *fptr;
    int i,ixo,iyo,status,j,hdutype;
    long naxis[2];
    dstrct *dd;
    mywcs *outwcs,*refwcs,*curwcs;
    char errmsg[BUFSIZ];

    /* Set some preliminary info */

    outstr->outname = out;
    outstr->outcname = outc;

    /* Work out the output image size. Initialise the first element of 
       the results arrays since we're doing all this relative to the first 
       image */

    lowerleft[0] = 100000000.0;
    lowerleft[1] = 100000000.0;
    upperleft[0] = 100000000.0;
    upperleft[1] = -100000000.0;
    lowerright[0] = -10000000.0;
    lowerright[1] = 10000000.0;
    upperright[0] = -10000000.0;
    upperright[1] = -10000000.0;

    /* Loop for all images relative to the first one */

    refwcs = NULL;
    for (i = 0; i < nfptrs; i++) {
	dd = fileptrs + i;

	/* Open the image and loop for each hdu */

	status = 0;
	(void)fits_open_file(&fptr,dd->infile,READONLY,&status);
	for (j = 1; j <= dd->nhdus; j++) {
	    (void)fits_movabs_hdu(fptr,j+1,&hdutype,&status);
 	    (void)fits_get_img_size(fptr,2,naxis,&status);
	    curwcs = getwcsinfo(fptr);
	    if (i == 0 && j == 1) {
		refwcs = getwcsinfo(fptr);
		refwcs->istan = 1;
	    }
	    outloc(curwcs,(double)1.0,(double)1.0,refwcs,&xout,&yout);
	    lowerleft[0] = min(lowerleft[0],xout);
	    lowerleft[1] = min(lowerleft[1],yout);
	    outloc(curwcs,(double)1.0,(double)naxis[1],refwcs,&xout,&yout);
	    upperleft[0] = min(upperleft[0],xout);
	    upperleft[1] = max(upperleft[1],yout);
	    outloc(curwcs,(double)naxis[0],(double)1.0,refwcs,&xout,&yout);
	    lowerright[0] = max(lowerright[0],xout);
	    lowerright[1] = min(lowerright[1],yout);
	    outloc(curwcs,(double)naxis[0],(double)naxis[1],refwcs,&xout,&yout);
	    upperright[0] = max(upperright[0],xout);
	    upperright[1] = max(upperright[1],yout);
	    freespace(curwcs);
	}
	closefits(fptr);
    }

    /* Ok, what are the limits? */

    ixo = cir_nint(max(lowerright[0]-lowerleft[0],upperright[0]-upperleft[0])) 
	+ 5;
    iyo = cir_nint(max(upperright[1]-lowerright[1],upperleft[1]-lowerleft[1])) 
	+ 5;
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
    outstr->outdata = cir_calloc(ixo*iyo,sizeof(float));
    outstr->outcdata = cir_calloc(ixo*iyo,sizeof(short int));
    outstr->work = cir_calloc(ixo*iyo,sizeof(float));
    (void)sprintf(errmsg,"%s[1]",outc);
    (void)fits_update_key(outstr->optr,TSTRING,"CIR_CPM",errmsg,NULL,&status);

    /* Update the reference point for the WCS */

    outwcs = cir_malloc(sizeof(mywcs));
    outwcs->crpix[0] = refwcs->crpix[0] - xmin + 1.0;
    outwcs->crpix[1] = refwcs->crpix[1] - ymin + 1.0;
    outwcs->crval[0] = refwcs->crval[0];
    outwcs->crval[1] = refwcs->crval[1];
    for (i = 0; i < 4; i++)
	outwcs->cd[i] = refwcs->cd[i];
    outwcs->pv21 = 1.0;
    outwcs->pv23 = -0.33333333333;
    outwcs->pv25 = 0.0;
    outwcs->nx = ixo;
    outwcs->ny = iyo;
    outwcs->istan = 1;
    outwcs->tand = tan(outwcs->crval[1]);
    outwcs->secd = 1.0/cos(outwcs->crval[1]);
    outstr->outwcs = outwcs;

    /* Now see what the RA,Dec of the centre of the map is and shift the
       reference point to that place */

    xout = 0.5*(double)ixo;
    yout = 0.5*(double)iyo;
    xytord(outwcs,xout,yout,&ra,&dec);
    outwcs->crpix[0] = xout;
    outwcs->crpix[1] = yout;
    outwcs->crval[0] = ra*DEGRAD;
    outwcs->crval[1] = dec*DEGRAD;
    outwcs->tand = tan(outwcs->crval[1]);
    outwcs->secd = 1.0/cos(outwcs->crval[1]);
    freespace(refwcs);
}

static void writewcs(fitsfile *fptr, mywcs *o) {
    char key[16],sval[16],proj[4];
    int status,i,j,n;
    float val;
    double dval;

    /* Write out CTYPE */

    if (o->istan)
	(void)strcpy(proj,"TAN");
    else 
	(void)strcpy(proj,"ZPN");
    status = 0;
    for (i = 1; i <= 2; i++) {
	(void)sprintf(key,"CTYPE%d",i);
	if (i == 1) 
	    (void)sprintf(sval,"RA---%s",proj);
	else
	    (void)sprintf(sval,"DEC--%s",proj);
	(void)fits_update_key(fptr,TSTRING,key,sval,NULL,&status);
    }

    /* CRVAL, CRPIX */

    for (i = 1; i <= 2; i++) {
	(void)sprintf(key,"CRVAL%d",i);
	dval = (o->crval[i-1])/DEGRAD;
	(void)fits_update_key(fptr,TDOUBLE,key,&dval,NULL,&status);
    }
    for (i = 1; i <= 2; i++) {
	(void)sprintf(key,"CRPIX%d",i);
	dval = (o->crpix[i-1]);
	(void)fits_update_key(fptr,TDOUBLE,key,&dval,NULL,&status);
    }

    /* CD matrix */

    n = 0;
    for (j = 1; j <= 2; j++) {
	for (i = 1; i <= 2; i++) {
	    (void)sprintf(key,"CD%d_%d",j,i);
 	    dval = (o->cd[n])/DEGRAD;
	    (void)fits_update_key(fptr,TDOUBLE,key,&dval,NULL,&status);
	    n++;
	}
    }    /* PV vector */

    if (! o->istan) {
	val = (float)(o->pv21);
	(void)fits_update_key(fptr,TFLOAT,"PV2_1",&val,NULL,&status);
	val = (float)(o->pv23);
	(void)fits_update_key(fptr,TFLOAT,"PV2_3",&val,NULL,&status);
	val = (float)(o->pv25);
	(void)fits_update_key(fptr,TFLOAT,"PV2_5",&val,NULL,&status);
    } else {
	(void)fits_delete_key(fptr,"PV2_1",&status);
	status = 0;
	(void)fits_delete_key(fptr,"PV2_3",&status);
	status = 0;
	(void)fits_delete_key(fptr,"PV2_5",&status);
	status = 0;
    }
}

static mywcs *getwcsinfo(fitsfile *fptr) {
    mywcs *w;
    int i,status;
    char ctype[16];

    /* Copy over the relevant info */

    w = cir_malloc(sizeof(mywcs));
    status = 0;
    (void)fits_read_key(fptr,TSTRING,"CTYPE1",ctype,NULL,&status);
    if (strstr(ctype,"TAN") == NULL) 
	w->istan = 0;
    else
	w->istan = 1;
    (void)fits_read_key(fptr,TDOUBLE,"CRVAL1",&(w->crval[0]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CRVAL2",&(w->crval[1]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CRPIX1",&(w->crpix[0]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CRPIX2",&(w->crpix[1]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD1_1",&(w->cd[0]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD1_2",&(w->cd[1]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD2_1",&(w->cd[2]),NULL,&status);
    (void)fits_read_key(fptr,TDOUBLE,"CD2_2",&(w->cd[3]),NULL,&status);
    if (! w->istan) {
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
    } else {
	w->pv21 = 1.0;
	w->pv23 = -0.3333333;
        w->pv25 = 0.0;
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

static int read_global(dstrct *dd, char *expkey, int skyflag, char *errmsg) {
    int i,nhdu,status,hdutype,npts,nalloc,anynul,offset,j;
    long naxis[2];
    fitsfile *fptr;
    float *data,percorr;

    /* Initialise the output in case of problems */

    dd->global_sky = 0.0;
    dd->global_noise = 0.0;
    dd->tpa = 0.0;
    dd->tpd = 0.0;
    dd->exptime = 0.0;

    /* Now open the file. We have theoretically already tested whether we
       can or cannot do this, so we don't have to be quite so neurotic about
       error checking here */

    status = 0;
    (void)fits_open_file(&fptr,dd->infile,READONLY,&status);

    /* Read the exposure time keyword. If it doesn't exist, then signal an
       error and get out of here */

    (void)fits_read_key(fptr,TFLOAT,expkey,&(dd->exptime),NULL,&status);
    if (status != 0) {
	(void)sprintf(errmsg,
		      "MOSAIC: File %s doesn't have exposure keyword %s",
		      dd->infile,expkey);
	closefits(fptr);
	return(CIR_FATAL);
    }

    /* Get the number of HDUs */

    (void)fits_get_num_hdus(fptr,&nhdu,&status);
    nhdu--;
    dd->nhdus = nhdu;

    /* Loop for all HDUs and read the data in */

    data = NULL;
    nalloc = 0;
    for (i = 1; i <= nhdu; i++) {
	(void)fits_movabs_hdu(fptr,i+1,&hdutype,&status);
	(void)fits_get_img_size(fptr,2,naxis,&status);
	npts = naxis[0]*naxis[1];
	if (i == 1) {
	    offset = 0;
	    nalloc = npts;
	    data = cir_malloc(nalloc*sizeof(float));
	} else {
	    offset = nalloc;
	    nalloc += npts;
	    data = cir_realloc(data,nalloc*sizeof(float));
	}
	(void)fits_read_img(fptr,TFLOAT,1,npts,NULL,data+offset,&anynul,
			    &status);
	
	/* Read the header keyword with the background offset */

	(void)fits_read_key(fptr,TFLOAT,"PERCORR",&percorr,NULL,&status);
	if (status != 0) {
	    status = 0;
	} else {
	    percorr = (float)pow(10.0,0.4*percorr);
	    for (j = 0; j < npts; j++) 
		data[offset+j] *= percorr;
	}

	/* Get the tanget points from the first extension only */

	if (i == 1) {
    	    (void)fits_read_key(fptr,TFLOAT,"CRVAL1",&(dd->tpa),NULL,&status);
    	    (void)fits_read_key(fptr,TFLOAT,"CRVAL2",&(dd->tpd),NULL,&status);
	}
    }
   
    /* Get the global sky value now */

    (void)cir_qmedsig(data,NULL,nalloc,3.0,3,-1000.0,65535.0,
		      &(dd->global_sky),&(dd->global_noise),errmsg);

    /* Free up workspace and get out of here */

    closefits(fptr);
    freespace(data);
    return(CIR_OK);
}

static float distort_corr(double x, double y, mywcs *w) {
    double xi,xn,r,rprime,drprime_bydr,dc,tanr;

    xi = (x - w->crpix[0])*(w->cd[0]) + (y - w->crpix[1])*(w->cd[1]);
    xn = (x - w->crpix[0])*(w->cd[2]) + (y - w->crpix[1])*(w->cd[3]);
    r = sqrt(xi*xi + xn*xn);
    if (r == 0.0) {
	dc = 1.0;
    } else {
	if (! w->istan) {
	    rprime = r*w->pv21 + pow(r,(double)3.0)*w->pv23 + 
		pow(r,(double)5.0)*w->pv25;
	    drprime_bydr = w->pv21 + 3.0*pow(r,(double)2.0)*w->pv23 + 
		5.0*pow(r,(double)4.0)*w->pv25;
	    tanr = tan(rprime);
	    dc = tanr*drprime_bydr*(1.0 + tanr*tanr)/r;
	} else {
	    tanr = tan(r);
	    dc = (1.0 + tanr*tanr)*tanr/r;
	}
    }
    return(dc);
}

static void tidy() {
    int i;
    dstrct *dd;

    for (i = 0; i < nfptrs; i++) {
	dd = fileptrs + i;
	freespace(dd->wcs_extn);
	freespace(dd->data_extn);
	freespace(dd->conf_extn);
    }
    freespace(fileptrs);

}

/*

$Log: cir_mosaic.c,v $
Revision 1.10  2013/02/07 16:05:22  jim
Fixed bug in outloc caused by CRPIX values that are exact integers

Revision 1.9  2012-10-01 19:03:43  jim
fixed zeropoint bug

Revision 1.8  2012-09-04 10:43:59  jim
Fixed small memory leak caused by not closing a FITS file

Revision 1.7  2012-08-20 08:01:24  jim
Fixed flux distortion

Revision 1.6  2012-08-15 16:18:57  jim
Fixed missing scalefactor

Revision 1.5  2012-08-13 09:53:36  jim
Fixed routine outloc so that zpn is calculated correctly

Revision 1.4  2011-01-12 13:15:51  jim
Changed some mallocs to callocs

Revision 1.3  2010-10-29 11:19:31  jim
Modified provenance routines to fix overruning FITS key values

Revision 1.2  2010-09-21 10:52:34  jim
Added verbose output. Also fixed bug affecting the reading of the confidence
map when only 1 is present.

Revision 1.1  2010-09-06 09:03:01  jim
New entry

Revision 1.7  2010/06/03 08:43:21  jim
Fixed so that if the mapped output pixel is out of range it get's skipped

Revision 1.6  2010/05/20 10:41:22  jim
Fixed bug in call to add dribble keyword

Revision 1.5  2010/05/20 08:59:49  jim
Adds confidence map name into output mosaic header

Revision 1.4  2010/05/20 08:49:19  jim
Adds DRIBBLE keyword and provenance

Revision 1.3  2010/05/14 12:43:23  jim
Adds confidence map name into mosaic header

Revision 1.2  2010/05/13 12:24:59  jim
Small change to output wcs specification

Revision 1.1  2010/05/05 08:38:55  jim
initial entry


*/
