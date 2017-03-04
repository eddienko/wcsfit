/*

$Id: cir_platesol.c,v 1.2 2010/08/02 16:09:28 jim Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/* #include <values.h> */
#include <limits.h>

#include <wcsfit.h>
#include <tools.h>

#define RAD2DEG 180.0/M_PI
#define TWOPI 2*M_PI
#define INITALLOC 64

static double *xptr = NULL,*yptr = NULL,*xiptr = NULL,*etaptr = NULL;
static double *xptr2 = NULL,*yptr2 = NULL,*wptr3 = NULL,*wptr4 = NULL;
static double *ra = NULL,*dec = NULL;
static unsigned char *isbad = NULL,*wptr2 = NULL;
static float *wptr = NULL;
static struct wcsprm *wcs = NULL;

static void tidy();

static int cir_plate6(double *, double *, double *, double *,
                      unsigned char *, int, double *, double *, double *,
                      double *, double *, double *);
static int cir_plate4(double *, double *, double *, double *, 
                      unsigned char *, int, double *, double *, double *,
                      double *, double *, double *);
/* static int dsolve(double *, double *, int); */

/*+
 *  Name:
 *      cir_platesol
 *
 *  Purpose:
 *      Do a plate solution for an image.
 *
 *  Description:
 *      Given an image and a list of astrometric coordinates matched to
 *      objects on that image, then do a plate solution. Write the WCS
 *      to the header of the image.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      image = char * (Given)
 *          The name of the input FITS image for the WCS fit.
 *      posfile = char * (Given)
 *          The name of the text file with the position info (both cartesian
 *          and equatorial) for any astrometric standards that happen to 
 *          be on the image.
 *      nconst = int (Given)
 *          The number of plate constants to fit. This is currently either
 *          4 or 6.
 *      pass = int (Given)
 *          The pass number for the WCS.  This signifies a level of refinement
 *          in the WCS and is only used to write to the header.
 *      fixtanpt = int (Given)
 *          If set, then the tangent point will not be moved from what is 
 *          already set in the current header.
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      There is a bug in dsolve which is why it currently isn't used.
 *              
 *  Dependencies:
 *      cfitsio, wcslib, cir_wcssubs.c, cur_stats.c
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *      Mike Irwin (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_platesol(char *image, char *posfile, int nconst, int pass, 
			int fixtanpt, char *errmsg) {
    int i,nalloc,npts,nrej,status,niter,ngood,nc2,retval;
    FILE *fd;
    float averr;
    double dxi,deta,daverr,xi,eta,r1,r2,d1,d2;
    double a,b,c,d,e,f,xifit,etafit,crpix1,crpix2,cd[4];
    double newcrval1,newcrval2;
    char msg[BUFSIZ],v1[16],v2[16],v3[16],v4[16],v5[16],v6[16];
    fitsfile *iptr;

    /* Open the fits image and read the current WCS info */

    retval = cir_wcsopen(image,&wcs,msg);
    if (retval != CIR_OK) {
        (void)sprintf(errmsg,"PLATESOL: Couldn't read file input WCS %s -- %s\n",
		      image,msg);
        tidy();
        return(CIR_FATAL);
    }

    /* Get some workspace for the positional arrays */

    nalloc = INITALLOC;
    xptr = cir_malloc(nalloc*sizeof(*xptr));
    yptr = cir_malloc(nalloc*sizeof(*yptr));
    xptr2 = cir_malloc(nalloc*sizeof(*xptr2));
    yptr2 = cir_malloc(nalloc*sizeof(*yptr2));
    ra = cir_malloc(nalloc*sizeof(*ra));
    dec = cir_malloc(nalloc*sizeof(*dec));

    /* Open the position file and read the relevant data. */

    if ((fd = fopen(posfile,"r")) == NULL) {
        (void)sprintf(errmsg,"PLATESOL: No input position file");
        tidy();
        return(CIR_FATAL);
    }
    npts = 0;
    while (fscanf(fd,"%s %s %s %s %s %s",v1,v2,v3,v4,v5,v6) != EOF) {
        xptr[npts] = atof(v1);
        yptr[npts] = atof(v2);
        xptr2[npts] = atof(v3);
        yptr2[npts] = atof(v4);
        ra[npts] = atof(v5);
        dec[npts] = atof(v6);
        npts++;
        if (npts == nalloc) {
            nalloc += INITALLOC;
	    xptr = cir_realloc(xptr,nalloc*sizeof(*xptr));
	    yptr = cir_realloc(yptr,nalloc*sizeof(*yptr));
	    xptr2 = cir_realloc(xptr2,nalloc*sizeof(*xptr2));
	    yptr2 = cir_realloc(yptr2,nalloc*sizeof(*yptr2));
	    ra = cir_realloc(ra,nalloc*sizeof(*ra));
	    dec = cir_realloc(dec,nalloc*sizeof(*dec));
        }
    }
    fclose(fd);
    nc2 = nconst/2;
    if (npts < nc2) {
        (void)sprintf(errmsg,"PLATESOL: Not enough position standards: %d",npts);
        tidy();
        return(CIR_FATAL);
    }

    /* Now if we want to move the tangent point do that now. */

    if (fixtanpt) {
        wptr3 = cir_malloc(npts*sizeof(*wptr3));
	wptr4 = cir_malloc(npts*sizeof(*wptr4));
        for (i = 0; i < npts; i += 1) {
            cir_xytoradec(wcs,xptr[i],yptr[i],&r1,&d1);
            cir_xytoradec(wcs,xptr2[i],yptr2[i],&r2,&d2);
            wptr3[i] = r2 - r1;
            wptr4[i] = d2 - d1;
        }
        (void)cir_dmed(wptr3,NULL,npts,&r1,msg);
        (void)cir_dmed(wptr4,NULL,npts,&d1,msg);
        freespace(wptr3);
        freespace(wptr4);
        newcrval1 = wcs->crval[0] + r1;
        newcrval2 = wcs->crval[1] + d1;
        status = 0;
        (void)fits_open_file(&iptr,image,READWRITE,&status);
        (void)fits_update_key(iptr,TDOUBLE,"CRVAL1",&newcrval1,NULL,&status);
        (void)fits_update_key(iptr,TDOUBLE,"CRVAL2",&newcrval2,NULL,&status);
        (void)fits_close_file(iptr,&status);
        if (status != 0) {
            fits_get_errstatus(status,msg);
            (void)sprintf(errmsg,"PLATESOL: Tangent point update failed %s",msg);
            tidy();
            return(CIR_FATAL);
        }
        freewcs(wcs);
	retval = cir_wcsopen(image,&wcs,msg);
    }

    /* Now calculate xi and eta for each object */

    isbad = cir_calloc(npts,sizeof(*isbad));
    xiptr = cir_malloc(npts*sizeof(*xiptr));
    etaptr = cir_malloc(npts*sizeof(*etaptr));
    for (i = 0; i < npts; i++) {
	cir_radectoxieta(wcs,ra[i],dec[i],&xi,&eta);
	xiptr[i] = xi;
	etaptr[i] = eta;
    }

    /* Get some workspace for the rejection cycles */

    wptr = cir_malloc(2*npts*sizeof(*wptr));
    wptr2 = cir_malloc(2*npts*sizeof(*wptr2));
     
    /* Right, now loop for maximum number of iterations or until
       convergence */    
   
    niter = 0;
    while (1) {

        /* Do a plate solution */

        switch (nconst) {
	case 6:
	    status = cir_plate6(xptr,yptr,xiptr,etaptr,isbad,npts,&a,&b,&c,
				&e,&d,&f);
            break;
        case 4:
	    status = cir_plate4(xptr,yptr,xiptr,etaptr,isbad,npts,&a,&b,&c,
				&e,&d,&f);
            break;
	default:
	    status = cir_plate6(xptr,yptr,xiptr,etaptr,isbad,npts,&a,&b,&c,
				&e,&d,&f);
            break;
        }
        if (status != CIR_OK) {
            (void)sprintf(errmsg,"PLATESOL: Solution failed");
            tidy();
            return(status);
        }

        /* Now look at the residuals and see if any should be rejected */

        for (i = 0; i < npts; i++) {
            xifit = xptr[i]*a + yptr[i]*b + c;
            etafit = xptr[i]*d + yptr[i]*e + f;
            dxi = fabs(xifit - xiptr[i]);
            deta = fabs(etafit - etaptr[i]);
            wptr[i*2] = dxi;
            wptr[i*2+1] = deta;
            wptr2[i*2] = isbad[i];
            wptr2[i*2+1] = isbad[i];
        }
        (void)cir_med(wptr,wptr2,2*npts,&averr,msg);
        averr *= 1.48;
        if (niter == 3)
            break;
        nrej = 0;
        ngood = 0;
        for (i = 0; i < npts; i++) {
 	    if (!isbad[i] && (wptr[i*2] > 3.0*averr || wptr[i*2+1] > 3.0*averr))
                nrej++;
            if (!isbad[i]) 
  	        ngood++;
        }
        ngood -= nrej;
        if (nrej == 0 || ngood < nconst)
            break;
        for (i = 0; i < npts; i++) {
 	    if (!isbad[i] && (wptr[i*2] > 3.0*averr || wptr[i*2+1] > 3.0*averr))
                isbad[i] = 1;
        }
        niter++;
    }

    /* If all objects were clipped out, then get out of here */

    if (ngood == 0) {
        (void)sprintf(errmsg,"PLATESOL: All objects clipped out");
        tidy();
        return(CIR_FATAL);
    }

    /* Convert values to degrees now */

    crpix1 = (e*c - b*f)/(d*b - e*a);
    crpix2 = (a*f - d*c)/(d*b - e*a);
    cd[0] = a*RAD2DEG;
    cd[1] = b*RAD2DEG;
    cd[2] = d*RAD2DEG;
    cd[3] = e*RAD2DEG;

    /* Number of good points fit and average error in arcsec*/

    ngood = 0;
    for (i = 0; i < npts; i++) 
        if (! isbad[i])
            ngood++;
    averr *= RAD2DEG*3600.0;
    daverr = (double)averr;
    
    /* Update the header now */

    status = 0;
    (void)fits_open_file(&iptr,image,READWRITE,&status);
    (void)fits_update_key(iptr,TDOUBLE,"CRPIX1",&crpix1,NULL,&status);
    (void)fits_update_key(iptr,TDOUBLE,"CRPIX2",&crpix2,NULL,&status);
    (void)fits_update_key(iptr,TDOUBLE,"CD1_1",&cd[0],NULL,&status);
    (void)fits_update_key(iptr,TDOUBLE,"CD1_2",&cd[1],NULL,&status);
    (void)fits_update_key(iptr,TDOUBLE,"CD2_1",&cd[2],NULL,&status);
    (void)fits_update_key(iptr,TDOUBLE,"CD2_2",&cd[3],NULL,&status);
    (void)fits_update_key(iptr,TINT,"NUMBRMS",&ngood,NULL,&status);
    (void)fits_update_key(iptr,TDOUBLE,"STDCRMS",&daverr,NULL,&status);
    (void)fits_update_key(iptr,TINT,"WCSPASS",&pass,NULL,&status);
    (void)fits_close_file(iptr,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"PLATESOL: Header update failed %s",msg);
        tidy();
        return(CIR_FATAL);
    }

    /* Free workspace */

    tidy();

    /* Get out of here */

    return(CIR_OK);
}

static void tidy() {

    freespace(wptr);
    freespace(wptr2);
    freespace(wptr3);
    freespace(wptr4);
    freespace(xptr);
    freespace(yptr);
    freespace(xptr2);
    freespace(yptr2);
    freespace(ra);
    freespace(dec);
    freespace(isbad);
    freespace(xiptr);
    freespace(etaptr);
    freewcs(wcs);
}

static int cir_plate6(double *xpos, double *ypos, double *xi, double *eta,
		      unsigned char *flag, int npts, double *a, double *b, 
                      double *c, double *d, double *e, double *f) {
/*     double xmatx[4],vect[2]; */
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

    /* Do solution for X */

/*      xmatx[0] = sx1sq; 1 */
/*      xmatx[1] = sx1y1; 2 */
/*      xmatx[2] = sx1y1; 4 */
/*      xmatx[3] = sy1sq; 5 */
/*      vect[0] = sx1x2;  3 */
/*      vect[1] = sy1x2;  6 */
/*      if (dsolve(xmatx,vect,2) != CIR_OK) */
/*          return(CIR_FATAL); */
/*      *a = vect[0]; */
/*      *b = vect[1]; */
/*      *a = (sx1y1*sy1x2 - sx1x2*sy1sq)/(sx1y1*sx1y1 - sx1sq*sy1sq); */
/*      *b = (sx1x2*sx1y1 - sx1sq*sy1x2)/(sx1y1*sx1y1 - sx1sq*sy1sq); */
/*      *c = -xposmean*(*a) - yposmean*(*b) + ximean; */
    *a = (sx1y1*sy1x2 - sx1x2*sy1sq)/(sx1y1*sx1y1 - sx1sq*sy1sq);
    *b = (sx1x2*sx1y1 - sx1sq*sy1x2)/(sx1y1*sx1y1 - sx1sq*sy1sq);
    *c = -xposmean*(*a) - yposmean*(*b) + ximean;

    /* Now the solution for Y */

/*      xmatx[0] = sy1sq; 1 */
/*      xmatx[1] = sx1y1; 2 */
/*      xmatx[2] = sx1y1; 4 */
/*      xmatx[3] = sx1sq; 5 */
/*      vect[0] = sy1y2; 3 */
/*      vect[1] = sx1y2; 6 */
/*      if (dsolve(xmatx,vect,2) != CIR_OK) */
/*          return(CIR_FATAL); */
/*      *d = vect[0]; */
/*      *e = vect[1]; */
/*      *d = (sx1y1*sx1y2 - sy1y2*sx1sq)/(sx1y1*sx1y2 - sy1sq*sx1sq); */
/*      *e = (sy1y2*sx1y1 - sy1sq*sx1y2)/(sx1y1*sx1y2 - sy1sq*sx1sq); */
/*      *f = -xposmean*(*e) - yposmean*(*d) + etamean; */
    *d = (sx1y1*sx1y2 - sy1y2*sx1sq)/(sx1y1*sx1y1 - sy1sq*sx1sq);
    *e = (sy1y2*sx1y1 - sy1sq*sx1y2)/(sx1y1*sx1y1 - sy1sq*sx1sq);
    *f = -xposmean*(*e) - yposmean*(*d) + etamean;

    /* Get outta here */

    return (CIR_OK);
}

static int cir_plate4(double *xpos, double *ypos, double *xi, double *eta,
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

/* static int dsolve(double *xmatx, double *vect, int n) { */
/*     int i,iu,k,l,j,jl,ir; */
/*     double big,rmax,temp,pivot; */

/*     iu = n - 1; */
/*     for (i = 0; i < iu; i++) { */

        /* Find largest remaining term in this column for pivot */

/*         big = 0.0; */
/*         for (k = i; k < n; k++) { */
/*             rmax = fabs(xmatx[k + i*n]); */
/*             if (rmax > big) { */
/*                 big = rmax; */
/*                 l = k; */
/*             } */
/* 	} */

        /* If this is zero, then return an error */

/*         if (big == 0.0)  */
/*             return(CIR_FATAL); */

        /* If this isn't the current row then swap rows back */

/*         if (i != l) { */
/*   	   for (j = 0; j < n; j++) { */
/*                temp = xmatx[i + j*n]; */
/*                xmatx[i + j*n] = xmatx[l + j*n]; */
/*                xmatx[l + j*n] = temp; */
/*            } */
/*            temp = vect[i]; */
/*            vect[i] = vect[l]; */
/*            vect[l] = temp; */
/*         } */

        /* Pivotal reduction */

/*         pivot = xmatx[i + i*n]; */
/*         jl = i + 1; */
/*         for (j = jl; j < n; j++) { */
/*             temp = xmatx[j + i*n]/pivot; */
/*             vect[j] -= temp*vect[i]; */
/*             for (k = i; k < n; k++)  */
/*                 xmatx[j + k*n] -= temp*xmatx[i + k*n]; */
/*         } */
/*     } */

    /* Back substitution for solution */

/*     for (i = 0; i < n; i++) { */
/*         ir = n - i - 1; */
/*         temp = vect[ir]; */
/*         if (ir != (n - 1)) { */
/*             for (j = 1; j < i; j++) { */
/*                 k = n - j; */
/*                 temp -= xmatx[ir + k*n]*vect[k]; */
/* 	    } */
/*         } */
/*         vect[ir] = temp/xmatx[ir + ir*n]; */
/*     } */
/*     return(CIR_OK); */
/* } */
/*

$Log: cir_platesol.c,v $
Revision 1.2  2010/08/02 16:09:28  jim
Uses <limits.h> now instead of <values.h>

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.13  2008/02/11 12:23:47  jim
Fixed sign error in plate4

Revision 1.12  2007/06/07 14:06:31  jim
Fixed little bug in one of the error messages

Revision 1.11  2007/04/04 05:22:19  jim
Modified to trap for condition of all the objects being clipped out during
the iteration loop.

Revision 1.10  2005/08/09 11:04:52  jim
*** empty log message ***

Revision 1.9  2004/09/07 14:18:54  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.8  2004/08/19 11:34:24  jim
Added new cir_memory routines for memeory allocation

Revision 1.7  2004/08/02 11:50:02  jim
Modified to use new version of cir_wcssubs routines

Revision 1.6  2003/10/01 20:35:18  jim
Modified so that the iterations don't continue if the number of objects
is less than the number of plate constants.

Revision 1.5  2003/09/11 14:16:19  jim
Fixed bug in which calculated the number of good entries rather than the number of bad ones

Revision 1.4  2003/09/11 11:58:09  jim
Modified to use revised statistics routines

Revision 1.3  2003/02/03 09:32:36  jim
Added history logging

Revision 1.2  2002/12/16 10:25:06  jim
Added prologues

Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS


*/
