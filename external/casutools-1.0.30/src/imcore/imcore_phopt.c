/*

$Id: imcore_phopt.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "floatmath.h"
#include "util.h"
#include "imcore.h"
#include "ap.h"

/* Function Prototypes */

static void dchole (double a[IMNUM+1][IMNUM+1], double b[IMNUM+1], int n);

/* Does multiple profile fitting to determine intensities */

void phopt(ap_t *ap, float parm[IMNUM][NPAR], int nbit, int naper, 
	   float apertures[], float cflux[], float badpix[], int nrcore,
           float avconf[]) {
    double aa[IMNUM+1][IMNUM+1],bb[IMNUM+1];
    float d,arg,*map,rcirc;
    float cn,parrad,xmin,xmax,xi,yi,ymin,ymax;
    float t,xj,yj,cnsq,tj,xk,yk,tk,ff;
    int i,ii,j,kk,ix1,ix2,iy1,iy2,nx,ny,k,iaper;
    unsigned char *mflag,mf;
    short int *conf;

    /* Set up some local variables */

    map = ap->data;
    conf = ap->conf;
    mflag = ap->mflag;
    nx = ap->lsiz;
    ny = ap->csiz;   
   
    /* Loop for each of the apertures */

    for (iaper = 0; iaper < naper; iaper++) {
	rcirc = apertures[iaper];
        parrad = rcirc + 0.5;
        cn = 1.0/(M_PI*rcirc*rcirc);       /* Profile normalising constant */
        cnsq = cn*cn;
    
	/* Set up covariance matrix - analytic special case for cores */

	for (i = 0; i < nbit; i++) {
	    aa[i][i] = cn;                 /* Overlaps totally area=pi*r**2 */
	    if (nbit > 1) {
		xi = parm[i][1];
		yi = parm[i][2];
		for (j = i+1; j < nbit; j++) {
		    d = sqrtf((xi-parm[j][1])*(xi-parm[j][1])
			+ (yi-parm[j][2])*(yi-parm[j][2]));
		    if (d >= 2.0*rcirc) {
			aa[j][i] = 0.0;
			aa[i][j] = aa[j][i];
		    } else {
			arg = d/(2.0*rcirc);
			aa[j][i] = cnsq*2.0*rcirc*rcirc*
			    (acosf(arg)-arg*(sqrtf(1.0-arg*arg)));
			aa[i][j] = aa[j][i];
		    }
		}
	    }
	}

	/* Clear accumulators */

	for (i = 0; i < nbit; i++) 
	    bb[i] = 0.0;

	/* Generate image-blend outer boundaries */

	xmin = 1.0e6;
	xmax = -1.0e6;
	ymin = 1.0e6;
	ymax = -1.0e6;
	for (i = 0; i < nbit; i++) {
	    xi = parm[i][1];
	    yi = parm[i][2];
	    xmin = MIN(xmin, xi);
	    xmax = MAX(xmax, xi);
	    ymin = MIN(ymin, yi);
	    ymax = MAX(ymax, yi);
	}
	ix1 = MAX(0,(int)(xmin-parrad)-1);
	ix2 = MIN(nx-1,(int)(xmax+parrad));
	iy1 = MAX(0,(int)(ymin-parrad)-1);
	iy2 = MIN(ny-1,(int)(ymax+parrad));

	/* Now go through pixel region */

	for (ii = iy1; ii <= iy2; ii++) {
	    kk = ii*nx;
	    for (i = ix1; i <= ix2; i++) {
		mf = mflag[kk+i];
		if (mf == MF_ZEROCONF || mf == MF_STUPID_VALUE) {
		    for (j = 0; j < nbit; j++) {
			xj = i - parm[j][1] + 1.0;
			yj = ii - parm[j][2] + 1.0;
			tj = fraction(xj,yj,rcirc);
			aa[j][j] -= tj*tj*cnsq;
			for (k = j + 1; k < nbit; k++) {
			    xk = i - parm[k][1] + 1.0;
			    yk = ii - parm[k][2] + 1.0;
			    tk = fraction(xk,yk,rcirc);
			    aa[k][j] -= tk*tj*cnsq;
			    aa[j][k] = aa[k][j];
			}
                        if (iaper == nrcore)
			    badpix[j] += tj;
		    }
		} else if (mf == MF_CLEANPIX || mf == MF_OBJPIX ||
			   mf == MF_SATURATED) {
		    t = map[kk+i];
		    for (j = 0; j < nbit; j++) {
			xj = i - parm[j][1] + 1.0;
			yj = ii - parm[j][2] + 1.0;
			ff = fraction(xj,yj,rcirc);
			bb[j] += ff*t;
			if (iaper == nrcore) 
			    avconf[j] += ff*(float)conf[kk+i];
		    }
		}
	    }
	}

	/* Trivial solution for single object */

	if (nbit == 1) {
	    cflux[iaper] = bb[0];

	/* Solve for profile intensities */

	} else {
	    for (i = 0; i < nbit; i++) 
		aa[i][i] = MAX(aa[i][i],cnsq);
  	    dchole (aa,bb,nbit);
	    for (i = 0; i < nbit; i++) 
	        cflux[i*naper+iaper] = cn*bb[i];
	}
    }
}

/* CHOLEsky decomposition of +ve definite symmetric matrix to solve Ax = b */

static void dchole (double a[IMNUM+1][IMNUM+1], double b[IMNUM+1], int n) {
  double sum, l[IMNUM+1][IMNUM+1], y[IMNUM+1];
  double aveigv, offset;
  int i, j, k;

restart:
    l[0][0] = sqrt(a[0][0]);

    for(k = 1; k < n; k++) {
        for(j = 0; j <= k-1; j++) {
	    sum = a[j][k];
	    if(j != 0) 
		for(i = 0; i <= j-1; i++) 
		    sum -= l[i][k]*l[i][j];
	    l[j][k] = sum/l[j][j];
        }
	sum = a[k][k];
	for(i = 0; i <= k-1; i++) 
	    sum -= l[i][k]*l[i][k];
	if(sum <= 0.0) {
/* 	    fprintf(stderr, "dchole: warning: matrix ill-conditioned\n"); */
	    aveigv = a[0][0];
	    for(i = 1; i < n; i++) 
		aveigv += a[i][i];
	    /* max eigenvalue < trace */
	    offset = 0.1*aveigv/((double) n);
	    for(i = 0; i < n; i++) 
		a[i][i] += offset;
/* 	    fprintf(stderr, "dchole: Offset added to diagonal = %f\n", offset); */
	    goto restart;
	}
	l[k][k] = sqrt(sum);
    }

    /* solve Ly = b */

    y[0] = b[0]/l[0][0];
    for(i = 1; i < n; i++) {
	sum = b[i];
	for(k = 0; k <= i-1; k++) 
	    sum -= l[k][i]*y[k];
	y[i] = sum/l[i][i];
    }

    /* solve L(T)x = y */

    b[n-1] = y[n-1]/l[n-1][n-1];
    for(i = n-2; i >= 0; i--) {
	sum = y[i];
	for(k = i+1; k < n; k++) 
	    sum -= l[i][k]*b[k];
	b[i] = sum/l[i][i];
    }
}

/* returns fraction of pixel bounded by 0 -  r_out
 * x,y coordinates relative to centre
 * Uses linear approximation ok if pixel located >>1 away from centre */

float fraction (float x, float y, float r_out) {
    float r,t,x_a,x_b,frac,tanao2,cosa,tanp2a,sqrt2o2;

    r = sqrtf(x*x + y*y);
    sqrt2o2 = 0.5*M_SQRT2;

    /* is it worth bothering? */

    if(r > r_out+sqrt2o2) 
	return(0.0);

    /* is it trivially all in? */

    if(r < r_out-sqrt2o2) 
	return(1.0);

    /* bugger - have to do some work then ... ok first ...
     * use 8-fold symmetry to convert to 0-45 degree range */

    x = fabsf(x);
    y = fabsf(y);
    if(y > x) {
	t = x;
	x = y;
	y = t;
    }

    /* If the angles are too close to cardinal points, then fudge something */

    if (x > 0.0 && y > 0.0) {
        tanao2 = 0.5*y/x;
        tanp2a = x/y;
        cosa = x/sqrt(x*x + y*y);
    } else {
        tanao2 = 0.00005;
        tanp2a = 10000.0;
        cosa = 1.0;
    }

    /* only outer radius - compute linear intersections top and bot of pixel */

    x_a = x - tanao2 + (r_out - r)/cosa;
    if(x_a < x+0.5) {

	/* intersects */

	x_b = x + tanao2 + (r_out - r)/cosa;

	/* three cases to consider */

	if(x_a < x-0.5)
	    frac = 0.5*MAX(0.0,x_b-(x-0.5))*MAX(0.0,x_b-(x-0.5))*tanp2a;
	else {
	    if(x_b > x+0.5)
		frac = 1.0 - 0.5*(x+0.5-x_a)*(x+0.5-x_a)*tanp2a;
	    else
		frac = 0.5-(x-x_a)+0.5*(x_b-x_a);
	}
    } else  /* missed entirely */
	frac = 1.0;

    return(frac);
}

/* calculates equivalent radius (along semi-major axis) for ellipse
 * Ratio is b/a */

/* static void checkpt (float x, float y, float x0, float y0, float ratio, float theta, float *radsq) { */
/*   float xx, yy, xnew, ynew; */

/*   xx = x-x0; */
/*   yy = y-y0; */
/*   xnew = xx*cos(theta) + yy*sin(theta); */
/*   ynew = -xx*sin(theta) + yy*cos(theta); */
/*   *radsq = xnew*xnew + (ynew/ratio)*(ynew/ratio); */
/* } */

/*

$Log: imcore_phopt.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.8  2009/12/17 11:32:06  jim
Modified to calculate the total confidence in the aperture

Revision 1.7  2009/01/22 14:10:19  jim
Modified how flagged pixels are dealt with...

Revision 1.6  2009/01/20 09:32:31  jim
Fixed bug with pixel flagging

Revision 1.5  2008/04/15 19:07:07  jim
Fixed code that accounts for bad pixels

Revision 1.4  2007/06/04 10:34:03  jim
Modified to add list driven routines

Revision 1.3  2005/08/26 04:46:20  jim
Modified to add new radii and error estimates

Revision 1.2  2004/09/07 14:18:57  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:54:59  jim
New version for rewrite of imcore


*/
