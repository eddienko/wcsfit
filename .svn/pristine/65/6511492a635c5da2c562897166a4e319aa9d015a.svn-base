/*

$Id: imcore_radii.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/

#include <stdio.h>
#include <math.h>

#include "imcore.h"
#include "imcore_radii.h"
#include "floatmath.h"
#include "util.h"
#include "ap.h"

extern float imcore_halflight(float rcores[], float cflux[], float halflight,
			      float peak, int naper) {
    float delr,halfrad;
    int i,gotone;

    /* Work out the half-light value from either isophotal flux or the 
       flux at Rcore. The find out roughly where the curve of growth
       exceeds this */

    gotone = 0;
    for (i = 0; i < naper; i++) {
	if (cflux[i] > halflight) {
	    gotone = 1;
	    break;
	}
    }
    if (! gotone) 
	i = naper - 1;
    
    /* Now work out what the radius of half light is */

    if (i == 0) {
	delr = (cflux[i] - halflight)/MAX(1.0,cflux[i]-peak);
	halfrad = rcores[0]*(1.0 - delr) + delr*sqrt(1.0/M_PI);
    } else {
        delr = (cflux[i] - halflight)/MAX(1.0,(cflux[i] - cflux[i-1]));
        halfrad = rcores[i-1]*delr + rcores[i]*(1.0-delr);
    }
    return(halfrad);
}   

extern float imcore_exprad(float thresh, float peak, float areal0, 
			   float rcores[], int naper) {
    float pk,r_t,rad;

    /* Work out the radius... */

    pk = MAX(1.5*thresh,peak);
    r_t = sqrtf(areal0/M_PI);
    rad = 5.0*r_t/logf(pk/thresh);
    rad = MAX(r_t,MIN(5.0*r_t,MIN(rad,rcores[naper-1])));
    return(rad);
}
    
extern float imcore_kronrad(float areal0, float rcores[], float cflux[], 
			    int naper) {
    int i,imax;
    float r_t,rad,wt,sum;

    /* Work out the radius... */

    r_t = sqrtf(areal0/M_PI);
    rad = 0.5*rcores[0]*cflux[0];
    sum = cflux[0];
    imax = MIN(naper,7);
    for (i = 1; i < imax; i++) {
	wt = MAX(0.0,cflux[i]-cflux[i-1]);
	rad += 0.5*(rcores[i] + rcores[i-1])*wt;
	sum += wt;
    }
    rad /= sum;
    rad = MAX(r_t,MIN(5.0*r_t,MIN(2.0*rad,rcores[naper-1])));    
    return(rad);
}

extern float imcore_petrad (float areal0, float rcores[], float cflux[], 
			    int naper) {
    int j;
    float eta,r_t,etaold,r1,r2,r3,r4,r5,r_petr;

    /* Work out petrosian radius */

    r_t = sqrtf(areal0/M_PI);
    eta = 1.0;
    etaold = eta;
    j = 1;
    while (eta > 0.2 && j < naper) {
	etaold = eta;
	r1 = rcores[j]*rcores[j]/(rcores[j-1]*rcores[j-1]) - 1.0;
	r2 = cflux[j]/cflux[j-1] - 1.0;
        eta = r2/r1;
	j++;
    }
    if (j == naper) {
	r_petr = rcores[naper-1];
    } else {
	r1 = rcores[j]*rcores[j];
	r2 = rcores[j-1]*rcores[j-1];
	r3 = rcores[j-2]*rcores[j-2];
        r4 = (etaold - 0.2)/(etaold - eta);
	r5 = (0.2 - eta)/(etaold - eta);
	r_petr = r4*sqrt(0.5*(r1 + r2)) + r5*sqrt(0.5*(r2 + r3));
    }
    r_petr = MAX(r_t,MIN(5.0*r_t,MIN(2.0*r_petr,rcores[naper-1])));
    return(r_petr);
}

void imcore_flux(ap_t *ap, float parm[IMNUM][NPAR], int nbit, float apers[], 
		 float fluxes[], int nr, float rcores[], float rfluxes[]) {
    float *map,t,xj,yj,sumiso,sumcf,delr;
    unsigned char *mflag,mf;
    long nx,ny;
    int xmin,xmax,ymin,ymax,ix1,ix2,iy1,iy2,i,j,kk,n;

    /* Set up some local variables */
 
    map = ap->data;
    mflag = ap->mflag;
    nx = ap->lsiz;
    ny = ap->csiz;

    /* Section for nbit == 1 */

    if (nbit == 1) {

	/* Generate image-blend outer boundaries */

	xmin = parm[0][1] - apers[0] - 0.5;
	xmax = parm[0][1] + apers[0] + 0.5;
	ymin = parm[0][2] - apers[0] - 0.5;
	ymax = parm[0][2] + apers[0] + 0.5;
	ix1 = MAX(0,(int)xmin-1);
	ix2 = MIN(nx-1,(int)xmax);
	iy1 = MAX(0,(int)ymin-1);
	iy2 = MIN(ny-1,(int)ymax);

	/* Now go through pixel region and add up the contributions inside
	   the aperture */

	fluxes[0] = 0.0;
	for(j = iy1; j <= iy2; j++) {
	    kk = j*nx;
	    for(i = ix1; i <= ix2; i++) {
		mf = mflag[kk+i];
		if (mf == MF_CLEANPIX || mf == MF_OBJPIX || 
		    mf == MF_SATURATED) {
		    t = map[kk+i];   
		    xj = (float)i - parm[0][1] + 1.0;
		    yj = (float)j - parm[0][2] + 1.0;
		    fluxes[0] += fraction(xj,yj,apers[0])*t;
		} 
	    }
	}
	if (fluxes[0] <= 0) 
	    fluxes[0] = parm[0][0];
	
    /* Section for blended images */
    
    } else {
    
	/* Interpolate circular aperture fluxes */
    
	sumiso = 0.0;
	sumcf = 0.0;
	for (j = 0; j < nbit; j++) {
	    sumiso += parm[j][0];
	    n = 1;
	    while (rcores[n] < apers[j] && n < nr-1)
		n++;
	    delr = (rcores[n] - apers[j])/(rcores[n] - rcores[n-1]);
	    fluxes[j] = rfluxes[j*nr+n]*(1.0 - delr) + rfluxes[j*nr+n-1]*delr;
	    sumcf += fluxes[j];
	}
/* 	fprintf(stderr,"***%g %g %g %g\n",sumiso,sumcf,fluxes[0],fluxes[1]); */

	/* Constrain the result so that the ratios are the same as for the
	   isophotal fluxes */

	for (j = 0; j < nbit; j++) {
	    fluxes[j] = sumcf*parm[j][0]/MAX(1.0,sumiso);
	    if (fluxes[j] < 0.0)
		fluxes[j] = parm[j][0];
/* 	    fprintf(stderr,"***%d %g %g\n",j,parm[j][0],fluxes[j]); */
	}
    }
}

/*

$Log: imcore_radii.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.11  2010/05/05 08:24:34  jim
Fixed imcore_halflight to cope with the situation where the half light level
is never found

Revision 1.10  2010/02/10 11:55:24  jim
Modified the calculation of half-light radius in case of seriously under-
sampled images

Revision 1.9  2009/12/17 11:34:33  jim
Added imcore_halflight

Revision 1.8  2009/01/30 11:02:51  jim
make sure wcslib include file is added correctly

Revision 1.7  2009/01/29 10:44:06  jim
new version of imcore_flux. Improved to clamp the ratio of fluxes for deblended
objects to that of the isophotal fluxes. Rewritten also so that loads of
unnecessary calculation is avoided for single isolated objects

Revision 1.6  2009/01/22 14:13:39  jim
Fixed the flagging of pixels in imcore_flux and fixed a bug in kron radius
determination

Revision 1.5  2008/10/21 13:05:41  jim
Modified imcore_flux to deal with possible instability

Revision 1.4  2007/06/04 10:34:03  jim
Modified to add list driven routines

Revision 1.3  2005/08/26 04:46:20  jim
Modified to add new radii and error estimates

Revision 1.2  2004/09/07 14:18:58  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:55:00  jim
New version for rewrite of imcore


*/
