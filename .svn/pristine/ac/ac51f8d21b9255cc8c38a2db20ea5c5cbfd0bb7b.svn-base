/*

$Id: moments.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include "imcore.h"
#include "util.h"

extern void moments(ap_t *ap, float results[]) {
    int i,np;
    float x,y,xoff,yoff,xsum,ysum,xsumsq,ysumsq,tsum,xysum,t,tmax,twelfth;
    float xbar,ybar,sxx,syy,sxy,xintmin,w,wsum,xsum_w,ysum_w;
    plstruct *plarray;

    /* Initialise a few things */

    xintmin = ap->xintmin;
    plarray = ap->plarray;
    np = ap->npl_pix;
    xoff = (float)plarray[0].x;
    yoff = (float)plarray[0].y;
    xsum = 0.0;
    ysum = 0.0;
    xsum_w = 0.0;
    ysum_w = 0.0;
    wsum = 0.0;
    xsumsq = 0.0;
    ysumsq = 0.0;
    tsum = 0.0;
    xysum = 0.0;
    tmax = plarray[0].z;
    twelfth = 1.0/12.0;

    /* Do a moments analysis on an object */

    for (i = 0; i < np; i++) {
        x = (float)plarray[i].x - xoff;
        y = (float)plarray[i].y - yoff;
        t = plarray[i].z;
	w = plarray[i].zsm;
	if (t < 0.0)
	    continue;
        xsum += t*x;
        ysum += t*y;
	tsum += t;
	xsum_w += w*t*x;
	ysum_w += w*t*y;
	wsum += w*t;
	tmax = MAX(tmax,plarray[i].z);
	xsumsq += (x*x + twelfth)*t;
	ysumsq += (y*y + twelfth)*t;
	xysum += x*y*t;
    }

    /* Check that the total intensity is enough and if it is, then do
       the final results */

    if (tsum >= xintmin) {
        xbar = xsum/tsum;
	ybar = ysum/tsum;
	sxx = MAX(0.0,(xsumsq/tsum-xbar*xbar));
	syy = MAX(0.0,(ysumsq/tsum-ybar*ybar));
	sxy = xysum/tsum - xbar*ybar;
	xbar = xsum_w/wsum;
	ybar = ysum_w/wsum;
        xbar += xoff;
        ybar += yoff;
	xbar = MAX(1.0,MIN(xbar,ap->lsiz));
	ybar = MAX(1.0,MIN(ybar,ap->csiz));
        results[0] = 1.0;
	results[1] = xbar;
	results[2] = ybar;
	results[3] = tsum;
	results[4] = sxx;
	results[5] = sxy;
	results[6] = syy;
	results[7] = tmax;
    } else {
	results[0] = -1.0;
    }
}

/*

$Log: moments.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.5  2006/08/21 09:06:54  jim
Modified centring algorithm

Revision 1.4  2005/05/17 08:32:44  jim
small bug fix

Revision 1.3  2005/05/09 08:45:07  jim
Fixed so that x,y coordinates can't take on stupid values

Revision 1.2  2005/05/03 08:36:51  jim
Changed declaration of cirdr.h

Revision 1.1  2004/04/02 10:55:00  jim
New version for rewrite of imcore


*/
