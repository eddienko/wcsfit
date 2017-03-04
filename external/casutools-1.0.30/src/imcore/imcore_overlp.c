/*

$Id: imcore_overlp.c,v 1.3 2010/10/05 09:16:56 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ap.h"
#include "util.h"
#include "imcore.h"
#include "floatmath.h"

#define NITER 6

static void moments_thr(ap_t *, float [], int []);
static void sort_on_zsm_rev(int, plstruct *);
static void update_ov(float [], float, float, float, float);
static void check_term(ap_t *, int *, float [IMNUM][NPAR+1], int [IMNUM][2], 
		       int *);

static float oldthr;
static float curthr;
static float nexthr;
static float lasthr;
static float xbar_start;
static float ybar_start;

void overlp(ap_t *ap, float parm[IMNUM][NPAR], int *nbit, float xbar,
            float ybar, float total, int npix, float tmax) {
    plstruct *pl;
    int npl,ipix,ipixo2,npl2,nbitprev,nobj,toomany,i,isnew,k,kk,j,ibitx[IMNUM];
    int ibity[IMNUM],ibitl[IMNUM],iwas,iupdate[IMNUM],npl3,iter,lastone,ic,jj;
    int ii,conv,ipks[IMNUM][2];
    float fconst,offset,tmul,smul,xintmn,itmaxlim,algthr,radmax,xb,yb,radius2;
    float results[IMNUM][NPAR+1],distmax,dx,dy,parmnew[IMNUM][NPAR],sumint;
    float xlevol,radold,slope,xx,xlevel,radius,xdat[NAREAL+1],xcor[NAREAL+1];
    float dlbydr,wt,dist,xeff,polycf[3],ttt,radthr,delb,deli,ratio;
    float bitx[IMNUM],bitl[IMNUM],sxx,syy;
    ap_t ap2;
    unsigned char *mflag;

    /* Initialise a few variables */

    pl = ap->plarray;
    npl = ap->npl_pix;
    ipix = ap->ipnop;
    oldthr = ap->thresh;
    fconst = ap->fconst;
    offset = ap->areal_offset;
    xbar_start = xbar;
    ybar_start = ybar;

    /* Initialise some constants that you might need later */

    tmul = 1.2589678;             /* 1/4 mag deblending contour increment */
    smul = 2.5;                   /* starting contour increment */
    ipixo2 = MAX(2,(ipix + 1)/2); /* ipix is the minimum image pixel size */
    xintmn = oldthr*ipixo2;       /* minimum intensity for fragments */
    itmaxlim = 0.9*tmax;          /* upper limit deblending 90% of peak */
    lasthr = itmaxlim;
    algthr = logf(oldthr);        /* Convenient bit of shorthand */
    radmax = sqrtf(((float)npix)/M_PI); /* Isophotal size of input data array */

    /* Get a maximum of IDBLIM points brighter than the new detection threshold
       by reverse sorting the input array. If there are still more than IDBLIM
       above the threshold, then revise the thresold until there aren't. Then
       use the variable npl2 to stop the rest of the routine accessing any of
       the fainter pixels in the list. This is to restrict processing time
       for large extended objects */
 
    curthr = smul*oldthr;
    sort_on_zsm_rev(npl,pl);
    while (1) {
        npl2 = 0;
  	while (pl[npl2].zsm > curthr && npl2 < npl-1)
	    npl2++;    
        if (npl2 > IDBLIM) 
            curthr += oldthr;
	else
	    break;
    }

    /* If there are fewer pixels above the new threshold than the minimum 
       specified in the input parameters, then you have no reason to be here */

    if (npl2 < ipix) {
        *nbit = 1;
        return;
    }

    /* Get a new ap structure */

    ap2.lsiz = ap->lsiz;
    ap2.csiz = ap->csiz;
    ap2.multiply = 1;
    ap2.ipnop = ipixo2;
    ap2.areal_offset = offset;
    ap2.fconst = fconst;
    mflag = calloc((ap2.lsiz)*(ap2.csiz),sizeof(*mflag));
    ap2.mflag = mflag;
    ap2.conf = ap->conf;

    /* Main analysis loop at new thresholds */

    *nbit = 0;
    nbitprev = 0;
    apinit(&ap2);
    while (1) {
        nexthr = MAX(curthr+oldthr,curthr*tmul);
        
        /* Locate objects in this cluster */

        ap2.thresh = curthr;
        apclust(&ap2,npl2,pl);
        check_term(&ap2,&nobj,results,ipks,&toomany);
	apreinit(&ap2);
        if (nobj == 0)
	    break;

        /* For each image check the number of points above the next threshold
           and flag. Do a moments analysis of each object */

	for (i = 0; i < nobj; i++) {

	    /* Ok, has this object already been detected?  If so, then
	       load the new results into the parmnew array. We'll check
               whether it's worthwhile updating the master parameters
	       list later */

            isnew = 1;
	    xb = results[i][1];
            yb = results[i][2];
	    sxx = MAX(1.0,results[i][4]);
	    syy = MAX(1.0,results[i][6]);
            for (k = 0; k < nbitprev; k++) {
		dx = xb - parm[k][1];
		dy = yb - parm[k][2];
		radius2 = dx*dx/sxx + dy*dy/syy;
		if ((ibitx[k] == ipks[i][0] && ibity[k] == ipks[i][1]) || 
		    radius2 < 1.0) {
		    isnew = 0;
		    for (kk = 0; kk < NPAR; kk++)
		        parmnew[k][kk] = results[i][kk];
		    ibitl[k] = (int)results[i][NPAR];
		    break;
		}
	    }

            /* If this is a new one and it's above the minimum threshold 
	       then check to make sure you don't already have too many. 
	       If you do, then flag this and break out to the next iteration.
	       If there is room for another, then store the moments analysis
	       profile */

	    if (isnew && results[i][0] > xintmn) {
		if (*nbit >= IMNUM) {
		    *nbit = IMNUM;
		    toomany = 1;
		    break;
		}
		ibitx[*nbit] = ipks[i][0];
		ibity[*nbit] = ipks[i][1];
                for (kk = 0; kk < NPAR; kk++)
		    parm[*nbit][kk] = results[i][kk];
		ibitl[*nbit] = (int)results[i][NPAR];
		(*nbit)++;
	    }
	} /* End of object loop */

        /* If too many objects were found, then skip the next bit...otherwise
	   go through and update parameters if necessary. This block of
	   code is a bit of a place holder waiting for something better to
	   be worked out...*/

	if (! toomany) {
            if (*nbit > nbitprev && nbitprev > 0) {
		for (i = 0; i < nbitprev; i++) 
		    iupdate[i] = 0;
		for (j = nbitprev; j < *nbit; j++) {
		    distmax = 0.0;
		    iwas = 0;
		    for (i = 0; i < nbitprev; i++) {
			if (parmnew[i][0] > 0.0) {
			    dx = parmnew[i][1] - parm[i][1];
			    dy = parmnew[i][2] - parm[i][2];
			    radius2 = dx*dx + dy*dy;
			    if (radius2 > distmax) {
				iwas = i;
				distmax = radius2;
			    }
			}
		    }
		    iupdate[iwas] = 1;
		}
		for (i = 0; i < nbitprev; i++)
		    if (iupdate[i] == 1 && parmnew[i][0] > 0.0) 
		        for (j = 0; j < NPAR; j++) 
			    parm[i][j] = parmnew[i][j];
	    }
	    
            /* Reset the update flag and prepare for next iteration*/

	    for (i = 0; i <= *nbit; i++)
		parmnew[i][0] = -1.0;
	    nbitprev = *nbit;
	}

        /* Where do we cut in the list now? */


        npl3 = 0;
  	while (pl[npl3].zsm > nexthr && npl3 < npl2-1)
	    npl3++;    
	npl2 = npl3;

	/* Now, do we need to move onto the next threshold? */

	if (npl2 == 0 || toomany || nexthr >= itmaxlim) 
	    break;

        /* If so, then reset some variables and continue */

	curthr = nexthr;

    } /* End of main analysis loop */

    /* Free up some workspace */

    free(ap2.mflag);
    apclose(&ap2);

    /* If there is only one then we can exit now */

    if (*nbit == 1)
	return;

    /* Find out which images terminated properly and remove those that didn't */

    j = -1;
    for (k = 0; k < *nbit; k++) {

	/* Commented this out as checking for terminations seems to miss some
	   if the total flux above the threshold for an object is negative */

/*         if (ibitl[k] == 1 && parm[k][0] > xintmn) { */ 

	if (parm[k][0] > xintmn) {
            j++;
	    if (j != k)
		for (i = 0; i < NPAR; i++)
		    parm[j][i] = parm[k][i];
	}
    }
    *nbit = j + 1;
    for (j = 0; j < *nbit; j++) {
	bitx[j] = 0.0;
	bitl[j] = 0.0;
    }

    /* For each image find true areal profile levels and iterate to find 
       local continuum */

    iter = 0;
    sumint = 0.0;
    lastone = 0;
    while (iter < NITER) {
	iter++;
	
	/* Loop for each of the objects and create a level vs radius array */

	for (k = 0; k < *nbit; k++) {
	    if (parm[k][0] < 0.0)
		continue;
	    xlevol = logf(parm[k][7] + parm[k][3] - bitl[k]); /* Pk + newthresh - cont */
	    xlevel = xlevol;
	    radold = 0.0;
	    radius = radold;
	    slope = 1.0;
	    ic = 0;
	    for (i = 1; i <= NAREAL; i++) {
		jj = NPAR - i;
		ii = NAREAL - i;
		xx = (float)ii + offset;
		if (parm[k][jj] > 0.5) {
		    if (ii == 0) 
			xlevel = logf(parm[k][3] - bitl[k] + 0.5);
		    else 
			xlevel = logf(powf(2.0,xx) - oldthr + parm[k][3] -
				      bitl[k] - 0.5);
		    radius = sqrt(parm[k][jj]/M_PI);
		    xdat[ic] = xlevel;
		    xcor[ic++] = radius;
		    dlbydr = (xlevol - xlevel)/MAX(0.01,radius-radold);
		    wt = MIN(1.0,MAX((radius-radold)*5.0,0.1));
		    slope = (1.0 - 0.5*wt)*slope + 0.5*wt*MIN(5.0,dlbydr);
		    radold = radius;
		    xlevol = xlevel;
		}
	    }

            /* If this is not the last iteration then work out the effect
	       on the local continuum from each image */

	    if (! lastone) {
		for (i = 0; i < *nbit; i++) {
		    if (parm[i][0] >= 0.0 && i != k) {
			dx = parm[k][1] - parm[i][1];
			dy = parm[k][2] - parm[i][2];
			dist = sqrtf(dx*dx + dy*dy);
			xeff = xlevel - MAX(0.0,MIN(50.0,slope*(dist-radius)));
			bitx[i] += expf(xeff);
		    }
		}

	    /* If this is the last iteration loop, then update the parameters
	       before exiting*/
	    
	    } else {
		if (ic > 2) {
		    polynm(xdat,xcor,ic,polycf,3,0);
		    ttt = polycf[1] + 2.0*polycf[2]*radius;
		} else
		    ttt = 0.0;
		slope = MAX(0.1,MAX(-ttt,slope));
		radthr = radius + (xlevel - algthr)/slope;
		if (radthr > radmax) {
		    slope = 1.0;
		    radthr = radmax;
		}
		
		/* Pixel area */
	       
		delb = parm[k][8]*(parm[k][3] - bitl[k]);
		parm[k][8] = M_PI*radthr*radthr;

		/* Peak height */

		parm[k][7] += (parm[k][3] - bitl[k]);
		
		/* Intensity */

		deli = 2.0*M_PI*((parm[k][3] - bitl[k])*(1.0 + slope*radius) -
				 oldthr*(1.0 + slope*radthr))/(slope*slope);
		parm[k][0] += delb + MAX(0.0,deli);
		for (i = 0; i < 7; i++)
		    parm[k][i+9] = -1.0;
		if (parm[k][0] > xintmn)
		    sumint += parm[k][0];
	    }
	}

        /* If this is not the last iteration then check and see how the
	   continuum estimates are converging. If they appear to be converging
	   then let the next iteration be the last one. */

	if (! lastone) {
	    conv = 1;
	    for (i = 0; i < *nbit; i++) {
		if (parm[i][0] >= 0.0) {
		    if (abs(bitx[i] - bitl[i]) > 3.0)
			conv = 0;
		    bitl[i] = bitx[i];
		    bitx[i] = 0;
		    bitl[i] = MIN(bitl[i],NINT(parm[i][3]-oldthr));
		}
	    }
	    lastone = (conv || (iter == NITER-1));
	} else {
	    break;
	}
    }

    /* Find the scaling if needed */

    if (sumint == 0.0) {
	*nbit = 1;
	return;
    } else 
	ratio = total/sumint;
    for (i = 0; i < *nbit; i++)
	parm[i][0] = ratio*parm[i][0];
}

        
static void sort_on_zsm_rev(int npts, plstruct *pts) {
    int i,j,ii,jj,ifin;
    plstruct tmp;

    jj = 4;
    while (jj < npts)
	jj = 2*jj;
    jj = MIN(npts,(3*jj)/4-1);
    while (jj > 1) {
	jj = jj/2;
	ifin = npts - jj;
	for (ii = 0; ii < ifin; ii++) {
	    i = ii;
	    j = i + jj;
	    if (pts[i].zsm > pts[j].zsm)
		continue;
	    tmp = pts[j];
	    do {
		pts[j] = pts[i];
		j = i;
		i = i - jj;
		if (i < 0) 
		    break;
	    } while (pts[i].zsm <= tmp.zsm);
	    pts[j] = tmp;
	}
    }
}

static void moments_thr(ap_t *ap, float results[NPAR+1], int ipk[2]) {
    int i,np,nnext;
    float x,y,xoff,yoff,xsum,ysum,xsumsq,ysumsq,tsum,xysum,t,tmax,twelfth;
    float xbar,ybar,sxx,syy,sxy,fconst,offset,xsum_w,ysum_w,wsum,w;
    plstruct *plarray;

    /* Copy some stuff to local variables */

    fconst = ap->fconst;
    offset = ap->areal_offset;
    plarray = ap->plarray;
    np = ap->npl_pix;

    /* Initialise a few things */

    xoff = xbar_start;
    yoff = ybar_start;
    xsum = 0.0;
    ysum = 0.0;
    xsum_w = 0.0;
    ysum_w = 0.0;
    wsum = 0.0;
    xsumsq = 0.0;
    ysumsq = 0.0;
    tsum = 0.0;
    xysum = 0.0;
    tmax = plarray[0].z - curthr;
    ipk[0] = plarray[0].x;
    ipk[1] = plarray[0].y;
    twelfth = 1.0/12.0;
    for (i = 8; i < NPAR; i++)
	results[i] = 0.0;

    /* Do a moments analysis on an object */

    nnext = 0;
    for (i = 0; i < np; i++) {
        x = (float)plarray[i].x - xoff;
        y = (float)plarray[i].y - yoff;
        t = plarray[i].z - curthr;
	w = plarray[i].zsm - curthr;
        if (w > nexthr)
	    nnext++;
        xsum += t*x;
        ysum += t*y;
	tsum += t;
	xsum_w += w*t*x;
	ysum_w += w*t*y;
	wsum += w*t;
	xsumsq += (x*x + twelfth)*t;
	ysumsq += (y*y + twelfth)*t;
	xysum += x*y*t;
	update_ov(results+8,t,oldthr,fconst,offset);
        if (t > tmax) {
	    ipk[0] = plarray[i].x;
	    ipk[1] = plarray[i].y;
	    tmax = t;
	}
    }

    /* Check that the total intensity is enough and if it is, then do
       the final results. Use negative total counts to signal an error */

    if (tsum > 0.0) {
        results[0] = tsum;
    } else {
	results[0] = -1.0;
        tsum = 1.0;
    }
    xbar = xsum/tsum;
    ybar = ysum/tsum;
    sxx = MAX(0.0,(xsumsq/tsum-xbar*xbar));
    syy = MAX(0.0,(ysumsq/tsum-ybar*ybar));
    sxy = xysum/tsum - xbar*ybar;
    wsum = MAX(1.0,wsum);
    xbar = xsum_w/wsum;
    ybar = ysum_w/wsum;
    xbar += xoff;
    ybar += yoff;
    xbar = MAX(1.0,MIN(xbar,ap->lsiz));
    ybar = MAX(1.0,MIN(ybar,ap->csiz));

    /* Store the results now */

    results[1] = xbar;
    results[2] = ybar;
    results[3] = curthr;
    results[4] = sxx;
    results[5] = sxy;
    results[6] = syy;
    results[7] = tmax;
    results[NPAR] = ((nnext > ap->ipnop && nexthr < lasthr) ? 0 : 1);
}
	
static void update_ov(float iap[NAREAL], float t, float thresh, float fconst, 
		      float offset) {
    int nup,i;

    /* Get out of here if the intensity is too small */

    if (t <= 0.0) 
	return;

    /* Otherwise update the relevant profile counts */

    nup = MAX(1,MIN(NAREAL,(int)(logf(t+thresh)*fconst-offset)+1));
    for (i = 0; i < nup; i++)
	iap[i] += 1.0;
}

static void check_term(ap_t *ap, int *nobj, float parm[IMNUM][NPAR+1], 
		       int peaks[IMNUM][2], int *toomany) {
    int ip,i,ipks[2];
    float momresults[NPAR+1];

    /* Search through all possible parents */

    *nobj = 0;
    *toomany = 0;
    for (ip = 1; ip <= ap->maxip; ip++) {
        if (ap->parent[ip].pnop != -1) {
/*             if (ap->parent[ip].pnop == ap->parent[ip].growing) { */
  
                /* That's a termination: */
 
                if ((ap->parent[ip].pnop >= ap->ipnop &&
		     ap->parent[ip].touch == 0)) {
		    extract_data(ap,ip);
 	 	    moments_thr(ap,momresults,ipks);
                    if (momresults[0] > 0.0) {
			if (*nobj == IMNUM-1) {
			    *toomany = 1;
			    break;
			}
			for (i = 0; i <= NPAR; i++) 
			    parm[*nobj][i] = momresults[i];
			for (i = 0; i < 2; i++)
			    peaks[*nobj][i] = ipks[i];
			(*nobj)++;
		    }
 	        }
	        restack(ap,ip);
/* 	    } else { */

/* 	        /\* This parent still active: *\/ */
 
/* 	        ap->parent[ip].growing = ap->parent[ip].pnop; */
/* 	    } */
        }
    }
}

/*

$Log: imcore_overlp.c,v $
Revision 1.3  2010/10/05 09:16:56  jim
Modified moments_thr to avoid possibly uninitialised values

Revision 1.2  2010-08-06 08:26:04  jim
Fixed slightly dodgy calloc for mflag

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.11  2009/01/22 14:03:10  jim
ap2.mflag is no longer a reference to the original

Revision 1.10  2008/04/15 19:05:24  jim
added reference to ap2.conf

Revision 1.9  2007/12/18 15:26:33  jim
Fixed bug in the bit that updates parm with parmnew so that an update isn't
done if the 'new' flux is -1 (i.e. no new information on the object is
available)

Revision 1.8  2006/08/21 09:06:54  jim
Modified centring algorithm

Revision 1.7  2006/06/29 13:29:45  jim
Modifications to smoothing kernel

Revision 1.6  2006/03/24 13:34:15  jim
a little fix to try and get around problems with noisy images that go
negative during at the high threshold levels

Revision 1.5  2005/05/17 08:32:44  jim
small bug fix

Revision 1.4  2005/05/09 08:45:07  jim
Fixed so that x,y coordinates can't take on stupid values

Revision 1.3  2004/09/07 14:18:57  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.2  2004/04/05 11:25:42  jim
Small modifications and bug fixes

Revision 1.1  2004/04/02 10:54:59  jim
New version for rewrite of imcore


*/
