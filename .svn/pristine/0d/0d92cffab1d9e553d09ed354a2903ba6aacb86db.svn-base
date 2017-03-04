/*

$Id: cir_filt1d.c,v 1.1 2012/12/08 07:31:00 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tools.h>

static void padext(float x[], int n);
static void median(float xbuf[], int npt, int nfilt);
static void linear(float *xbuf, int npt, int nfilt);
static void sortm(float ia[], int ib[], int n);
static void quicksort(float x[], int point[], int l, int nfilt);

static float BLANKVAL;

extern void cir_filt1d(float *ybuf, int mpt, int mfilt, int lfilt, 
		       float blankval) {
    float *wbuf;
    int i,irc;

    /* Allocate temporary storage */

    BLANKVAL = blankval;
    wbuf = cir_malloc(mpt*sizeof(float));

    /* Copy data over so long as it doesn't have the flag value */

    irc = 0;
    for (i = 0; i < mpt; i++){
        if (ybuf[i] != BLANKVAL){
            wbuf[irc] = ybuf[i];
            irc++;
        }
    }

    /* If they were all bad, then do nothing more */

    if (irc == 0) {
        freespace(wbuf);
        return;
    }

    /* Otherwise apply median filter */

    median(wbuf,irc,mfilt);

    /* Now copy it back */

    irc = 0;
    for (i = 0; i < mpt; i++) {
        if (ybuf[i] != BLANKVAL) {
            ybuf[i] = wbuf[irc];
	    irc++;
	}
    }
    padext(ybuf,mpt);
    linear(ybuf,mpt,lfilt);

    /* Tidy up */
    
    freespace(wbuf);
}
 
static void padext(float x[], int n) {
    int i,j,ilow,ihih=0,ic;
    float xlow,xhih,slope,t1,t2;

    i = 0;
    while(x[i] == BLANKVAL && i < n-1)
	i++;
    ilow = i;
    for (i = ilow+1; i < n; i++) {
	if (x[i] == BLANKVAL) {
	    ic = 1;
	    while (i+ic < n-1 && x[i+ic] == BLANKVAL)
		ic++;
	    if (i+ic < n-1) {
		xlow = x[i-1];
		xhih = x[i+ic];
		for (j = 0; j < ic; j++) {
		    t2 = ((float) j+1)/((float) ic+1);
		    t1 = 1.0 - t2;
		    x[i+j] = t1*xlow+t2*xhih;
		}
	    }
	} else {
	    ihih = i;
	}
    }

    /* Linear extrapolation of ends */

    if (ilow > 0) {
	slope = x[ilow+1] - x[ilow];
        for (i = 0; i < ilow; i++) 
	    x[i] = x[ilow] - slope*(float)(ilow-i);
    }
    if (ihih < n-1) {
	slope = x[ihih] - x[ihih-1];
        for (i = ihih+1; i < n; i++) 
	    x[i] = x[ihih] + slope*(float)(i-ihih);
    }
}

/* performs median filtering on array xbuf */

static void median(float xbuf[], int npt, int nfilt) {
    float *ybuf,*array;
    float xmns,xmnf;
    int *point;
    int nfo2p1,i,il,ilow,j,jl,jh,nelem,l=0;
    char errmsg[BUFSIZ];

    /* Make sure the filter number is odd and that there are enough
       data points */

    if (nfilt == 0)
	return;
    if ((nfilt/2)*2 == nfilt) 
	nfilt++;
    if (npt <= nfilt) 
	return;
    nfo2p1 = nfilt/2;

    /* Get some workspace */

    nelem = npt + nfilt;
    ybuf = cir_malloc(nelem*sizeof(float));
    array = cir_malloc(nfilt*sizeof(float));
    point = cir_malloc(nfilt*sizeof(int));

    /* Set first and last edges equal */

    il = nfilt/2;
    ilow = max(3, nfilt/4);
    ilow = (ilow/2)*2 + 1;
    for (i = 0; i < ilow; i++) 
	array[i] = xbuf[i];
    (void)cir_med(array,NULL,ilow,&xmns,errmsg);
    for (i = 0; i < ilow; i++) 
	array[i] = xbuf[npt-1-i];
    (void)cir_med(array,NULL,ilow,&xmnf,errmsg);

    /* Reflect edges before filtering */

    for (i = 0; i < il; i++) {
	ybuf[i] = 2.0*xmns - xbuf[il+ilow-1-i];
	ybuf[npt+i+il] = 2.0*xmnf - xbuf[npt-i-ilow-1];
    }
    for (i = 0; i < npt; i++) 
	ybuf[i+il] = xbuf[i];

    /* Do median filtering on rest */

    for (i = 0; i < nfilt; i++) {
	array[i] = ybuf[i];
	point[i] = i+1;
    }
    sortm(array,point,nfilt);
    xbuf[0] = array[nfo2p1];
    jl = nfilt;
    jh = nfilt+npt-1;
    for (j = jl; j < jh; j++) {
	for (i = 0; i < nfilt; i++) {
	    if (point[i] != 1) {
		point[i]--;
		continue;
	    }
	    point[i] = nfilt;
	    array[i] = ybuf[j];
	    l = i;
	}
	quicksort(array,point,l,nfilt);
	xbuf[j-jl+1] = array[nfo2p1];
    }

    /* Free temporary arrays */

    freespace(point);
    freespace(array);
    freespace(ybuf);
}

static void linear(float *xbuf, int npt, int nfilt) {
    int i,il,ilow,xmns,xmnf,nelem;
    float sum,*ybuf,fnfilt;

    /* Make sure you have an odd number of pixels in the filter window */

    if (nfilt == 0)
	return;
    if (! nfilt % 2) 
	nfilt++;
    if (npt <= nfilt)
	return;

    /* Set first and last edges equal */

    il = nfilt/2;
    ilow = max(3,nfilt/4);
    ilow = (ilow/2)*2 + 1;
    sum = 0.0;
    for (i = 0; i < ilow; i++)
	sum += xbuf[i];
    xmns = sum/((float)ilow);
    sum = 0.0;
    for (i = 0; i < ilow; i++)
	sum += xbuf[npt-1-i];
    xmnf = sum/(float)ilow;

    /* Allocate ybuf array to the maximum number of elements required */

    nelem = npt + nfilt;
    ybuf = cir_malloc(nelem*sizeof(float));

    /* Reflect edges before filtering */

    for (i = 0; i < il; i++) {
        ybuf[i] = 2.0*xmns - xbuf[il+ilow-1-i];
        ybuf[npt+i+il] = 2.0*xmnf - xbuf[npt-i-ilow-1];
    } 
    for (i = 0; i < npt; i++)
        ybuf[i+il] = xbuf[i];

    /* Do the filtering now */

    fnfilt = (float)nfilt;
    sum = 0.0;
    for (i = 0; i < nfilt; i++) 
	sum += ybuf[i];
    xbuf[0] = sum/fnfilt;
    for (i = 1; i < npt; i++) {
	sum += (ybuf[i+nfilt-1] - ybuf[i-1]);
	xbuf[i] = sum/fnfilt;
    }
    freespace(ybuf);
    return;
}	    

static void sortm(float ia[], int ib[], int n) {
    int i,j,ii,jj,ifin,iu;
    float it;

    jj = 4;
    while(jj < n) 
	jj = 2 * jj;
    jj = min(n,(3*jj)/4 - 1);
    while (jj > 1) {
	jj = jj/2;
	ifin = n - jj;
	for (ii = 0; ii < ifin; ii++) {
	    i = ii;
	    j = i + jj;
	    if(ia[i] <= ia[j]) 
		continue;
	    it = ia[j];
	    iu = ib[j];
	    do {
		ia[j] = ia[i];
		ib[j] = ib[i];
		j = i;
		i = i - jj;
		if (i < 0)
		    break;
	    } while(ia[i] > it);
	    ia[j] = it;
	    ib[j] = iu;
	}
    }
}


static void quicksort(float x[], int point[], int l, int nfilt) {
    float test,temp;
    int i,it,j,npt,ii;

    test = x[l];
    j = nfilt;
    for (i = 0; i < nfilt; i++) {
	if (i != l && test <= x[i]) {
	    j = i;
	    break;
	}
    }
    if (j - 1 == l) return;
  
    if (j < l) {
	temp = x[l];
	it = point[l];
	npt = l - j;
	for (i = 0; i < npt; i++) {
	    ii = l - i - 1;
	    x[ii+1] = x[ii];
	    point[ii+1] = point[ii];
	}
	x[j] = temp;
	point[j] = it;
    } else if (j > l) {
	temp = x[l];
	it = point[l];
	j--;
	npt = j - l;
	if (npt != 0) {
	    for (i = 0; i < npt; i++) {
		ii = l + i + 1;
		x[ii-1] = x[ii];
		point[ii-1] = point[ii];
	    }
	}
	x[j] = temp;
	point[j] = it;
    }
}

/*

$Log: cir_filt1d.c,v $
Revision 1.1  2012/12/08 07:31:00  jim
New entry


*/
