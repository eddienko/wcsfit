/*

$Id: imcore_filter.c,v 1.3 2012/07/20 09:34:44 jim Exp $

*/

#include <stdlib.h>

#include "util.h"
#include "imcore.h"

/* Function prototypes */

static void sortm (float ia[], int ib[], int n);
static void quicksort (float x[], int point[], int l, int nfilt);


/* does bilinear median and linear filtering on background values */
 
extern void bfilt (float **xbuf, int nx, int ny) {
   float *ybuf, *save;
   int mfilt = 5, j, k;
 
/* Allocate temporary storage */
   ybuf = (float *) malloc(MAX(nx,ny) * sizeof(float));
   save = (float *) malloc((nx+1) * ny * sizeof(float));
/*    if(!ybuf || !save) bombout(1, "malloc"); */
 
/* median filter across */
   for(k = 0; k < ny; k++) {
     for(j = 0; j < nx; j++) {
       save[(nx+1)*k+j] = xbuf[k][j];
       ybuf[j] = xbuf[k][j];
     }
     filt1d(ybuf, nx, mfilt);
     for(j = 0; j < nx; j++) xbuf[k][j] = ybuf[j];
   }
 
/* and now down */
   for(k = 0; k < nx; k++) {
     for(j = 0; j < ny; j++) ybuf[j] = xbuf[j][k];
     filt1d(ybuf, ny, mfilt);
     for(j = 0; j < ny; j++)
/* make sure median filtered values are not large than original */
       if(save[(nx+1)*j+k] > -1000.0)
       xbuf[j][k] = MIN(save[(nx+1)*j+k], ybuf[j]);
   }

/* now repeat with linear filters across */
   for(k = 0; k < ny; k++) {
     for(j = 0; j < nx; j++) ybuf[j] = xbuf[k][j];
     hanning(ybuf, nx);
     for(j = 0; j < nx; j++) xbuf[k][j] = ybuf[j];
   }

/* and now down */
   for(k = 0; k < nx; k++) {
     for(j = 0; j < ny; j++) ybuf[j] = xbuf[j][k];
     hanning(ybuf, ny);
     for(j = 0; j < ny; j++) xbuf[j][k] = ybuf[j];
   }

/* Free temporary storage */
   free((void *) ybuf);
   free((void *) save);
}/* --------------------------------------- ------------------------------ */
 
/* does median filtering allowing for unmeasured entries */
 
void filt1d (float ybuf[], int mpt, int mfilt) {
   float *wbuf;
   int i, irc;
/* Allocate temporary storage */
   wbuf = (float *) malloc(mpt * sizeof(float));
/*    if(!wbuf) bombout(1, "malloc"); */
 
   irc = 0;
   for(i = 0; i < mpt; i++){
     if(ybuf[i] > -1000.0){
       wbuf[irc] = ybuf[i];
       irc++;
     }
   }
   if(irc == 0) {
     free((void *) wbuf);
     return;
   }
   median(wbuf, irc, mfilt);
   irc = 0;
   for(i = 0; i < mpt; i++){
     if(ybuf[i] > -1000.0){
       ybuf[i] = wbuf[irc];
       irc++;
     }
   }
   padext(ybuf, mpt);
/* Free temporary storage */
   free((void *) wbuf);
}


/* pads out array with missing points and linearly extrapolates the ends */
 
void padext (float x[], int n) {
   int i, j, ilow, ihih=0, ic;
   float xlow, xhih, slope, t1 ,t2;
/* elements <= 0.0 are treated as missing */
   i = 0;
   while(x[i] <= -1000.0) i++;
   ilow = i;
   for(i = ilow+1; i < n; i++){
     if(x[i] <= -1000.0) {
       ic = 1;
       while(i+ic < n-1 && x[i+ic] <= -1000.0) ic++;
       if(i+ic < n-1){
         xlow = x[i-1];
         xhih = x[i+ic];
         for(j = 0; j < ic; j++){
           t2 = ((float) j+1)/((float) ic+1);
           t1 = 1.0 - t2;
           x[i+j] = t1*xlow+t2*xhih;
         }
       }
     } else {
       ihih = i;
     }
   }
/* linear extrapolation of ends */
   if(ilow > 0){
     slope = x[ilow+1]-x[ilow];
     for(i = 0; i < ilow; i++) x[i] = x[ilow]-slope*(ilow-i);
   }
   if(ihih < n-1) {
     slope = x[ihih]-x[ihih-1];
     for(i = ihih+1; i < n; i++) x[i] = x[ihih]+slope*(i-ihih);
   }
}

/* performs linear filtering on array xbuf */

void hanning (float xbuf[], int npt) {
  float *ybuf;
  float sum = 0.0, xmns, xmnf;
  int nfilt = 3, i, il, ilow, nelem;

  if(npt <= nfilt)
    return;

  /* set first and last edges equal */
  il   = nfilt/2;
  ilow = MAX(3,nfilt/4);
  ilow = (ilow/2)*2 + 1;

  for(i = 0; i < ilow; i++)
    sum += xbuf[i];

  xmns = sum/((float) ilow);

  sum=0.0;
  for(i = 0; i < ilow; i++)
    sum += xbuf[npt-1-i];

  xmnf = sum/((float) ilow);

  /* allocate ybuf array */
  nelem = npt + nfilt;  /* Max. number of elements req'd */

  ybuf = (float *) malloc(nelem * sizeof(float));
/*   if(!ybuf) */
/*     bombout(1, "malloc"); */

  /* reflect edges before filtering */
  for(i = 0; i < il; i++) {
    ybuf[i] = 2.0 * xmns - xbuf[il+ilow-1-i];
    ybuf[npt+i+il] = 2.0 * xmnf - xbuf[npt-i-ilow-1];
  }

  for(i = 0; i < npt; i++)
    ybuf[i+il] = xbuf[i];

  /* do linear filtering on rest */
  for(i = 0; i < npt; i++)
    xbuf[i] = 0.25 * (ybuf[i] + 2.0 * ybuf[i+1] + ybuf[i+2]);  /* 1-2-1 Hanning weighting */

  free((void *) ybuf);
}

/* performs median filtering on array xbuf */

void median (float xbuf[], int npt, int nfilt) {
  float *ybuf, *array;
  float xmns, xmnf;
  int *point;
  int nfo2p1, i, il, ilow, j, jl, jh, nelem, l=0;

  if((nfilt/2)*2 == nfilt) nfilt++;
  if(npt <= nfilt) return;
  nfo2p1 = nfilt/2;

  /* allocate ybuf, array, point */
  nelem = npt + nfilt;  /* Max. number of elements req'd */
  ybuf = (float *) malloc(nelem * sizeof(float));
/*   if(!ybuf) */
/*     bombout(1, "malloc"); */
  array = (float *) malloc(nfilt * sizeof(float));
  point = (int *) malloc(nfilt * sizeof(int));
/*   if(!array || !point) */
/*     bombout(1, "malloc"); */

  /* set first and last edges equal */
  il   = nfilt/2;
  ilow = MAX(3, nfilt/4);
  ilow = (ilow/2)*2 + 1;

  for(i = 0; i < ilow; i++) array[i] = xbuf[i];
  sortm(array, point, ilow);
  xmns = array[ilow/2];

  for(i = 0; i < ilow; i++) array[i] = xbuf[npt-1-i];
  sortm(array, point, ilow);
  xmnf = array[ilow/2];

  /* reflect edges before filtering */
  for(i = 0; i < il; i++) {
    ybuf[i] = 2.0 * xmns - xbuf[il+ilow-1-i];
    ybuf[npt+i+il] = 2.0 * xmnf - xbuf[npt-i-ilow-1];
  }
  for(i = 0; i < npt; i++) ybuf[i+il] = xbuf[i];

  /* do median filtering on rest */
  for(i = 0; i < nfilt; i++) {
    array[i] = ybuf[i];
    point[i] = i+1;
  }

  sortm(array, point, nfilt);

  xbuf[0] = array[nfo2p1];
  jl = nfilt;
  jh = nfilt+npt-1;
  for(j = jl; j < jh; j++) {

    for(i = 0; i < nfilt; i++) {
      if(point[i] != 1) {
	point[i]--;
	continue;
      }
      point[i] = nfilt;
      array[i] = ybuf[j];
      l = i;
    }
    quicksort(array, point, l, nfilt);
    xbuf[j-jl+1] = array[nfo2p1];
  }

  /* Free temporary arrays */
  free((void *) point);
  free((void *) array);
  free((void *) ybuf);
}

static void sortm (float ia[], int ib[], int n) {
  int i,j, ii, jj, ifin, iu;
  float it;

  jj = 4;
  while(jj < n) jj = 2 * jj;
  jj = MIN(n,(3 * jj)/4 - 1);
  while(jj > 1) {
    jj = jj/2;
    ifin = n - jj;
    for(ii = 0; ii < ifin; ii++) {
      i = ii;
      j = i + jj;
      if(ia[i] <= ia[j]) continue;
      it = ia[j];
      iu = ib[j];
      do {
	ia[j] = ia[i];
	ib[j] = ib[i];
	j = i;
	i = i - jj;
	if (i < 0) break;
      } while(ia[i] > it);
      ia[j] = it;
      ib[j] = iu;
    }
  }
}

static void quicksort (float x[], int point[], int l, int nfilt) {
  float test, temp;
  int i, it, j, npt, ii;

  test = x[l];
  j = nfilt;
  for(i = 0; i < nfilt; i++) {
    if(i != l && test <= x[i]) {
      j = i;
      break;
    }
  }
  if(j - 1 == l) return;
  
  if(j < l) {
    temp = x[l];
    it = point[l];
    npt = l - j;
    for(i = 0; i < npt; i++) {
      ii = l - i - 1;
      x[ii+1] = x[ii];
      point[ii+1] = point[ii];
    }
    x[j] = temp;
    point[j] = it;
  }
  else if(j > l) {
    temp = x[l];
    it = point[l];
    j--;
    npt = j - l;
    if(npt != 0) {
      for(i = 0; i < npt; i++) {
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

$Log: imcore_filter.c,v $
Revision 1.3  2012/07/20 09:34:44  jim
fixed minor bug in sort routine

Revision 1.2  2010-10-05 09:16:35  jim
Swapped the order of an if test in padext to avoid possibly encountering
uninitialised values

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.4  2008/02/11 12:26:15  jim
Fixed dodgy memory allocation

Revision 1.3  2006/07/27 12:54:47  jim
Fixed problem where null value was being recorded as 0.0 rather than
-1000.0

Revision 1.2  2004/09/07 14:18:57  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:54:59  jim
New version for rewrite of imcore


*/
