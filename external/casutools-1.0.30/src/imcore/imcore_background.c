/*

$Id: imcore_background.c,v 1.3 2014/07/31 12:45:16 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "floatmath.h"
#include "util.h"
#include "errcodes.h"
#include "imcore.h"

static int **hist = NULL;
static int *nnp = NULL;
static int npvx;
static int npvy;
static void tidy();
static void sortit (float [], int);

int imcore_background(ap_t *ap, int nbsize, float nullval, int verbose,
		      char *errstr) {
    float fracx,fracy,skymed,sigma,skymedc,sigmac,avsky,fnbsize,dely,delx;
    float t1,t2,dsky,*map,**bvals,*work;
    int ifracx,ifracy,nbsizx,nbsizy,nbx,nby,npixstripe,l,i,ll;
    int isquare,ilev,j,iclip,mcpix,iloop,irej,nbsizo2,kk,k,iby,ibyp1,ibx,ibxp1;
    int *shist;
    unsigned char *mflag;
    long nx,ny;

    /* Set up some variables */

    map = ap->data;
    mflag = ap->mflag;
    nx = ap->lsiz;
    ny = ap->csiz;

    /* check to see if nbsize is close to exact divisor */

    fracx = ((float)nx)/((float)nbsize);
    fracy = ((float)ny)/((float)nbsize);
    ifracx = (int)(fracx + 0.1);
    ifracy = (int)(fracy + 0.1);
    nbsizx = nx/ifracx;
    nbsizy = ny/ifracy;
    nbsize = MAX(NINT(0.9*nbsize), MIN(nbsize, MIN(nbsizx,nbsizy)));
    nbsize = MIN(nx,MIN(ny,nbsize)); /* trap for small maps */

    /* Divide the map into partitions */

    nbx = nx/nbsize;
    nby = ny/nbsize;
    npixstripe = nbsize*nx;
    npvx = nbx;
    npvy = nby;

    /* Get histogram workspace if you can */

    hist = malloc(nbx*sizeof(int *));
    if (hist == NULL) {
        sprintf(errstr,"unable to allocate workspace in imcore_background");
	tidy();
	return(ERRCODE_MALLOC);
    }    
    for (l = 0; l < nbx; l++) {
	hist[l] = malloc(MAXHIST*sizeof(int));
	if (hist[l] == NULL) {
	    sprintf(errstr,"unable to allocate workspace in imcore_background");
	    tidy();
	    return(ERRCODE_MALLOC);
	}
    }

    /* Same for background values array */

    bvals = ap->backmap.bvals;
    if (bvals == NULL) {
        bvals = malloc(nby*sizeof(float *));
        if (bvals == NULL) {
            sprintf(errstr,"unable to allocate workspace in imcore_background");
 	    tidy();
 	    return(ERRCODE_MALLOC);
        }    
        for (l = 0; l < nby; l++) {
	    bvals[l] = malloc(nbx*sizeof(float));
	    if (bvals[l] == NULL) {
	        sprintf(errstr,"unable to allocate workspace in imcore_background");
	        tidy();
	        return(ERRCODE_MALLOC);
	    }
	}
    }

    /* Store some of this away for use later */

    ap->backmap.nbx = nbx;
    ap->backmap.nby = nby;
    ap->backmap.nbsize = nbsize;
    ap->backmap.bvals = bvals;

    /* Finally a counter array */

    nnp = malloc(nbx*sizeof(int));
    if (nnp == NULL) {
        sprintf(errstr,"unable to allocate workspace in imcore_background");
	tidy();
	return(ERRCODE_MALLOC);
    }    

    /* Loop for each row of background squares. Start by initialising
       the accumulators and histograms */

    for(l = 0; l < nby; l++) {
	memset((char *)nnp,0,nbx*sizeof(*nnp));
        for (i = 0; i < nbx; i++) 
   	    memset((char *)hist[i],0,MAXHIST*sizeof(int));

        /* Skim through the data in this stripe. Find out which square each
	   belongs to and add it it to the relevant histogram */

	ll = l*npixstripe;
        for (i = 0; i < npixstripe; i++) {
	    if (map[ll+i] != nullval && mflag[ll+i] != MF_ZEROCONF &&
		mflag[ll+i] != MF_STUPID_VALUE) {
		isquare = (int)((float)(i % nx)/(float)nbsize);
		isquare = MIN(nbx-1,MAX(0,isquare));
		ilev = MIN(MAXHISTVAL,MAX(MINHISTVAL,NINT(map[i+ll])));
		hist[isquare][ilev-MINHISTVAL] += 1;
		nnp[isquare] += 1;
	    }
	}
  
        /* but only do background estimation if enough pixels ----------- */

	for(j = 0; j < nbx; j++) {
	    if(nnp[j] > 0.25*nbsize*nbsize){
		shist = hist[j];
		imcore_medsig(shist,MAXHIST,MINHISTVAL-1,nnp[j],&skymed,&sigma);

                /* do an iterative 3-sigma upper clip to give a more robust
		   estimator */

		iclip = MAXHISTVAL;
		mcpix = nnp[j];
		skymedc = skymed;
		sigmac = sigma;
		for(iloop = 0; iloop < 3; iloop++) {
		    irej = 0;
 	            for(i = NINT(skymedc+3.0*sigmac); i <= iclip; i++)
                        irej += shist[i-MINHISTVAL];
		    if (irej == 0)
			break;
		    iclip = NINT(skymedc+3.0*sigmac) - 1;
		    mcpix = mcpix - irej;
		    imcore_medsig(shist,MAXHIST,MINHISTVAL-1,mcpix,&skymedc,
				  &sigmac);
		}
		bvals[l][j] = skymedc;
	    } else {
		bvals[l][j] = -1000.0;
	    }
	}
    }

    /* Write background values out */

    if (verbose) {
	printf("\nBackground values:\n");
        for(l = nby-1; l >= 0; l--)
	    for(j = 0; j < nbx; j++) printf("%8.1f", bvals[l][j]);
	printf("\n");
    }

    /* filter raw background values */

    bfilt(bvals,nbx,nby);
    if (verbose) {
	printf("\nFiltered background values:\n");
        for(l = nby-1; l >= 0; l--)
	    for(j = 0; j < nbx; j++) 
		printf("%8.1f", bvals[l][j]);
	printf("\n\n");
    }

    /* compute average sky level */

    work = malloc(nbx*nby*sizeof(*work));
    k = 0;
    for(l = 0; l < nby; l++)
        for(j = 0; j < nbx; j++) 
	    work[k++] = bvals[l][j];
    sortit(work,nbx*nby);
    avsky = work[(nbx*nby)/2];
    free(work);

    /* ok now correct map for background variations and put avsky back on */

    nbsizo2 = nbsize/2;
    fnbsize = 1.0/((float)nbsize);
    for(k = 0; k < ny; k++) {
	kk = k*nx;

	/* Nearest background pixel vertically */

        iby = (k + 1 + nbsizo2)/nbsize;
	ibyp1 = iby + 1;
	iby = MIN(nby,MAX(1,iby));
	ibyp1 = MIN(nby,ibyp1);
	dely = (k + 1 - nbsize*iby + nbsizo2)*fnbsize;

	for(j = 0; j < nx; j++) {
	    if (map[kk+j] == nullval) 
		continue;

            /* nearest background pixel across */

	    ibx = (j + 1 + nbsizo2)/nbsize;
	    ibxp1 = ibx + 1;
	    ibx = MIN(nbx,MAX(1,ibx));
	    ibxp1 = MIN(nbx,ibxp1);
	    delx = (j + 1 - nbsize*ibx + nbsizo2)*fnbsize;

            /* bilinear interpolation to find background */
	    
	    t1 = (1.0 - dely)*bvals[iby-1][ibx-1] + dely*bvals[ibyp1-1][ibx-1];
	    t2 = (1.0 - dely)*bvals[iby-1][ibxp1-1] + dely*bvals[ibyp1-1][ibxp1-1];
	    dsky = avsky - (1.0 - delx)*t1 - delx*t2;
            map[kk+j] += dsky;
	}
    }

    /* Free some workspace */

    tidy();
    return(ERRCODE_OK);
}

int imcore_backstats(ap_t *ap, float nullval, int satonly, float *skymed, 
		     float *skysig, float *sat, char *errstr) {
    int ilev,iclip,iloop,i,*ihist,isat,iter;
    long mpix,npix,k,mcpix,irej,lpix,nx,ny;
    float sigsq,skymedc,sigmac,*map,sata,fac,skyref;
    unsigned char *mflag;

    /* Get some info from the ap structure */

    map = ap->data;
    nx = ap->lsiz;
    ny = ap->csiz;
    mflag = ap->mflag;

    /* Check to make sure there are some non-zero values here */

    ilev = 1;
    for (i = 0; i < nx*ny; i++) {
	if (map[i] != nullval && mflag[i] != MF_ZEROCONF &&
	    mflag[i] != MF_STUPID_VALUE) {
	    ilev = 0;
	    break;
	}
    }
    if (ilev == 1) {
	*skymed = 0.0;
	*skysig = 0.0;
	*sat = 0.0;
	sprintf(errstr,"all values appear to be zero");
	return(ERRCODE_FILE_DATA);
    }

    /* First, get some workspace for the background histogram */

    ihist = calloc(MAXHIST,sizeof(*ihist));
    if (ihist == NULL) {
        sprintf(errstr,"unable to allocate workspace in imcore_backstats");
        return(ERRCODE_MALLOC);
    }

    /* Loop for up to 10 iterations. For each iteration we multiply the 
       input data by a successively higher power of 2 in order to 
       try and deal with data that has very small noise estimates */

    fac = 0.5;
    skyref = 0.0;
    for (iter = 0; iter <=9; iter++) {
	if (iter == 1)
	    skyref = skymedc;
	fac *= 2.0;
	for (k = 0; k < MAXHIST; k++)
	    ihist[k] = 0;

	/* Now form the histogram of all pixel intensities */

	mpix = 0;
	isat = 0;
	npix = nx*ny;
	for (k = 0; k < npix; k++) {
	    if (map[k] != nullval && mflag[k] != MF_ZEROCONF &&
		mflag[k] != MF_STUPID_VALUE) {
		ilev = MIN(MAXHISTVAL,MAX(MINHISTVAL,NINT(fac*(map[k]-skyref))));
		ihist[ilev - MINHISTVAL] += 1;
		isat = MAX(isat,ilev);
		mpix++;
	    }
	}
	sata = MIN(MAXHISTVAL,MAX(MINSATURATE,0.9*((float)isat))/fac);
	lpix = ihist[isat - MINHISTVAL];
	while (lpix < mpix/1000 && isat > MINHISTVAL) {
	    isat--;
	    lpix += ihist[isat - MINHISTVAL];
	}
	*sat = ((float)isat)/fac + skyref;
	*sat = MIN(MAXHISTVAL,MAX(MINSATURATE,MAX(0.95*(*sat),sata)));

	/* If all you want is the saturation level, then get out of here...*/

	if (satonly) {
	    free(ihist);
	    return(ERRCODE_OK);
	}

	/* Now find the median and sigma */

	imcore_medsig(ihist,MAXHIST,MINHISTVAL-1,mpix,skymed,skysig);
	sigsq = (*skysig)*(*skysig);

	/* Do an iterative 3-sigma upper clip to give a more robust 
	   estimator */

	iclip = MAXHISTVAL;
	mcpix = mpix;
	skymedc = *skymed;
	sigmac = *skysig;
	for (iloop = 0; iloop < 3; iloop++) {
	    irej = 0;
	    for (i = NINT(skymedc+3.0*sigmac); i <= iclip; i++)
		irej += ihist[i - MINHISTVAL];
	    if (irej == 0)
		break;
	    iclip = NINT(skymedc+3.0*sigmac)-1;
	    mcpix = mcpix-irej;
	    imcore_medsig(ihist,MAXHIST,MINHISTVAL-1,mcpix,&skymedc,&sigmac);
	}
	if (sigmac > 2.5)
	    break;
    }

    /* Set the final answer */

    *skymed = skymedc/fac + skyref;
    *skysig = sigmac/fac;
    free(ihist);
    return(ERRCODE_OK);
}   

void imcore_backest(ap_t *ap, float x, float y, float *skylev, float *skyrms) {
    int i,j,nbx,nby,nbsize,nbsizo2,iby,ibyp1,ibx,ibxp1;
    float **bvals,fnbsize,dely,delx,t1,t2;

    /* Define some local variables */

    nbx = ap->backmap.nbx;
    nby = ap->backmap.nby;
    nbsize = ap->backmap.nbsize;
    bvals = ap->backmap.bvals;

    /* If background tracking hasn't been done, then just give the global 
       result */

    if (nbsize <= 0) {
        *skylev = ap->background;
        *skyrms = ap->sigma;
        return;
    }

    /* Get closest pixel to the input location */

    i = NINT(x);
    j = NINT(y);

    /* Now, work out where in the map to do the interpolation */

    nbsizo2 = nbsize/2;
    fnbsize = 1.0/((float)nbsize);
    iby = (j + nbsizo2)/nbsize;
    ibyp1 = iby + 1;
    iby = MIN(nby,MAX(1,iby));
    ibyp1 = MIN(nby,ibyp1);
    dely = (j  - nbsize*iby + nbsizo2)*fnbsize;
    ibx = (i + nbsizo2)/nbsize;
    ibxp1 = ibx + 1;
    ibx = MIN(nbx,MAX(1,ibx));
    ibxp1 = MIN(nbx,ibxp1);
    delx = (i - nbsize*ibx + nbsizo2)*fnbsize;

    /* Now do a linear interpolation to find the background. Calculate MAD of
       the four adjacent background cells as an estimate of the RMS */

    t1 = (1.0 - dely)*bvals[iby-1][ibx-1] + dely*bvals[ibyp1-1][ibx-1];
    t2 = (1.0 - dely)*bvals[iby-1][ibxp1-1] + dely*bvals[ibyp1-1][ibxp1-1];
    *skylev = (1.0 - delx)*t1 + delx*t2;
    *skyrms = 0.25*(fabsf(bvals[iby-1][ibx-1] - *skylev) +
		    fabsf(bvals[ibyp1-1][ibx-1] - *skylev) +
		    fabsf(bvals[iby-1][ibxp1-1] - *skylev) +
		    fabsf(bvals[ibyp1-1][ibxp1-1] - *skylev));
}

void imcore_medsig(int *hist, int nh, int ist, int itarg, float *med, 
		   float *sig) {
  int isum, medata;
  float ffrac,sigmed;
 
  /* median */ 

  isum = 0;
  medata = ist;
  while (isum <= (itarg+1)/2 && (medata-MINHISTVAL) < MAXHIST) {
    medata++;
    isum += hist[medata-MINHISTVAL];
  }
  if (hist[medata-MINHISTVAL] == 0) {
      ffrac = 0.0;
  } else {
      ffrac = (float)(isum - (itarg+1)/2)/(float)hist[medata-MINHISTVAL];
  }
  *med = (float)medata - ffrac + 0.5;
 
  /* sigma */

  isum = 0;
  medata = ist;
  while (isum <= (itarg+3)/4 && (medata-MINHISTVAL) < MAXHIST) {
    medata++;
    isum += hist[medata-MINHISTVAL];
  }
  if (hist[medata-MINHISTVAL] == 0) {
      ffrac = 0.0;
  } else {
      ffrac = (float)(isum - (itarg+3)/4)/(float)hist[medata-MINHISTVAL];
  }
  sigmed = (float)medata - ffrac + 0.5;
  *sig = 1.48*(*med - sigmed);
  *sig = MAX(0.5,*sig);
}

static void sortit (float ia[], int n) {
    int i, j, ii, jj, ifin;
    float it;
 
    jj = 4;
    while (jj < n) 
	jj = 2 * jj;
    jj = MIN(n,(3 * jj)/4 - 1);
    while (jj > 1) {
        jj = jj/2;
        ifin = n - jj;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + jj;
            if (ia[i] <= ia[j]) 
		continue;
            it = ia[j];
            do {
                ia[j] = ia[i];
                j = i;
                i = i - jj;
                if (i < 0) 
		    break;
            } while (ia[i] > it);
            ia[j] = it;
        }
    }
    return;
}


static void tidy() {
  int i;
 
  freespace(nnp);
  if (hist != NULL) {
      for (i = 0; i < npvx; i++) 
          freespace(hist[i]);
  }
  freespace(hist);
  return;
}

/*

$Log: imcore_background.c,v $
Revision 1.3  2014/07/31 12:45:16  jim
Modified so that nbsize <= 0 gives a constant background

Revision 1.2  2010/12/10 12:26:49  jim
Fixed bug in backstats to cater for images with high background and low
noise

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.18  2009/01/22 13:47:38  jim
Changed pixel flag conditions when doing stats

Revision 1.17  2008/06/20 11:09:34  jim
Added in missing error message

Revision 1.16  2008/04/15 19:01:02  jim
Changed how the routines look for flagged pixels

Revision 1.15  2007/12/18 15:22:53  jim
Modified the background level iteration parameters

Revision 1.14  2007/07/31 12:04:33  jim
allocates bval now only if the pointer is NULL

Revision 1.13  2007/05/10 09:53:55  jim
Modified backstats to check there are some non-zero values

Revision 1.12  2006/07/28 19:26:15  jim
Modified clause of while loops

Revision 1.11  2006/07/27 12:54:47  jim
Fixed problem where null value was being recorded as 0.0 rather than
-1000.0

Revision 1.10  2006/07/24 11:41:25  jim
Fixed some problems with the background estimation

Revision 1.5  2006/03/08 11:31:49  jim
Modified minimum background sigma

Revision 1.4  2005/11/28 14:22:43  jim
Fixed possible bug in tidy

Revision 1.3  2005/03/07 01:14:53  jim
Fixed rejection loop in background routine to quit if no pixels are rejected

Revision 1.2  2004/09/07 14:18:57  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:54:58  jim
New version for rewrite of imcore


*/
