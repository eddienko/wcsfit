#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tools.h>
#include <nebuliser.h>

static float sum;
static float sumw;
static int   naver;
static float nextw;
static float lastw;
static float nextval;
static float lastval;
static short int nextc;
static short int lastc;

static void docols_2(float *, unsigned char *, float *, unsigned char *,
		     int, int, int, int);
static void dorows_2(float *, unsigned char *, float *, unsigned char *,
		     int, int, int, int);
static void wraparound(float *, unsigned char *, int, int, int, float **,
		       unsigned char **, int *);
static void medavg(float *, unsigned char *, int *, int, int, int, float *,
		   unsigned char *);
static void quickie(float *, unsigned char *, int *, int, int);
static void sortm(float *, unsigned char *, int *, int);
static void plugholes(float *, unsigned char *, int);

extern void twodfilt(float *data, unsigned char *bpm, int nx, int ny, 
		     int medfilt, int linfilt, int niter, int axis, 
		     int twod, int takeout_sky, int inorm, int wantback, 
		     float signeg, float sigpos, float **backmap) {
    int i,nn,iter,nter,nmask;
    float *buffer,*orig,*orig_sm,*work,medsky,sigsky,rescale,lthr,hthr;
    float diff;
    unsigned char *bpmcopy;
    char msg[BUFSIZ];
    
    /* Get some workspace. One holds a copy of the original data. The
       others are for work */

    nn = nx*ny;
    buffer = cir_malloc(3*nn*sizeof(*buffer));
    orig = buffer;
    orig_sm = orig + nn;
    work = orig_sm + nn;
    memmove((char *)orig,(char *)data,nn*sizeof(*data));
    memmove((char *)orig_sm,(char *)data,nn*sizeof(*data));
    memmove((char *)work,(char *)data,nn*sizeof(*data));

    /* Copy the bad pixel mask, so that the pre-existing bad pixels are 
       now flagged with 1. */

    bpmcopy = cir_calloc(nn,sizeof(*bpmcopy));
    for (i = 0; i < nn; i++) 
	bpmcopy[i] = (bpm[i] ? 1 : 0);

    /* Do gentle smooth on the original data */

    if (niter > 1 && medfilt > 10) {
	cir_bfilt2(orig_sm,bpmcopy,nx,ny,5,MEDIANCALC,axis);
	cir_bfilt2(orig_sm,bpmcopy,nx,ny,3,MEANCALC,axis);
    }

    /* Now do an iteration loop */

    for (iter = 1; iter <= niter; iter++) {
	if (iter > 1)
            memmove((char *)data,(char *)orig,nn*sizeof(*data));

	/* Filter the original input data, using the latest interation
	   on the pixel mask */

	if (! twod) 
	    cir_bfilt2(data,bpmcopy,nx,ny,medfilt,MEDIANCALC,axis);
	else
	    cir_bfilt_2d(data,bpmcopy,nx,ny,medfilt,MEDIANCALC);
	if (iter == niter) 
	    break;

        /* Look at the difference between the smoothed map and the (possibly
	   gently smoothed) original data */

	for (i = 0; i < nn; i++) 
	    work[i] = orig_sm[i] - data[i];

	/* What is the median level and RMS of the residual map? We may need
	   to iterate on this */

        (void)cir_qmedsig(work,bpmcopy,nn,3.0,3,-1000.0,65535.0,&medsky,&sigsky,
			  msg);
	rescale = 2.0;
	nter = 0;
	while (sigsky < 2.5 && nter < 16) {
	    nter++;
	    for (i = 0; i < nn; i++)
		work[i] *= rescale;
            (void)cir_qmedsig(work,bpmcopy,nn,3.0,3,-1000.0,65535.0,&medsky,
			      &sigsky,msg);
	}
	if (nter > 0) {
	    rescale = (float)pow(2.0,(double)nter);
	    for (i = 0; i < nn; i++)
		work[i] /= rescale;
	    medsky /= rescale;
	    sigsky /= rescale;
	}
	lthr = -signeg*sigsky;
	hthr = sigpos*sigsky;

	/* Clip out discordant points */

	nmask = 0;
	for (i = 0; i < nn; i++) {
	    if (bpmcopy[i] == 1)
		continue;
	    diff = work[i] - medsky;
	    if (diff > hthr || diff < lthr) {
		bpmcopy[i] = 2;
		nmask++;
	    } else {
		bpmcopy[i] = 0;
	    }
	}
    }
      
    /* Now do the linear filter */

    if (! twod) 
	cir_bfilt2(data,bpm,nx,ny,linfilt,MEANCALC,axis);
    else
	cir_bfilt_2d(data,bpm,nx,ny,linfilt,MEANCALC);
    
    /* Get the sky level if you want to keep it */

    if (! takeout_sky)
        (void)cir_qmedsig(orig,bpmcopy,nn,3.0,3,-1000.0,65535.0,&medsky,
			  &sigsky,msg);
    else 
	medsky = 0.0;

    /* Do we want a background map */

    if (wantback) {
	*backmap = cir_malloc(nn*sizeof(**backmap));
	for (i = 0; i < nn; i++) 
	    (*backmap)[i] = data[i];
    } else {
	*backmap = NULL;
    }

    /* How do we want to normalise? */

    if (inorm == 0) {
	for (i = 0; i < nn; i++) 
	    data[i] = orig[i] - data[i] + medsky;
    } else {
	for (i = 0; i < nn; i++) 
	    data[i] = orig[i]/max(1.0,data[i]);
    }

    /* Tidy and exit */

    freespace(buffer);
    freespace(bpmcopy);
}

extern void cir_bfilt2(float *data, unsigned char *bpm, int nx, int ny, 
		       int filt, int stat, int axis) {

    float *dbuf;
    unsigned char *bbuf;
    int nbuf;

    /* Get some workspace */

    nbuf = max(nx,ny);
    dbuf = cir_malloc(nbuf*sizeof(*dbuf));
    bbuf = cir_malloc(nbuf*sizeof(*bbuf));

    /* Order the reset correction so that the first smoothing is done
       across the axis of the anomaly */

    if (axis == 1) {
        dorows_2(data,bpm,dbuf,bbuf,nx,ny,filt,stat);
        docols_2(data,bpm,dbuf,bbuf,nx,ny,filt,stat);
    } else {
        docols_2(data,bpm,dbuf,bbuf,nx,ny,filt,stat);
        dorows_2(data,bpm,dbuf,bbuf,nx,ny,filt,stat);
    }

    /* Ditch workspace */

    free(dbuf);
    free(bbuf);
}

extern void cir_bfilt_2d(float *data, unsigned char *bpmcopy, int nx, int ny,
			 int filt, int stat) {
    float *dbuf,*outmap,value,*om;
    unsigned char *outbpm,*ob;
    int nbuf,j,i,nf2,nalloc,jj,ii,ind,ind1,j1old,j2old,i1old,i2old;
    char msg[BUFSIZ];

    /* Filter halfwidth */

    nf2 = filt/2;

    /* Get some workspace */

    nalloc = (2*filt+1)*(2*filt+1);
    dbuf = cir_malloc(nalloc*sizeof(*dbuf));
    outmap = cir_malloc(nx*ny*sizeof(*outmap));
    outbpm = cir_malloc(nx*ny*sizeof(*outbpm));

    /* Loop for each input pixel */

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
	    nbuf = 0;
	    for (jj = j - nf2; jj <= j + nf2; jj++) {
		if (jj < 0 || jj >= ny)
		    continue;
		ind1 = jj*nx;
		for (ii = i - nf2; ii <= i + nf2; ii++) {
		    if (ii < 0 || ii >= nx) 
			continue;
		    ind = ind1 + ii;
		    if (bpmcopy[ind])
			continue;
		    dbuf[nbuf++] = data[ind];
		}
	    }
	    
	    /* If we don't have enough, try and increase the window size.
	       This will only affect the edges */

	    if (nbuf < filt/4) {
		j1old = j - nf2;
		j2old = j + nf2;
		i1old = i - nf2;
		i2old = i + nf2;
		for (jj = j - filt; jj <= j + filt; jj++) {
		    if (jj < 0 || jj >= ny || (jj >= j1old && jj <= j2old))
			continue;
		    ind1 = jj*nx;
		    for (ii = i - filt; ii <= i + filt; ii++) {
			if (ii < 0 || ii >= nx || (ii >= i1old && ii <= i2old)) 
			    continue;
			ind = ind1 + ii;
			if (bpmcopy[ind])
			    continue;
			dbuf[nbuf++] = data[ind];
		    }
		}
	    }

	    /* Right, assuming we have enough entries, then get a median */

	    ind = j*nx + i;
	    if (nbuf > filt/4) {
  	        if (stat == MEDIANCALC) 
		    (void)cir_med(dbuf,NULL,nbuf,&value,msg);
	        else
		    (void)cir_mean(dbuf,NULL,nbuf,&value,msg);
		outmap[ind] = value;
		outbpm[ind] = 0;
	    } else {
		outmap[ind] = -1000.0;
		outbpm[ind] = 1;
	    }
	}
    }

    /* Right, fill in the holes and then transfer the filtered data to
       the input array */

    for (j = 0; j < ny; j++) {
	om = outmap + j*nx;
	ob = outbpm + j*nx;
	plugholes(om,ob,nx);
	for (i = 0; i < nx; i++) 
	    data[j*nx+i] = om[i];
    }

    /* Tidy up and get out of here */

    freespace(outmap);
    freespace(outbpm);
    freespace(dbuf);
    return;
}

static void docols_2(float *data, unsigned char *bpm, float *dbuf, 
		     unsigned char *bbuf, int nx, int ny, int filter, 
		     int stat) {

    int j,k,indx,nn;
    unsigned char *goodval,*b;
    float *t;

    if (filter <= 0)
        return;

    goodval = cir_malloc(ny*sizeof(*goodval));
    t = cir_malloc(ny*sizeof(*t));
    b = cir_malloc(ny*sizeof(*b));
    for (k = 0; k < nx; k++) {
        memset((char *)goodval,0,ny);
	nn = 0;
	for (j = 0; j < ny; j++) {
	    indx = j*nx + k;
	    if (bpm[indx] == 0) {
		dbuf[nn] = data[indx];
		bbuf[nn++] = 0;
	    }
	} 
	dostat(dbuf,bbuf,goodval,nn,filter,stat);
	nn = 0;
	for (j = 0; j < ny; j++) {
	    indx = j*nx + k;
	    if (bpm[indx] == 0) {
		t[j] = dbuf[nn++];
		b[j] = 0;
	    } else {
		t[j] = -999.0;
		b[j] = 1;
	    }
	}
	plugholes(t,b,ny);
	nn = 0;
	for (j = 0; j < ny; j++) {
	    indx = j*nx + k;
	    data[indx] = t[j];
	}
    } 
    free(goodval);
    free(t);
    free(b);
}

static void dorows_2(float *data, unsigned char *bpm, float *dbuf, 
		     unsigned char *bbuf, int nx, int ny, int filter, 
		     int stat) {

    int j,k,indx,nn;
    unsigned char *goodval,*b;
    float *t;

    if (filter <= 0)
        return;

    goodval = cir_malloc(nx*sizeof(*goodval));
    t = cir_malloc(nx*sizeof(*t));
    b = cir_malloc(nx*sizeof(*b));
    for (k = 0; k < ny; k++) {
        memset((char *)goodval,0,nx);
	nn = 0;
	for (j = 0; j < nx; j++) {
	    indx = k*nx + j;
	    if (bpm[indx])
		continue;
	    dbuf[nn] = data[indx];
	    bbuf[nn++] = 0;
	}
	dostat(dbuf,bbuf,goodval,nn,filter,stat);
	nn = 0;
	for (j = 0; j < nx; j++) {
	    indx = k*nx + j;
	    if (bpm[indx] == 0) {
		t[j] = dbuf[nn++];
		b[j] = 0;
	    } else {
		t[j] = -999.0;
		b[j] = 1;
	    }
	}
	plugholes(t,b,nx);
        for (j = 0; j < nx; j++) {
	    indx = k*nx + j;
	    data[indx] = t[j];
	}
    }
    free(goodval);
    free(t);
    free(b);
}

extern void dostat(float *data, unsigned char *bpm, unsigned char *goodval,
		   int npts, int nfilt, int whichstat) {
    int nbuf,jl,jh,j,*ipoint,ifree,i;
    unsigned char *ybbuf,*barray,bval;
    float *ybuf,*darray,val;
    
    /* Check to make sure there are some points in the array... */

    if (npts < nfilt || npts < 10)
	return;

    /* Check to make sure the filter size is odd */

    if ((nfilt/2)*2 == nfilt)
        nfilt++;

    /* Do the wrap around and load the data into an oversized array*/

    wraparound(data,bpm,npts,nfilt,whichstat,&ybuf,&ybbuf,&nbuf);

    /* Start doing the filtering...*/

    darray = cir_malloc(nfilt*sizeof(*darray));
    barray = cir_malloc(nfilt*sizeof(*barray));
    ipoint = cir_malloc(nfilt*sizeof(*ipoint));
    memmove((char *)darray,(char *)ybuf,nfilt*sizeof(*ybuf));
    memmove((char *)barray,(char *)ybbuf,nfilt*sizeof(*ybbuf));
    for (j = 0; j < nfilt; j++)
        ipoint[j] = j;
    ifree = 0;
    medavg(darray,barray,ipoint,nfilt,whichstat,-1,&val,&bval);
    if (! bval) 
        data[0] = val;
    goodval[0] = bval;
    jl = nfilt;
    jh = nfilt + npts - 2;
    for (j = jl; j <= jh; j++) {
        for (i = 0; i < nfilt; i++) {
            if (ipoint[i] == 0) {
                ifree = i;
                ipoint[i] = nfilt - 1;
                lastval = darray[ifree];
                lastw = 0.0;
                lastc = 0;
                if (barray[ifree] == 0) {
                    lastw = 1.0;
                    lastc = 1;
                }                
                darray[ifree] = ybuf[j];
                barray[ifree] = ybbuf[j];
                nextval = darray[ifree];
                nextw = 0.0;
                nextc = 0;
                if (barray[ifree] == 0) {
                    nextw = 1.0;
                    nextc = 1;
                }                
            } else
                ipoint[i]--;
        }
        medavg(darray,barray,ipoint,nfilt,whichstat,ifree,&val,&bval);
        if (! bval) 
            data[j-jl+1] = val;
        goodval[j-jl+1] = bval;
    }

    /* Ditch workspace */

    free(darray);
    free(barray);
    free(ipoint);
    free(ybuf);
    free(ybbuf);
}

static void wraparound(float *data, unsigned char *bpm, int npts, int nfilt, 
		       int whichstat, float **ybuf, unsigned char **ybbuf, 
		       int *nbuf) {

    float *darray,xmns,xmnf;
    int i1,ilow,i,*ipoint;
    unsigned char *barray,bxmns,bxmnf;

    /* Do some padding at the edges */

    i1 = nfilt/2;
    ilow = max(3,nfilt/4);
    ilow = (ilow/2)*2 + 1;

    /* Get some workspace */

    darray = cir_malloc(nfilt*sizeof(*darray));
    barray = cir_malloc(nfilt*sizeof(*barray));
    ipoint = cir_calloc(nfilt,sizeof(*ipoint));
    *nbuf = npts + 2*i1;
    *ybuf = cir_malloc(*nbuf*sizeof(float));
    *ybbuf = cir_malloc(*nbuf*sizeof(unsigned char));

    /* Do the wrap around.*/

    memmove((char *)darray,(char *)data,ilow*sizeof(*data));
    memmove((char *)barray,(char *)bpm,ilow*sizeof(*bpm));
    medavg(darray,barray,ipoint,ilow,whichstat,-1,&xmns,&bxmns);
    memmove((char *)darray,(char *)(data+npts-ilow),ilow*sizeof(*data));
    memmove((char *)barray,(char *)(bpm+npts-ilow),ilow*sizeof(*bpm));
    medavg(darray,barray,ipoint,ilow,whichstat,-1,&xmnf,&bxmnf);
    for (i=0; i < i1; i++) {
        if (! bxmns) {
            (*ybuf)[i] = 2.0*xmns - data[i1+ilow-i-1];
            (*ybbuf)[i] = bpm[i1+ilow-i-1];
        } else {
            (*ybuf)[i] = data[i1+ilow-i-1];
            (*ybbuf)[i] = 1;
        }
        if (! bxmnf) {
            (*ybuf)[npts+i1+i] = 2.0*xmnf - data[npts-i-ilow-1];
            (*ybbuf)[npts+i1+i] = bpm[npts-i-ilow-1];
        } else {
            (*ybuf)[npts+i1+i] = data[npts-i-ilow-1];
            (*ybbuf)[npts+i1+i] = 1;
        }
    }

    /* Now place the full line into the buffer */

    memmove((char *)(*ybuf+i1),data,npts*sizeof(*data));
    memmove((char *)(*ybbuf+i1),bpm,npts*sizeof(*bpm));

    /* Free workspace */

    free(darray);
    free(barray);
    free(ipoint);
}

static void medavg(float *array, unsigned char *bpm, int *ipoint, int npix, 
                   int whichstat, int newl, float *outval, 
                   unsigned char *outbp) { 

    float *buf = NULL;
    int m,i;

    /* If there is a new element and that new element is bad, then
       the buffer is already sorted from last time (just one element
       sorter */
    
    m = 0;
    if (whichstat == MEDIANCALC) {
        if (newl == -1) 
            sortm(array,bpm,ipoint,npix);
        else  
            quickie(array,bpm,ipoint,newl,npix);

        /* Get some workspace */

        buf = cir_malloc(npix*sizeof(*buf));

        /* Now put everything that's good in the buffer */

        m = 0;
        for (i = 0; i < npix; i++) {
            if (bpm[i] == 0) {
                buf[m] = array[i];
                m++;
            }
        }
    } else if (whichstat == MEANCALC) {
        if (newl == -1) {
            sum = 0.0;
            sumw = 0.0;
            naver = 0;
            for (i = 0; i < npix; i++) {
                if (bpm[i] == 0) {
                    sum += array[i];
                    sumw += 1.0;
                    naver += 1;
                }
            }
            m = naver;
        } else {
            sum += (nextw*nextval - lastw*lastval);
            sumw += (nextw - lastw);
            naver += (nextc - lastc);
            m = naver;
        }
    }
        
    /* If they were all bad, then send a null result back */

    if (m == 0) {
        *outval = 0.0;
        *outbp = 1;
        if (whichstat == MEDIANCALC)
            free(buf);

    /* Otherwise calculate the relevant stat */

    } else {
        if (whichstat == MEDIANCALC) {
            *outval = buf[m/2];
            free(buf);
        } else if (whichstat == MEANCALC)
            *outval = sum/sumw;
        *outbp = 0;
    }
}

static void quickie(float *array, unsigned char *iarray, int *iarray2, 
		    int lll, int narray) {

    float test;
    int i,j,npt,it2;
    unsigned char it;

    test = array[lll];
    it = iarray[lll];
    it2 = iarray2[lll];
    j = -1;
    for (i = 0; i < narray; i++) {
        if (i != lll && test <= array[i]) {
            j = i;
            break;
        }
    }
    if (j == -1) 
        j = narray;
    if (j - 1 == lll)
        return;

    if (j - lll < 0) {
        npt = lll - j;
        for (i = 0; i < npt; i++) {
            array[lll-i] = array[lll-i-1];
            iarray[lll-i] = iarray[lll-i-1];
            iarray2[lll-i] = iarray2[lll-i-1];
        }
        array[j] = test;
        iarray[j] = it;
        iarray2[j] = it2;
    } else {
        j--;
        npt = j - lll;
        if (npt != 0) {
            for (i = 0; i < npt; i++) {
                array[lll+i] = array[lll+i+1];
                iarray[lll+i] = iarray[lll+i+1];
                iarray2[lll+i] = iarray2[lll+i+1];
            }
        }
        array[j] = test;
        iarray[j] = it;
        iarray2[j] = it2;
    }
}
            
static void sortm(float *a1, unsigned char *a2, int *a3, int n) {
    int iii,ii,i,ifin,j,b3;
    unsigned char b2;
    float b1;

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
            if (a1[i] > a1[j]) {
                b1 = a1[j];
                b2 = a2[j];
                b3 = a3[j];
                while (1) {
                    a1[j] = a1[i];
                    a2[j] = a2[i];
                    a3[j] = a3[i];
                    j = i;
                    i = i - iii;
                    if (i < 0 || a1[i] <= b1) 
                        break;
                }
                a1[j] = b1;
                a2[j] = b2;
                a3[j] = b3;
            }
        }
    }
}

static void plugholes(float *data, unsigned char *bpm, int nx) {
    int i,ifirst,ilast,i1,i2,j;
    float nc,d1,d2,t1,t2,slope;

    /* First of all, find the first good value in the array */

    i = 0;
    while (i < nx && bpm[i] != 0)
        i++;
    ifirst = i;

    /* If all the values in the array are bad, then do nothing */

    if (ifirst == nx)
        return;

    /* Find the last good value in the array */

    i = nx - 1;
    while (i >= 0 && bpm[i] != 0) 
        i--;
    ilast = i;

    /* Right, now start from the first good value and fill in any holes in the
       middle part of the array */

    i = ifirst;
    while (i <= ilast) {
        if (bpm[i] == 0) {
            i++;
            continue;
        }
        i1 = i - 1;
        while (bpm[i] != 0) 
            i++;
        i2 = i;
        nc = (float)(i2 - i1 + 1);
        d1 = data[i1];
        d2 = data[i2];
        for (j = i1+1; j <= i2-1; j++) {
            t1 = 1.0 - (float)(j - i1)/nc;
            t2 = 1.0 - t1;
            data[j] = t1*d1 + t2*d2;
        }
    }

    /* Now the left bit... */

    if (ifirst > 0) {
        slope = data[ifirst+1] - data[ifirst];
        for (j = 0; j < ifirst; j++)
            data[j] = slope*(float)(j - ifirst) + data[ifirst];
    }

    /* Now the right bit... */

    if (ilast < nx - 1) {
        slope = data[ilast] - data[ilast-1];
        for (j = ilast; j < nx; j++) 
            data[j] = slope*(float)(j - ilast) + data[ilast];
    }
}
