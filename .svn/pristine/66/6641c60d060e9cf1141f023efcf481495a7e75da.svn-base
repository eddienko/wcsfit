#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>
#include <stdint.h>

#include "usnob.h"
#include "common.h"
#include "cvtunit.h"
#include "floatmath.h"
#include "util.h"

/* Linux seems to call uint32_t u_int32_t... */
/* #ifdef linux */
/* typedef u_int32_t uint32_t; */
/* #endif */

#define READ_INT32(h, b, p, o) {		\
  if(h->swap) {					\
    *(o)       = (b)[(p) + 3];			\
    *((o) + 1) = (b)[(p) + 2];			\
    *((o) + 2) = (b)[(p) + 1];			\
    *((o) + 3) = (b)[(p)];			\
  }						\
  else {					\
    *(o)       = (b)[(p)];			\
    *((o) + 1) = (b)[(p) + 1];			\
    *((o) + 2) = (b)[(p) + 2];			\
    *((o) + 3) = (b)[(p) + 3];			\
  }						\
}

#define READ_INT32_ARRAY(h, b, p, o) {				\
  READ_INT32(h, b, p, (unsigned char *) &((o)[0]));		\
  READ_INT32(h, b, (p) + 4, (unsigned char *) &((o)[1]));	\
  READ_INT32(h, b, (p) + 8, (unsigned char *) &((o)[2]));	\
  READ_INT32(h, b, (p) + 12, (unsigned char *) &((o)[3]));	\
  READ_INT32(h, b, (p) + 16, (unsigned char *) &((o)[4]));	\
}

static int usnob_buffer (struct usnob_handle *h, char *errstr);
static unsigned char isbigendian (void);

int usnob_open (struct usnob_handle *h, char *path, char *errstr) {
  int len;
  char *p;

  /* Determine byte-ordering of system */
  h->swap = isbigendian();

  /* Decode filename */
  p = strstr(path, "%%%%");
  if(!p) {
    report_err(errstr, "no token '%%' found in path: %s", path);
    goto error;
  }

  len = strlen(path);
  if(len > FNBUFSIZE) {
    report_err(errstr, "pathname too long");
    goto error;
  }

  memcpy(h->fname, path, len+1);
  h->fnptr = (int) (p - path);

  /* Allocate buffers */
  h->buf = (struct usnob_row *) malloc(USNOB_BLKNREC * sizeof(struct usnob_row));
  h->cbuf = (char *) malloc(USNOB_BLKNREC * USNOB_RECSIZE * sizeof(struct usnob_row));
  if(!h->buf || !h->cbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  h->f = (struct usnob_file *) NULL;
  h->nf = 0;

  return(0);

 error:
  return(1);
}

static unsigned char isbigendian (void) {
  short one = 1;
  unsigned char *p;

  p = (unsigned char *) &one;

  return(*p == 1 ? 0 : 1);
}

int usnob_close (struct usnob_handle *h, char *errstr) {
  int i;

  free((void *) h->buf);
  h->buf = (struct usnob_row *) NULL;
  free((void *) h->cbuf);
  h->cbuf = (char *) NULL;

  if(h->f) {
    for(i = 0; i < h->nf; i++)
      if(h->f[i].fp)
	fclose(h->f[i].fp);

    free((void *) h->f);
    h->f = (struct usnob_file *) NULL;
  }

  return(0);
}

static int usnob_findfiles (struct usnob_handle *h, float dec_low, float dec_high,
			   char *errstr) {
  int32_t spd_low, spd_high, spd;
  int i, rv;
  struct usnob_file *pf;
  struct stat sb;
  char tmpnam[5];

  /* Make sure the previous lot are cleaned up */
  if(h->f) {
    for(i = 0; i < h->nf; i++)
      if(h->f[i].fp)
	fclose(h->f[i].fp);

    free((void *) h->f);
    h->f = (struct usnob_file *) NULL;
  }

  /* Calculate SPD limits */
  spd_low = (int32_t) floor((dec_low * RAD_TO_DEG + 90) * 10);
  spd_high = (int32_t) ceil((dec_high * RAD_TO_DEG + 90) * 10);

  if(spd_low < 0)
    spd_low = 0;
  if(spd_high > 1800)
    spd_high = 1800;

  /* Allocate array for files */
  h->nf = (spd_high - spd_low) + (spd_high < 1800 ? 1 : 0);

  h->f = (struct usnob_file *) malloc(h->nf * sizeof(struct usnob_file));
  if(!h->f) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  for(spd = spd_low, pf = h->f; spd <= spd_high && spd < 1800; spd++, pf++) {
    /* Form filename */
    memcpy(pf->path, h->fname, FNBUFSIZE);
    snprintf(tmpnam, sizeof(tmpnam), "%04d", (int) spd);
    memcpy(&(pf->path[h->fnptr]), tmpnam, 4);

#if 0
    fprintf(stderr, "%s\n", pf->path);
#endif

    /* Find out how big the catalogue is */
    rv = stat(pf->path, &sb);
    if(rv == -1) {
      report_syserr(errstr, "could not stat %s", pf->path);
      goto error;
    }
    
    if((sb.st_size % USNOB_RECSIZE) != 0) {
      report_err(errstr, "file %s record size not %ld bytes", pf->path, USNOB_RECSIZE);
      goto error;
    }
    
    pf->nrows = sb.st_size / USNOB_RECSIZE;

    /* Store parameters */
    pf->fp = (FILE *) NULL;
    pf->min_dec = DEG_TO_RAD * (((float) spd) / 10.0 - 90.0);
    pf->max_dec = DEG_TO_RAD * (((float) (spd + 1)) / 10.0 - 90.0);
  }

  return(0);

 error:
  if(h->f)
    free((void *) h->f);
  
  return(1);
}

long usnob_find (struct usnob_handle *h, float ra, float dec, float width, char *errstr) {
  float l_ra_low, l_ra_high, ra_low = 0.0, ra_high = 0.0, o_ra_low = 0.0, o_ra_high = 0.0;
  long start, finish, first_index, last_index, ntot = 0L;
  unsigned char i, ra_wrap = 0;

  struct usnob_file *pf;
  int rv, ifile;
  size_t srv;
  char buf[USNOB_RECSIZE];
  float rec_ra;

#define FETCH_RA(rnum) {							\
  int recno;									\
  uint32_t tmp;									\
										\
  recno = (rnum)-1;								\
										\
  if(pf->rec != recno) {								\
    rv = fseek(pf->fp, recno * USNOB_RECSIZE, SEEK_SET);				\
    if(rv == -1) {								\
      report_syserr(errstr, "seek: record %ld", recno+1);			\
      goto error;								\
    }										\
										\
    pf->rec = recno;								\
  }										\
										\
  srv = fread((void *) buf, 1, USNOB_RECSIZE, pf->fp);				\
  if(srv != USNOB_RECSIZE) {							\
    if(ferror(pf->fp)) {								\
      report_syserr(errstr, "read: record %ld", pf->rec+1);			\
      goto error;								\
    }										\
    else {									\
      report_err(errstr, "read: unexpected EOF before record %ld", pf->rec+1);	\
      goto error;								\
    }										\
  }										\
										\
  pf->rec++;									\
										\
  READ_INT32(h, buf, USNOB_BUF_RA, (unsigned char *) &tmp);			\
  rec_ra = AS_TO_RAD * ((float) tmp) / 100.0;					\
}

  /* Calculate RA, Dec limits */
  radeclimits(ra, dec, width, &l_ra_low, &l_ra_high, &(h->dec_low), &(h->dec_high));

  h->ra = ra;
  h->dec = dec;
  h->width = width;

  /* Find the appropriate files */
  if(usnob_findfiles(h, h->dec_low, h->dec_high, errstr))
    goto error;

  /* For each file... */
  for(ifile = 0, pf = h->f; ifile < h->nf; ifile++, pf++) {
    /* Open the catalogue */
    pf->fp = fopen(pf->path, "rb");
    if(!pf->fp) {
      report_syserr(errstr, "open: %s", pf->path);
      goto error;
    }
    
    pf->rec = 0;

    /* Find our starting and finishing offsets, using binary chop on the RA column. */
    ra_low = l_ra_low;
    ra_high = l_ra_high;

    if(ra_low < 0) {
      o_ra_low  = M_PI * 2 + ra_low;
      o_ra_high = M_PI * 2;
      
      ra_low  = 0;
      ra_wrap = 1;
    }
    
    if(ra_high > M_PI * 2) {
      o_ra_low  = ra_low;
      o_ra_high = M_PI * 2;
      
      ra_low  = 0;
      ra_high = ra_high - M_PI * 2;
      
      ra_wrap = 1;
    }
    
    for(i = 0; ; i++) {
      start  = 1;
      finish = pf->nrows + 1;
      first_index = pf->nrows / 2;
      
      while(finish - start >= 2) {
	/* Fetch RA of row 'first_index' */
	FETCH_RA(first_index);
	
	if(rec_ra < ra_low) {
	  start = first_index;
	  first_index = (first_index + finish) / 2;
	}
	else {
	  finish = first_index;
	  first_index  = (first_index + start) / 2;
	}
      }
      
      /* 'first_index' now gives the start RA */
      
      start  = first_index;
      finish = pf->nrows + 1;
      last_index = (finish - start) / 2;
      
      while(finish - start >= 2) {
	/* Fetch RA of row 'last_index' */
	FETCH_RA(last_index);

	if(rec_ra < ra_high) {
	  start = last_index;
	  last_index = (last_index + finish) / 2;
	}
	else {
	  finish = last_index;
	  last_index  = (last_index + start) / 2;
	}
      }
      
      /* 'last_index' now gives the finish RA */
      
      if(last_index < first_index)
	/* Trap for case where first_index at end of catalogue */
	last_index = first_index;

      pf->first_idx[i] = first_index;
      pf->last_idx[i] = last_index;
      
      ntot += last_index - first_index;
      
      if(ra_wrap == 0)
	break;
      
      ra_low  = o_ra_low;
      ra_high = o_ra_high;
      ra_wrap = 0;
    }

    pf->nidx = i + 1;

    /* Close file */
    fclose(pf->fp);
    pf->fp = (FILE *) NULL;
  }

#undef FETCH_RA

  pf = h->f;

  h->cf = 0;
  h->ridx = 0;
  h->nbuf = 0;

  h->frow = pf->first_idx[0];
  h->nrem = pf->last_idx[0] - pf->first_idx[0];

  /* Open the first file */
  pf->fp = fopen(pf->path, "rb");
  if(!pf->fp) {
    report_syserr(errstr, "open: %s", pf->path);
    goto error;
  }

  pf->rec = 0;

  /* Fill buffer */
  if(usnob_buffer(h, errstr))
    goto error;

  return(ntot);

 error:
  return(-1);
}

static int usnob_buffer (struct usnob_handle *h, char *errstr) {
  long r;
  struct usnob_file *pf;

  int rv, recno;
  size_t srv;
  char *bp;
  struct usnob_row *rp;
  uint32_t tmp_ra, tmp_spd, tmp_pm, tmp_pmerr, tmp_poserr, tmp_mag[5], tmp_magerr[5];
  uint32_t tmp_index[5];

  pf = h->f + h->cf;

  h->nbuf = (h->nrem > USNOB_BLKNREC ? USNOB_BLKNREC : h->nrem);

  /* Seek to the correct part of the file */
  recno = h->frow-1;

  if(pf->rec != recno) {
    rv = fseek(pf->fp, recno * USNOB_RECSIZE, SEEK_SET);
    if(rv == -1) {
      report_syserr(errstr, "seek: record %ld", recno+1);
      goto error;
    }

    pf->rec = recno;
  }

  /* Read 'nbuf' records into the buffer */
  srv = fread((void *) h->cbuf, 1, h->nbuf * USNOB_RECSIZE, pf->fp);
  if(srv != h->nbuf * USNOB_RECSIZE) {
    if(ferror(pf->fp)) {
      report_syserr(errstr, "read: record %ld", pf->rec+1);
      goto error;
    }
    else {
      report_err(errstr, "read: unexpected EOF before record %ld", pf->rec+1);
      goto error;
    }
  }

  pf->rec += h->nbuf;

  for(r = 0, bp = h->cbuf, rp = h->buf; r < h->nbuf; r++, bp += USNOB_RECSIZE, rp++) {
    READ_INT32(h, bp, USNOB_BUF_RA, (unsigned char *) &tmp_ra);
    READ_INT32(h, bp, USNOB_BUF_SPD, (unsigned char *) &tmp_spd);
    READ_INT32(h, bp, USNOB_BUF_PM, (unsigned char *) &tmp_pm);
    READ_INT32(h, bp, USNOB_BUF_PMERR, (unsigned char *) &tmp_pmerr);
    READ_INT32(h, bp, USNOB_BUF_POSERR, (unsigned char *) &tmp_poserr);
    READ_INT32_ARRAY(h, bp, USNOB_BUF_MAG, tmp_mag);
    READ_INT32_ARRAY(h, bp, USNOB_BUF_MAGERR, tmp_magerr);
    READ_INT32_ARRAY(h, bp, USNOB_BUF_INDEX, tmp_index);

    rp->flags = 0;

    rp->ra  = AS_TO_RAD * ((float) tmp_ra) / 100.0;
    rp->dec = AS_TO_RAD * ((float) tmp_spd) / 100.0 - M_PI_2;

    rp->pmra = (tmp_pm % 10000) * 0.002 - 10;
    tmp_pm /= 10000;
    rp->pmdec = (tmp_pm % 10000) * 0.002 - 10;
    tmp_pm /= 10000;
    rp->pmprob = (tmp_pm % 10) * 0.1;
    tmp_pm /= 10;
    if(tmp_pm % 10)
      rp->flags |= USNOB_FLAG_MOTION;
    
    rp->pmerrra = (tmp_pmerr % 1000) * 0.001;
    tmp_pmerr /= 1000;
    rp->pmerrdec = (tmp_pmerr % 1000) * 0.001;
    tmp_pmerr /= 1000;
    rp->sigra = (tmp_pmerr % 10) * 0.1;
    tmp_pmerr /= 10;
    rp->sigdec = (tmp_pmerr % 10) * 0.1;
    tmp_pmerr /= 10;
    rp->ndet = (tmp_pmerr % 10);
    tmp_pmerr /= 10;
    if(tmp_pmerr % 10)
      rp->flags |= USNOB_FLAG_SPIKE;

    rp->raerr = (tmp_poserr % 1000) * 0.001;
    tmp_poserr /= 1000;
    rp->decerr = (tmp_poserr % 1000) * 0.001;
    tmp_poserr /= 1000;
    rp->epoch = (tmp_poserr % 1000) * 0.1 + 1950.0;
    tmp_poserr /= 1000;
    if(tmp_poserr % 10)
      rp->flags |= USNOB_FLAG_YSCORR;

#define GET_MAG(t, o) {				\
  o.m = ((t) % 10000) * 0.01;			\
  (t) /= 10000;					\
  o.field = ((t) % 1000);			\
  (t) /= 1000;					\
  o.survey = ((t) % 10);			\
  (t) /= 10;					\
  o.star = ((t) % 100);				\
  (t) /= 100;					\
};

    GET_MAG(tmp_mag[0], rp->b1);
    GET_MAG(tmp_mag[1], rp->r1);
    GET_MAG(tmp_mag[2], rp->b2);
    GET_MAG(tmp_mag[3], rp->r2);
    GET_MAG(tmp_mag[4], rp->n);

#undef GET_MAG

#define GET_MAGERR(t, o) {			\
  o.xires = ((t) % 10000) * 0.01;		\
  (t) /= 10000;					\
  o.xnres = ((t) % 10000) * 0.01;		\
  (t) /= 10000;					\
  o.calsrc = ((t) % 10);			\
  (t) /= 10;					\
};

    GET_MAGERR(tmp_magerr[0], rp->b1);
    GET_MAGERR(tmp_magerr[1], rp->r1);
    GET_MAGERR(tmp_magerr[2], rp->b2);
    GET_MAGERR(tmp_magerr[3], rp->r2);
    GET_MAGERR(tmp_magerr[4], rp->n);

#undef GET_MAGERR

    rp->b1.index = tmp_index[0];
    rp->r1.index = tmp_index[1];
    rp->b2.index = tmp_index[2];
    rp->r2.index = tmp_index[3];
    rp->n.index = tmp_index[4];

    rp->ptr = recno + r + 1;
  }

  h->frow += h->nbuf;
  h->nrem -= h->nbuf;
  h->brow = 0;  /* NOTE: this counts from zero! */

  return(0);

 error:

  return(1);
}

int usnob_read (struct usnob_handle *h, struct usnob_row *row, char *errstr) {
  float rra, rdec, xi, xn;
  struct usnob_file *pf;

  if(h->cf >= h->nf)
    return(0);

  pf = h->f + h->cf;

 fetch_row:
  if(h->brow >= h->nbuf) {
    /* Ran out of rows in buffer, need to fetch more */
    if(h->nrem < 1) {
      /* Done the current segment, find the next */
      h->ridx++;

      if(h->ridx >= pf->nidx) {
	/* Finished that file, get the next */
	fclose(pf->fp);
	pf->fp = (FILE *) NULL;

	h->cf++;

	if(h->cf >= h->nf)
	  return(0);
      
	pf++;
	h->nbuf = 0;

	/* Open it */
	pf->fp = fopen(pf->path, "rb");
	if(!pf->fp) {
	  report_syserr(errstr, "open: %s", pf->path);
	  goto error;
	}
	
	pf->rec = 0;
	h->ridx = 0;
      }

      h->frow = pf->first_idx[h->ridx];
      h->nrem = pf->last_idx[h->ridx] - pf->first_idx[h->ridx];
    }

    if(usnob_buffer(h, errstr))
      goto error;
  }

  rra  = (h->buf + h->brow)->ra;
  rdec = (h->buf + h->brow)->dec;

  /* Reject rows not in the correct Dec range */
  if(!(rdec > h->dec_low && rdec < h->dec_high)) {
    h->brow++;
    goto fetch_row;
  }

  /* Convert to standard coordinates and check if in range */
  STANDC(h->ra, h->dec, rra, rdec, xi, xn);
  if(!(xi > -h->width && xi < h->width && xn > -h->width && xn < h->width)) {
    h->brow++;
    goto fetch_row;
  }

  /* Found one */
  memcpy((void *) row, (void *) (h->buf + h->brow), sizeof(struct usnob_row));
  row->xi = xi;
  row->xn = xn;
  h->brow++;

  return(1);

 error:
  return(-1);
}
