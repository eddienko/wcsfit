#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>

#include "psc.h"
#include "common.h"
#include "cvtunit.h"
#include "floatmath.h"
#include "util.h"

#define DEC_COLUMN 1

static int psc_findfiles (struct psc_handle *h, float ra_low, float ra_high, char *errstr);
static int psc_buffer (struct psc_handle *h, char *errstr);

static struct {
  char *name;
  float ra_low;
  float ra_high;
} psc_filelist[] = {
  { "0",    0.000,  1.000 },
  { "1",    1.000,  2.000 },
  { "2",    2.000,  3.000 },
  { "3",    3.000,  4.000 },
  { "4",    4.000,  5.000 },
  { "5a",   5.000,  5.500 },
  { "5b",   5.500,  6.000 },
  { "6a",   6.000,  6.250 },
  { "6b",   6.250,  6.500 },
  { "6c",   6.500,  6.750 },
  { "6d",   6.750,  7.000 },
  { "7a",   7.000,  7.250 },
  { "7b",   7.250,  7.500 },
  { "7c",   7.500,  7.750 },
  { "7d",   7.750,  8.000 },
  { "8a",   8.000,  8.500 },
  { "8b",   8.500,  9.000 },
  { "9",    9.000, 10.000 },
  { "10",  10.000, 11.000 },
  { "11",  11.000, 12.000 },
  { "12",  12.000, 13.000 },
  { "13",  13.000, 14.000 },
  { "14",  14.000, 15.000 },
  { "15",  15.000, 16.000 },
  { "16a", 16.000, 16.500 },
  { "16b", 16.500, 17.000 },
  { "17a", 17.000, 17.125 },
  { "17b", 17.125, 17.250 },
  { "17c", 17.250, 17.375 },
  { "17d", 17.375, 17.500 },
  { "17e", 17.500, 17.625 },
  { "17f", 17.625, 17.750 },
  { "17g", 17.750, 17.875 },
  { "17h", 17.875, 18.000 },
  { "18a", 18.000, 18.250 },
  { "18b", 18.250, 18.500 },
  { "18c", 18.500, 18.750 },
  { "18d", 18.750, 19.000 },
  { "19a", 19.000, 19.250 },
  { "19b", 19.250, 19.500 },
  { "19c", 19.500, 19.750 },
  { "19d", 19.750, 20.000 },
  { "20a", 20.000, 20.250 },
  { "20b", 20.250, 20.500 },
  { "20c", 20.500, 20.750 },
  { "20d", 20.750, 21.000 },
  { "21",  21.000, 22.000 },
  { "22",  22.000, 23.000 },
  { "23",  23.000, 24.000 }
};

int psc_open (struct psc_handle *h, char *path, char *errstr) {
  fitsfile *fits;
  int status = 0, len1, len2;

  char fname[FLEN_FILENAME], *p;
  char *colnames[PSC_NUM_COLUMNS] = 
  { "RA", "Dec", "maj", "min", "pa",
    "j_m", "j_msig", "j_msigcom", "j_qual", "j_rd", "j_bl", "j_cc",
    "h_m", "h_msig", "h_msigcom", "h_qual", "h_rd", "h_bl", "h_cc",
    "k_m", "k_msig", "k_msigcom", "k_qual", "k_rd", "k_bl", "k_cc",
    "ptskey" };
  int c;

  /* Form name of first catalogue to extract parameters */
  p = strstr(path, "%%");
  if(!p) {
    report_err(errstr, "no token '%%' found in path: %s", path);
    goto error;
  }

  len1 = p - path;
  len2 = strlen(p + 2);
  if(len1 + len2 + 4 > FLEN_FILENAME) {
    report_err(errstr, "pathname too long");
    goto error;
  }

  memcpy(fname, path, len1);
  fname[len1] = '0';
  memcpy(fname + len1 + 1, p + 2, len2);
  fname[len1 + len2 + 1] = '\0';

  /* Keep a copy of the bits for later */
  memcpy(h->fpref, path, len1);
  h->fpref[len1] = '\0';
  h->preflen = len1;

  memcpy(h->fpost, p + 2, len2);
  h->fpost[len2] = '\0';

  /* Open the first catalogue and move to the first extension */
  ffopen(&fits, fname, READONLY, &status);
  if(status) {
    fitsio_err(errstr, status, "ffopen: %s", path);
    goto error;
  }

  ffmahd(fits, 2, (int *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "could not move to HDU 2");
    goto error;
  }

  for(c = 0; c < PSC_NUM_COLUMNS; c++) {
    /* Get column number */
    ffgcno(fits, CASEINSEN, colnames[c], &(h->gcols[c]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", colnames[c]);
      goto error;
    }
  }

  /* Determine optimum blocksize and allocate buffer */
  ffgrsz(fits, &(h->totbuf), &status);
  if(status) {
    fitsio_err(errstr, status, "ffgrsz");
    goto error;
  }

  h->buf = (struct psc_row *) malloc(h->totbuf * sizeof(struct psc_row));
  h->cbuf = (float *) malloc(h->totbuf * sizeof(float));
  h->bbuf = (unsigned char *) malloc(h->totbuf * sizeof(unsigned char));
  h->lbuf = (long *) malloc(h->totbuf * sizeof(long));
  h->pbuf = (char **) malloc(h->totbuf * sizeof(char *));
  if(!h->buf || !h->cbuf || !h->bbuf || !h->lbuf || !h->pbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Close file */
  ffclos(fits, &status);
  if(status) {
    fitsio_err(errstr, status, "ffclos");
    goto error;
  }

  h->f = (struct psc_file *) NULL;
  h->nf = 0;
  
  return(0);

 error:
  return(1);
}

int psc_close (struct psc_handle *h, char *errstr) {
  int status = 0, i;

  free((void *) h->buf);
  h->buf = (struct psc_row *) NULL;
  free((void *) h->cbuf);
  h->cbuf = (float *) NULL;
  free((void *) h->bbuf);
  h->bbuf = (unsigned char *) NULL;
  free((void *) h->lbuf);
  h->lbuf = (long *) NULL;
  free((void *) h->pbuf);
  h->pbuf = (char **) NULL;

  if(h->f) {
    for(i = 0; i < h->nf; i++) {
      if(h->f[i].fits) {
	ffclos(h->f[i].fits, &status);
	if(status) {
	  fitsio_err(errstr, status, "ffclos");
	  goto error;
	}
      }
      if(h->f[i].path)
	free((void *) h->f[i].path);
    }

    free((void *) h->f);
    h->f = (struct psc_file *) NULL;
  }

  return(0);

 error:

  return(1);
}

static int psc_findfiles (struct psc_handle *h, float ra_low, float ra_high, char *errstr) {
  float o_ra_low = 0.0, o_ra_high = 0.0;
  unsigned char wrapped = 0;
  int status = 0, i, ilim;
  struct psc_file *pf;
  char fnbuf[FLEN_FILENAME];

  /* Make sure the previous lot are cleaned up */
  if(h->f) {
    for(i = 0, pf = h->f; i < h->nf; i++, pf++) {
      ffclos(pf->fits, &status);
      if(status) {
	fitsio_err(errstr, status, "ffclos");
	goto error;
      }
    }

    free((void *) h->f);
    h->f = (struct psc_file *) NULL;
  }

  /* Calculate RA limits. */
  ra_low  *= RAD_TO_HR;
  ra_high *= RAD_TO_HR;

  if(ra_low < 0) {
    o_ra_low  = 24.0 + ra_low;
    o_ra_high = 24.0;
    
    ra_low = 0;
    wrapped = 1;
  }

  if(ra_high > 24.0) {
    o_ra_low  = 0;
    o_ra_high = ra_high - 24.0;
    
    ra_high = 24.0;
    wrapped = 1;
  }

  ilim = ARRAYLEN(psc_filelist);

  /* Allocate initial array big enough for everything */
  h->f = (struct psc_file *) malloc(ilim * sizeof(struct psc_file));
  if(!h->f) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  h->nf = 0;
  pf = h->f;

  memcpy(fnbuf, h->fpref, h->preflen);

  for(i = 0; i < ilim; i++)
    if(((psc_filelist[i].ra_low >= ra_low && psc_filelist[i].ra_low <= ra_high) ||
	(psc_filelist[i].ra_high > ra_low && psc_filelist[i].ra_high < ra_high) ||
	(ra_low >= psc_filelist[i].ra_low && ra_low < psc_filelist[i].ra_high) ||
	(ra_high >= psc_filelist[i].ra_low && ra_high < psc_filelist[i].ra_high)) ||
       (wrapped &&
	((psc_filelist[i].ra_low >= o_ra_low && psc_filelist[i].ra_low <= o_ra_high) ||
	 (psc_filelist[i].ra_high > o_ra_low && psc_filelist[i].ra_high < o_ra_high) ||
	 (o_ra_low >= psc_filelist[i].ra_low && o_ra_low < psc_filelist[i].ra_high) ||
	 (o_ra_high >= psc_filelist[i].ra_low && o_ra_high < psc_filelist[i].ra_high)))) {
      /* Found a candidate, form it's name */
      snprintf(fnbuf + h->preflen, sizeof(fnbuf) - h->preflen, "%s%s",
	       psc_filelist[i].name, h->fpost);

      /* Open it and move to the first extension */
      ffopen(&(pf->fits), fnbuf, READONLY, &status);
      if(status) {
	fitsio_err(errstr, status, "ffopen: %s", fnbuf);
	goto error;
      }
      
      ffmahd(pf->fits, 2, (int *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "could not move to HDU 2");
	goto error;
      }

      /* Get the number of rows in the table */
      ffgnrw(pf->fits, &(pf->nrows), &status);
      if(status) {
	fitsio_err(errstr, status, "Could not get table dimensions");
	goto error;
      }

      /* Get the range of RA covered by this catalogue. */
      (void) ffgky(pf->fits, TFLOAT, "RA-MIN", &(pf->min_ra), (char *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgky: RA-MIN");
	goto error;
      }
      (void) ffgky(pf->fits, TFLOAT, "RA-MAX", &(pf->max_ra), (char *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgky: RA-MAX");
	goto error;
      }

      /* Done, close file */
      ffclos(pf->fits, &status);
      if(status) {
	fitsio_err(errstr, status, "ffclos");
	goto error;
      }

      /* Keep a copy of the name */
      pf->path = strdup(fnbuf);
      if(!pf->path) {
	report_syserr(errstr, "malloc");
	goto error;
      }

      pf->fits = (fitsfile *) NULL;

      h->nf++;
      pf++;
    }
 
  /* Return any spare array */
  h->f = (struct psc_file *) realloc((void *) h->f, h->nf * sizeof(struct psc_file));
  if(!h->f) {
    report_syserr(errstr, "realloc");
    goto error;
  }
 
  return(0);

 error:
  if(h->f)
    free((void *) h->f);

  return(1);
}

long psc_find (struct psc_handle *h, float ra, float dec, float width, char *errstr) {
  float ra_low, ra_high, dec_low, dec_high;
  long start, finish, first_index, last_index, ntot = 0L;
  unsigned char ifile;

  int status = 0;
  struct psc_file *pf;
  float rec_dec;

  /* Calculate RA, Dec limits */
  radeclimits(ra, dec, width, &ra_low, &ra_high, &dec_low, &dec_high);

  h->ra = ra;
  h->dec = dec;
  h->width = width;

  /* Find the appropriate FITS files */
  if(psc_findfiles(h, ra_low, ra_high, errstr))
    goto error;

  /* For each file... */
  for(ifile = 0, pf = h->f; ifile < h->nf; ifile++, pf++) {
    /* Open it and move to the first extension */
    ffopen(&(pf->fits), pf->path, READONLY, &status);
    if(status) {
      fitsio_err(errstr, status, "ffopen: %s", pf->path);
      goto error;
    }
    
    ffmahd(pf->fits, 2, (int *) NULL, &status);
    if(status) {
      fitsio_err(errstr, status, "could not move to HDU 2");
      goto error;
    }

    /* Determine RA range */
    pf->ra_low  = ra_low;
    pf->ra_high = ra_high;
      
    /* Clip at range of this table */
    if(pf->ra_low < 0.0)
      pf->ra_low += 2 * M_PI;
    if(pf->ra_high > 2 * M_PI)
      pf->ra_high -= 2 * M_PI;
    
    if(pf->ra_low < pf->min_ra || pf->ra_low > pf->max_ra)
      pf->ra_low = pf->min_ra;
    if(pf->ra_high > pf->max_ra || pf->ra_high < pf->min_ra)
      pf->ra_high = pf->max_ra;

    start  = 1;
    finish = pf->nrows + 1;
    first_index = pf->nrows / 2;
    
    while(finish - start >= 2) {
      /* Fetch Dec of row 'first_index' */
      ffgcve(pf->fits, h->gcols[DEC_COLUMN], first_index, 1, 1, 0, &rec_dec,
	     (int *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgcve");
	goto error;
      }
      
      if(rec_dec < dec_low) {
	start = first_index;
	first_index = (first_index + finish) / 2;
      }
      else {
	finish = first_index;
	first_index  = (first_index + start) / 2;
      }
    }
    
    /* 'first_index' now gives the start Dec */
    
    start  = first_index;
    finish = pf->nrows + 1;
    last_index = (finish - start) / 2;
    
    while(finish - start >= 2) {
      /* Fetch Dec of row 'last_index' */
      ffgcve(pf->fits, h->gcols[DEC_COLUMN], last_index, 1, 1, 0, &rec_dec,
	     (int *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgcve");
	goto error;
      }
      
      if(rec_dec < dec_high) {
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

    pf->first_idx = first_index;
    pf->last_idx = last_index;
    
    ntot += last_index - first_index;
    
    /* Close file */
    ffclos(pf->fits, &status);
    if(status) {
      fitsio_err(errstr, status, "ffclos");
      goto error;
    }
    pf->fits = (fitsfile *) NULL;
  }

  pf = h->f;

  h->cf = 0;
  h->nbuf = 0;

  h->frow = pf->first_idx;
  h->nrem = pf->last_idx - pf->first_idx;

  /* Open first file */
  ffopen(&(pf->fits), pf->path, READONLY, &status);
  if(status) {
    fitsio_err(errstr, status, "ffopen: %s", pf->path);
    goto error;
  }
  
  ffmahd(pf->fits, 2, (int *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "could not move to HDU 2");
    goto error;
  }

  /* Fill buffer */
  if(psc_buffer(h, errstr))
    goto error;

  return(ntot);

 error:
  return(-1);
}

static int psc_buffer (struct psc_handle *h, char *errstr) {
  int status = 0;
  long r;
  struct psc_file *pf;

  pf = h->f + h->cf;

  h->nbuf = (h->nrem > h->totbuf ? h->totbuf : h->nrem);

#define READ_COL(num, name) {						\
  ffgcve(pf->fits, h->gcols[(num)], h->frow, 1, h->nbuf, 0, h->cbuf,	\
	 (int *) NULL, &status);					\
  if(status) {								\
    fitsio_err(errstr, status, "ffgcve");				\
    goto error;								\
  }									\
									\
  for(r = 0; r < h->nbuf; r++)						\
    (h->buf + r)->name = *(h->cbuf + r);				\
}

#define READ_COL_CHAR(num, name) {					\
  for(r = 0; r < h->nbuf; r++)						\
    h->pbuf[r] = &(h->buf[r].name[0]);					\
									\
  ffgcvs(pf->fits, h->gcols[(num)], h->frow, 1, h->nbuf, "", h->pbuf,	\
	 (int *) NULL, &status);					\
  if(status) {								\
    fitsio_err(errstr, status, "ffgcvs");				\
    goto error;								\
  }									\
}

#define READ_COL_BYT(num, name) {					\
  ffgcvb(pf->fits, h->gcols[(num)], h->frow, 1, h->nbuf, 0, h->bbuf,	\
	 (int *) NULL, &status);					\
  if(status) {								\
    fitsio_err(errstr, status, "ffgcvb");				\
    goto error;								\
  }									\
									\
  for(r = 0; r < h->nbuf; r++)						\
    (h->buf + r)->name = *(h->bbuf + r);				\
}

#define READ_COL_LNG(num, name) {					\
  ffgcvj(pf->fits, h->gcols[(num)], h->frow, 1, h->nbuf, 0, h->lbuf,	\
	 (int *) NULL, &status);					\
  if(status) {								\
    fitsio_err(errstr, status, "ffgcvj");				\
    goto error;								\
  }									\
									\
  for(r = 0; r < h->nbuf; r++)						\
    (h->buf + r)->name = *(h->lbuf + r);				\
}

  READ_COL(0, ra);
  READ_COL(1, dec);
  READ_COL(2, maj);
  READ_COL(3, min);
  READ_COL(4, pa);
  READ_COL(5, j.m);
  READ_COL(6, j.msig);
  READ_COL(7, j.msigcom);
  READ_COL_CHAR(8, j.qual);
  READ_COL_BYT(9, j.rd);
  READ_COL_BYT(10, j.bl);
  READ_COL_CHAR(11, j.cc);
  READ_COL(12, h.m);
  READ_COL(13, h.msig);
  READ_COL(14, h.msigcom);
  READ_COL_CHAR(15, h.qual);
  READ_COL_BYT(16, h.rd);
  READ_COL_BYT(17, h.bl);
  READ_COL_CHAR(18, h.cc);
  READ_COL(19, k.m);
  READ_COL(20, k.msig);
  READ_COL(21, k.msigcom);
  READ_COL_CHAR(22, k.qual);
  READ_COL_BYT(23, k.rd);
  READ_COL_BYT(24, k.bl);
  READ_COL_CHAR(25, k.cc);
  READ_COL_LNG(26, ptskey);

#undef READ_COL
#undef READ_COL_LNG

  for(r = 0; r < h->nbuf; r++)
    h->buf[r].ptr = h->frow + r;

  h->frow += h->nbuf;
  h->nrem -= h->nbuf;
  h->brow = 0;  /* NOTE: this counts from zero! */

  return(0);

 error:

  return(1);
}

int psc_read (struct psc_handle *h, struct psc_row *row, char *errstr) {
  int status = 0;
  float rra, rdec, xi, xn;
  struct psc_file *pf;

  if(h->cf >= h->nf)
    return(0);

  pf = h->f + h->cf;

 fetch_row:
  if(h->brow >= h->nbuf) {
    /* Ran out of rows in buffer, need to fetch more */
    if(h->nrem < 1) {
      /* Finished that file, get the next */
      ffclos(pf->fits, &status);
      if(status) {
	fitsio_err(errstr, status, "ffclos");
	goto error;
      }
      pf->fits = (fitsfile *) NULL;
      
      h->cf++;
      
      if(h->cf >= h->nf)
	return(0);
      
      pf++;
      h->nbuf = 0;
      
      /* Open it */
      ffopen(&(pf->fits), pf->path, READONLY, &status);
      if(status) {
	fitsio_err(errstr, status, "ffopen: %s", pf->path);
	goto error;
      }
      
      ffmahd(pf->fits, 2, (int *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "could not move to HDU 2");
	goto error;
      }
      
      h->frow = pf->first_idx;
      h->nrem = pf->last_idx - pf->first_idx;
    }

    if(psc_buffer(h, errstr))
      goto error;
  }

  rra  = (h->buf + h->brow)->ra;
  rdec = (h->buf + h->brow)->dec;

  /* Reject rows not in the correct RA range */
  if(!(rra > pf->ra_low && rra < pf->ra_high)) {
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
  memcpy((void *) row, (void *) (h->buf + h->brow), sizeof(struct psc_row));
  row->xi = xi;
  row->xn = xn;
  h->brow++;

  return(1);

 error:
  return(-1);
}
