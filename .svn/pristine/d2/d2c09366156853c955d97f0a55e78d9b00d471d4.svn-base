#ifndef __PSC_H__
#define __PSC_H__

#include <fitsio.h>

/* Number of FITS file columns we use */
#define PSC_NUM_COLUMNS 27

struct psc_band {
  float m;
  float msig;
  float msigcom;
  char qual[2];
  unsigned char rd;
  unsigned char bl;
  char cc[2];
};

struct psc_row {
  float xi;
  float xn;
  long ptr;

  float ra;
  float dec;
  float maj;
  float min;
  float pa;

  struct psc_band j;
  struct psc_band h;
  struct psc_band k;

  long ptskey;
};

struct psc_file {
  fitsfile *fits;
  char *path;
  long nrows;

  float min_ra;
  float max_ra;

  long first_idx;
  long last_idx;

  float ra_low;
  float ra_high;
};

struct psc_handle {
  char fpref[FLEN_FILENAME];
  int preflen;
  char fpost[FLEN_FILENAME];
  int gcols[PSC_NUM_COLUMNS];

  struct psc_file *f;
  int nf;
  int cf;  /* Current file */

  struct psc_row *buf;
  float *cbuf;
  unsigned char *bbuf;
  long *lbuf;
  char **pbuf;
  long totbuf;
  long nbuf;

  long frow;  /* Current row in file */
  long brow;  /* Current row in buffer */
  long nrem;

  float ra;
  float dec;
  float width;
};

int psc_open (struct psc_handle *h, char *path, char *errstr);
int psc_close (struct psc_handle *h, char *errstr);
long psc_find (struct psc_handle *h, float ra, float dec, float width, char *errstr);
int psc_read (struct psc_handle *h, struct psc_row *row, char *errstr);

#endif  /* __PSC_H__ */
