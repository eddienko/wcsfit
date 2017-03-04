#ifndef __USNOB_H__
#define __USNOB_H__

/* Internal defines */
#define USNOB_RECSIZE    80
#define USNOB_BLKNREC  1024
#define USNOB_BUF_RA      0
#define USNOB_BUF_SPD     4
#define USNOB_BUF_PM      8
#define USNOB_BUF_PMERR  12
#define USNOB_BUF_POSERR 16
#define USNOB_BUF_MAG    20
#define USNOB_BUF_MAGERR 40
#define USNOB_BUF_INDEX  60

#define FNBUFSIZE 1024

struct usnob_band {
  float m;
  short field;
  unsigned char survey;
#define USNOB_SURVEY_POSS_I_O          0
#define USNOB_SURVEY_POSS_I_E          1
#define USNOB_SURVEY_POSS_II_J         2
#define USNOB_SURVEY_POSS_II_F         3
#define USNOB_SURVEY_SERC_J            4
#define USNOB_SURVEY_ESO_R             5
#define USNOB_SURVEY_AAO_R             6
#define USNOB_SURVEY_POSS_II_N         7
#define USNOB_SURVEY_SERC_I            8
#define USNOB_SURVEY_POSS_II_N_SERC_I  9
  unsigned char star;

  float xires;
  float xnres;
  unsigned char calsrc;

  int index;
};

struct usnob_row {
  float xi;
  float xn;
  long ptr;

  float ra;
  float dec;

  float pmra;
  float pmdec;
  float pmprob;

  float pmerrra;
  float pmerrdec;
  float sigra;
  float sigdec;
  unsigned char ndet;

  float raerr;
  float decerr;
  float epoch;

  unsigned char flags;
#define USNOB_FLAG_NONE       0
#define USNOB_FLAG_MOTION     1
#define USNOB_FLAG_SPIKE      2
#define USNOB_FLAG_YSCORR     4

  struct usnob_band b1;
  struct usnob_band r1;
  struct usnob_band b2;
  struct usnob_band r2;
  struct usnob_band n;
};

struct usnob_file {
  FILE *fp;
  char path[FNBUFSIZE];
  long nrows;

  float min_dec;
  float max_dec;

  long rec;   /* File record pointer */

  unsigned char nidx;
  long first_idx[2];
  long last_idx[2];

  float ra_low;
  float ra_high;
};

struct usnob_handle {
  char fname[FNBUFSIZE];
  int fnptr;
  unsigned char swap;

  struct usnob_file *f;
  int nf;
  int cf;

  struct usnob_row *buf;
  char *cbuf;
  long nbuf;

  unsigned char ridx;

  long frow;  /* Current row in file */
  long brow;  /* Current row in buffer */
  long nrem;

  float dec_low;
  float dec_high;

  float ra;
  float dec;
  float width;
};

int usnob_open (struct usnob_handle *h, char *path, char *errstr);
int usnob_close (struct usnob_handle *h, char *errstr);
long usnob_find (struct usnob_handle *h, float ra, float dec, float width, char *errstr);
int usnob_read (struct usnob_handle *h, struct usnob_row *row, char *errstr);

#endif  /* __USNOB_H__ */
