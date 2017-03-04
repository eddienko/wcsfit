#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "floatmath.h"
#include "util.h"

#define MIN8(a) \
  MIN(a[0], MIN(a[1], MIN(a[2], MIN(a[3], MIN(a[4], MIN(a[5], MIN(a[6], a[7])))))))
#define MAX8(a) \
  MAX(a[0], MAX(a[1], MAX(a[2], MAX(a[3], MAX(a[4], MAX(a[5], MAX(a[6], a[7])))))))

void radeclimits (float ra, float dec, float width,
		  float *ra_low, float *ra_high, float *dec_low, float *dec_high) {
  float ralim[8], declim[8];

  /* Get RA, Dec limits */
  XIXN(-width, -width, ra, dec, ralim[0], declim[0]);
  XIXN(-width,  width, ra, dec, ralim[1], declim[1]);
  XIXN( width, -width, ra, dec, ralim[2], declim[2]);
  XIXN( width,  width, ra, dec, ralim[3], declim[3]);
  XIXN(   0.0, -width, ra, dec, ralim[4], declim[4]);
  XIXN(   0.0,  width, ra, dec, ralim[5], declim[5]);
  XIXN(-width,    0.0, ra, dec, ralim[6], declim[6]);
  XIXN( width,    0.0, ra, dec, ralim[7], declim[7]);

  /* Check if the poles are within the range */
  if(dec + width > M_PI_2 && dec - width < -M_PI_2) {
    /* Contains both poles, RA and Dec range of all */
    *ra_low  = 0.0;
    *ra_high = 2 * M_PI;

    *dec_low  = -M_PI_2;
    *dec_high = M_PI_2;
  }
  else if(dec + width > M_PI_2) {
    /* Contains the North pole, RA range of all */
    *ra_low  = 0.0;
    *ra_high = 2 * M_PI;

    *dec_low  = MIN8(declim);
    *dec_high = M_PI_2;
  }
  else if(dec - width < -M_PI_2) {
    /* Contains the South pole, RA range of all */
    *ra_low  = 0.0;
    *ra_high = 2 * M_PI;

    *dec_low  = -M_PI_2;
    *dec_high = MAX8(declim);
  }
  else {
    *ra_low   = MIN8(ralim);
    *ra_high  = MAX8(ralim);
    *dec_low  = MIN8(declim);
    *dec_high = MAX8(declim);
  }
}

