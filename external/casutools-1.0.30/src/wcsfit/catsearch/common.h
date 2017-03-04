#ifndef __COMMON_H__
#define __COMMON_H__

#include "floatmath.h"

/* --Macros-- */

/* RA, Dec -> standard coordinates */

#define STANDC(tpa, tpd, a, d, xi, xn) {	\
  float sd, cd, tdec, c, denom;			\
						\
  sd = sinf(tpd);				\
  cd = cosf(tpd);				\
  tdec = tanf(d);				\
  c = cosf((a) - (tpa));			\
  denom = sd * tdec + cd * c;			\
  (xi) = sinf((a) - (tpa)) / denom;		\
  (xn) = (cd * tdec - sd * c) / denom;		\
}

/* Standard coordinates -> Ra, Dec */

#define XIXN(xi, xn, tpa, tpd, a, d) {				\
  float tand, secd, aa;						\
								\
  tand = tanf(tpd);						\
  secd = 1.0 / cosf(tpd);					\
  aa = atanf((xi) * secd / (1.0 - (xn) * tand));		\
								\
  (a) = aa + tpa;						\
  if((xi) == 0.0)						\
    (d) = xn + tpd;						\
  else								\
    (d) = atanf(((xn) + tand) * sinf(aa) / ((xi) * secd));	\
								\
  if(fabsf((d) - (tpd)) > M_PI_2) {				\
    (d) = -(d);							\
    (a) += M_PI;						\
  }								\
}

/* --Prototypes-- */

void radeclimits (float ra, float dec, float width,
		  float *ra_low, float *ra_high, float *dec_low, float *dec_high);

#endif  /* __COMMON_H__ */
