/*

$Id: floatmath.h,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/

#ifndef __FLOATMATH_H__
#define __FLOATMATH_H__

#include <math.h>

/* Oft used constants */
 
#define M_PI            3.14159265358979323846  /* pi */
#define M_SQRT2         1.41421356237309504880  /* sqrt(2) */
#define M_LN2           0.69314718055994530942  /* log_e 2 */
#define M_PI_2          1.57079632679489661923  /* pi/2 */
#define M_2_SQRTPI      1.12837916709551257390  /* 2/sqrt(pi) */
 
#if defined(__SUNPRO_C)

/* The Sun Workshop compiler suite has the 'float' math functions
 * in sunmath.h.
 */

#include <sunmath.h>

#else

/* Emulate 'float' versions of the math.h functions we use on systems which
 * don't have them */

/* sqrtf */
#define sqrtf(a) ((float) sqrt((double) (a)))

/* fabsf */
#define fabsf(a) ((float) fabs((double) (a)))

/* logf */
#define logf(a) ((float) log((double) (a)))

/* log10f */
#define log10f(a) ((float) log10((double) (a)))

/* expf */
#define expf(a) ((float) exp((double) (a)))

/* sinf */
#define sinf(a) ((float) sin((double) (a)))

/* cosf */
#define cosf(a) ((float) cos((double) (a)))

/* tanf */
#define tanf(a) ((float) tan((double) (a)))

/* asinf */
#define asinf(a) ((float) asin((double) (a)))

/* acosf */
#define acosf(a) ((float) acos((double) (a)))

/* atanf */
#define atanf(a) ((float) atan((double) (a)))

/* atan2f */
#define atan2f(x, y) ((float) atan2((double) (x), (double) (y)))

/* powf */
#define powf(x, y) ((float) pow((double) (x), (double) (y)))

#endif

#endif  /* __FLOATMATH_H__ */

/*

$Log: floatmath.h,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.2  2004/09/07 14:18:56  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:54:58  jim
New version for rewrite of imcore

Revision 1.2  2003/09/01 13:00:01  jim
New version

Revision 1.1  2002/09/10 08:23:48  jim
Initial entry into CVS


*/
