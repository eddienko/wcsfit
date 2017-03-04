/*

$Id: util.h,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/

#ifndef __UTIL_H__
#define __UTIL_H__

/* --Macros-- */

#undef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#undef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define ARRAYLEN(a) (sizeof(a)/sizeof((a)[0]))

#define NINT(a) ((int) ((a)+((a) < 0 ? -0.5 : 0.5)))

#define ARRAYELEM2(a, s, x, y) ((a)[(x)*(s) + (y)])

#endif
/*

$Log: util.h,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.1  2004/04/02 10:55:01  jim
New version for rewrite of imcore

Revision 1.1  2002/09/10 08:24:05  jim
Initial entry into CVS


*/
