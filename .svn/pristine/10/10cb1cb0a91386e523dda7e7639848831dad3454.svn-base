/*

$Id: tools.h,v 1.3 2010/11/01 10:24:39 jim Exp $

*/

#include <fitsio.h>
#include <misc.h>
#include <cir_memory.h>
#include <casu_wcs.h>

/* Various definitions probably stolen from cirdr */

#define CIR_OK    0
#define CIR_WARN  1
#define CIR_FATAL 2

/* Codes needed for various file creation schemes...*/

#define DIRECTCOPY 1
#define MAKE2D     2
#define NEWBITPIX  3
#define NEWNAXES   4
#define NEWDATASZ  5

/* A number of extra header records just to make sure you have enough space
   so that you don't have to rewrite the whole file */

#define MOREKEYS 36

/* Degrees to radians */

#define DEGRAD  0.017453292519943

/* Useful macros */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

/* Macros that deal with tidying up */

#define freespace(_p) if (_p != NULL) {free(_p); _p = NULL;}
#define closefits(_p) {int status = 0; if (_p != NULL) {(void)fits_close_file(_p,&status); _p = NULL;}}
#define deletefits(_p) {int status = 0; if (_p != NULL) {(void)fits_delete_file(_p,&status); _p = NULL;}}
#define closefile(_p) if (_p != NULL) {fclose(_p); _p = NULL;}
#define freewcs(_p) if (_p != NULL) {cir_wcsclose(_p); _p = NULL;}

/* Final routine declarations */

extern int cir_catcoord(char *infile, char *incat, char *errmsg);

/*

$Log: tools.h,v $
Revision 1.3  2010/11/01 10:24:39  jim
Modified macros deletefits and closefits

Revision 1.2  2010-09-06 09:05:10  jim
Tidied some docs


*/
