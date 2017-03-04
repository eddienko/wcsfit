/*

$Id: casu_wcs.h,v 1.4 2013/05/30 14:36:23 jim Exp $

*/

#include <wcs.h>
#include <wcshdr.h>

extern void cir_xytoradec(struct wcsprm *wcs, double x, double y, 
			  double *ra, double *dec);
extern void cir_radectoxy(struct wcsprm *wcs, double ra, double dec,
			  double *x, double *y);
extern void cir_wcsclose(struct wcsprm *wcs);
extern int cir_wcsopen(char *infile, struct wcsprm **wcs, char *errmsg);
extern void cir_radectoxieta(struct wcsprm *wcs, double ra, double dec,
			     double *xi, double *eta);
extern void cir_xytoxy_list(struct wcsprm *wcs1, struct wcsprm *wcs2, int nc,
			    double *x1, double *y1, double *ra, double *dec, 
			    double *x2, double *y2);
/*

$Log: casu_wcs.h,v $
Revision 1.4  2013/05/30 14:36:23  jim
Removed extraneous declaration

Revision 1.3  2010-11-22 10:42:42  jim
Added cir_xytoxy_list

Revision 1.2  2010-09-06 09:05:10  jim
Tidied some docs


*/
