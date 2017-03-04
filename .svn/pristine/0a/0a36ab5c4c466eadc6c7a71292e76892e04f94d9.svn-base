/* 

$Id: wcsfit.h,v 1.2 2010/09/06 09:05:59 jim Exp $

*/

#include <fitsio.h>
#include <tools.h>


extern int wcsfit(char *infile, char *incat, char *catsrc, char *site, 
		  char *catpath, char *errmsg);
extern int cir_getstds(char *infile, char *outfile, char *catsrc, char *site,
		       char *catpath, float equinox, int fudge, int cache,
		       char *errmsg);
extern int cir_matchstds(char *infxy, char *infradec, float srad, int nx,
                         int ny, char *output, int *nm, char *errmsg);
extern int cir_platesol(char *image, char *posfile, int nconst, int pass, 
			int fixtanpt, char *errmsg);

/*

$Log: wcsfit.h,v $
Revision 1.2  2010/09/06 09:05:59  jim
Tidied some docs


*/
