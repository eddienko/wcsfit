/*

$Id: misc.h,v 1.4 2012/12/08 07:28:48 jim Exp $

*/

#include <fitsio.h>

/* Miscellaneous routines taken from parts of cirdr. */

extern int cir_nint(float value);
extern int cir_open_output(char *newfile, char *template, fitsfile **iptr,
			   int recipe, int newbitpix, int newaxis,
			   long *newaxes, char *errmsg);
extern int cir_crfile(char *newfile, fitsfile **iptr, char *errmsg);
extern int cir_copy_hdu(fitsfile *iptr, fitsfile *tptr, int recipe, 
			int newbitpix, int newnaxis, long *newnaxes, 
			char *errmsg);
extern int cir_copy_cards(fitsfile *iptr, fitsfile *tptr);
extern int cir_qmedsig(float *data, unsigned char *bpm, int npts, 
		       float thresh, int niter, float lowv, float highv,
		       float *median,float *sigma, char *errmsg);
extern float histexam(int *histo, int nhist, int level);
extern int cir_med(float *data, unsigned char *bpm, int npts, float *value, 
		   char *errmsg);
extern int cir_dmed(double *data, unsigned char *bpm, int npts, double *value,
		    char *errmsg);
extern int cir_mean(float *data, unsigned char *bpm, int npts, float *value, 
		    char *errmsg);
extern int cir_dmean(double *data, unsigned char *bpm, int npts, double *value,
		     char *errmsg);
extern int cir_medmad(float *data, unsigned char *bpm, int npts, float *medval,
		      float *madval, char *errmsg);
extern float kselect(float *a, int n, int k);
extern double dkselect(double *a, int n, int k);
extern int cir_sumbpm(unsigned char *bpm, int npts, int *sumb, char *errmsg);
extern int fndmatch(float x, float y, float *xlist, float *ylist, int nlist, 
		    float err);
extern int casu_stamp(fitsfile *fptr, char *mainprog);

/*

$Log: misc.h,v $
Revision 1.4  2012/12/08 07:28:48  jim
Added cir_medmad

Revision 1.3  2010-09-06 09:04:58  jim
added cir_mean


*/
