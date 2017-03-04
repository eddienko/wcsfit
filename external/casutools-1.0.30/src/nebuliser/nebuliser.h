/*

$Id: nebuliser.h,v 1.2 2010/09/20 09:04:04 jim Exp $

*/

#define MEANCALC    1
#define MEDIANCALC  2

extern int cir_2dfilt(char *infile, char *confmap, int nfiltmed, int nfiltlin,
		      char *backmap, int niter, int axis, int twod, 
		      int takeout_sky, int inorm, float signeg, float sigpos, 
		      char *errmsg);
extern void twodfilt(float *data, unsigned char *bpm, int nx, int ny, 
		     int medfilt, int linfilt, int niter, int axis, 
		     int twod, int takeout_sky, int inorm, int wantback, 
		     float signeg, float sigpos, float **backmap);
extern void cir_bfilt2(float *data, unsigned char *bpm, int nx, int ny, 
		       int filt, int stat, int axis);
extern void cir_bfilt_2d(float *data, unsigned char *bpmcopy, int nx, int ny,
			 int filt, int stat);
extern void dostat(float *data, unsigned char *bpm, unsigned char *goodval,
		   int npts, int nfilt, int whichstat);

/*

$Log: nebuliser.h,v $
Revision 1.2  2010/09/20 09:04:04  jim
Added missing declaration to cir_2dfilt

Revision 1.1  2010-09-06 09:03:50  jim
New entry


*/
