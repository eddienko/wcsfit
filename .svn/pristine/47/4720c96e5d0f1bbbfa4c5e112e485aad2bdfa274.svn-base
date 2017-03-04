/*

$Id: cir_classify.h,v 1.2 2011/01/12 13:12:56 jim Exp $

*/
#define MAXHIST2  111
#define STEP      0.05
#define NSAMPLE   150
#define MAXLOOP   5
#define BLIMDEF   15.0;
#define FLIMDEF   11.0;
#define CMINDEF   7.5
#define CMAXDEF   15.0
#define NAREAL    8

#define COREMAG(A,B,C) 2.5*log10((double)(max(C,A-B)))

/* Make the data arrays and header values global */

long nrows;
float thresh,skylevel,skynoise,rcore,exptime;

/* Derived values */

int poor;
float sigell,fitell,elllim,sigellf,fitellf,sigpa,fitpa;
float blim,flim,cmin,cmax;
float fit1,fit2,fit3,fit4,fit5,sigma1,sigma2,sigma3,sigma4,sigma5;
float fit6,fit7,sigma6,sigma7;
float fit_final,sigma_final;
float *lower1,*lower2,*lower3,*upper1,*upper2,*upper3,*lower,*upper,*uppere;
float avsig1,avsig2,avsig3,wt1,wt2,wt3;

/* Classification values */

int nstar,ngal,njunk,ncmp;

/* Values for the data quality and aperture corrections */

float avsat,corlim,cormin,apcpkht,apcor,apcor1,apcor2,apcor3,apcor4,apcor5;
float apcor6,apcor7;

/* Subroutine prototypes */

static void anhist(float *, int, float *, float *);
static void boundaries(float *, float *, float *, float, float, float, float, 
		       int, float, float, float *, float *, float *, float *);
static void boundpk(float *, float *, float, float, float *, float *, 
		    float *, float *);
static void classify();
static void classstats(float *, float *, int, float, float *, float *);
static void classstats_ap0(float *, float *);
static void classstats_ap67(float *, float *, float *, float *);
static void classstats_el();
static void classstats_ellf(float);
static void classstats_pa();
static void classstats_final();
static void medstat(float *, int, float *, float *);
static void sort1(float *, int);
static void sort2(float *, float *, int);
/*

$Log: cir_classify.h,v $
Revision 1.2  2011/01/12 13:12:56  jim
Added pa results, classstats_ap67 and classstats_pa

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.6  2009/11/15 18:42:54  jim
Added argument to classstats

Revision 1.5  2005/11/25 15:41:02  jim
fixed stupid typo

Revision 1.4  2005/05/13 09:24:55  jim
Modified some declarations

Revision 1.3  2004/08/26 07:42:12  jim
Modified to bring it up to date with MJI's latest version

Revision 1.2  2002/12/12 11:43:52  jim
Changed declaration of supporting routines to static

Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS


*/
