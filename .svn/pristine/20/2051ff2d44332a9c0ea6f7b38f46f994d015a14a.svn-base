/*

$Id: ap.h,v 1.2 2010/09/06 08:57:12 jim Exp $

*/
#ifndef APDEF
#define APDEF

typedef struct {
    int x;
    int y;
    float z;
    float zsm;
    int iobj;
} plstruct;

typedef struct {
    short int areal[8];	/* height above thresh of areal-prof cuts */
    int lsiz;		/* size of a line */
    int csiz;           /* size of a column */
    int maxip;		/* max no. of parents ever used. */
    int maxbl;		/* size of pixel-storage block stack */
    int maxpa;		/* size of parent-stack. */
    int ipnop;		/* parent-number-of-pixels, min size of image */
    int nimages;	/* count of images */
    int ipstack;	/* parent-name stack pointer */
    int ibstack;	/* pixel-block name stack pointer */
    float thresh;       /* threshold for image detection */
    float background;   /* background value */
    float sigma;	/* median background sigma */
    int multiply;       /* smoothing multiplication */
    float xintmin;      /* minimum intensity for consideration */
    int mulpix;         /* minimum size for considering multiple images */
    float areal_offset; /* offset in areal profile levels */
    float fconst;       /* Normalisation constant for areal profiles */
    float saturation;   /* saturation level from background analysis */
    int icrowd;         /* true if deblending routine is to be used */

    int *blink;		/* block-link array */
    int *bstack;	/* stack of pixel names */
    struct {		/* Image control block array */
	int first;	/* link to first data block */
	int last;	/* current last block	*/
	int pnop;	/* Parent no. pixels (-1 = inactive) */
	int growing;
	int touch;	/* 0 = does not touch an edge */
	int pnbp;       /* Parent no of bad pixels */
    } *parent;

    short int *pstack;	/* stack of parent names */
    plstruct *plessey;  /* x,y,i storage array */
    short int *lastline;/* Parents on last line */

    float *data;        /* Pointer to original image data */
    short int *conf;    /* Pointer to original confidence map */
    unsigned char *mflag; /* Pointer to mflag array for tracking merges */
    unsigned char *opm; /* Object pixel mask */
    float rcore;        /* Core radius for aperture photometry */
    float filtfwhm;     /* FWHM of smoothing kernel for detection algorithm */
    plstruct *plarray;  /* Plessey structure workspace for passing data to 
			   various processing routines */
    int npl;            /* Size of the above */
    int npl_pix;        /* Number of pixels in the above structure */
    int nbsize;         /* size of smoothing cells for background estimation */
    
    struct {
	int nbx;        /* X dimension of background map */
	int nby;        /* Y dimension of background map */
	int nbsize;     /* Size of a side of background map cell */
        float **bvals;  /* Pointer to background map */
    } backmap;
} ap_t;

typedef struct {
    float x;	        /* x position 				*/
    float y;	        /* y position				*/
    float total;        /* total integrated intensity		*/
    int area;           /* image area in pixels			*/
    float peak;	        /* peak image intensity above sky	*/
    float xx;           /* 2nd moment x				*/
    float xy;	        /* 2nd moment cross term		*/
    float yy;	        /* 2nd moment y				*/
    float ecc;          /* Eccentricity                         */
    int areal[8];       /* areal profile of image		*/
} apmCat_t;
#endif

#define PNOPDEF     5
#define GRIDDEF     64
#define THRESHDEF   -2.0
#define MAXBL       500000
#define NAREAL      8
/*

$Log: ap.h,v $
Revision 1.2  2010/09/06 08:57:12  jim
Added nbsize to structure

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.3  2008/02/11 12:30:04  jim
Upped the number of pixels in MAXBL

Revision 1.2  2006/07/31 13:21:09  jim
Modified imcore now allows for a smoothing kernel with variable FWHM

Revision 1.1  2004/04/02 10:54:56  jim
New version for rewrite of imcore

Revision 1.3  2002/07/08 11:51:31  jim
Increased value of MAXBL again

Revision 1.2  2002/07/04 14:43:45  jim
Raised the value of MAXBL to allow for large objects

Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS


*/

