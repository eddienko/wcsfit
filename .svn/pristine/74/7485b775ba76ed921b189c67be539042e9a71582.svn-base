/*

$Id: imcore.h,v 1.2 2010/11/01 10:22:55 jim Exp $

*/

/* Required information */

#include <fitsio.h>
#include "ap.h"
#include "imcore_list.h"

/* imcore specific macros */

#define MINHISTVAL -1000 /* Minimum value to histogram for background stats */
#define MAXHISTVAL 65535 /* Maximum value to histogram for background stats */
#define MAXHIST (MAXHISTVAL-MINHISTVAL+1) /* maximum size of histogram array */
#define HIST_ELEM(a, i) ((a)[(i) - MINHISTVAL]);

/* Catalogue generation parameters */

#define MINSATURATE 20000 /* Minimum background saturation level */
#define IMNUM 200         /* Maximum number of images to be deblended */
#define NPAR 16           /* Number of parameters in a basic results array */
#define NAREAL 8          /* Number of areal profiles */
#define IDBLIM 10000      /* Maximum number of pixels to use in deblending */

#define STUPID_VALUE -1000.0 /* Minimum value of a pixel */

/* MFLAG values used for tracking the quality of individual pixels */

#define MF_CLEANPIX     0
#define MF_OBJPIX       1
#define MF_SATURATED    2
#define MF_ZEROCONF     3
#define MF_STUPID_VALUE 4
#define MF_3SIG         5
#define MF_POSSIBLEOBJ  6

/* Flags for types of catalogues to be created */

#define CAT_INTWFC      1   /* Original 32 column catalogue used for INT WFC */
#define CAT_WFCAM       2   /* 80 column catalogue for WFCAM */
#define CAT_BASIC       3   /* Very basic with positions, moments and areals */
#define CAT_OBJMASK     4   /* No table, just an object mask */
#define CAT_VIRCAM      6   /* 80 column catalogue for VIRCAM */

/* Useful constants */

#define M_PI            3.14159265358979323846  /* pi */
#define M_SQRT2         1.41421356237309504880  /* sqrt(2) */
#define M_LN2           0.69314718055994530942  /* log_e 2 */

/* Tidy-up macros */

#define freespace(_p) if (_p != NULL) {free(_p); _p = NULL;}
#define closefits(_p) {int status = 0; if (_p != NULL) {(void)fits_close_file(_p,&status); _p = NULL;}}
#define closefile(_p) if (_p != NULL) {fclose(_p); _p = NULL;}

/* External Function Prototypes. First, main processing routines */

extern int imcore_conf(char *, char *, int, float, int, float, int, float, 
		       char *, char *, int, int, char *);
extern int imcore_colour(char *, char *, int, float, int, float, int, float, 
			 char *, char *, int, int, float, char *);
extern int imcore_list(char *, char *, char *, char *, float, float, int, 
		       char *, char *, int, int, char *);
extern int imcore_opm(char *, char *, int, float, int, float, char *, int, int,
		      char *);
extern int imcore_opm_dataonly(float *, short int *, int, int, int, float, int, 
			       float, int, int, unsigned char **, int *, 
			       char *);
extern int imcore_rdbuf_mef(char *, int, void **, long *, long *, int, char *);
extern int imcore_rdbuf_colour(char *, int, void **, long *, long *, int, 
			       char *);
extern int imcore_rdbuf_conf(char *, short int **, float **, long, long, int, 
			     char *);

extern int imcore_background(ap_t *, int, float, int, char *);
extern int imcore_backstats(ap_t *, float, int, float *, float *, float *, 
			    char *);
extern void imcore_backest(ap_t *, float, float, float *, float *);
extern void imcore_medsig(int *, int, int, int, float *, float *);


extern int extend(ap_t *, float, float, float, float, float, 
		  float, float, float, float *);
extern void overlp(ap_t *ap, float [IMNUM][NPAR], int *, float, float, float,
		   int, float);
extern void phopt(ap_t *, float [IMNUM][NPAR], int, int, float [], 
		  float [], float [], int, float []);
extern float fraction(float, float, float);
extern void seeing(ap_t *, int, float *, float *, float **, float *, float *);

/* Filter routines */

extern void bfilt(float **, int, int);
extern void filt1d(float [], int, int);
extern void padext(float [], int);
extern void hanning(float [], int);
extern void median(float[], int, int);
extern void polynm(float [], float [], int, float [], int, int);
extern void solve (double a[25][25], double b[25], int m);

/* Routines that generate and read the catalogues */

extern int tabinit(ap_t *,char *, char *, char *);
extern int tabclose(ap_t *, char *);
extern int do_seeing(ap_t *, char *);
extern int process_results(ap_t *, char *);
extern int process_results_list(ap_t *, objstruct *[], int, char *);
extern int readtab_list(char *, objstruct **, int *, int, char *);
extern int get_ncol(int);
extern int tabinit_1(char *, char *, char *);
extern int do_seeing_1(ap_t *, char *);
extern int process_results_1(ap_t *, char *);
extern int process_results_list_1(ap_t *, objstruct *[], int, char *);
extern int tabclose_1(ap_t *, char *);
extern int readtab_list_1(char *, objstruct **, int *, char *);
extern int get_ncol_1(void);
extern int tabinit_2(char *, char *, char *);
extern int do_seeing_2(ap_t *, char *);
extern int process_results_2(ap_t *, char *);
extern int process_results_list_2(ap_t *, objstruct *[], int, char *);
extern int tabclose_2(ap_t *, char *);
extern int readtab_list_2(char *, objstruct **, int *, char *);
extern int get_ncol_2(void);
extern int tabinit_3(char *, char *, char *);
extern int do_seeing_3(ap_t *, char *);
extern int process_results_3(ap_t *, char *);
extern int process_results_list_3(ap_t *, objstruct *[], int, char *);
extern int tabclose_3(ap_t *, char *);
extern int readtab_list_3(char *, objstruct **, int *, char *);
extern int get_ncol_3(void);
extern int tabinit_4(ap_t *, char *, char *, char *);
extern int do_seeing_4(ap_t *, char *);
extern int process_results_4(ap_t *, char *);
extern int tabclose_4(ap_t *, char *);
extern int get_ncol_4(void);
extern int tabinit_6(char *, char *, char *);
extern int do_seeing_6(ap_t *, char *);
extern int process_results_6(ap_t *, char *);
extern int process_results_list_6(ap_t *, objstruct *[], int, char *);
extern int tabclose_6(ap_t *, char *);
extern int readtab_list_6(char *, objstruct **, int *, char *);
extern int get_ncol_6(void);
extern int tabinit_gen(char *, char *, int, char *[], char *[], char *[], 
		       char *);
extern int tabinit_genmef(char *, char *, int, char *[], char *[], char *[], 
			  char *);
extern int do_seeing_gen(ap_t *, int, int, int[], char *);
extern int tabclose_gen(ap_t *, char *);

/* AP routines */

extern void apclose(ap_t *);
extern void apfu(ap_t *);
extern void apinit(ap_t *);
extern void apreinit(ap_t *);
extern void apline(ap_t *, float [], float [], float [], float [], 
		   int, unsigned char *);
extern void apclust(ap_t *, int, plstruct *);
extern void moments (ap_t *, float []);
extern void areals(ap_t *, int [NAREAL]);
extern void restack(ap_t *, int);
extern void terminate(ap_t *, unsigned char *);
extern void extract_data(ap_t *, int);

/* Global variables */
 
int verbose;               /* Verbose flag */
fitsfile *tptr;            /* CFITSIO pointer for the output table file */
int cattype;               /* Flag for the type of table */
FILE *ellfp;               /* File descriptor for .ell file */
float gain;                /* Gain of the chip */
int dribble;               /* Set if map has been interpolated in some way */

/*

$Log: imcore.h,v $
Revision 1.2  2010/11/01 10:22:55  jim
*** empty log message ***

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.22  2010/02/11 21:55:33  jim
changed a few routine declarations

Revision 1.21  2010/02/10 11:54:43  jim
Gave dribble global variable an initial value

Revision 1.20  2010/01/19 12:01:47  jim
removed reference to colour stuff

Revision 1.19  2009/12/17 11:33:24  jim
Modified API for phopt. Also added dribble global variable

Revision 1.18  2009/09/24 13:49:52  jim
Modified for create_table_6 issues

Revision 1.17  2009/09/14 12:53:26  jim
Added imcore_opm_dataonly

Revision 1.16  2009/01/29 12:13:42  jim
removed imcore_version into a separate file

Revision 1.15  2009/01/29 10:49:40  jim
Added imcore_version global variable

Revision 1.14  2009/01/22 13:41:52  jim
Added to MF_ flag values

Revision 1.13  2008/04/15 18:59:58  jim
Redefined data flagging definitions. Added prototypes for colour routines

Revision 1.12  2007/07/31 12:02:55  jim
added imcore_opm

Revision 1.11  2007/06/04 10:34:02  jim
Modified to add list driven routines

Revision 1.10  2006/07/31 13:21:09  jim
Modified imcore now allows for a smoothing kernel with variable FWHM

Revision 1.9  2006/07/24 11:41:25  jim
Fixed some problems with the background estimation

Revision 1.8  2006/06/26 15:13:51  jim
Background stats routines improved to deal with small noise estimates

Revision 1.7  2006/06/22 15:07:05  jim
Fixed to flag pixels that very low values

Revision 1.6  2006/06/05 11:21:03  jim
Fixed apline so that it takes confidence maps into account better

Revision 1.5  2005/12/08 10:31:59  jim
modified definition of closefits

Revision 1.4  2005/11/25 15:37:23  jim
Redefined some of the MF_ constants so that they are not negative

Revision 1.3  2005/08/26 04:46:20  jim
Modified to add new radii and error estimates

Revision 1.2  2004/09/07 14:18:57  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:54:58  jim
New version for rewrite of imcore



*/
