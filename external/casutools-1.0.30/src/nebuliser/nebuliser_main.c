/* 

$Id: nebuliser_main.c,v 1.7 2012/05/09 09:45:52 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <fitsio.h>

#include <nebuliser.h>
#include <tools.h>

/* Command line arguments and defaults */

enum {BACKMAP_ARG,
      NITER_ARG,
      AXIS_ARG,
      TWOD_ARG,
      NOTWOD_ARG,
      TAKEOUTSKY_ARG,
      NOTAKEOUTSKY_ARG,
      INORM_ARG,
      NOINORM_ARG,
      SIGNEG_ARG,
      SIGPOS_ARG,
      VERB_ARG,
      NOVERB_ARG};

static struct option myoptions [] = {
    {"backmap",required_argument,NULL,BACKMAP_ARG},
    {"niter",required_argument,NULL,NITER_ARG},
    {"axis",required_argument,NULL,AXIS_ARG},
    {"twod",no_argument,NULL,TWOD_ARG},
    {"notwod",no_argument,NULL,NOTWOD_ARG},
    {"takeout_sky",no_argument,NULL,TAKEOUTSKY_ARG},
    {"notakeout_sky",no_argument,NULL,NOTAKEOUTSKY_ARG},
    {"inorm",no_argument,NULL,INORM_ARG},
    {"noinorm",no_argument,NULL,NOINORM_ARG},
    {"signeg",required_argument,NULL,SIGNEG_ARG},
    {"sigpos",required_argument,NULL,SIGPOS_ARG},
    {"verbose",no_argument,NULL,VERB_ARG},
    {"noverbose",no_argument,NULL,NOVERB_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: nebuliser infile confmap outfile medfiltsize linfiltsize\n"
    "[--backmap=%s] [--niter=%d] [--axis=%d] [--(no)twod (%s}]\n"
    "[--(no)takeoutsky (%s)] [--(no)inorm (%s)] [--signeg=%g] [--sigpos=%g]\n"
    "[--(no)verbose (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define BACKMAP_DEF ""
#define NITER_DEF 3
#define AXIS_DEF 2
#define TWOD_DEF 0
#define TAKEOUTSKY_DEF 0
#define INORM_DEF 0
#define SIGNEG_DEF 10.0
#define SIGPOS_DEF 3.0
#define VERB_DEF 0

static fitsfile *fptr = NULL;
static fitsfile *cptr = NULL;
static fitsfile *optr = NULL;
static void copy_file(char *inextn, fitsfile *optr);
static int istab(char *inextn, char *outfile, char *backmap);

static void tidy();

int main (int argc, char *argv[]) {
    char *infile,*confmap,*outfile,confextn[BUFSIZ],outextn[BUFSIZ];
    char backextn[BUFSIZ],errmsg[BUFSIZ],inextn[BUFSIZ];
    int medfiltsize,linfiltsize,c,status,nhdui,simple,nhduc,nerr,i,j,hdutype;
    int retval,option_index,yep=1,hdutype2,nim;
    long naxisi[2],naxisc[2];

    /* Set up defaults for command line switches */

    char backmap[BUFSIZ];
    (void)strcpy(backmap,BACKMAP_DEF);
    int niter = NITER_DEF;
    int axis = AXIS_DEF;
    int twod = TWOD_DEF;
    int takeout_sky = TAKEOUTSKY_DEF;
    int inorm = INORM_DEF;
    float signeg = SIGNEG_DEF;
    float sigpos = SIGPOS_DEF;
    int verbose = VERB_DEF;

    /* First get the command line arguments */

    if (argc < 6) {
	fprintf(stderr,USAGE,backmap,niter,axis,YESNO(twod),YESNO(takeout_sky),
		YESNO(inorm),signeg,sigpos,YESNO(verbose));
	exit(1);
    } else {
	infile = argv[1];
	confmap = argv[2];
	outfile = argv[3];
	medfiltsize = atoi(argv[4]);
	linfiltsize = atoi(argv[5]);
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case BACKMAP_ARG:
	    (void)strcpy(backmap,optarg);
	    break;
	case NITER_ARG:
	    niter = atoi(optarg);
	    break;
	case AXIS_ARG:
	    axis = atoi(optarg);
	    if (axis != 1 && axis != 2) {
		fprintf(stderr,"axis = 0 or 1 only\n");
		exit(1);
	    }
	    break;
	case TWOD_ARG:
	    twod = 1;
	    break;
	case NOTWOD_ARG:
	    twod = 0;
	    break;
	case TAKEOUTSKY_ARG:
	    takeout_sky = 1;
	    break;
	case NOTAKEOUTSKY_ARG:
	    takeout_sky = 0;
	    break;
	case INORM_ARG:
	    inorm = 1;
	    break;
	case NOINORM_ARG:
	    inorm = 0;
	    break;
	case SIGNEG_ARG:
	    signeg = (float)atof(optarg);
	    break;
	case SIGPOS_ARG:
	    sigpos = (float)atof(optarg);
	    break;
	case VERB_ARG:
	    verbose = 1;
	    break;
	case NOVERB_ARG:
	    verbose = 0;
	    break;
	default:
 	    fprintf(stderr,USAGE,BACKMAP_DEF,NITER_DEF,AXIS_DEF,
		    YESNO(TWOD_DEF),YESNO(TAKEOUTSKY_DEF),YESNO(INORM_DEF),
		    SIGNEG_DEF,SIGPOS_DEF,YESNO(VERB_DEF));
	    exit(1);
	}
    }

    /* Do some sanity checks. First make sure we can open the input files and
       see that they match up in terms of the number of extensions. */

    status = 0;
    (void)fits_open_file(&fptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(fptr,&nhdui,&status);
    nim = 0;
    for (i = 0; i < nhdui; i++) {
	(void)fits_movabs_hdu(fptr,i+1,&hdutype,&status);
	if (hdutype == IMAGE_HDU)
	    nim++;
    }
    simple = (nhdui == 1 ? 1 : 0);
    if (! simple)
	nhdui--;
    closefits(fptr);
    if (status != 0) {
	fprintf(stderr,"NEBULISER: Can't open file %s\n",infile);
	tidy();
	exit(1);
    }
    if (strcmp(confmap,"noconf") != 0) {
	(void)fits_open_file(&fptr,confmap,READONLY,&status);
	(void)fits_get_num_hdus(fptr,&nhduc,&status);
	closefits(fptr);
	if (! simple)
	    nhduc--;
	if (status != 0) {
	    fprintf(stderr,"NEBULISER: Can't open confidence map %s\n",confmap);
	    tidy();
	    exit(1);
	}
	if (nhdui != nhduc) {
	    fprintf(stderr,"NEBULISER: Image %s has %d extensions\n",infile,
		    nhdui);
	    fprintf(stderr,"           Conf %s has %d extensions\n",confmap,
		    nhduc);
	    fprintf(stderr,"           These must match\n");
	    tidy();
	    exit(1);

	}
    }

    /* Scroll through the image extensions and see if the image dimensions
       match up (so long as we're using a confidence map...). If the input
       file is a mixture of images and tables, then the confidence map 
       is expected to have the same structture. Tables can then just be
       skipped. */
 
    (void)fits_open_file(&fptr,infile,READONLY,&status);
    if (strcmp(confmap,"noconf") != 0) {
	nerr = 0;
	(void)fits_open_file(&cptr,confmap,READONLY,&status);
	for (i = 1; i <= nhdui; i++) {
	    j = (simple ? 1 : i + 1);
	    (void)fits_movabs_hdu(fptr,j,&hdutype,&status);
	    (void)fits_movabs_hdu(cptr,j,&hdutype2,&status);
	    if (hdutype != hdutype2) {
		fprintf(stderr,
			"NEBULISER: Image %s[%d] and Conf %s[%d] type mismatch\n");
		nerr++;
		continue;
	    }
	    if (hdutype != IMAGE_HDU)
		continue;
	    (void)fits_get_img_size(fptr,2,naxisi,&status);
	    (void)fits_get_img_size(cptr,2,naxisc,&status);
	    if (naxisi[0] != naxisc[0] || naxisi[1] != naxisc[1]) {
		fprintf(stderr,"NEBULISER: Image %s[%d] dims [%d:%d]\n",infile,
			j-1,naxisi[0],naxisi[1]);
		fprintf(stderr,"           Conf  %s[%d] dims [%d:%d]\n",confmap,
			j-1,naxisc[0],naxisc[1]);
		fprintf(stderr,"           Data array size mismatch\n");
		nerr++;
	    }
	}
	closefits(cptr);
	if (nerr != 0) {
	    tidy();
	    exit(1);
	}
    }
        
    /* If the output file exists, then delete it */

    if (access(outfile, F_OK) == 0)
	remove(outfile);
    closefits(fptr);

    /* Ok, let's get on with it. */

    if (verbose) 
	fprintf(stdout,"File: %s, %d image extensions will be filtered\n",
		infile,nim);
    for (i = 1; i <= nhdui; i++) {
	if (simple) {
	    (void)strcpy(inextn,infile);
	    (void)strcpy(confextn,confmap);
	    (void)strcpy(outextn,outfile);
	    if (strlen(backmap) > 0) 
		(void)strcpy(backextn,backmap);
	    else
		backextn[0] = '\0';
	} else {
	    if (strcmp(confmap,"noconf") != 0) 
		(void)sprintf(confextn,"%s[%d]",confmap,i);
	    else
		(void)strcpy(confextn,confmap);
	    (void)sprintf(inextn,"%s[%d]",infile,i);
	    (void)sprintf(outextn,"%s[%d]",outfile,i);
	    if (strlen(backmap) > 0) 
		(void)sprintf(backextn,"%s[%d]",backmap,i);
	    else
		backextn[0] = '\0';
	}
	if (istab(inextn,outfile,backmap)) 
	    continue;
	if (verbose) 
	    fprintf(stdout,"    filtering %s to %s\n",inextn,outextn);
	retval = cir_open_output(outextn,inextn,&optr,NEWBITPIX,FLOAT_IMG,0,
			      NULL,errmsg);
	if (retval != CIR_OK) {
	    fprintf(stderr,"NEBULISER: Error opening output image %s\n",
		    errmsg);
	    tidy();
	    exit(1);
	}
	copy_file(inextn,optr);
	closefits(optr);
	retval = cir_2dfilt(outextn,confextn,medfiltsize,linfiltsize,
			    backextn,niter,axis,twod,takeout_sky,inorm,
			    signeg,sigpos,errmsg);
	if (retval != CIR_OK) {
	    fprintf(stderr,"NEBULISER: Error %s\n",errmsg);
	    tidy();
	    exit(1);
	}
	status = 0;
	(void)fits_open_file(&fptr,outextn,READWRITE,&status);
	(void)fits_update_key(fptr,TLOGICAL,"NEBULSED",&yep,
			      "Nebuliser has been used on this image",&status);
	closefits(fptr);

	/* Stamp the primary */

	if (i == 1) {
	    status = 0;
	    (void)fits_open_file(&fptr,outfile,READWRITE,&status);
	    retval = casu_stamp(fptr,"nebuliser");
	    closefits(fptr);
	}
    }
    if (verbose) 
	fprintf(stdout,"Done\n");

    /* Get out of here */

    tidy();
    exit(0);
}

static void copy_file(char *inextn, fitsfile *optr) {
    int status,anynul;
    float *data;
    long naxes[2],npts;
    fitsfile *iptr;

    /* Get the data from the input file */

    status = 0;
    (void)fits_open_file(&iptr,inextn,READONLY,&status);
    (void)fits_get_img_size(iptr,2,naxes,&status);
    npts = naxes[0]*naxes[1];
    data = cir_malloc(npts*sizeof(float));
    (void)fits_read_img(iptr,TFLOAT,1,npts,NULL,data,&anynul,&status);
    (void)fits_close_file(iptr,&status);

    /* Now write it out */

    (void)fits_write_img(optr,TFLOAT,1,npts,data,&status);
    (void)fits_flush_file(optr,&status);
    freespace(data);
}


static int istab(char *inextn, char *outfile, char *backmap) {
    int status,nh,nhdu,hdutype;
    fitsfile *iptr,*optr;

    /* Open the extension and see if this is a table. If it isn't then
       close up and get out of here */

    status = 0;
    (void)fits_open_file(&iptr,inextn,READONLY,&status);
    (void)fits_get_hdu_type(iptr,&hdutype,&status);
    if (hdutype == IMAGE_HDU) {
	closefits(iptr);
	return(0);
    }

    /* If it's a table, then just copy the hdu over to the output file. If
       there is a background map, then do exactly the same. If the output
       file doesn't already exist, then the table is in the first extension
       and the primary is a dummy. Copy the dummy to the output file */

    if (access(outfile,F_OK) == 0) {
        (void)fits_open_file(&optr,outfile,READWRITE,&status);
	(void)fits_get_num_hdus(optr,&nhdu,&status);
	(void)fits_movabs_hdu(optr,nhdu,&hdutype,&status);
    } else {
	(void)fits_create_file(&optr,outfile,&status);
	(void)fits_get_hdu_num(iptr,&nh);
	(void)fits_movabs_hdu(iptr,1,&hdutype,&status);
	(void)fits_copy_hdu(iptr,optr,0,&status);
	(void)fits_movabs_hdu(iptr,nh,&hdutype,&status);
    }
    (void)fits_copy_hdu(iptr,optr,0,&status);
    (void)fits_close_file(optr,&status);
    if (access(backmap,F_OK) == 0) {
        (void)fits_open_file(&optr,backmap,READWRITE,&status);
	(void)fits_get_num_hdus(optr,&nhdu,&status);
	(void)fits_movabs_hdu(optr,nhdu,&hdutype,&status);
    } else {
	(void)fits_create_file(&optr,backmap,&status);
	(void)fits_get_hdu_num(iptr,&nh);
	(void)fits_movabs_hdu(iptr,1,&hdutype,&status);
	(void)fits_copy_hdu(iptr,optr,0,&status);
	(void)fits_movabs_hdu(iptr,nh,&hdutype,&status);
    }
    (void)fits_copy_hdu(iptr,optr,0,&status);
    (void)fits_close_file(optr,&status);

    /* Get out of here */

    (void)fits_close_file(iptr,&status);
    return(1);
}
	    

static void tidy() {
    int status;

    closefits(fptr);
    closefits(cptr);
    closefits(optr);
}

/*

$Log: nebuliser_main.c,v $
Revision 1.7  2012/05/09 09:45:52  jim
Modified so that if there is a table in the input file it skips over it

Revision 1.6  2011-06-09 11:38:58  jim
Adds keyword boolean NEBULSED to header

Revision 1.5  2010-11-02 10:38:32  jim
Added static routine to copy the input data across to the output file

Revision 1.4  2010-11-01 10:27:39  jim
Changed readwrite access to readonly for input files

Revision 1.3  2010-09-21 10:50:19  jim
Changed the way that the output file is created to avoid the problem of
an input compressed file. Also insist that output must be a floating point
image. Added verbose switch

Revision 1.2  2010-09-20 09:03:35  jim
Moved call to fits_open_file to avoid bug when noconf option is used

Revision 1.1  2010-09-06 09:03:50  jim
New entry


*/
