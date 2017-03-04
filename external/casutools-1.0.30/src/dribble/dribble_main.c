/*

$Id: dribble_main.c,v 1.1 2010/11/22 10:43:47 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <fitsio.h>

#include <tools.h>

/* Command line arguments and defaults */

enum {ISIZE_ARG,
      CONS3_ARG,
      NOCONS3_ARG,
      VERB_ARG,
      NOVERB_ARG};

static struct option myoptions [] = {
    {"isize",required_argument,NULL,ISIZE_ARG},
    {"cons3",no_argument,NULL,CONS3_ARG},
    {"nocons3",no_argument,NULL,NOCONS3_ARG},
    {"verbose",no_argument,NULL,VERB_ARG},
    {"noverbose",no_argument,NULL,NOVERB_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: dribble infile confmap outfile outconf\n"
    "[--isize=%d] [--(no)cons3 (%s)] [--(no)verbose (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define ISIZE_DEF 2048
#define CONS3_DEF 1
#define VERB_DEF 0

static int dribble(char *inextn, char *confextn, char *outextn, char *outcextn,
		   int iconf, int isize, int cons3, char *errmsg);

int main (int argc, char *argv[]) {
    char *infile,*confmap,*outfile,*outconf,inextn[BUFSIZ],confextn[BUFSIZ];
    char outextn[BUFSIZ],outcextn[BUFSIZ],errmsg[BUFSIZ];
    int status,nhdui,simple,nhduc,nerr,i,j,hdutype,option_index,retval,iconf,c;
    long naxisi[2],naxisc[2];
    fitsfile *iptr=NULL,*cptr=NULL;

    /* Set up defaults for command line switches */

    int isize = ISIZE_DEF;
    int cons3 = CONS3_DEF;
    int verbose = VERB_DEF;

    /* First get the command line arguments */

    if (argc < 5) {
	fprintf(stderr,USAGE,isize,YESNO(cons3),YESNO(verbose));
	exit(1);
    } else {
	infile = argv[1];
	confmap = argv[2];
	outfile = argv[3];
	outconf = argv[4];
    }

    /* Are we using a confidence map? */

    iconf = strcmp(confmap,"noconf");

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case ISIZE_ARG:
	    isize = atoi(optarg);
	    break;
	case CONS3_ARG:
	    cons3 = 1;
	    break;
	case NOCONS3_ARG:
	    cons3 = 0;
	    break;
	case VERB_ARG:
	    verbose = 1;
	    break;
	case NOVERB_ARG:
	    verbose = 0;
	    break;
	default:
	    fprintf(stderr,USAGE,ISIZE_DEF,YESNO(CONS3_DEF),YESNO(VERB_DEF));
	    exit(1);
	}
    }

    /* Do some sanity checks. First make sure we can open the input files and
       see that they match up in terms of the number of extensions. */

    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(iptr,&nhdui,&status);
    simple = (nhdui == 1 ? 1 : 0);
    if (! simple)
	nhdui--;
    closefits(iptr);
    if (status != 0) {
	fprintf(stderr,"DRIBBLE: Can't open file %s\n",infile);
	exit(1);
    }
    if (iconf) {
	(void)fits_open_file(&iptr,confmap,READONLY,&status);
	(void)fits_get_num_hdus(iptr,&nhduc,&status);
	closefits(iptr);
	if (! simple)
	    nhduc--;
	if (status != 0) {
	    fprintf(stderr,"DRIBBLE: Can't open confidence map %s\n",confmap);
	    exit(1);
	}
	if (nhdui != nhduc) {
	    fprintf(stderr,"DRIBBLE: Image %s has %d extensions\n",infile,
		    nhdui);
	    fprintf(stderr,"         Conf  %s has %d extensions\n",confmap,
		    nhduc);
	    fprintf(stderr,"         These must match\n");
	    exit(1);

	}
    }

    /* Scroll through the image extensions and see if the image dimensions
       match up */

    if (iconf) {
	(void)fits_open_file(&iptr,infile,READONLY,&status);
	(void)fits_open_file(&cptr,confmap,READONLY,&status);
	nerr = 0;
	for (i = 1; i <= nhdui; i++) {
	    j = (simple ? 1 : i + 1);
	    (void)fits_movabs_hdu(iptr,j,&hdutype,&status);
	    if (hdutype != IMAGE_HDU) {
		fprintf(stderr,"DRIBBLE: Image %s[%d] is not a FITS image\n",
			infile,i);
		nerr++;
		continue;
	    }
	    (void)fits_movabs_hdu(cptr,j,&hdutype,&status);
	    if (hdutype != IMAGE_HDU) {
		fprintf(stderr,"DRIBBLE: Conf %s[%d] is not a FITS image\n",
			infile,i);
		nerr++;
		continue;
	    }
	    (void)fits_get_img_size(iptr,2,naxisi,&status);
	    (void)fits_get_img_size(cptr,2,naxisc,&status);
	    if (naxisi[0] != naxisc[0] || naxisi[1] != naxisc[1]) {
		fprintf(stderr,"DRIBBLE: Image %s[%d] dims [%d:%d]\n",infile,
			j-1,naxisi[0],naxisi[1]);
		fprintf(stderr,"         Conf   %s[%d] dims [%d:%d]\n",confmap,
			j-1,naxisc[0],naxisc[1]);
		fprintf(stderr,"         Data array size mismatch\n");
		nerr++;
	    }
	} 
        closefits(iptr);
        closefits(cptr);
    }
    if (nerr != 0) 
	exit(1);

    /* If the output file exists, then delete it */

    if (access(outfile, F_OK) == 0)
	remove(outfile);

    /* Friendly message */

    if (verbose) {
	if (iconf) 
 	    fprintf(stdout,
	 	    "Files: %s,%s, %d image extensions\n-->  : %s,%s\n",
 		    infile,confmap,nhdui,outfile,outconf);
	else
 	    fprintf(stdout,
	 	    "File: %s, %d image extensions\n--> : %s\n",
 		    infile,nhdui,outfile);
    }
	    
    /* Ok, let's get on with it. Loop for each extension and form the extended
       input/output file names */

    if (! iconf) {
	confextn[0] = '\0';
	outcextn[0] = '\0';
    }
    for (i = 1; i <= nhdui; i++) {
	if (simple) {
	    (void)strcpy(inextn,infile);
	    (void)strcpy(outextn,outfile);
	    if (iconf) {
		(void)strcpy(confextn,confmap);
		(void)strcpy(outcextn,outconf);
	    }
	} else {
	    (void)sprintf(inextn,"%s[%d]",infile,i);
	    (void)sprintf(outextn,"%s[%d]",outfile,i);
	    if (iconf) {
		(void)sprintf(confextn,"%s[%d]",confmap,i);
		(void)sprintf(outcextn,"%s[%d]",outconf,i);
	    }
	}
	if (verbose) 
	    fprintf(stdout,"    dribbling %s to %s\n",inextn,outextn);

	/* Do the dribble */

	retval = dribble(inextn,confextn,outextn,outcextn,iconf,isize,
			 cons3,errmsg);
	if (retval != CIR_OK) {
	    fprintf(stderr,"DRIBBLE: Error %s\n",errmsg);
	    exit(1);
	}

	/* Stamp the primary */

	if (i == 1) {
	    status = 0;
	    (void)fits_open_file(&iptr,outfile,READWRITE,&status);
	    retval = casu_stamp(iptr,"dribble");
	    closefits(iptr);
	    (void)fits_open_file(&iptr,outconf,READWRITE,&status);
	    retval = casu_stamp(iptr,"dribble");
	    closefits(iptr);
	}
    }
    if (verbose) 
	fprintf(stdout,"Done\n");

    /* Get out of here */

    exit(0);
}

static int dribble(char *inextn, char *confextn, char *outextn, char *outcextn,
		   int iconf, int isize, int cons3, char *errmsg) {
    int status,anynul,retval,inter,i,j,iarg,j1,jj,i1,ii,id;
    long naxis[2],npts;
    short int *sdata=NULL,*cdata,*ocdata,ic;
    char errstr[BUFSIZ];
    float *fdata=NULL,*idata,*odata,xx;
    fitsfile *iptr=NULL,*cptr=NULL,*optr=NULL,*ocptr=NULL;

    /* First open the image and confidence map extensions. (We've already 
       tested that we can do this and that the images match) */

    status = 0;
    (void)fits_open_file(&iptr,inextn,READONLY,&status);
    if (iconf) 
	(void)fits_open_file(&cptr,confextn,READONLY,&status);
    (void)fits_get_img_size(iptr,2,naxis,&status);

    /* Create the output files */

    retval = cir_open_output(outextn,inextn,&optr,NEWBITPIX,FLOAT_IMG,0,
			     NULL,errstr);
    if (retval != CIR_OK) {
	(void)sprintf(errmsg,"Error opening output image %s\n",errstr);
	closefits(iptr);
	closefits(cptr);
	return(CIR_FATAL);
    }
    if (iconf) {
	retval = cir_open_output(outcextn,confextn,&ocptr,DIRECTCOPY,SHORT_IMG,
				 0,NULL,errstr);
	if (retval != CIR_OK) {
	    (void)sprintf(errmsg,"Error opening output confmap %s\n",errstr);
	    closefits(iptr);
	    closefits(cptr);
	    closefits(optr);
	    return(CIR_FATAL);
	}
    }

    /* What was the original size of the interleaving? */

    inter = (int)((float)naxis[0]/(float)isize + 0.5);

    /* Get some workspace for the data arrays */

    npts = naxis[0]*naxis[1];
    fdata = cir_calloc(2*npts,sizeof(*fdata));
    idata = fdata;
    odata = fdata + npts;
    sdata = cir_calloc(2*npts,sizeof(*sdata));
    cdata = sdata;
    ocdata = sdata + npts;

    /* Read the input data and confidence. Check the confidence for
       silly values */

    (void)fits_read_img(iptr,TFLOAT,1,npts,NULL,idata,&anynul,&status);
    if (iconf) {
	(void)fits_read_img(cptr,TSHORT,1,npts,NULL,cdata,&anynul,&status);
	for (i = 0; i < npts; i++)
	    cdata[i] = max(0,min(110,cdata[i]));
    } else {
	for (i = 0; i < npts; i++)
	    cdata[i] = 100;
    }

    /* If this hasn't really been interleaved, then just copy it over
       and be done with it */
    
    if (inter == 1) {
	memmove(odata,idata,npts*sizeof(*idata));
	if (iconf)
	    memmove(ocdata,cdata,npts*sizeof(*cdata));
    } else {

	/* Right, now go through the grid and do the dribbling */

	for (j = 0; j < naxis[1]; j++) {
	    for (i = 0; i < naxis[0]; i++) {
		iarg = j*naxis[0] + i;
		xx = idata[iarg];
		ic = cdata[iarg];
		for (j1 = -1; j1 <= 1; j1++) {
		    jj = j + j1;
		    for (i1 = -1; i1 <= 1; i1++) {
			ii = i + i1;
			if (ii > -1 && ii < naxis[0] && jj > -1 && 
			    jj < naxis[1]) {
			    iarg = jj*naxis[0] + ii;
			    if (inter == 3 && cons3) {
				odata[iarg] += xx*(float)ic;
				ocdata[iarg] += ic;
			    } else {
				id = i1*i1 + j1*j1;
				switch (id) {
				case 0:
				    odata[iarg] += 4.0*xx;
				    ocdata[iarg] += 4*ic;
				    break;
				case 1:
				    odata[iarg] += 2.0*xx;
				    ocdata[iarg] += 2*ic;
				    break;
				case 2:
				    odata[iarg] += xx;
				    ocdata[iarg] += ic;
				}
			    }
			}
		    }
		}
	    }
	}

	/* Now renormalise */

	for (i = 0; i < npts; i++) {
	    if (inter == 3 && cons3) {
		odata[i] /= (float)ocdata[i];
		ocdata[i] = (ocdata[i] + 4)/9;
	    } else {
		odata[i] /= 16.0;
		ocdata[i] = (ocdata[i] + 8)/16.0;
	    }
	}
    }

    /* Write the data out */

    (void)fits_write_img(optr,TFLOAT,1,npts,odata,&status);
    (void)sprintf(errstr,"%dx%d",inter,inter);
    (void)fits_update_key(optr,TSTRING,"DRIBBLE",errstr,
			  "Data have been dribbled",&status);
    if (iconf) {
	(void)fits_write_img(ocptr,TSHORT,1,npts,ocdata,&status);
	(void)fits_update_key(optr,TSTRING,"DRIBBLE",errstr,
			      "Data have been dribbled",&status);
    }
    
    /* Tidy and exit */

    closefits(iptr);
    closefits(optr);
    if (iconf) {
	closefits(cptr);
	closefits(ocptr);
    }
    freespace(fdata);
    freespace(sdata);
    return(CIR_OK);
}

/* 

$Log: dribble_main.c,v $
Revision 1.1  2010/11/22 10:43:47  jim
new entry


*/
