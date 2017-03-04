/*

$Id: wcsfit_main.c,v 1.6 2014/06/10 08:57:16 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <fitsio.h>

#include <wcsfit.h>
#include <tools.h>

enum {CATSRC_ARG,
      SITE_ARG,
      CATPATH_ARG};

static struct option myoptions [] = {
    {"catsrc",required_argument,NULL,CATSRC_ARG},
    {"site",required_argument,NULL,SITE_ARG},
    {"catpath",required_argument,NULL,CATPATH_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: wcsfit infile incat \n"
    "[--catsrc=%s] [--site=%s] [--catpath=%s]\n";

#define CATSRC_DEF "viz2mass"
#define SITE_DEF "casu"
#define CATPATH_DEF ""

int main (int argc, char *argv[]) {
    int retval,nhdu_im,nhdu_cat,status,c,option_index,hdutype,nerr,i,issimple;
    char errmsg[BUFSIZ],*inimage,*incat,imextn[BUFSIZ],catextn[BUFSIZ];
    fitsfile *fptr;

    /* Set up defaults for command line switches */

    char catsrc[16];
    (void)strcpy(catsrc,CATSRC_DEF);
    char site[16];
    (void)strcpy(site,SITE_DEF);
    char catpath[BUFSIZ];
    (void)strcpy(catpath,CATPATH_DEF);

    /* Get command line arguments */

    if (argc < 3) {
	fprintf(stderr,USAGE,catsrc,site,catpath);
	return(1);
    } else {
	inimage = argv[1];
	incat = argv[2];
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case CATSRC_ARG:
	    if (! strcmp(optarg,"viz2mass") || 
		! strcmp(optarg,"vizusnob") ||
		! strcmp(optarg,"vizppmxl") ||
		! strcmp(optarg,"local2mass") ||
		! strcmp(optarg,"localfits") ||
		! strcmp(optarg,"localusnob")){
  	        strcpy(catsrc,optarg);
	    } else {
		fprintf(stderr,"--catsrc = \"viz2mass\",\"vizusnob\",\"vizppmxl\",\"local2mass\" or \"localusnob\" only\n");
		fprintf(stderr,"You entered: %s\n",optarg);
		return(1);
	    }
	    break;
	case SITE_ARG:
	    if (! strcmp(optarg,"casu") || ! strcmp(optarg,"cds") || 
		! strcmp(optarg,"ukirt")) {
  	        strcpy(site,optarg);
	    } else {
		fprintf(stderr,"--site = \"casu\", \"cds\", or \"ukirt\" only\n");
		return(1);
	    }
	    break;
	case CATPATH_ARG:
  	    strcpy(catpath,optarg);
	    break;
	default:
	    fprintf(stderr,USAGE,CATSRC_DEF,SITE_DEF,CATPATH_DEF);
	    return(1);
	}
    }

    /* Open the image MEF. Get the number of extensions and make sure that they
       are all images */

    status = 0;
    (void)fits_open_file(&fptr,inimage,READWRITE,&status);
    if (status != 0) {
	fprintf(stderr,"WCSFIT: No readwrite access to %s\n",inimage);
	return(1);
    }
    (void)fits_get_num_hdus(fptr,&nhdu_im,&status);
    issimple = 0;
    if (nhdu_im == 1) 
	issimple = 1;
    else
	nhdu_im--;
    nerr = 0;
    for (i = 1; i <= nhdu_im; i++) {
	if (! issimple)
	    (void)fits_movabs_hdu(fptr,i+1,&hdutype,&status);
	else
	    (void)fits_movabs_hdu(fptr,i,&hdutype,&status);
	if (hdutype != IMAGE_HDU) {
	    fprintf(stderr,"%s[%d] is not a image extension\n",inimage,i);
	    nerr++;
	}
    }
    (void)fits_close_file(fptr,&status);
    if (nerr != 0) 
	return(1);

    /* Now do the same for the catalogue. Only make sure these are all 
       binary tables */

    status = 0;
    (void)fits_open_file(&fptr,incat,READWRITE,&status);
    if (status != 0) {
	fprintf(stderr,"WCSFIT: No readwrite access to %s\n",incat);
	return(1);
    }
    (void)fits_get_num_hdus(fptr,&nhdu_cat,&status);
    nhdu_cat--;
    nerr = 0;
    for (i = 1; i <= nhdu_cat; i++) {
	(void)fits_movabs_hdu(fptr,i+1,&hdutype,&status);
	if (hdutype != BINARY_TBL) {
	    fprintf(stderr,"%s[%d] is not a binary table extension\n",incat,i);
	    nerr++;
	}
    }
    (void)fits_close_file(fptr,&status);
    if (nerr != 0) 
	return(1);

    /* Finally compare the number of extensions to make sure that the
       catalogue and the image at least match up */

    if (nhdu_im != nhdu_cat) {
        fprintf(stderr,"WCSFIT: %s has %d extensions and %s has %d extensions\n",
		inimage,nhdu_im,incat,nhdu_cat);
	fprintf(stderr,"        These must match!\n");
	return(1);
    }

    /* Ok, now loop for each extension and do the job */

    
    for (i = 1; i <= nhdu_im; i++) {
	if (issimple) 
	    (void)sprintf(imextn,"%s",inimage);
	else 
	    (void)sprintf(imextn,"%s[%d]",inimage,i);
	(void)sprintf(catextn,"%s[%d]",incat,i);
	retval = wcsfit(imextn,catextn,catsrc,site,catpath,errmsg);
	if (retval == CIR_FATAL) {
	    fprintf(stderr,"%s",errmsg);
	    return(1);
	} else if (retval == CIR_WARN) {
            fprintf(stderr,"%s",errmsg);
        }
    }

    /* Get out of here */

    return(0);
}
	    
/*

$Log: wcsfit_main.c,v $
Revision 1.6  2014/06/10 08:57:16  jim
Modified so that if the input catalogue has no rows you get a non fatal
error

Revision 1.5  2012/08/13 10:48:12  jim
Added PPMXL to list of astrometric catalogues

Revision 1.4  2011-01-19 12:47:09  jim
Fixed options menu

Revision 1.3  2010-09-06 09:05:59  jim
Tidied some docs


*/
