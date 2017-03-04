/* 

$Id: imcore_main.c,v 1.5 2012/08/13 09:53:11 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

#include "errcodes.h"
#include "imcore.h"
#include <tools.h>

/* Command line arguments and defaults */

enum {CROWD_ARG,
      NOCROWD_ARG,
      RCORE_ARG,
      NBSIZE_ARG,
      FILTFWHM_ARG,
      ELL_ARG,
      NOELL_ARG,
      VERBOSE_ARG,
      NOVERBOSE_ARG,
      CATTYPE_ARG};

static struct option myoptions [] = {
    {"crowd",no_argument,NULL,CROWD_ARG},
    {"nocrowd",no_argument,NULL,NOCROWD_ARG},
    {"rcore",required_argument,NULL,RCORE_ARG},
    {"nbsize",required_argument,NULL,NBSIZE_ARG},
    {"filtfwhm",required_argument,NULL,FILTFWHM_ARG},
    {"ell",no_argument,NULL,ELL_ARG},
    {"noell",no_argument,NULL,NOELL_ARG},
    {"verbose",no_argument,NULL,VERBOSE_ARG},
    {"noverbose",no_argument,NULL,NOVERBOSE_ARG},
    {"cattype",required_argument,NULL,CATTYPE_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: imcore infile confmap outfile ipix thresh\n"
    "[--(no)crowd (%s)] [--rcore=%g] [--nbsize=%d] [--filtfwhm=%g]\n"
    "[--(no)ell (%s)] [--(no)verbose (%s)] [--cattype=%d]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define CROWD_DEF 1
#define RCORE_DEF 5.0
#define NBSIZE_DEF 64
#define FILTFWHM_DEF 3.0
#define ELL_DEF 1
#define VERBOSE_DEF 0
#define CATTYPE_DEF 6

extern int cir_classify(char *infile, char *expkey, float minsize, 
			char *errmsg);    

int main (int argc, char *argv[]) {
    char errmsg[BUFSIZ],*infile,*conf,*outfile,*p,inf_ext[BUFSIZ];
    char conf_ext[BUFSIZ],ell_ext[BUFSIZ],out_ext[BUFSIZ];
    int retval,ipix,c,option_index,status,i,ist,ifn,nhdu,first;
    float thresh;
    fitsfile *fptr;

    /* Set up defaults for command line switches */

    int crowd = CROWD_DEF;
    float rcore = RCORE_DEF;
    int nbsize = NBSIZE_DEF;
    float filtfwhm = FILTFWHM_DEF;
    int ell = ELL_DEF;
    int verb = VERBOSE_DEF;
    int cattyp = CATTYPE_DEF;

    /* Read command line */

    if (argc < 6) {
	fprintf(stderr,USAGE,YESNO(CROWD_DEF),RCORE_DEF,NBSIZE_DEF,
		FILTFWHM_DEF,YESNO(ELL_DEF),YESNO(VERBOSE_DEF),CATTYPE_DEF);
	exit(1);
    } else {
	infile = argv[1];
	conf = argv[2];
	outfile = argv[3];
	ipix = atoi(argv[4]);
	thresh = (float)atof(argv[5]);
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case CROWD_ARG:
	    crowd = 1;
	    break;
	case NOCROWD_ARG:
	    crowd = 0;
	    break;
	case RCORE_ARG:
	    rcore = (float)atof(optarg);
	    break;
	case NBSIZE_ARG:
	    nbsize = atoi(optarg);
	    break;
	case FILTFWHM_ARG:
	    filtfwhm = (float)atof(optarg);
	    break;
	case ELL_ARG:
	    ell = 1;
	    break;
	case NOELL_ARG:
	    ell = 0;
	    break;
	case VERBOSE_ARG:
	    verb = 1;
	    break;
	case NOVERBOSE_ARG:
	    verb = 0;
	    break;
	case CATTYPE_ARG:
	    cattyp = atoi(optarg);
	    break;
	default:
	    fprintf(stderr,USAGE,YESNO(CROWD_DEF),RCORE_DEF,NBSIZE_DEF,
		    FILTFWHM_DEF,YESNO(ELL_DEF),YESNO(VERBOSE_DEF),CATTYPE_DEF);
	    exit(1);
	}
    }

    /* Get rid of the output file if it already exists */

    if (access(outfile,F_OK) == 0)
	remove(outfile);

    /* Find out how many image extensions there are */

    status = 0;
    (void)fits_open_file(&fptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(fptr,&nhdu,&status);
    (void)fits_close_file(fptr,&status);
    if (status != 0) {
        fprintf(stderr,"Unable to open input file %s\n",infile);
        exit(1);
    }
    if (nhdu == 1) {
	ist = 0;
        ifn = 0;
    } else {
        ist = 1;
        ifn = nhdu - 1;
    }
    
    /* Loop for all extensions */

    first = 1;
    for (i = ist; i <= ifn; i++) {
        (void)sprintf(inf_ext,"%s[%d]",infile,i);
	if (strcmp(conf,"noconf")) 
       	    (void)sprintf(conf_ext,"%s[%d]",conf,i);
        else
	    (void)strcpy(conf_ext,"noconf");
	if (ell) {
	    (void)strcpy(ell_ext,outfile);
	    p = strrchr(ell_ext,'.');
	    if (p != NULL) {
		if (i == 0) 
		    (void)strcpy(p,".ell");
		else
		    (void)sprintf(p,"_%d.ell",i);    
	    } else {
		if (i == 0) 
		    (void)strcat(ell_ext,".ell");
		else
		    (void)sprintf(ell_ext,"%s_%d.ell",outfile,i);
	    }
	} else {
	    ell_ext[0] = '\0';
	}
        if (access(ell_ext,F_OK) == 0)
	    remove(ell_ext);

	/* Run imcore */

        retval = imcore_conf(inf_ext,conf_ext,ipix,thresh,crowd,rcore,nbsize,
			     filtfwhm,outfile,ell_ext,verb,cattyp,errmsg);
        if (retval != ERRCODE_OK)
  	    fprintf(stderr,"Error running imcore_conf\n%s\n",errmsg);

	/* Now run the classifier */

	if (i == 0)
	    (void)sprintf(out_ext,"%s[%d]",outfile,1);
	else
	    (void)sprintf(out_ext,"%s[%d]",outfile,i);
	retval = cir_classify(out_ext,"",16,errmsg);
        if (retval != 0)
  	    fprintf(stderr,"Error running classify\n%s\n",errmsg);

	/* Now stamp the primary */

	if (first) {
	    status = 0;
	    (void)fits_open_file(&fptr,outfile,READWRITE,&status);
	    retval = casu_stamp(fptr,"imcore");
	    if (retval != 0) 
		fprintf(stderr,"Error stamping output table primary %s\n",
		        outfile);
	    (void)fits_close_file(fptr,&status);
	    first = 0;
	}
    }
    return(0);
}

/*

$Log: imcore_main.c,v $
Revision 1.5  2012/08/13 09:53:11  jim
Fixed so that it doesn't barf if the output image doesn't have a file
name extension

Revision 1.4  2011-01-19 12:46:32  jim
Fixed options menu in header

Revision 1.3  2010-09-06 08:59:02  jim
Tidied up some documentation


*/
