#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>

#include "errcodes.h"
#include "imcore.h"

/* Command line arguments and defaults */

enum {RCORE_ARG,
      NBSIZE_ARG,
      ELL_ARG,
      NOELL_ARG,
      VERBOSE_ARG,
      NOVERBOSE_ARG,
      TRANS_ARG,
      CATTYPE_ARG};

static struct option myoptions [] = {
    {"rcore",required_argument,NULL,RCORE_ARG},
    {"nbsize",required_argument,NULL,NBSIZE_ARG},
    {"ell",no_argument,NULL,ELL_ARG},
    {"noell",no_argument,NULL,NOELL_ARG},
    {"verbose",no_argument,NULL,VERBOSE_ARG},
    {"noverbose",no_argument,NULL,NOVERBOSE_ARG},
    {"trans",required_argument,NULL,TRANS_ARG},
    {"cattype",required_argument,NULL,CATTYPE_ARG}};

static char *USAGE = "Usage: imcore_list infile confmap listfile outfile thresh\n"
    "[--rcore=%g] [--nbsize=%d] [--trans=%s] [--(no)ell (%s)]\n"
    "[--(no)verbose (%s)] [--cattype=%d]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define RCORE_DEF 5.0
#define NBSIZE_DEF 64
#define ELL_DEF 1
#define VERBOSE_DEF 0
#define TRANS_DEF ""
#define CATTYPE_DEF 1
    
int main (int argc, char *argv[]) {
    char errmsg[BUFSIZ],*infile,*conf,*outfile,*p,inf_ext[BUFSIZ];
    char conf_ext[BUFSIZ],ell_ext[BUFSIZ],*listfile,list_ext[BUFSIZ];
    int retval,c,option_index,status,i,ist,ifn,nhdu,isfits;
    float thresh;
    fitsfile *fptr;

    /* Set up defaults for command line switches */

    float rcore = RCORE_DEF;
    int nbsize = NBSIZE_DEF;
    int ell = ELL_DEF;
    int verb = VERBOSE_DEF;
    int cattyp = CATTYPE_DEF;
    char trans[BUFSIZ];
    (void)strcpy(trans,TRANS_DEF);

    /* Read command line */

    if (argc < 6) {
	fprintf(stderr,USAGE,RCORE_DEF,NBSIZE_DEF,TRANS_DEF,
		YESNO(ELL_DEF),YESNO(VERBOSE_DEF),CATTYPE_DEF);
	exit(1);
    } else {
	infile = argv[1];
	conf = argv[2];
	listfile = argv[3];
	outfile = argv[4];
	thresh = (float)atof(argv[5]);
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case RCORE_ARG:
	    rcore = (float)atof(optarg);
	    break;
	case NBSIZE_ARG:
	    nbsize = atoi(optarg);
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
	case TRANS_ARG:
	    (void)strcpy(trans,optarg);
	    break;
	default:
	    fprintf(stderr,USAGE,RCORE_DEF,NBSIZE_DEF,TRANS_DEF,
		    YESNO(ELL_DEF),YESNO(VERBOSE_DEF),CATTYPE_DEF);
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
    
    /* Is the list file a FITS table? */

    status = 0;
    (void)fits_open_file(&fptr,listfile,READONLY,&status);
    isfits = (status == 0);
    (void)fits_close_file(fptr,&status);
      
    /* Loop for all extensions */

    for (i = ist; i <= ifn; i++) {
        (void)sprintf(inf_ext,"%s[%d]",infile,i);
	if (strcmp(conf,"noconf"))
 	    (void)sprintf(conf_ext,"%s[%d]",conf,i);
	else
	    (void)strcpy(conf_ext,"noconf");
	if (isfits) {
	    if (i == 0) 
		(void)sprintf(list_ext,"%s[1]",listfile);
            else 
		(void)sprintf(list_ext,"%s[%d]",listfile,i);
        } else {
	    (void)strcpy(list_ext,listfile);
        }
	if (ell) {
	    (void)strcpy(ell_ext,outfile);
	    p = strrchr(ell_ext,'.');
	    if (i == 0) 
		(void)strcpy(p,".ell");
	    else
   	        (void)sprintf(p,"_%d.ell",i);    
	} else {
	    ell_ext[0] = '\0';
	}
        if (access(ell_ext,F_OK) == 0)
	    remove(ell_ext);

	/* Run imcore */

        retval = imcore_list(inf_ext,conf_ext,list_ext,trans,thresh,rcore,
			     nbsize,outfile,ell_ext,verb,cattyp,errmsg);
        if (retval != ERRCODE_OK)
  	    fprintf(stderr,"Error running imcore_list\n%s\n",errmsg);
    }
    return(0);
}



