/*

$Id: mosaic_main.c,v 1.3 2010/09/21 10:51:36 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <fitsio.h>

#include "mosaic.h"
#include <tools.h>

/* Command line arguments and defaults */

enum {INTERP_ARG,
      SKYFLAG_ARG,
      SKYWISH_ARG,
      EXPKEY_ARG,
      CONFLIM_ARG,
      VERB_ARG,
      NOVERB_ARG};

static struct option myoptions [] = {
    {"interp",required_argument,NULL,INTERP_ARG},
    {"skyflag",required_argument,NULL,SKYFLAG_ARG},
    {"skywish",required_argument,NULL,SKYWISH_ARG},
    {"expkey",required_argument,NULL,EXPKEY_ARG},
    {"conflim",required_argument,NULL,CONFLIM_ARG},
    {"verbose",no_argument,NULL,VERB_ARG},
    {"noverbose",no_argument,NULL,NOVERB_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: mosaic infiles confmaps outfile outconf\n"
    "[--interp=%d] [--skyflag=%d] [--skywish=%g] [--expkey=%s]\n"
    "[--conflim=%d] [--(no)verbose (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define INTERP_DEF 1
#define SKYFLAG_DEF 2
#define SKYWISH_DEF 0.0
#define EXPKEY_DEF "EXPTIME"
#define CONFLIM_DEF 25
#define VERB_DEF 0

static int getlist(char *inlist, char ***files, int *nfiles);
static void tidy();

static char **infiles = NULL;
static int ninfiles = 0;
static char **inconfs = NULL;
static int ninconfs = 0;

int main (int argc, char *argv[]) {
    int c,option_index,retval,status;
    char errmsg[BUFSIZ],*inlist,*inconflist,*outfile,*outconf;
    fitsfile *fptr;

    /* Set up defaults for optional command line arguments */

    int interp = INTERP_DEF;
    int skyflag = SKYFLAG_DEF;
    float skywish = SKYWISH_DEF;
    char expkey[81];
    (void)strcpy(expkey,EXPKEY_DEF);
    int conflim = CONFLIM_DEF;
    int verbose = VERB_DEF;

    /* Set up required command line arguments */

    if (argc < 5) {
	fprintf(stderr,USAGE,interp,skyflag,skywish,expkey,conflim,
		YESNO(verbose));
	exit(1);
    } else {
	inlist = argv[1];
	inconflist = argv[2];
	outfile = argv[3];
	outconf = argv[4];
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case INTERP_ARG:
	    interp = atoi(optarg);
	    if (interp < 0 || interp > 1) {
		fprintf(stderr,"MOSAIC: --interp only takes values 0,1\n");
		exit(1);
	    }
	    break;
	case SKYFLAG_ARG:
	    skyflag = atoi(optarg);
	    if (skyflag < 0 || skyflag > 2) {
		fprintf(stderr,"MOSAIC: --skyflag only takes values 0,1,2\n");
		exit(1);
	    }
	    break;
	case SKYWISH_ARG:
	    skywish = (float)atof(optarg);
	    break;
	case EXPKEY_ARG:
	    (void)strcpy(expkey,optarg);
	    break;
	case CONFLIM_ARG:
	    conflim = atoi(optarg);
	    break;
	case VERB_ARG:
	    verbose = 1;
	    break;
	case NOVERB_ARG:
	    verbose = 0;
	    break;
	default:
	    fprintf(stderr,USAGE,INTERP_DEF,SKYFLAG_DEF,SKYWISH_ARG,EXPKEY_ARG,
		    CONFLIM_ARG,YESNO(VERB_DEF));
	    exit(1);
	}
    }

    /* Get rid of the output files if they already exists */

    if (access(outfile,F_OK) == 0)
	remove(outfile);
    if (access(outconf,F_OK) == 0)
	remove(outconf);

    /* Get the input file list */

    if (getlist(inlist,&infiles,&ninfiles) != CIR_OK) {
	tidy();
	exit(1);
    }
    if (getlist(inconflist,&inconfs,&ninconfs) != CIR_OK) {
	tidy();
	exit(1);
    }

    /* Ok, do a quick sanity check here. There must be some input images
       (obviously!) There can be either one confidence map for all the input 
       files or a confidence map for each. The main routine cir_mosaic
       does a lot of checking about the sizes of the images match with
       those of the confidence maps, so we don't need to repeat these here */

    if (ninfiles == 0) {
	fprintf(stderr,"MOSAIC: No input files\n");
	tidy();
	exit(1);
    }
    if (ninconfs != 1 && ninconfs != ninfiles) {
	fprintf(stderr,"MOSAIC: Nimages = %d, Nconfs = %d\n",ninfiles,ninconfs);
	fprintf(stderr,"        Nconfs must be 1 or equal to Nimages\n");
	tidy();
	exit(1);
    }

    /* Do the work now...*/

    if (verbose) {
	setvbuf(stdout,NULL,_IONBF,0);
	fprintf(stdout,"%d files will be mosaicked to %s\n",ninfiles,outfile);
    }
    retval = cir_mosaic(infiles,inconfs,ninfiles,ninconfs,interp,skyflag,
			skywish,expkey,conflim,outfile,outconf,verbose,errmsg);
    if (retval != CIR_OK) {
	fprintf(stderr,"MOSAIC: Error %s\n",errmsg);
	tidy();
	exit(1);
    }

    /* Stamp the primaries of the output files */

    status = 0;
    (void)fits_open_file(&fptr,outfile,READWRITE,&status);
    retval = casu_stamp(fptr,"mosaic");
    (void)fits_close_file(fptr,&status);
    (void)fits_open_file(&fptr,outconf,READWRITE,&status);
    retval = casu_stamp(fptr,"mosaic");
    (void)fits_close_file(fptr,&status);
    if (verbose)
	fprintf(stdout,"Done\n");

    /* Get out of here */

    tidy();
    return(0);
}

static int getlist(char *inlist, char ***infiles, int *ninfiles) {
    FILE *fd;
    char fname[BUFSIZ],*v;
    int i,nerr,n;

    /* If there is no input string then send back some zero values */

    if (strlen(inlist) == 0) {
	*ninfiles = 0;
	*infiles = NULL;
	return(CIR_OK);
    }
    
    /* Check to see if the first character of the input list is an
       ampersand. If it is, then this is a list file. Check that this
       file exists and if it does try and open it */

    if (! strncmp(inlist,"@",1)) {
        if (access(inlist+1,R_OK) == 0) {
	    fd = fopen(inlist+1,"r");
	    if (fd == NULL) {
		fprintf(stderr,"MOSAIC: Unable to open list file %s\n",inlist+1);
		return(CIR_FATAL);
	    }
	} else {
	    fprintf(stderr,"MOSAIC: Unable to find list file %s\n",inlist+1);
	    return(CIR_FATAL);
	}
	
	/* Run through the list first to see how many files there are */

	*ninfiles = 0;
	while (fscanf(fd,"%s",fname) != EOF) 
	    (*ninfiles)++;
	rewind(fd);
	/* Allocate space for file names */

	*infiles = cir_malloc(*ninfiles*sizeof(char *));
	for (i = 0; i < *ninfiles; i++) 
	    (*infiles)[i] = cir_malloc(BUFSIZ);
	
	/* Now read the file */

	i = 0;
	nerr = 0;
	while (fscanf(fd,"%s",fname) != EOF) {
	    strcpy((*infiles)[i],fname);
	    if (access(fname,R_OK) != 0) {
		fprintf(stderr,"MOSAIC: Unable to find file %s\n",fname);
		nerr++;
	    }
	    i++;
	}
	fclose(fd);

    /* Otherwise assume that this is a list separated by commas. Count
       the number of commas to define the number of files */

    } else {
	*ninfiles = 0;
	n = strlen(inlist);
	for (i = 0; i < n; i++) {
	    if (inlist[i] == ',')
		(*ninfiles)++;
	}
	(*ninfiles)++;
	
	/* Allocate space for file names */

	*infiles = cir_malloc(*ninfiles*sizeof(char *));
	for (i = 0; i < *ninfiles; i++) 
	    (*infiles)[i] = cir_malloc(BUFSIZ);
	
	/* Now grab all the file nmes */

	nerr = 0;
	for (i = 0; i < *ninfiles; i++) {
	    if (i == 0) 
		v = strtok(inlist,",");
	    else
		v = strtok(NULL,",");
	    strcpy((*infiles)[i],v);
	    if (access(v,R_OK) != 0) {
		fprintf(stderr,"MOSAIC: Unable to find file %s\n",v);
		nerr++;
	    }
	}
    }

    /* Were there any errors? */

    if (nerr != 0) 
	return(CIR_FATAL);
    else 
	return(CIR_OK);
}

static void tidy() {
    int i;
    for (i = 0; i < ninfiles; i++)
	freespace(infiles[i]);
    freespace(infiles);
    for (i = 0; i < ninconfs; i++)
	freespace(inconfs[i]);
    freespace(inconfs);
}

/*

$Log: mosaic_main.c,v $
Revision 1.3  2010/09/21 10:51:36  jim
Added verbose switch

Revision 1.2  2010-09-16 12:23:22  jim
fixed typo

Revision 1.1  2010-09-06 09:03:01  jim
New entry


*/
