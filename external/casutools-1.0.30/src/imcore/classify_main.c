/*

$Id: classify_main.c,v 1.2 2011/01/19 12:46:32 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <fitsio.h>
#include <math.h>

#include <tools.h>

enum {VERBOSE_ARG,
      NOVERBOSE_ARG};

static struct option myoptions [] = {
    {"verbose",no_argument,NULL,VERBOSE_ARG},
    {"noverbose",no_argument,NULL,NOVERBOSE_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: classify infile\n"
    "[--(no)verbose (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define VERBOSE_DEF 0

extern int cir_classify(char *infile, char *expkey, float minsize, 
			char *errmsg);

int main (int argc, char *argv[]) {
    int status,c,option_index,nhdu,ist,ifn,first,i,retval;
    fitsfile *fptr;
    char *infile,inf_ext[BUFSIZ],errmsg[BUFSIZ];

    /* Set up defaults for command line switches */

    int verbose = VERBOSE_DEF;

    /* Read command line */

    if (argc < 2) {
	fprintf(stderr,USAGE,YESNO(VERBOSE_DEF));
	exit(1);
    } else {
	infile = argv[1];
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case VERBOSE_ARG:
	    verbose = 1;
	    break;
	case NOVERBOSE_ARG:
	    verbose = 0;
	    break;
	default:
	    fprintf(stderr,USAGE,YESNO(VERBOSE_DEF));
	    exit(1);
	}
    }

    /* If doing verbose output set stdout to have no buffering */

    if (verbose)
	setvbuf(stdout,NULL,_IONBF,0);

    /* Find out how many image extensions there are */

    status = 0;
    (void)fits_open_file(&fptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(fptr,&nhdu,&status);
    (void)fits_close_file(fptr,&status);
    if (status != 0) {
        fprintf(stderr,"Unable to open input file %s\n",infile);
        exit(1);
    }
    ist = 1;
    ifn = nhdu - 1;
    
    /* Loop for all extensions */

    first = 1;
    for (i = ist; i <= ifn; i++) {
        (void)sprintf(inf_ext,"%s[%d]",infile,i);
	if (verbose)
	    fprintf(stdout,"Working on %s\n",inf_ext);
	retval = cir_classify(inf_ext,"",16,errmsg);
        if (retval != 0)
  	    fprintf(stderr,"Error running classify\n%s\n",errmsg);

	/* Now stamp the primary */

	if (first) {
	    status = 0;
	    (void)fits_open_file(&fptr,infile,READWRITE,&status);
	    retval = casu_stamp(fptr,"classify");
	    if (retval != 0) 
		fprintf(stderr,"Error stamping output table primary %s\n",
		        infile);
	    (void)fits_close_file(fptr,&status);
	    first = 0;
	}
    }
    return(0);
}

/*

$Log: classify_main.c,v $
Revision 1.2  2011/01/19 12:46:32  jim
Fixed options menu in header

Revision 1.1  2010-11-23 09:52:12  jim
New entry


*/
