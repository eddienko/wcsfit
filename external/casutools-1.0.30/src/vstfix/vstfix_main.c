/*

$Id: vstfix_main.c,v 1.1 2013/05/30 14:28:26 jim Exp $
 
*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>

#include <vstfix.h>
#include <tools.h>

/* Command line arguments and defaults */

enum {CONF_ARG,
      VERBOSE_ARG,
      NOVERBOSE_ARG};

static struct option myoptions[] = {
    {"conf",required_argument,NULL,CONF_ARG},
    {"verbose",no_argument,NULL,VERBOSE_ARG},
    {"noverbose",no_argument,NULL,NOVERBOSE_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: vstfix incat corr_file outcat\n"
    "[--conf=%s] [--(no)verbose (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")
#define CONF_DEF "auto"
#define VERBOSE_DEF 0

int main (int argc, char *argv[]) {
    char msg[BUFSIZ],*incat,*corr_file,*outcat;
    int c,option_index,retval;

    /* Set up defaults for command line switches */

    char confmap[BUFSIZ];
    (void)strcpy(confmap,CONF_DEF);
    int verbose = VERBOSE_DEF;

    /* Check that we have enough information */

    if (argc < 4) {
	fprintf(stderr,USAGE,confmap,YESNO(verbose));
        exit(1);
    } else {
	incat = argv[1];
	corr_file = argv[2];
	outcat = argv[3];
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case CONF_ARG:
	    (void)strcpy(confmap,optarg);
	    break;
	case VERBOSE_ARG:
	    verbose = 1;
	    break;
	case NOVERBOSE_ARG:
	    verbose = 0;
	    break;
	default:
	    fprintf(stderr,USAGE,confmap,YESNO(verbose));
	    exit(1);
	}
    }

    /* Call the main processing routine */

    retval = cir_vstfix(incat,corr_file,confmap,outcat,msg);
    if (retval != CIR_OK) {
	fprintf(stderr,"VSTFIX: Error %s\n",msg);
	exit(1);
    }
    exit(0);
}

/*

$Log: vstfix_main.c,v $
Revision 1.1  2013/05/30 14:28:26  jim
New Entry


*/
