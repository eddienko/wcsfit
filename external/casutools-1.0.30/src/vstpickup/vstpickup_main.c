/*

$Id: vstpickup_main.c,v 1.2 2013/05/30 14:32:01 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>

#include <vstpickup.h>
#include <tools.h>

/* Command line arguments and defaults */

enum {FLAT_ARG,
      CHECKONLY_ARG,
      NOCHECKONLY_ARG};

static struct option myoptions [] = {
    {"flat",required_argument,NULL,FLAT_ARG},
    {"checkonly",no_argument,NULL,CHECKONLY_ARG},
    {"nocheckonly",no_argument,NULL,NOCHECKONLY_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: vstpickup: infile\n"
    "[--flat=%s] [--(no)checkonly (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define FLAT_DEF ""
#define CHECKONLY_DEF 0

int main (int argc, char *argv[]) {
    int option_index,c,haspickup,retval;
    char msg[BUFSIZ],*infile;

    /* Set up defaults for command line switches */

    char flatfile[BUFSIZ];
    (void)strcpy(flatfile,FLAT_DEF);
    int checkonly = CHECKONLY_DEF;

    /* First get the command line arguments */

    if (argc < 2) {
	fprintf(stderr,USAGE,flatfile,YESNO(checkonly));
	exit(1);
    } else {
	infile = argv[1];
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case FLAT_ARG:
	    (void)strcpy(flatfile,optarg);
	    break;
	case CHECKONLY_ARG:
	    checkonly = 1;
	    break;
	case NOCHECKONLY_ARG:
	    checkonly = 0;
	    break;
	default:
	    fprintf(stderr,USAGE,flatfile,YESNO(checkonly));
	    exit(1);
	}
    }

    /* Call the main processing routine */

    retval = cir_vstpickup(infile,flatfile,checkonly,&haspickup,msg);
    if (retval != CIR_OK) {
	fprintf(stderr,"VSTPICKUP: Error %s\n",msg);
	exit(1);
    }

    /* If this was a check-only, then send out a little message */

    if (checkonly) 
	(void)fprintf(stdout,"VSTPICKUP: File %s had %d rows of pickup\n",
		      infile,haspickup);
    exit(0);
}

/*

$Log: vstpickup_main.c,v $
Revision 1.2  2013/05/30 14:32:01  jim
Fixed little problem in one of the error messages

Revision 1.1  2012-12-08 07:31:00  jim
New entry


*/
