/*

$Id: imstack_main.c,v 1.7 2012/12/08 07:29:26 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <fitsio.h>

#include "imstack.h"
#include <tools.h>

/* Command line arguments and defaults */

enum {LTHR_ARG,
      HTHR_ARG,
      METHOD_ARG,
      NPLATE_ARG,
      EXPKEY_ARG,
      SEEWT_ARG,
      NOSEEWT_ARG,
      MAGZPTSCL_ARG,
      NOMAGZPTSCL_ARG};

static struct option myoptions [] = {
    {"lthr",required_argument,NULL,LTHR_ARG},
    {"hthr",required_argument,NULL,HTHR_ARG},
    {"method",required_argument,NULL,METHOD_ARG},
    {"nplate",required_argument,NULL,NPLATE_ARG},
    {"expkey",required_argument,NULL,EXPKEY_ARG},
    {"seewt",no_argument,NULL,SEEWT_ARG},
    {"noseewt",no_argument,NULL,NOSEEWT_ARG},
    {"magzptscl",no_argument,NULL,MAGZPTSCL_ARG},
    {"nomagzptscl",no_argument,NULL,NOMAGZPTSCL_ARG},
    {0,0,0,0}};

static char *USAGE = "Usage: imstack infiles confmaps cats outfile outconf\n"
    "[--lthr=%g] [--hthr=%g] [--method=%d] [--nplate=%d] [--expkey=%s]\n"
    "[--(no)seewt (%s)] [--(no)magzptscl (%s)]\n";

#define YESNO(a) (a == 0 ? "no" : "yes")

#define LTHR_DEF 5.0
#define HTHR_DEF 5.0
#define METHOD_DEF 1
#define NPLATE_DEF 6
#define EXPKEY_DEF "EXPTIME"
#define SEEWT_DEF 0
#define MAGZPTSCL_DEF 0

static void dummy_image(char *infile, char *inconf, int i, char *outextn,
			char *outcextn);
static int getlist(char *inlist, char ***files, int *nfiles);
static int allzeros(char *conf, int extn);
static void tidy();

static char **infiles = NULL;
static int ninfiles = 0;
static char **inconfs = NULL;
static int ninconfs = 0;
static char **incats = NULL;
static int nincats = 0;
static char **inextn = NULL;
static char **incextn = NULL;
static char **incatextn = NULL;

int main (int argc, char *argv[]) {
    int c,i,j,nerr,nextn,status,nhdu,retval,option_index,confok,n_in,n_cf,n_ct;
    int first,simple,iex,hdutype;
    char *inlist,*inconflist,*catlist,*outfile,*outconf,outextn[BUFSIZ];
    char outcextn[BUFSIZ],errmsg[BUFSIZ];
    float *magzpts,magzptref;
    fitsfile *fptr;

    /* Set up defaults for command line switches */

    float lthr = LTHR_DEF;
    float hthr = HTHR_DEF;
    int method = METHOD_DEF;
    int nplate = NPLATE_DEF;
    char expkey[81];
    (void)strcpy(expkey,EXPKEY_DEF);
    int seewt = SEEWT_DEF;
    int magzptscl = MAGZPTSCL_DEF;

    /* First get the command line arguments */

    if (argc < 6) {
        fprintf(stderr,USAGE,lthr,hthr,method,nplate,expkey,YESNO(seewt),
		YESNO(magzptscl));
	exit(1);
    } else {
	inlist = argv[1];
	inconflist = argv[2];
	catlist = argv[3];
	outfile = argv[4];
	outconf = argv[5];
    }

    /* Get any optional information */

    while ((c = getopt_long(argc,argv,"",myoptions,&option_index)) != -1) {
	switch (c) {
	case LTHR_ARG:
	    lthr = (float)atof(optarg);
	    break;
	case HTHR_ARG:
	    hthr = (float)atof(optarg);
	    break;
	case METHOD_ARG:
	    method = atoi(optarg);
	    if (method != 0 && method != 1) {
		fprintf(stderr,"method = 0 or 1 only\n");
		exit(1);
	    }
	    break;
	case NPLATE_ARG:
	    nplate = atoi(optarg);
	    if (nplate != 4 && nplate != 6) {
		fprintf(stderr,"nplate = 4 or 6 only\n");
		exit(1);
	    }
	    break;
	case EXPKEY_ARG:
	    (void)strcpy(expkey,optarg);
	    break;
	case SEEWT_ARG:
	    seewt = 1;
	    break;
	case NOSEEWT_ARG:
	    seewt = 0;
	    break;
	case MAGZPTSCL_ARG:
	    magzptscl = 1;
	    break;
	case NOMAGZPTSCL_ARG:
	    magzptscl = 0;
	    break;
	default:
	    fprintf(stderr,USAGE,LTHR_DEF,HTHR_DEF,METHOD_DEF,NPLATE_DEF,
		    EXPKEY_DEF,YESNO(SEEWT_DEF),YESNO(MAGZPTSCL_DEF));
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
    if (getlist(catlist,&incats,&nincats) != CIR_OK) {
	tidy();
	exit(1);
    }

    /* Ok, do a quick sanity check here. If there are catalogues included
       then there must be one for each input file. There can be either
       one confidence map for all the input files or a confidence map
       for each. Obviously there must be at least 1 input file! */

    if (ninfiles == 0) {
	fprintf(stderr,"IMSTACK: No input files\n");
	tidy();
	exit(1);
    }
    if (nincats > 0 && nincats != ninfiles) {
	fprintf(stderr,"IMSTACK: Nimages = %d, Ncats = %d\n",ninfiles,nincats);
	fprintf(stderr,"         Ncats must be zero or equal to Nimages\n");
	tidy();
	exit(1);
    }
    if (ninconfs != 1 && ninconfs != ninfiles) {
	fprintf(stderr,"IMSTACK: Nimages = %d, Nconfs = %d\n",ninfiles,ninconfs);
	fprintf(stderr,"         Nconfs must be 1 or equal to Nimages\n");
	tidy();
	exit(1);
    }

    /* Check the number of extensions in each file compared to the
       number of extensions in the first file in the list */

    nerr = 0;
    nextn = -1;
    simple = 0;
    for (i = 0; i < ninfiles; i++) {
	status = 0;
	(void)fits_open_file(&fptr,infiles[i],READONLY,&status);
	(void)fits_get_num_hdus(fptr,&nhdu,&status);
	(void)fits_close_file(fptr,&status);
	if (status != 0) {
	    fprintf(stderr,"IMSTACK: %s unreadable\n",infiles[i]);
	    nerr++;
	    continue;
	}
	if (nextn == -1) {
	    if (nhdu == 1) {
		nextn = 1;
		simple = 1;
	    } else {
		nextn = nhdu - 1;
	    }
	    continue;
	} else {
	    if (! simple) 
		nhdu--;
	}
	if (nhdu != nextn) {
	    fprintf(stderr,"IMSTACK: %s has %d extensions\n",infiles[i],nhdu);
	    fprintf(stderr,"         should have %d\n",nextn);
	    nerr++;
	}
    }
    for (i = 0; i < ninconfs; i++) {
	status = 0;
	(void)fits_open_file(&fptr,inconfs[i],READONLY,&status);
	(void)fits_get_num_hdus(fptr,&nhdu,&status);
	(void)fits_close_file(fptr,&status);
	if (status != 0) {
	    fprintf(stderr,"IMSTACK: %s unreadable\n",inconfs[i]);
	    nerr++;
	    continue;
	}
	if (! simple) 
	    nhdu--;
	if (nhdu != nextn) {
	    fprintf(stderr,"IMSTACK: %s has %d extensions\n",inconfs[i],nhdu);
	    fprintf(stderr,"         should have %d\n",nextn);
	    nerr++;
	}
    }
    for (i = 0; i < nincats; i++) {
	status = 0;
	(void)fits_open_file(&fptr,incats[i],READONLY,&status);
	(void)fits_get_num_hdus(fptr,&nhdu,&status);
	(void)fits_close_file(fptr,&status);
	if (status != 0) {
	    fprintf(stderr,"IMSTACK: %s unreadable\n",incats[i]);
	    nerr++;
	    continue;
	}
	nhdu--;
	if (nhdu != nextn) {
	    fprintf(stderr,"IMSTACK: %s has %d extensions\n",incats[i],nhdu);
	    fprintf(stderr,"         should have %d\n",nextn);
	    nerr++;
	}
    }
    if (nerr != 0) {
	tidy();
	exit(1);
    }

    /* If scaling by the magnitude zeropoints, work out the reference 
       zeropoint from the first image in the set */

    if (magzptscl) {
	status = 0;
	(void)fits_open_file(&fptr,infiles[0],READONLY,&status);
	magzpts = cir_malloc(nextn*sizeof(float));
	for (i = 1; i <= nextn; i++) {
	    if (! simple) {
		(void)fits_movabs_hdu(fptr,i+1,&hdutype,&status);
		iex = 0;
	    } else {
		iex = 1;
	    }
	    (void)fits_read_key(fptr,TFLOAT,"MAGZPT",magzpts+i-1,NULL,&status);
	    if (status != 0) {
		fprintf(stderr,"IMSTACK: %s[%d] header has no MAGZPT",
			infiles[0],iex);
		status = 0;
		nerr++;
	    }
	}
	closefits(fptr);
	if (nerr != 0) {
	    freespace(magzpts);
	    tidy();
	    exit(1);
	}
	(void)cir_med(magzpts,NULL,nextn,&magzptref,errmsg);
	freespace(magzpts);
    } else {
	magzptref = -1.0;
    }

    /* Get some workspace for the filename + extension specifiers */

    inextn = cir_malloc(ninfiles*sizeof(char *));
    for (i = 0; i < ninfiles; i++) 
	inextn[i] = cir_malloc(BUFSIZ);
    incextn = cir_malloc(ninconfs*sizeof(char *));
    for (i = 0; i < ninconfs; i++) 
	incextn[i] = cir_malloc(BUFSIZ);
    if (nincats != 0) {
	incatextn = cir_malloc(nincats*sizeof(char *));
	for (i = 0; i < nincats; i++) 
	    incatextn[i] = cir_malloc(BUFSIZ);
    } else {
	incatextn = NULL;
    }

    /* Loop now for each extension */

    first = 1;
    for (i = 1; i <= nextn; i++) {
	if (simple) {
	    (void)sprintf(outextn,"%s",outfile);
	    (void)sprintf(outcextn,"%s",outconf);
	    iex = 0;
	} else {
	    (void)sprintf(outextn,"%s[%d]",outfile,i);
	    (void)sprintf(outcextn,"%s[%d]",outconf,i);
	    iex = i;
	}

	/* Check the current extension of the confidence map if there is 
	   only one */

	if (ninconfs == 1) {
	    retval = allzeros(inconfs[0],iex);
	    confok = (retval == CIR_OK);
	    if (! confok) {
		fprintf(stderr,"IMSTACK: Extension %d confidence map zero\n",
			iex);
		dummy_image(infiles[0],inconfs[0],iex,outextn,outcextn);
		continue;
	    }
	}
	       
	/* Create the list of input images */

	n_in = 0;
	n_cf = 0;
	n_ct = 0;
	for (j = 0; j < ninfiles; j++) {
	    if (ninconfs != 1) {
		retval = allzeros(inconfs[j],iex);
		confok = (retval == CIR_OK);
	    }
	    if (confok) {
 	        (void)sprintf(inextn[n_in],"%s[%d]",infiles[j],iex);
		n_in++;
		if (nincats > 0) {
		    (void)sprintf(incatextn[n_ct],"%s[%d]",incats[j],i);
		    n_ct++;
		}
		if (ninconfs == 1 && j == 0) {
		    (void)sprintf(incextn[n_cf],"%s[%d]",inconfs[0],iex);
		    n_cf++;
		} else if (ninconfs != 1) {
		    (void)sprintf(incextn[n_cf],"%s[%d]",inconfs[j],iex);
		    n_cf++;
		}
	    }
	}
	if (n_in == 0) {
	    fprintf(stderr,"IMSTACK: Extension %d confidence all zero\n",i);
	    dummy_image(infiles[0],inconfs[0],iex,outextn,outcextn);
	    continue;
	}
		    
	/* Call the stacking routine */

	retval = cir_imstack_cat(inextn,incextn,incatextn,n_in,n_cf,n_ct,
				 lthr,hthr,method,nplate,expkey,seewt,
				 magzptref,outextn,outcextn,errmsg);
	if (retval != CIR_OK) {
	    fprintf(stderr,"IMSTACK: Error %s\n",errmsg);
	    tidy();
	    exit(1);
	}
	
	/* Stamp the primary */

	if (first) {
	    status = 0;
	    (void)fits_open_file(&fptr,outfile,READWRITE,&status);
	    retval = casu_stamp(fptr,"imstack");
	    (void)fits_close_file(fptr,&status);
	    (void)fits_open_file(&fptr,outconf,READWRITE,&status);
	    retval = casu_stamp(fptr,"imstack");
	    (void)fits_close_file(fptr,&status);
	    first = 0;
	}
    }

    /* Get the out of here */

    tidy();
    exit(0);
}

static void dummy_image(char *infile, char *inconf, int i, char *outextn, 
			char *outcextn) {
    int status,*cdata,true=1;
    short int *ocdata;
    fitsfile *optr,*ocptr;
    char inextn[BUFSIZ],cextn[BUFSIZ],errmsg[BUFSIZ];
    long naxis[2],npts;
    float *odata;

    /* Create input file names with extension */

    (void)sprintf(inextn,"%s[%d]",infile,i);
    (void)sprintf(cextn,"%s[%d]",inconf,i);

    /* Create the output images */

    cir_open_output(outextn,inextn,&optr,DIRECTCOPY,FLOAT_IMG,2,NULL,errmsg);
    cir_open_output(outcextn,cextn,&ocptr,DIRECTCOPY,SHORT_IMG,2,NULL,errmsg);

    /* Get the data size for the output file. This will match the confidence
       map too */

    status = 0;
    (void)fits_get_img_size(optr,2,naxis,&status);
    npts = naxis[0]*naxis[1];
    odata = cir_calloc(npts,sizeof(*odata));
    (void)fits_write_img(optr,TFLOAT,1,npts,odata,&status);
    (void)fits_update_key(optr,TLOGICAL,"IMADUMMY",&true,
			  "Dummy image from input with all zero confidence",
	                  &status);
    freespace(odata);
    closefits(optr);
    ocdata = cir_calloc(npts,sizeof(*ocdata));
    (void)fits_write_img(ocptr,TSHORT,1,npts,ocdata,&status);
    (void)fits_update_key(ocptr,TLOGICAL,"IMADUMMY",&true,
			  "Dummy image from input with all zero confidence",
			  &status);
    freespace(ocdata);
    closefits(ocptr);
}

static int allzeros(char *conf, int extn) {
    fitsfile *fptr;
    int i,status,hdutype,anynul,retval;
    long naxis[2],npts;
    short int *cdata;

    /* Open the file and move to the correct extension */

    status = 0;
    (void)fits_open_file(&fptr,conf,READONLY,&status);
    (void)fits_movabs_hdu(fptr,extn+1,&hdutype,&status);

    /* Get the array size, get some workspace and read it in */

    (void)fits_get_img_size(fptr,2,naxis,&status);
    npts = naxis[0]*naxis[1];
    cdata = cir_malloc(npts*sizeof(*cdata));
    (void)fits_read_img(fptr,TSHORT,1,npts,NULL,cdata,&anynul,&status);
    closefits(fptr);

    /* Loop for the array. As soon as you get a non-zero value then 
       break out */

    retval = CIR_FATAL;
    for (i = 0; i < npts; i++) {
	if (cdata[i] != 0) {
	    retval = CIR_OK;
	    break;
	}
    }
    freespace(cdata);
    return(retval);
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
		fprintf(stderr,"IMSTACK: Unable to open list file %s\n",inlist+1);
		return(CIR_FATAL);
	    }
	} else {
	    fprintf(stderr,"IMSTACK: Unable to find list file %s\n",inlist+1);
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
		fprintf(stderr,"IMSTACK: Unable to find file %s\n",fname);
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
		fprintf(stderr,"IMSTACK: Unable to find file %s\n",v);
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
    for (i = 0; i < ninfiles; i++)
	freespace(inextn[i]);
    freespace(inextn);
    for (i = 0; i < ninconfs; i++)
	freespace(inconfs[i]);
    freespace(inconfs);
    for (i = 0; i < ninconfs; i++)
	freespace(incextn[i]);
    freespace(incextn);
    for (i = 0; i < nincats; i++)
	freespace(incats[i]);
    freespace(incats);
    for (i = 0; i < nincats; i++)
	freespace(incatextn[i]);
    freespace(incatextn);
}

/*

$Log: imstack_main.c,v $
Revision 1.7  2012/12/08 07:29:26  jim
Modified to allow for magnitude zeropoint scaling

Revision 1.6  2011-06-27 10:57:05  jim
Modified so that the input files and the confidence maps can come from
simple fits files. NB: if one is a simple fits file, then they must all be
so. I'll change that one day...

Revision 1.5  2010-09-20 11:02:19  jim
Fixed missing memory free

Revision 1.4  2010-09-20 09:04:47  jim
Added routine dummy_image to deal with cases where all confidence values are
zero

Revision 1.3  2010-09-06 09:01:48  jim
Tidied up documentation


*/
