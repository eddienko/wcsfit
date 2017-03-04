#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>

#include "extract.h"
#include "usnob.h"
#include "psc.h"
#include "floatmath.h"
#include "util.h"

#define RAD2DEG (180.0/M_PI)
#define freespace(_p) if (_p != NULL) {free(_p); _p = NULL;}

#define NCOLUSNOB 6
static char *ttype_usnob[6] = {"ra","dec","Rmag1","Bmag1","Rmag2","Bmag2"};
static char *tform_usnob[6] = {"1E","1E","1E","1E","1E","1E"};
static char *tunit_usnob[6] = {"degrees","degrees","mags","mags","mags",
			       "mags"};
#define NCOL2MASS 8
static char *ttype_2mass[8] = {"ra","dec","Jmag","Hmag","Kmag","Jmag_err",
			       "Hmag_err","Kmag_err"};
static char *tform_2mass[8] = {"1E","1E","1E","1E","1E","1E","1E","1E"};
static char *tunit_2mass[8] = {"degrees","degrees","mags","mags","mags","mags",
			       "mags","mags"};

extern int extract_cat(char *catname, char *path, float ra_cent, 
		       float dec_cent, float half_width, long *n, 
		       char *outfile, char *errmsg) {
    char catpath[BUFSIZ],msg[BUFSIZ];
    struct usnob_handle uh;
    struct usnob_row ur;
    struct psc_handle th;
    struct psc_row tr;
    int status,i,colnum;
    long nrows;
    float *workspace = NULL;
    float *ra,*dec;
    float *rmag1 = NULL;
    float *rmag2 = NULL;
    float *bmag1 = NULL;
    float *bmag2 = NULL;
    float *jmag = NULL;
    float *kmag = NULL;
    float *hmag = NULL;
    float *jmagerr = NULL;
    float *kmagerr = NULL;
    float *hmagerr = NULL;
    fitsfile *iptr;

    /* Put RA, Dec and halfwidth into radians */

    ra_cent /= RAD2DEG;
    dec_cent /= RAD2DEG;
    half_width /= RAD2DEG;

    /* Which catalogue are you interested in? */

    if (! strcmp(catname,"usnob")) {

        /* Open USNOB handle structure */

        (void)sprintf(catpath,"%s/b%%%%%%%%.cat",path);
	status = usnob_open(&uh,catpath,errmsg);
	if (status != 0)
	    return(1);

        /* Are there any rows? If so then get some workspace for them */

        *n = usnob_find(&uh,ra_cent,dec_cent,half_width,errmsg);
        if (*n < 0) {
	    return(1);
	} else if (*n == 0) {
	    ra = NULL;
	    dec = NULL;
        } else {
	    workspace = malloc(NCOLUSNOB*(*n)*sizeof(*workspace));
	    ra = workspace;
            dec = ra + *n;
            rmag1 = dec + *n;
  	    rmag2 = rmag1 + *n;
            bmag1 = rmag2 + *n;
            bmag2 = bmag1 + *n;
	    if (! workspace) {
		sprintf(errmsg,"EXTRACT: Memory allocation failed\n");
		status = usnob_close(&uh,errmsg);
		return(1);
	    }
	}

        /* Now read the damn things */

        for (i = 0; i < *n; i++) {
	    status = usnob_read(&uh,&ur,errmsg);
	    if (status == -1) {
		(void)usnob_close(&uh,errmsg);
		freespace(workspace);
		return(1);
	    } else if (status == 0)
		break;
	    ra[i] = ((float)ur.ra)*RAD2DEG;
	    dec[i] =((float)ur.dec)*RAD2DEG;
	    rmag1[i] = (float)ur.r1.m;
	    rmag2[i] = (float)ur.r2.m;
	    bmag1[i] = (float)ur.b1.m;
	    bmag2[i] = (float)ur.b2.m;
	}

        /* Redefine the number of rows to those actually read */

        *n = i;

	/* Tidy up and get out of here */

        status = usnob_close(&uh,errmsg);

        /* If now rows were returned, then leave now... (NB: This isn't 
           necessarily an error) */

        if (*n == 0) {
	    freespace(workspace);
 	    return(0);
        }
	nrows = (long)*n;
    
        /* Create the FITS file */

	status = 0;
	(void)fits_create_file(&iptr,outfile,&status);
	(void)fits_create_tbl(iptr,BINARY_TBL,nrows,NCOLUSNOB,ttype_usnob,
			      tform_usnob,tunit_usnob,NULL,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"ra",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,ra,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"dec",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,dec,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Rmag1",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,rmag1,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Rmag2",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,rmag2,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Bmag1",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,bmag1,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Bmag2",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,bmag2,&status);
	(void)fits_close_file(iptr,&status);
	if (status != 0) {
	    (void)fits_get_errstatus(status,msg);
	    (void)sprintf(errmsg,"EXTRACT_CAT: Error writing FITS table: %s\n",
			  msg);
	    freespace(workspace);
	    return(1);
	}
	freespace(workspace);

    } else if (! strcmp(catname,"2mass")) {

        /* Open 2MASS handle structure */

        (void)sprintf(catpath,"%s/psc%%%%.fits",path);
	status = psc_open(&th,catpath,errmsg);
	if (status != 0)
	    return(1);

        /* Are there any rows? If so then get some workspace for them */

        *n = psc_find(&th,ra_cent,dec_cent,half_width,errmsg);
        if (*n < 0) {
	    return(1);
	} else if (*n == 0) {
	    ra = NULL;
	    dec = NULL;
        } else {
	    workspace = malloc(NCOL2MASS*(*n)*sizeof(*workspace));
	    ra = workspace;
            dec = ra + *n;
	    jmag = dec + *n;
	    hmag = jmag + *n;
	    kmag = hmag + *n;
	    jmagerr = kmag + *n;
	    hmagerr = jmagerr + *n;
            kmagerr = hmagerr + *n;
	    if (! workspace) {
		sprintf(errmsg,"EXTRACT: Memory allocation failed\n");
		status = usnob_close(&uh,errmsg);
		return(1);
	    }
	}

        /* Now read the damn things */

        for (i = 0; i < *n; i++) {
	    status = psc_read(&th,&tr,errmsg);
	    if (status == -1) {
		(void)psc_close(&th,errmsg);
		return(1);
	    } else if (status == 0)
		break;
	    ra[i] = ((float)tr.ra)*RAD2DEG;
	    dec[i] = ((float)tr.dec)*RAD2DEG;
	    jmag[i] = (float)tr.j.m;
	    hmag[i] = (float)tr.h.m;
	    kmag[i] = (float)tr.k.m;
	    jmagerr[i] = (float)tr.j.msig;
	    hmagerr[i] = (float)tr.h.msig;
	    kmagerr[i] = (float)tr.k.msig;
	}

        /* Redefine the number of rows to those actually read... */

        *n = i;

	/* Tidy up and get out of here */

        status = psc_close(&th,errmsg);

        /* If now rows were returned, then leave now... (NB: This isn't 
           necessarily an error) */

        if (*n == 0) {
	    freespace(ra);
   	    freespace(dec);
 	    return(0);
        }
    
        /* If now rows were returned, then leave now... (NB: This isn't 
           necessarily an error) */

        if (*n == 0) {
	    freespace(workspace);
 	    return(0);
        }
	nrows = (long)*n;
    
        /* Create the FITS file */

	status = 0;
	(void)fits_create_file(&iptr,outfile,&status);
	(void)fits_create_tbl(iptr,BINARY_TBL,nrows,NCOL2MASS,ttype_2mass,
			      tform_2mass,tunit_2mass,NULL,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"ra",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,ra,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"dec",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,dec,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Jmag",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,jmag,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Hmag",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,hmag,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Kmag",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,kmag,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Jmag_err",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,jmagerr,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Hmag_err",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,hmagerr,&status);
	(void)fits_get_colnum(iptr,CASEINSEN,"Kmag_err",&colnum,&status);
	(void)fits_write_col(iptr,TFLOAT,colnum,1,1,*n,kmagerr,&status);
	(void)fits_close_file(iptr,&status);
	if (status != 0) {
	    (void)fits_get_errstatus(status,msg);
	    (void)sprintf(errmsg,"EXTRACT_CAT: Error writing FITS table: %s\n",
			  msg);
	    freespace(workspace);
	    return(1);
	}
	freespace(workspace);

    } else {
	(void)sprintf(errmsg,"EXTRACT: Unrecognised catalogue name: %s\n",
		      catname);
	return(1);
    }
    

    return(0);
}
