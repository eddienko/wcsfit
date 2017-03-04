/*

$Id: imcore_rdbuf_mef.c,v 1.2 2010/08/02 16:04:51 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>

#include "errcodes.h"
#include "imcore.h"
#include "util.h"

extern int imcore_rdbuf_mef(char *infile, int datatype, void **map, long *nx, 
			    long *ny, int verbose, char *errstr) {
    int status,naxes,anynull;
    long naxis[2],npix;
    char errmsg[FLEN_STATUS],strval[32];
    fitsfile *fptr;

    /* Open file readonly */

    status = 0;
    ffopen(&fptr,infile,READONLY,&status);
    if (status != 0) {
        ffgerr(status,errmsg);
        sprintf(errstr,"ffopen of file %s failed: %s (status = %d)",
		infile,errmsg,status);
        return(ERRCODE_FILE_IO);
    }

    /* Find array dimensions */

    ffgidm(fptr,&naxes,&status);
    if (naxes != 2) {
        sprintf(errstr,"Image is not 2d");
        return(ERRCODE_FILE_DATA);
    }
    ffgisz(fptr,2,naxis,&status);

    /* Set up output axes */

    *nx = naxis[0];
    *ny = naxis[1];
    npix = naxis[0]*naxis[1];
    if (verbose) {
        printf("FITS header information:  Extension %s\n                          No. of columns =%8ld\n                          No. of rows  =%8ld\n",
	     infile,*nx,*ny);
    }

    /* ok now read all the data */

    if (verbose) 
        printf("\nReading file ......\n\n");

    /* Assume we are processing all the data and malloc the big array */

    *map = NULL;
    switch (datatype) {
    case TFLOAT:
        *map = (float *)malloc(npix*sizeof(float));
        break;
    case TSHORT:
        *map = (short int *)malloc(npix*sizeof(short int));
        break;
    default:
        sprintf(errstr,"routine doesn't cater for datatype = %d\n",datatype);
        closefits(fptr);
        return(ERRCODE_MALLOC);
        break;
    }
    if (! *map) {
        sprintf(errstr,"Couldn't get space image data");
        closefits(fptr);
        return(ERRCODE_MALLOC);
    }
   
    /* OK now read all the data */

    ffgpv(fptr,datatype,1L,npix,NULL,*map,&anynull,&status);

    /* Report any errors */

    if (status != 0) {
        ffgerr(status,errmsg);
        sprintf(errstr,"ffgpvj: %s (status = %d)",errmsg,status);
        freespace(map);
        closefits(fptr);
        return(ERRCODE_FILE_IO);
    }

    /* Read the gain here. This is a bit of a kludge to cover us until I have
       time to change the api to include the gain keyword */

    (void)fits_read_key(fptr,TFLOAT,"GAIN",&gain,NULL,&status);
    if (status != 0) {
	status = 0;
	gain = 5.0;
    }

    /* See if it's been interpolated at some stage */

    (void)fits_read_key(fptr,TSTRING,"DRIBBLE",strval,NULL,&status);
    dribble = (status == 0);

    /* Tidy and get out of here */

    closefits(fptr);    
    return(ERRCODE_OK);
}


extern int imcore_rdbuf_conf(char *infile, short int **map, float **mapsqrt,
			     long nx, long ny, int verbose, char *errstr) {
    long npix,nxc,nyc;
    int retval,i,old_dribble;

    /* If no confidence map is specified, then just create an array of 100s
       (100% confidence) */

    if (!strcmp(infile,"noconf")) {
        npix = nx*ny;
        *map = (short *)malloc(npix*sizeof(short int));
        if (! *map) {
            sprintf(errstr,"Couldn't allocate space for confidence map\n");
            return(ERRCODE_MALLOC);
	}
        for (i = 0; i < npix; i++) 
            (*map)[i] = 100;

	/* Create sqrt map */

        *mapsqrt = (float *)malloc(npix*sizeof(float));
	for (i = 0; i < npix; i++) 
            (*mapsqrt)[i] = sqrt(0.01*(float)((*map)[i]));

    /* Otherwise, open the confidence map and read it. Save and restore
       current status of dribble global variable in case the confidence
       map doesn't have dribble info in the header. */

    } else {
	old_dribble = dribble;
	npix = nx*ny;
        retval = imcore_rdbuf_mef(infile,TSHORT,(void *)map,&nxc,&nyc,verbose,
				  errstr);
        if (retval != ERRCODE_OK)
	    return(retval);
        if (nxc != nx || nyc != ny) {
            sprintf(errstr,"Confidence map and image dimensions do no match");
            freespace(map);
            return(ERRCODE_FILE_DATA);
	}
	dribble = old_dribble;
	
	/* Trap for stupid values */

	for (i = 0; i < npix; i++)
	    (*map)[i] = MAX(0,MIN(110,(*map)[i]));

	/* Create sqrt map */

        *mapsqrt = (float *)malloc(npix*sizeof(float));
	for (i = 0; i < npix; i++) 
            (*mapsqrt)[i] = sqrt(0.01*(float)((*map)[i]));

    }
    return(ERRCODE_OK);
}

/*

$Log: imcore_rdbuf_mef.c,v $
Revision 1.2  2010/08/02 16:04:51  jim
Fixed bug where sqrt map wasn't being provided when running without a
confidence map

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.6  2010/02/11 21:55:33  jim
changed a few routine declarations

Revision 1.5  2010/01/19 12:01:47  jim
removed reference to colour stuff

Revision 1.4  2009/12/17 11:34:07  jim
added check for dribble keyword

Revision 1.3  2008/04/15 19:07:48  jim
Added colour routine

Revision 1.2  2005/08/26 04:46:20  jim
Modified to add new radii and error estimates

Revision 1.1  2004/04/02 10:55:00  jim
New version for rewrite of imcore


*/
