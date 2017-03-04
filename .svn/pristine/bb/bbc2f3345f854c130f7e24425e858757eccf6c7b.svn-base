/*

$Id: cir_catcoord.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/
#include <stdio.h>
#include <stdlib.h>

#include <tools.h>
 
static void tidy();

static struct wcsprm *wcs = NULL;
static fitsfile *fptr = NULL;
static float *xptr = NULL;
static float *yptr = NULL;
static float *raptr = NULL;
static float *decptr = NULL;
static int status = 0;


/*+
 *  Name:
 *      cir_catcoord
 *
 *  Purpose:
 *      Calculate RA and Dec values for objects in a catalogue.
 *
 *  Description:
 *      The World Coordinates for objects in a catalogue are calculated
 *      from the header WCS.  The RA and Dec is written into the catalogue
 *      for each object.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infile = char * (Given)
 *          The name of the input FITS image. This image must have a world
 *          cooordinate system in its header.
 *      incat = char * (Given)
 *          The name of the input catalogue.  This should be a fits table with
 *          the extension number specified.  The required columns are listed in
 *          the Notes below.
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      The input catalogue must have the following columns:
 *          X_coordinate = float (pixels)
 *              The X coordinate of the object
 *          Y_coordinate = float (pixels)
 *              The Y coordinate of the object
 *      This routine writes the output world coordiantes in columns previously
 *      defined as 'RA' and 'DEC'. It's an error if they don't already exist.
 *              
 *  Dependencies:
 *      cfitsio, cir_wcssubs.c
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_catcoord(char *infile, char *incat, char *errmsg) {
    int colx,coly,anynul,i,racol,deccol,retval;
    long nrows;
    char msg[BUFSIZ];
    double xx,yy,rr,dd;

    /* First open the file and get the WCS info */

    retval = cir_wcsopen(infile,&wcs,msg);
    if (retval != CIR_OK) {
        (void)sprintf(errmsg,"CATCOORD: Couldn't get WCS in %s -- %s\n",
		      infile,msg);
        return(CIR_FATAL);
    }

    /* Open the catalogue file */

    (void)fits_open_file(&fptr,incat,READWRITE,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"CATCOORD: Couldn't open input catalogue %s -- %s\n",
		      incat,msg);
	tidy();
        return(CIR_FATAL);
    }

    /* Find the number of rows in the table and get the workspace for 
       x,y,ra and dec positions */

    (void)fits_get_num_rows(fptr,&nrows,&status);
    xptr = cir_malloc(nrows*sizeof(float));
    yptr = cir_malloc(nrows*sizeof(float));
    raptr = cir_malloc(nrows*sizeof(float));
    decptr = cir_malloc(nrows*sizeof(float));

    /* Get the information from each column */

    (void)fits_get_colnum(fptr,CASEINSEN,"X_coordinate",&colx,&status);
    (void)fits_read_col(fptr,TFLOAT,colx,1,1,nrows,NULL,xptr,&anynul,&status);
    (void)fits_get_colnum(fptr,CASEINSEN,"Y_coordinate",&coly,&status);
    (void)fits_read_col(fptr,TFLOAT,coly,1,1,nrows,NULL,yptr,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"CATCOORD: Error reading columns in %s -- %s\n",
		      incat,msg);
        tidy();
        return(CIR_FATAL);
    }

    /* Find RA and Dec columns */

    (void)fits_get_colnum(fptr,CASEINSEN,"RA",&racol,&status);
    (void)fits_get_colnum(fptr,CASEINSEN,"DEC",&deccol,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"CATCOORD: Error locating RA/DEC columns in %s -- %s\n",
		      incat,msg);
        tidy();
        return(CIR_FATAL);
    }


    /* Now loop for each object in the catalogue and define equitorial coord */

    for (i = 0; i < nrows; i++) {
        xx = (double)xptr[i];
        yy = (double)yptr[i];
        cir_xytoradec(wcs,xx,yy,&rr,&dd);
        raptr[i] = (float)rr;
        decptr[i] = (float)dd;
    }

    /* Now write the column data */

    (void)fits_write_col(fptr,TFLOAT,racol,1,1,nrows,raptr,&status);
    (void)fits_write_col(fptr,TFLOAT,deccol,1,1,nrows,decptr,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errmsg,"CATCOORD: Error writing column data in %s -- %s\n",
		      incat,msg);
        tidy();
        return(CIR_FATAL);
    }

    /* Get out of here */

    tidy();
    return(CIR_OK);
}

static void tidy() {
    int status;

    closefits(fptr);
    freewcs(wcs);
    freespace(xptr);
    freespace(yptr);
    freespace(raptr);
    freespace(decptr);
}
/*

$Log: cir_catcoord.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.7  2004/08/19 11:34:23  jim
Added new cir_memory routines for memeory allocation

Revision 1.6  2004/08/02 12:16:45  jim
Fixed missing reference to DEGRAD

Revision 1.5  2004/08/02 11:50:01  jim
Modified to use new version of cir_wcssubs routines

Revision 1.4  2004/06/11 08:47:47  jim
Added cvsid and cir_stamp call

Revision 1.3  2004/04/05 11:26:29  jim
Modified for new version of imcore

Revision 1.2  2002/12/16 10:25:06  jim
Added prologues

Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS


*/
