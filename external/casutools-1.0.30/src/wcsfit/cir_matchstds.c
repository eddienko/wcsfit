/*

$Id: cir_matchstds.c,v 1.3 2012/07/20 09:34:44 jim Exp $

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/* #include <values.h> */
#include <limits.h>
#include <string.h>

#include <wcsfit.h>
#include <tools.h>

typedef struct {
	float x;
	float y;
	float ra;
	float dec;
} xyrd_struct;

#define INITALLOC 512
#define NGRIDMAX 61

static void sort_y(xyrd_struct *, int);
static void tidy();

static fitsfile *iptr = NULL;
static int *matches = NULL;
static xyrd_struct *obj_xy = NULL;
static xyrd_struct *obj_rd = NULL;
static float *xptr = NULL;
static float *yptr = NULL;
static float *xoffs = NULL;
static float *yoffs = NULL;
static unsigned char *bpm = NULL;
static FILE *fdo = NULL;
static int status = 0;

/*+
 *  Name:
 *      cir_matchstds
 *
 *  Purpose:
 *      Match astrometric standards to objects located on an image.
 *
 *  Description:
 *      A list of astrometric standards as generated by cir_getstds is
 *      matched to a list of x,y coordinates of objects on an image as
 *      generated by a piece of cataloging software (cir_apm or cir_imcore).
 *      Output is a list of objects that matched.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infxy = char * (Given)
 *          The name of a text file that contains the x,y positions of
 *          objects on an image.
 *      infradec = char * (Given)
 *          The name of a text file that contains the RA, Dec, x and y of
 *          astrometric standards.  The RA and Dec come from the astrometric
 *          catalogue being used, whereas the x,y coordinates come from a
 *          transformation of the equitorial coordiantes using the WCS in
 *          the header of the image from which the objects in 'infxy' have been
 *          extracted.
 *      srad = float (Given)
 *          A search radius in pixels. This helps define the number of grid
 *          points that are used in the search.
 *      nx = int (Given)
 *          The X dimension of the image.
 *      ny = int (Given)
 *          The Y dimension of the image.
 *      output = char * (Given)
 *          The name of a text file to be used for writing the results. The
 *          actual content is given in the Notes below.
 *      nm = int * (Returned)
 *          The number of objects that matched.
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      Each row of the output file will contain the following:
 *          x,y:
 *              The actual cartesian coordinates of a matched object
 *          X,Y:
 *              The cartesian coordinates for the matched object that was
 *              given in the astrometric file (calculated from the image
 *              WCS).  Useful for reference.
 *          ra,dec:
 *              The RA and Dec of the matched object.
 *              
 *  Dependencies:
 *      cfitsio, cir_misc.c, cir_stats.c
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *      Mike Irwin (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_matchstds(char *infxy, char *infradec, float srad, int nx,
                         int ny, char *output, int *nm, char *errmsg) {
    int i,nxy,nrd,ngrid,ngrid2,ibest,nmatch,jm,j,k,ig,l,dont;
    int xcol,ycol,elcol,racol,deccol,anynul;
    long nrows;
    float xoff,yoff,x,y,xoffmed,yoffmed,sigx,sigy,x1,y1,x2,y2,r1,r2;
    float val,aveden,errlim,xoffbest,yoffbest;
    char msg[BUFSIZ];

    /* Grap the info from the object catalogue fits file */

    status = 0;
    (void)fits_open_file(&iptr,infxy,READONLY,&status);
    (void)fits_get_colnum(iptr,CASESEN,"X_coordinate",&xcol,&status);
    (void)fits_get_colnum(iptr,CASESEN,"Y_coordinate",&ycol,&status);
    (void)fits_get_colnum(iptr,CASESEN,"Ellipticity",&elcol,&status); 
    (void)fits_get_num_rows(iptr,&nrows,&status);
    if (status != 0) {
	fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"MATCHSTDS: Unable to open object catalogue\n%s\n",msg);
        tidy();
	return(CIR_FATAL);
    } else if (nrows == 0) {
        (void)sprintf(errmsg,"MATCHSTDS: No rows in object catalogue\n");
        tidy();
        return(CIR_FATAL);
    }
    obj_xy = cir_malloc(nrows*sizeof(xyrd_struct));
    nxy = 0;
    for (i = 1; i <= nrows; i++) {
	(void)fits_read_col(iptr,TFLOAT,elcol,i,1,1,NULL,&val,&anynul,
			    &status);
        if (val > 0.5)
            continue;        
	(void)fits_read_col(iptr,TFLOAT,xcol,i,1,1,NULL,&val,&anynul,
			    &status);
        (obj_xy+nxy)->x = val;
	(void)fits_read_col(iptr,TFLOAT,ycol,i,1,1,NULL,&val,&anynul,
			    &status);
        (obj_xy+nxy)->y = val;
        nxy++;
    }
    closefits(iptr);
    if (status != 0) {
	fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"MATCHSTDS: Unable to read object catalogue\n%s\n",msg);
        tidy();
	return(CIR_FATAL);
    }
    
    /* Now do the same for the RA and Dec standards catalogue */

    status = 0;
    (void)fits_open_file(&iptr,infradec,READONLY,&status);
    (void)fits_get_colnum(iptr,CASEINSEN,"X_coordinate",&xcol,&status);
    (void)fits_get_colnum(iptr,CASEINSEN,"Y_coordinate",&ycol,&status); 
    (void)fits_get_colnum(iptr,CASEINSEN,"ra",&racol,&status);
    (void)fits_get_colnum(iptr,CASEINSEN,"dec",&deccol,&status);
    (void)fits_get_num_rows(iptr,&nrows,&status);
    if (status != 0) {
	fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"MATCHSTDS: Unable to open position standard catalogue\n%s\n",msg);
        tidy();
	return(CIR_FATAL);
    } else if (nrows == 0) {
        (void)sprintf(errmsg,"MATCHSTDS: No rows in position standards catalogue\n");
        tidy();
        return(CIR_FATAL);
    }
    obj_rd = cir_malloc(nrows*sizeof(xyrd_struct));
    for (nrd = 0; nrd < nrows; nrd++) {
	(void)fits_read_col(iptr,TFLOAT,xcol,nrd+1,1,1,NULL,&val,&anynul,
			    &status);
        (obj_rd+nrd)->x = val;
	(void)fits_read_col(iptr,TFLOAT,ycol,nrd+1,1,1,NULL,&val,&anynul,
			    &status);
        (obj_rd+nrd)->y = val;
	(void)fits_read_col(iptr,TFLOAT,racol,nrd+1,1,1,NULL,&val,&anynul,
			    &status);
        (obj_rd+nrd)->ra = val;
	(void)fits_read_col(iptr,TFLOAT,deccol,nrd+1,1,1,NULL,&val,&anynul,
			    &status);
        (obj_rd+nrd)->dec = val;
    }
    nrd = nrows;
    closefits(iptr);
    if (status != 0) {
	fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"MATCHSTDS: Unable to read position standard catalogue\n%s\n",msg);
        tidy();
	return(CIR_FATAL);
    }
    
    /* Open the output file */
    
    fdo = fopen(output,"w+");

    /* Now get some workspace for the X,Y coordinates of the objects */

    xptr = cir_malloc(nxy*sizeof(*xptr));
    yptr = cir_malloc(nxy*sizeof(*yptr));

    /* Sort the XY positions into order of ascending Y and put them into
       the single arrays */

    sort_y(obj_xy,nxy);
    for (i = 0; i < nxy; i++) {
        xptr[i] = (obj_xy+i)->x; 
        yptr[i] = (obj_xy+i)->y; 
    }

    /* Calculate error limit */

    aveden = (float)max(nrd,nxy)/(float)(nx*ny);
    errlim = 1.0/sqrt(8.0*M_PI*aveden);
    errlim = min(errlim,50.0);
    ngrid = (int)(srad/errlim);
    ngrid = (ngrid/2)*2 + 1;    
    ngrid = max(5,min(NGRIDMAX,ngrid));
    ngrid2 = ngrid/2 + 1;

    /* Now search for the best solution */

    ibest = 0;
    xoffbest = 0.0;
    yoffbest = 0.0;
    for (ig = -ngrid2; ig <= ngrid2; ig++) {
	xoff = (float)ig*errlim*M_SQRT2;
	for (j = -ngrid2; j <= ngrid2; j++) {
	    yoff = (float)j*errlim*M_SQRT2;
	    nmatch = 0;
	    for (k = 0; k < nrd; k++) {                  
		x = (obj_rd+k)->x + xoff;
		y = (obj_rd+k)->y + yoff;
		jm = fndmatch(x,y,xptr,yptr,nxy,errlim);
		if (jm > -1) 
		    nmatch++;
	    }
	    if (nmatch > ibest) {
		ibest = nmatch;
		xoffbest = xoff;
		yoffbest = yoff;
	    }
	}
    }

    /* Ok, go through once more at the best grid point and find a 
       median x,y positional difference */

    xoffs = cir_malloc(nrd*sizeof(*xoffs));
    yoffs = cir_malloc(nrd*sizeof(*yoffs));
    matches = cir_malloc(nrd*sizeof(*matches));
    for (k = 0; k < nrd; k++)
        matches[k] = -1;
    nmatch = 0;
    for (k = 0; k < nrd; k++) {
	x = (obj_rd+k)->x + xoffbest;
	y = (obj_rd+k)->y + yoffbest;
	jm = fndmatch(x,y,xptr,yptr,nxy,errlim);
	if (jm > -1) {
            dont = 0;
            for (l = 0; l < nrd; l++) {
                if (matches[l] == jm) {
                    x1 = xptr[jm] - ((obj_rd+l)->x + xoffbest);
                    y1 = yptr[jm] - ((obj_rd+l)->y + yoffbest);
                    x2 = xptr[jm] - ((obj_rd+k)->x + xoffbest);
                    y2 = yptr[jm] - ((obj_rd+k)->y + yoffbest);
                    r1 = sqrt(x1*x1 + y1*y1);
                    r2 = sqrt(x2*x2 + y2*y2);
                    if (r2 < r1) 
                        matches[l] = -1;
                    else
                        dont = 1;
                    break;
		}
            }
            if (dont == 0)
                matches[k] = jm;
	}
    }
    for (k = 0; k < nrd; k++) {
        if (matches[k] != -1) {
            jm = matches[k];
            xoffs[nmatch] = (xptr[jm] - (obj_rd+k)->x);
            yoffs[nmatch] = (yptr[jm] - (obj_rd+k)->y);
	    nmatch++;
	}
    }
    bpm = cir_calloc(nmatch,sizeof(unsigned char));
    (void)cir_qmedsig(xoffs,bpm,nmatch,1.0,0,-(float)nx,(float)nx,&xoffmed,
		      &sigx,msg);
    (void)cir_qmedsig(yoffs,bpm,nmatch,1.0,0,-(float)ny,(float)ny,&yoffmed,
		      &sigy,msg);
    freespace(xoffs);
    freespace(yoffs);
    freespace(bpm);

    /* Now go through one final time with a reduced error box and 
       write out the final results */

    errlim = 3.0*max(sigx,sigy);
    for (k = 0; k < nrd; k++)
        matches[k] = -1;
    for (k = 0; k < nrd; k++) {
	x = (obj_rd+k)->x + xoffmed;
	y = (obj_rd+k)->y + yoffmed;
	jm = fndmatch(x,y,xptr,yptr,nxy,errlim);
	if (jm > -1) {
            dont = 0;
            for (l = 0; l < nrd; l++) {
                if (matches[l] == jm) {
                    x1 = xptr[jm] - ((obj_rd+l)->x + xoffmed);
                    y1 = yptr[jm] - ((obj_rd+l)->y + yoffmed);
                    x2 = xptr[jm] - ((obj_rd+k)->x + xoffmed);
                    y2 = yptr[jm] - ((obj_rd+k)->y + yoffmed);
                    r1 = sqrt(x1*x1 + y1*y1);
                    r2 = sqrt(x2*x2 + y2*y2);
                    if (r2 < r1) 
                        matches[l] = -1;
                    else
                        dont = 1;
		}
            }
            if (dont == 0)
                matches[k] = jm;
	}
    }
    *nm = 0;
    for (k = 0; k < nrd; k++) {
        jm = matches[k];
	if (jm > -1) {
            fprintf(fdo,"%7.2f %7.2f %7.2f %7.2f %12.8f %12.8f\n",
                xptr[jm],yptr[jm],(obj_rd+k)->x,(obj_rd+k)->y,(obj_rd+k)->ra,(obj_rd+k)->dec);
            (*nm)++;
        }
    }
    
    /* Tidy things up... */

    tidy();

    /* Get outta here */

    return(CIR_OK);
}


static void sort_y(xyrd_struct *os, int n) {
    int iii,ii,i,ifin,j;
    float b1,b2,b3,b4;

    iii = 4;
    while (iii < n)
        iii *= 2;
    iii = min(n,(3*iii)/4 - 1);

    while (iii > 1) {
        iii /= 2;
        ifin = n - iii;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + iii;
            if ((os+i)->y > (os+j)->y) {
                b1 = (os+j)->x;
                b2 = (os+j)->y;
                b3 = (os+j)->ra;
                b4 = (os+j)->dec;
                while (1) {
                    (os+j)->x = (os+i)->x;
                    (os+j)->y = (os+i)->y;
                    (os+j)->ra = (os+i)->ra;
                    (os+j)->dec = (os+i)->dec;
                    j = i;
                    i = i - iii;
                    if (i < 0 || (os+i)->y <= b2)
                        break;
                }
                (os+j)->x = b1;
                (os+j)->y = b2;
                (os+j)->ra = b3;
                (os+j)->dec = b4;
            }
        }
    }
}

static void tidy() {

    closefits(iptr);
    freespace(matches);
    freespace(obj_xy);
    freespace(obj_rd);
    freespace(xptr);
    freespace(yptr);
    freespace(xoffs);
    freespace(yoffs);
    closefile(fdo);
}
/*

$Log: cir_matchstds.c,v $
Revision 1.3  2012/07/20 09:34:44  jim
fixed minor bug in sort routine

Revision 1.2  2010-08-02 16:10:22  jim
Now uses <limits.h> instead of <values.h>

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.10  2009/11/13 12:41:50  jim
nrd was getting set to more than the number of rows.

Revision 1.9  2009/11/05 12:01:46  jim
Increased the value of NGRIDMAX

Revision 1.8  2004/09/07 14:18:53  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.7  2004/08/19 11:59:45  jim
Little fix to doc

Revision 1.6  2004/08/19 11:34:24  jim
Added new cir_memory routines for memeory allocation

Revision 1.5  2003/11/05 14:08:15  jim
Fixed bug where fits_close_file was failing to initialise the file pointer
to NULL and hence the closefits call in tidy was falling over

Revision 1.4  2003/10/05 20:55:43  jim
Modified to read FITS tables rather than text files

Revision 1.3  2003/09/11 11:58:09  jim
Modified to use revised statistics routines

Revision 1.2  2002/12/16 10:25:06  jim
Added prologues

Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS

*/
