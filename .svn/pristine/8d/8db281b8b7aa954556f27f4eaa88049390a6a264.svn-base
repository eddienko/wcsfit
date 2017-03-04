/*

$Id: wcsfit.c,v 1.4 2014/11/03 08:23:34 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>

#include <tools.h>
#include <wcsfit.h>

static void edittabs(char *incat, char *tmpfile, char *newcat);
static void update_wcskeys(char *infile, char *incat);
static void update_wcskeys_table(char *infile, char *incat);
static void update_wcskeys_image(char *infile, char *incat);
static void sort1(float *a, int n);

#define NMAGCOL 3
static char *magcol[NMAGCOL] = {"Kmag","Rmag1","R1mag"};

/*+
 *  Name:
 *      wcsfit
 *
 *  Purpose:
 *      Work out the WCS solution for an image and its associated object
 *      catalogue.
 *
 *  Description:
 *      An input FITS image with a rough WCS defined in the header is given
 *      along with an object catalogue. The rough WCS is used to work out
 *      the coverage of the image and astrometric standards in that region
 *      are extracted from a catalogue. These are matched to the objects in
 *      the given object catalogue. The matched objects are used to define
 *      a WCS whose parameters are written to the headers of both the image
 *      and the catalogue using the standard FITS keywords.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infile = char * (Given)
 *          An input FITS extension. This must have a WCS defined in the header.
 *      incat = char * (Given)
 *          An input source catalogue for the input image. 
 *      catsrc = char * (Given)
 *          This can be one of several values.  An internet connection is
 *          required for any starting with 'viz':
 *              localfits:
 *                  A locally held fits table will be searched. The full
 *                  path to this table must be given in 'catpath'.  The only
 *                  columns that the table must have are 'ra' and 'dec'.  Both
 *                  must be in degrees.
 *              viz2mass:
 *                  The 2MASS point source catalogue searched with VizieR
 *              vizlandolt:
 *                  The Landolt SA region standards catalogue search with VizieR
 *              vizusnob:
 *                  The USNO-B catalogue searched with Vizier.
 *              localusnob:
 *                  A local copy of the USNOB catalogue is searched in its 
 *                  native format. The path to the directory holding all the 
 *                  .cat and .acc files should be given in 'catpath'.
 *              local2mass:
 *                  A local copy of the 2MASS psc catalogue is searched in its 
 *                  native format. The path to the directory holding all the 
 *                  fits table files should be given in 'catpath'.
 *      site = char * (Given)
 *          The site for the VizieR search (e.g. catsrc = "vizusnob"). The copy
 *          of VizieR that is searched can be one of the following:
 *              casu:
 *                  CASU, Cambridge, UK
 *              cds:
 *                  CDS, Strasbourg, France
 *              ukirt:
 *                  UKIRT, Hawaii, USA
 *      catpath = char * (Given)
 *          The full path to a FITS table with a catalogue (catsrc=localfits)
 *          or the full path to a directory with local catalogue copies 
 *          (e.g. catsrc=vizusnob). 
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      None
 *
 *  Dependencies:
 *      cfitsio, wcslib, tools
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2010-2013 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int wcsfit(char *infile, char *incat, char *catsrc, char *site, 
		  char *catpath, char *errmsg) {
    int status,retval,maxsize,nmatch;
    long naxis[2],ncat;
    fitsfile *iptr,*cptr;
    char tmpfile[BUFSIZ],inroot[BUFSIZ],errstr[BUFSIZ],newcat[BUFSIZ];
    char useme[BUFSIZ],matchfile[BUFSIZ],tmpextn[BUFSIZ];
    float srad;

    /* Open the image and catalogue FITS files and make sure it's ok to
       read and write to them */

    status = 0;
    (void)fits_open_file(&iptr,infile,READWRITE,&status);
    (void)fits_get_img_size(iptr,2,naxis,&status);
    if (status != 0) {
        (void)sprintf(errmsg,"WCSFIT: Unable to open image %s\n",infile);
	closefits(iptr)
	return(CIR_FATAL);
    }
    closefits(iptr);
    status = 0;
    (void)fits_open_file(&cptr,incat,READWRITE,&status);
    if (status != 0) {
        (void)sprintf(errmsg,"WCSFIT: Unable to open catalogue %s\n",incat);
	closefits(iptr);
	closefits(cptr);
	return(CIR_FATAL);
    }
    (void)fits_get_num_rows(cptr,&ncat,&status);
    closefits(cptr);
    if (ncat < 1) {
        (void)sprintf(errmsg,"WCSFIT: No rows in input catalogue %s\n",incat);
        closefits(iptr);
        return(CIR_WARN);
    }

    /* Get some standards for the input image. Create a temporary file name 
       from the rootname of the FITS image */

    (void)fits_parse_rootname(infile,inroot,&status);
    (void)sprintf(tmpfile,"stds_%s",basename(inroot));
    if (access(tmpfile,F_OK) == 0) 
	remove(tmpfile);
    retval = cir_getstds(infile,tmpfile,catsrc,site,catpath,2000.0,5,1,
			 errstr);
    if (retval != CIR_OK) {
	(void)sprintf(errmsg,"WCSFIT: Unable to get standards -- %s\n",
		      errstr);
        if (access(tmpfile,F_OK) == 0) 
	    remove(tmpfile);
	return(CIR_FATAL);
    }

    /* Edit out some of the standards so that we remove the ones with high
       photometric error. Also cut down the catalogue if the total number of
       stars is much higher than what we have in the standards table */

    newcat[0] = '\0';
    edittabs(incat,tmpfile,newcat);

    /* Match the XY positions in the catalogues to the predicted XY positions
       in the catalogues */

    if (strlen(newcat) == 0) 
	strcpy(useme,incat);
    else 
	strcpy(useme,newcat);
    maxsize = max(naxis[0],naxis[1]);
    srad = (int)(0.25*(float)maxsize);
    (void)sprintf(matchfile,"match_%s",basename(inroot));
    if (access(matchfile,F_OK) == 0)
	remove(matchfile);
    (void)sprintf(tmpextn,"%s[1]",tmpfile);
    retval = cir_matchstds(useme,tmpextn,srad,naxis[0],naxis[1],matchfile,
			   &nmatch,errstr);
    if (strlen(newcat)) {
        (void)fits_parse_rootname(newcat,inroot,&status);
	remove(inroot);
    }
    remove(tmpfile);
    if (retval != CIR_OK) {
        (void)sprintf(errmsg,"WCSFIT: Matching routine failed -- %s\n",
		      errstr);
	return(CIR_FATAL);
    } else if (nmatch == 0) {
        (void)sprintf(errmsg,"WCSFIT: Matching routine found 0 matches\n");
	return(CIR_FATAL);
    }

    /* Fit the plate solution */

    retval = cir_platesol(infile,matchfile,6,2,1,errstr);
    remove(matchfile);
    if (retval != CIR_OK) {
	(void)sprintf(errmsg,"WCSFIT: Plate solution routine failed -- %s\n",
		      errstr);
	return(CIR_FATAL);
    }

    /* Update the catalogue file */

    cir_catcoord(infile,incat,errstr);

    /* Update header of input catalogue with WCS keywords */

    update_wcskeys(infile,incat);

    /* Get out of here */

    return(CIR_OK);
}

static void edittabs(char *incat, char *stdsfile, char *newcat) {
    char *jerr1 = "Jmag_err";
    char *jerr2 = "e_Jmag";
    char *kerr1 = "Kmag_err";
    char *kerr2 = "e_Kmag";
    char *jerr,*kerr,expr[64],catroot[BUFSIZ],colname[16];
    fitsfile *fptr,*nptr;
    int status,hdutype,jerrcol,kerrcol,afluxcol,level,ind,anynul,i,colnum;
    int maxstds = 2000;
    long nrows_stds,nrows_cat;
    float cut,*fluxes;

    /* Get the number of standards */

    status = 0;
    (void)fits_open_file(&fptr,stdsfile,READWRITE,&status);
    (void)fits_movabs_hdu(fptr,2,&hdutype,&status);
    (void)fits_get_num_rows(fptr,&nrows_stds,&status);

    /* Check to see if magnitude error columns are available for 2mass
       data */

    jerr = jerr1;
    (void)fits_get_colnum(fptr,CASESEN,jerr1,&jerrcol,&status);
    if (status != 0) {
	status = 0;
	jerr = jerr2;
	(void)fits_get_colnum(fptr,CASESEN,jerr2,&jerrcol,&status);
    }
    kerr = kerr1;
    if (status == 0) {
        (void)fits_get_colnum(fptr,CASESEN,kerr1,&kerrcol,&status);
	if (status != 0) {
	    status = 0;
	    kerr = kerr2;
	    (void)fits_get_colnum(fptr,CASESEN,kerr2,&kerrcol,&status);
	}
    }
    
    /* If they are, then select out the rows with good photometry */
    
    if (status == 0) {
	(void)sprintf(expr,"(%s < 0.2 && %s < 0.2)",jerr,kerr);
	(void)fits_select_rows(fptr,fptr,expr,&status);
	(void)fits_get_num_rows(fptr,&nrows_stds,&status);
	(void)fits_close_file(fptr,&status);
    } else {
	status = 0;
	fits_close_file(fptr,&status);
    }

    /* See if there are too many standard stars */

    if (nrows_stds > maxstds) {
        (void)fits_open_file(&fptr,stdsfile,READWRITE,&status);
        (void)fits_movabs_hdu(fptr,2,&hdutype,&status);
        for (i = 0; i < NMAGCOL; i++) {
            (void)fits_get_colnum(fptr,CASEINSEN,magcol[i],&colnum,&status);
            if (status == 0) {
                strcpy(colname,magcol[i]);
                break;
            } else {
                colnum = 0;
                status = 0;
            }
        }
        if (colnum != 0) {
            fluxes = cir_malloc(nrows_stds*sizeof(float));
            (void)fits_read_col(fptr,TFLOAT,colnum,1,1,nrows_stds,NULL,fluxes,
                                &anynul,&status);
            sort1(fluxes,nrows_stds);
            cut = fluxes[maxstds-1];
            (void)sprintf(expr,"(%s < %g)",colname,cut);
            (void)fits_select_rows(fptr,fptr,expr,&status);
            (void)fits_get_num_rows(fptr,&nrows_stds,&status);
            (void)fits_close_file(fptr,&status);
        } else {
            (void)fits_close_file(fptr,&status);
        }
    }

    /* See how many rows there are in the input catalogue. If it totally
       overwhelms the number in the standards list, then we need to cut
       things down a bit */
    
    (void)fits_open_file(&fptr,incat,READONLY,&status);
    (void)fits_get_num_rows(fptr,&nrows_cat,&status);
    if (nrows_cat > 500 && nrows_cat > 2.0*nrows_stds) {
    
        /* Look at the fluxes and see at what level we need to cut */
    
	(void)fits_get_colnum(fptr,CASESEN,"Aper_flux_3",&afluxcol,&status);
	fluxes = cir_malloc(nrows_cat*sizeof(float));
	(void)fits_read_col(fptr,TFLOAT,afluxcol,1,1,nrows_cat,NULL,fluxes,
			    &anynul,&status);
	(void)fits_close_file(fptr,&status);
	sort1(fluxes,nrows_cat);
	level = (2*nrows_stds > 500 ? 2*nrows_stds : 500);
	level = max(1,min(level,5000));
	ind = max(1,nrows_cat-level-1);
        cut = fluxes[ind];
	freespace(fluxes);
        (void)fits_parse_rootname(incat,catroot,&status);
	(void)sprintf(newcat,"newcat_%s",basename(catroot));
	if (access(newcat,F_OK) == 0)
	    remove(newcat);
	(void)fits_open_file(&fptr,catroot,READONLY,&status);
	(void)fits_create_file(&nptr,newcat,&status);
	(void)fits_copy_hdu(fptr,nptr,0,&status);
	(void)fits_close_file(fptr,&status);
	(void)fits_open_file(&fptr,incat,READONLY,&status);
	(void)fits_copy_hdu(fptr,nptr,0,&status);
	(void)fits_close_file(fptr,&status);
	(void)sprintf(expr,"(Aper_flux_3 > %g && Ellipticity < 0.5)",cut);
	(void)fits_select_rows(nptr,nptr,expr,&status);
	(void)fits_close_file(nptr,&status);
	(void)sprintf(newcat,"newcat_%s[1]",basename(catroot));
    } else {
	closefits(fptr);
    }
}

static void update_wcskeys(char *infile, char *incat) {
    int status;
    fitsfile *fptr;
    char *testme = "TCRVL*",card[81];

    /* First test to see whether you want image or table WCS keywords. This
       is a historical difference between wfcam and vircam images as the
       table keywords weren't invented by the time largescale processing
       of wfcam data started. */
    
    status = 0;
    (void)fits_open_file(&fptr,incat,READWRITE,&status);
    (void)fits_find_nextkey(fptr,&testme,1,NULL,0,card,&status);
    if (status == 0) {
        (void)fits_close_file(fptr,&status);
	update_wcskeys_table(infile,incat);
    } else {
	status = 0;
	(void)fits_close_file(fptr,&status);
	update_wcskeys_image(infile,incat);
    }
}

static void update_wcskeys_image(char *infile, char *incat) {
    fitsfile *iptr,*tptr;
    int status,i;
    char card[81];
    char *keys[20] = {"CRVAL1","CRVAL2","CRPIX1","CRPIX2","CD1_1",
		      "CD1_2","CD2_1","CD2_2","CUNIT1","CUNIT2",
		      "CTYPE1","CTYPE2","PV2_1","PV2_2","PV2_3",
		      "PV2_4","PV2_5","WCSPASS","NUMBRMS","STDCRMS"};

    /* Open both of the files */

    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_open_file(&tptr,incat,READWRITE,&status);

    /* Read the key value from the header of the image and transfer the 
       whole card over to the table header */

    for (i = 0; i < 20; i++) {
	status = 0;
	(void)fits_read_card(iptr,keys[i],card,&status);
	if (status != 0)
	    continue;
	(void)fits_update_card(tptr,keys[i],card,&status);
    }

    /* Close up and get out of here */

    closefits(iptr);
    closefits(tptr);
}

static void update_wcskeys_table(char *infile, char *incat) {
    fitsfile *iptr,*tptr;
    int status,cols[2],i,j;
    float val;
    double dval;
    char card[81],key1[9],key2[9];
    char *keys[3] = {"WCSPASS","NUMBRMS","STDCRMS"};

    /* Open both of the files */

    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_open_file(&tptr,incat,READWRITE,&status);	

    /* Determine where the x,y coordinate columns are */

    (void)fits_get_colnum(tptr,CASEINSEN,"X_Coordinate",cols,&status);
    (void)fits_get_colnum(tptr,CASEINSEN,"Y_Coordinate",cols+1,&status);
    if (status != 0) {
	closefits(iptr);
	closefits(tptr);
	return;
    }

    /* Update the CRVALn headers */

    status = 0;
    for (i = 1; i <= 2; i++) {
	(void)sprintf(key1,"CRVAL%d",i);
	(void)sprintf(key2,"TCRVL%d",cols[i-1]);
	(void)fits_read_key(iptr,TDOUBLE,key1,&dval,NULL,&status);
	(void)fits_update_key(tptr,TDOUBLE,key2,&dval,NULL,&status);
    }

    /* Update the CRPIXn headers */

    status = 0;
    for (i = 1; i <= 2; i++) {
	(void)sprintf(key1,"CRPIX%d",i);
	(void)sprintf(key2,"TCRPX%d",cols[i-1]);
	(void)fits_read_key(iptr,TFLOAT,key1,&val,NULL,&status);
	(void)fits_update_key(tptr,TFLOAT,key2,&val,NULL,&status);
    }
    
    /* Now the CD matrix */

    status = 0;
    for (j = 1; j <= 2; j++) {
	for (i = 1; i <= 2; i++) {
	    (void)sprintf(key1,"CD%d_%d",j,i);
	    (void)sprintf(key2,"TC%d_%d",cols[j-1],cols[i-1]);
  	    (void)fits_read_key(iptr,TDOUBLE,key1,&dval,NULL,&status);
	    (void)fits_update_key(tptr,TDOUBLE,key2,&dval,NULL,&status);
	}
    }

    /* The PV2_n vector */

    status = 0;
    for (i = 1; i <= 5; i++) {
	(void)sprintf(key1,"PV2_%d",i);
	(void)sprintf(key2,"TV%d_%d",cols[1],i);
	status = 0;
	(void)fits_read_key(iptr,TFLOAT,key1,&val,NULL,&status);
	if (status != 0)
	    continue;
	(void)fits_update_key(tptr,TFLOAT,key2,&val,NULL,&status);
    }

    /* CTYPEn vector */

    status = 0;
    for (i = 1; i <= 2; i++) {
	(void)sprintf(key1,"CTYPE%d",i);
	(void)sprintf(key2,"TCTYP%d",cols[i-1]);
	(void)fits_read_key(iptr,TSTRING,key1,card,NULL,&status);
	(void)fits_update_key(tptr,TSTRING,key2,card,NULL,&status);
    }

    /* CUNITn vector */
    
    status = 0;
    for (i = 1; i <= 2; i++) {
	(void)sprintf(key1,"CUNIT%d",i);
	(void)sprintf(key2,"TCUNI%d",cols[i-1]);
	(void)fits_read_key(iptr,TSTRING,key1,card,NULL,&status);
	if (status != 0) {
	    status = 0;
	    continue;
	}
	(void)fits_update_key(tptr,TSTRING,key2,card,NULL,&status);
    }

    /* Now just do a straight copy of the fit quality keywords */

    for (i = 0; i < 3; i++) {
	status = 0;
	(void)fits_read_card(iptr,keys[i],card,&status);
	if (status != 0)
	    continue;
	(void)fits_update_card(tptr,keys[i],card,&status);
    }

    /* Close up and get out of here */

    closefits(iptr);
    closefits(tptr);
}
			  
static void sort1(float *a, int n) {
    int iii,ii,i,ifin,j;
    float b;

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
            if (a[i] > a[j]) {
                b = a[j];
                while (1) {
                    a[j] = a[i];
                    j = i;
                    i = i - iii;
                    if (i < 0 || a[i] <= b) 
                        break;
                }
                a[j] = b;
            }
        }
    }
}

/*

$Log: wcsfit.c,v $
Revision 1.4  2014/11/03 08:23:34  jim
Now does a magnitude cut on the standard star catalogue if there are too many

Revision 1.3  2014/06/10 08:57:16  jim
Modified so that if the input catalogue has no rows you get a non fatal
error

Revision 1.2  2010/09/06 09:05:59  jim
Tidied some docs


*/
