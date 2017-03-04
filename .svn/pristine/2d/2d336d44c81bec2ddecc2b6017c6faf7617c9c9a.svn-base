/* 

$Id: create_table.c,v 1.4 2012/01/04 10:17:25 jim Exp $

*/

#include <stdio.h>
#include <string.h>
#include "imcore.h"
#include "errcodes.h"
#include "imcore_version.h"

#define BINNAME "APM-BINARYTABLE"

static float *ellipt = NULL;
static float *pkht = NULL;
static float **areal = NULL;
static float *work = NULL;
static fitsfile* iptr = NULL;

static void tabtag(void);
static void tidy();

extern int tabinit(ap_t *ap, char *infile, char *outtab, char *errstr) {
    int status;

    switch (cattype) {
    case CAT_INTWFC:
	status = tabinit_1(infile,outtab,errstr);
	break;
    case CAT_WFCAM:
	status = tabinit_2(infile,outtab,errstr);
	break;
    case CAT_BASIC:
	status = tabinit_3(infile,outtab,errstr);
	break;
    case CAT_OBJMASK:
	status = tabinit_4(ap,infile,outtab,errstr);
	break;
    case CAT_VIRCAM:
	status = tabinit_6(infile,outtab,errstr);
	break;
    default:
	status = ERRCODE_UNDEF_OPT;
	sprintf(errstr,"tabinit: option %d does not exist",cattype);
	return(status);
    }
    if (status == ERRCODE_OK)
	tabtag();
    return(status);
}

extern int tabclose(ap_t *ap, char *errstr) {
    int status;

    switch (cattype) {
    case CAT_INTWFC:
	status = tabclose_1(ap,errstr);
	break;
    case CAT_WFCAM:
	status = tabclose_2(ap,errstr);
	break;
    case CAT_BASIC:
	status = tabclose_3(ap,errstr);
	break;
    case CAT_OBJMASK:
	status = tabclose_4(ap,errstr);
	break;
    case CAT_VIRCAM:
	status = tabclose_6(ap,errstr);
	break;
    default:
	status = ERRCODE_UNDEF_OPT;
	sprintf(errstr,"tabclose: option %d does not exist",cattype);
	break;
    }
    return(status);
}

extern int do_seeing(ap_t *ap, char *errstr) {
    int status;

    switch (cattype) {
    case CAT_INTWFC:
	status = do_seeing_1(ap,errstr);
	break;
    case CAT_WFCAM:
	status = do_seeing_2(ap,errstr);
	break;
    case CAT_BASIC:
	status = do_seeing_3(ap,errstr);
	break;
    case CAT_OBJMASK:
	status = do_seeing_4(ap,errstr);
	break;
    case CAT_VIRCAM:
	status = do_seeing_6(ap,errstr);
	break;
    default:
	status = ERRCODE_UNDEF_OPT;
	sprintf(errstr,"do_seeing: option %d does not exist",cattype);
	break;
    }
    return(status);
}

extern int readtab_list(char *infile, objstruct **ob, int *nr, int ctype, 
			char *errstr) {
    int status;

    switch (ctype) {
    case CAT_INTWFC:
	status = readtab_list_1(infile,ob,nr,errstr);
	break;
    case CAT_WFCAM:
	status = readtab_list_2(infile,ob,nr,errstr);
	break;
    case CAT_BASIC:
	status = readtab_list_3(infile,ob,nr,errstr);
	break;
    case CAT_OBJMASK:
	status = ERRCODE_UNDEF_OPT;
	*nr = 0;
	*ob = NULL;
	(void)sprintf(errstr,
		      "readtab_list: requested table is an object mask");
	break;
    case CAT_VIRCAM:
	status = readtab_list_6(infile,ob,nr,errstr);
	break;
    default:
	status = ERRCODE_UNDEF_OPT;
	*nr = 0;
	*ob = NULL;
	(void)sprintf(errstr,"readtab_list: option %d does not exist",cattype);
	break;
    }
    return(status);
}
	
extern int process_results(ap_t *ap, char *errstr) {
    int status;
    
    switch (cattype) {
    case CAT_INTWFC:
	status = process_results_1(ap,errstr);
	break;
    case CAT_WFCAM:
	status = process_results_2(ap,errstr);
	break;
    case CAT_BASIC:
	status = process_results_3(ap,errstr);
	break;
    case CAT_OBJMASK:
	status = process_results_4(ap,errstr);
	break;
    case CAT_VIRCAM:
	status = process_results_6(ap,errstr);
	break;
    default:
	status = ERRCODE_UNDEF_OPT;
	sprintf(errstr,"process_results: option %d does not exist",cattype);
	break;
    }
    return(status);
}

extern int process_results_list(ap_t *ap, objstruct *oblist[IMNUM], int nbit, 
				char *errstr) {
    int status;

    switch (cattype) {
    case CAT_INTWFC:
	status = process_results_list_1(ap,oblist,nbit,errstr);
	break;
    case CAT_WFCAM:
	status = process_results_list_2(ap,oblist,nbit,errstr);
	break;
    case CAT_BASIC:
/* 	status = process_results_list_3(ap,oblist,nbit,errstr); */
	status = ERRCODE_OK; 
	break;
    case CAT_VIRCAM:
	status = process_results_list_6(ap,oblist,nbit,errstr);
	break;
    default:
	status = ERRCODE_UNDEF_OPT;
	sprintf(errstr,"process_results_list: option %d does not exist",
		cattype);
	break;
    }
    return(status);
}

extern int get_ncol(int testcat) {
    int nc;

    switch (testcat) {
    case CAT_INTWFC:
	nc = get_ncol_1();
	break;
    case CAT_WFCAM:
	nc = get_ncol_2();
	break;
    case CAT_BASIC:
	nc = get_ncol_3();
	break;
    case CAT_OBJMASK:
	nc = get_ncol_4();
	break;
    case CAT_VIRCAM:
	nc = get_ncol_6();
	break;
    default:
	nc = 0;
	break;
    }
    return(nc);
}

extern int tabinit_gen(char *infile, char *outtab, int ncols, char *ttype[],
    char *tunit[], char *tform[], char *errstr) {
    int i,status,nhdus,hdunum,hdutype,nkeys,nmore;
    long naxes[2],n;
    char card[FLEN_CARD],errmsg[FLEN_STATUS],val[FLEN_CARD],comment[FLEN_CARD];

    /* First decide whether this is a new table output file or whether we're
       just creating a new table extension. If this is a simple FITS image
       or if this is the first extension, then create a new file. Otherwise
       just tack on a new fits table. NB: this scenario assumes that there
       will be a catalogue generated for each image in the MEF */
 
    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(iptr,&nhdus,&status);
    if (nhdus > 1) {
        (void)fits_get_hdu_num(iptr,&hdunum);
        hdunum--;
    } else {
        hdunum = 1;
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to get extension information: %s (status = %d)",
                errmsg,status);
	closefits(iptr);
        return(ERRCODE_FILE_IO);
    }
    (void)fits_open_file(&tptr,outtab,READWRITE,&status);
    if (hdunum == 1 || nhdus == 1 || status != 0) {
	status = 0;
        (void)fits_create_file(&tptr,outtab,&status);
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to create/open FITS container: %s (status = %d)",
                errmsg,status);
	closefits(iptr);
        return(ERRCODE_FILE_IO);
    }
 
    /* Create a new table extension */
 
    (void)fits_create_tbl(tptr,BINARY_TBL,0,ncols,ttype,tform,tunit,BINNAME,
			  &status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to create new FITS table: %s (status = %d)",
                errmsg,status);
	closefits(iptr);
        return(ERRCODE_FILE_IO);
    }

    /* Add a date string */

    (void)fits_write_date(tptr,&status);
 
   /* Copy header cards from the PHU of the input file to the table HU. Make
       sure these aren't anything to do with the structure of the image data
       or compression */

    (void)fits_movabs_hdu(iptr,1,&hdutype,&status);
    (void)fits_get_hdrspace(iptr,&nkeys,&nmore,&status);
    for (i = 1; i <= nkeys; i++) {
        (void)fits_read_record(iptr,i,card,&status);
	if (! strncmp(card,"DATE    =",9) || ! strncmp(card,"MAGZ",4) ||
	    ! strncmp(card,"EXTINCT =",9) || ! strncmp(card,"VSA_",4) ||
	    ! strncmp(card,"IAMFUD",6))
	    continue;
        if (fits_get_keyclass(card) > TYP_SCAL_KEY)
            (void)fits_write_record(tptr,card,&status);
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to add image PHU header cards: %s (status = %d)",
                errmsg,status);
	closefits(iptr);
        return(ERRCODE_FILE_IO);
    }

    /* If this isn't a simple FITS input image, then copy the relevant cards
       from the image extension to the catalogue header */

    if (nhdus > 1) {
        (void)fits_movabs_hdu(iptr,hdunum+1,&hdutype,&status);
        (void)fits_write_comment(tptr,"Copy of header items in 2D map file............",
                                 &status);
        (void)fits_get_hdrspace(iptr,&nkeys,&nmore,&status);
	(void)fits_read_key(iptr,TSTRING,"EXTNAME",val,comment,&status);
	if (status == 0) 
	    (void)fits_update_key(tptr,TSTRING,"EXTNAME",val,comment,&status);
	else
	    status = 0;
        for (i = 1; i <= nkeys; i++) {
	    status = 0;
            (void)fits_read_record(iptr,i,card,&status);
	    if (strncmp(card,"EXTNAME",7) == 0)
		continue;
            if (fits_get_keyclass(card) > TYP_SCAL_KEY)
                (void)fits_write_record(tptr,card,&status);
        }
        if (status != 0) {
            fits_get_errstatus(status,errmsg);
            sprintf(errstr,"Unable to add image EHU header cards: %s (status = %d)",
                    errmsg,status);
  	    closefits(iptr);
            return(ERRCODE_FILE_IO);
        }
    }

    /* Write the original image size to the table header */

    (void)fits_get_img_size(iptr,2,naxes,&status);
    n = naxes[0];
    (void)fits_update_key(tptr,TINT,"NXOUT",&n,"X-axis size of image",
			  &status);
    n = naxes[1];
    (void)fits_update_key(tptr,TINT,"NYOUT",&n,"Y-axis size of image",
			  &status);
    if (status != 0) {
	fits_get_errstatus(status,errmsg);
	sprintf(errstr,"Unable to add image data size to table header: %s (status = %d)",
		errmsg,status);
	closefits(iptr);
	return(ERRCODE_FILE_IO);
    }

    /* Write the imcore version string */

    (void)sprintf(card,"IMCORE version: %s",imcore_version);
    (void)fits_write_history(tptr,card,&status);

    /* Get out of here */

    closefits(iptr);
    return(ERRCODE_OK);
}

extern int tabinit_genmef(char *infile, char *outtab, int ncols, char *ttype[],
			  char *tunit[], char *tform[], char *errstr) {
    int i,status,nhdus,hdunum,hdutype,nkeys,nmore,newfile;
    long naxes[2],n;
    char card[FLEN_CARD],errmsg[FLEN_STATUS],val[FLEN_CARD],comment[FLEN_CARD];

    /* First decide whether this is a new table output file or whether we're
       just creating a new table extension. If this is a simple FITS image
       or if this is the first extension, then create a new file. Otherwise
       just tack on a new fits table. NB: this scenario assumes that there
       will be a catalogue generated for each image in the MEF */
 
    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(iptr,&nhdus,&status);
    if (nhdus > 1) {
        (void)fits_get_hdu_num(iptr,&hdunum);
        hdunum--;
    } else {
        hdunum = 1;
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to get extension information: %s (status = %d)",
                errmsg,status);
	closefits(iptr);
        return(ERRCODE_FILE_IO);
    }

    /* Try to open an existing file. If it exists then we just need to
       tack on a new binary table extension. If it doesn't then we need
       to create a whole new file */

    (void)fits_open_file(&tptr,outtab,READWRITE,&status);
    if (status != 0) {
	status = 0;
        (void)fits_create_file(&tptr,outtab,&status);
	newfile = 1;

    /* Situation where we have a simple fits file or a MEF and we're just
       looking at the first extension. In that case we need to delete the 
       old file and start afresh */
           
    } else if (status == 0 && (nhdus == 1 || hdunum == 1)) {
	(void)fits_delete_file(tptr,&status);
        (void)fits_create_file(&tptr,outtab,&status);
	newfile = 1;
    } else {
	newfile = 0;
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to create/open FITS container: %s (status = %d)",
                errmsg,status);
	closefits(iptr);
        return(ERRCODE_FILE_IO);
    }

    /* If this is a new output file, then copy the image primary header to the
       primary header of the output file */

    if (newfile) {
	(void)fits_create_img(tptr,BYTE_IMG,0,NULL,&status);
	(void)fits_movabs_hdu(iptr,1,&hdutype,&status);
	(void)fits_get_hdrspace(iptr,&nkeys,&nmore,&status);
	for (i = 1; i <= nkeys; i++) {
	    (void)fits_read_record(iptr,i,card,&status);
	    if (! strncmp(card,"DATE    =",9) || ! strncmp(card,"MAGZ",4) ||
		! strncmp(card,"EXTINCT =",9) || ! strncmp(card,"VSA_",4) ||
		! strncmp(card,"IAMFUD",6))
		continue;
	    if (fits_get_keyclass(card) > TYP_SCAL_KEY)
		(void)fits_write_record(tptr,card,&status);
	}
	if (status != 0) {
	    fits_get_errstatus(status,errmsg);
	    sprintf(errstr,"Unable to add image PHU header cards: %s (status = %d)",
		    errmsg,status);
	    closefits(iptr);
	    return(ERRCODE_FILE_IO);
	}
    }

    /* Add a date string */

    (void)fits_write_date(tptr,&status);
 
    /* Create a new table extension */
 
    (void)fits_create_tbl(tptr,BINARY_TBL,0,ncols,ttype,tform,tunit,BINNAME,
                            &status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to create new FITS table: %s (status = %d)",
                errmsg,status);
	closefits(iptr);
        return(ERRCODE_FILE_IO);
    }

    /* If this isn't a simple FITS input image, then copy the relevant cards
       from the image extension to the catalogue header. If it is a simple
       fits file then just copy the one available header to the catalogue */

    if (nhdus > 1) {
        (void)fits_movabs_hdu(iptr,hdunum+1,&hdutype,&status);
        (void)fits_write_comment(tptr,"Copy of header items in 2D map file............",
                                 &status);
        (void)fits_get_hdrspace(iptr,&nkeys,&nmore,&status);
	(void)fits_read_key(iptr,TSTRING,"EXTNAME",val,comment,&status);
	if (status == 0) 
	    (void)fits_update_key(tptr,TSTRING,"EXTNAME",val,comment,&status);
	else
	    status = 0;
    }
    for (i = 1; i <= nkeys; i++) {
        (void)fits_read_record(iptr,i,card,&status);
	if (strncmp(card,"EXTNAME",7) == 0)
	    continue;
        if (fits_get_keyclass(card) > TYP_SCAL_KEY)
            (void)fits_write_record(tptr,card,&status);
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to add image EHU header cards: %s (status = %d)",
                errmsg,status);
        closefits(iptr);
        return(ERRCODE_FILE_IO);
    }

    /* Add a date string */

    (void)fits_write_date(tptr,&status);
 
    /* Write the original image size to the table header */

    (void)fits_get_img_size(iptr,2,naxes,&status);
    n = naxes[0];
    (void)fits_update_key(tptr,TINT,"NXOUT",&n,"X-axis size of image",
			  &status);
    n = naxes[1];
    (void)fits_update_key(tptr,TINT,"NYOUT",&n,"Y-axis size of image",
			  &status);
    if (status != 0) {
	fits_get_errstatus(status,errmsg);
	sprintf(errstr,"Unable to add image data size to table header: %s (status = %d)",
		errmsg,status);
	closefits(iptr);
	return(ERRCODE_FILE_IO);
    }

    /* Write the imcore version string */

    (void)sprintf(card,"IMCORE version: %s",imcore_version);
    (void)fits_write_history(tptr,card,&status);

    /* Get out of here */

    closefits(iptr);
    return(ERRCODE_OK);
}

extern int do_seeing_gen(ap_t *ap, int col_ellipt, int col_peakheight,
			 int col_areals[NAREAL], char *errstr) {
    int i,status,anynul;
    long nrows;
    float fwhm;
    char errmsg[FLEN_STATUS];

    /* Get some space to do the work... */

    status = 0;
    (void)fits_get_num_rows(tptr,&nrows,&status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to get nrows in table: %s (status = %d)",
		errmsg,status);
	return(ERRCODE_FILE_IO);
    }
    if (nrows >= 3) {
	ellipt = malloc(nrows*sizeof(*ellipt));
	pkht = malloc(nrows*sizeof(*pkht));
	work = malloc(nrows*sizeof(*work));
	areal = malloc(NAREAL*sizeof(float *));
	if (! ellipt || ! pkht || ! work || ! areal) {
	    sprintf(errstr,"Unable to get workspace in create_table_1\n");
	    tidy();
	    return(ERRCODE_MALLOC);
	}
	for (i = 0; i < NAREAL; i++) {
	    areal[i] = malloc(nrows*sizeof(float));
	    if (! areal[i]) {
		sprintf(errstr,"Unable to get workspace in create_table_1\n");
		tidy();
		return(ERRCODE_MALLOC);
	    }
	}

	/* Now read the info from the table */

	(void)fits_read_col(tptr,TFLOAT,col_ellipt,1,1,nrows,NULL,ellipt,
			    &anynul,&status);
	(void)fits_read_col(tptr,TFLOAT,col_peakheight,1,1,nrows,NULL,pkht,
			    &anynul,&status);
	for (i = 0; i < NAREAL; i++)
	    (void)fits_read_col(tptr,TFLOAT,col_areals[i],1,1,nrows,NULL,
				areal[i],&anynul,&status);
	if (status != 0) {
	    fits_get_errstatus(status,errmsg);
	    sprintf(errstr,"Unable to read columns in table: %s (status = %d)",
		    errmsg,status);
	    tidy();
	    return(ERRCODE_FILE_IO);
	}

	/* Do the seeing calculation */

	seeing(ap,nrows,ellipt,pkht,areal,work,&fwhm);
    } else 
	fwhm = 0.0;

    /* Now write info to header about sky and seeing */

    (void)fits_update_key_fixflt(tptr,"SKYLEVEL",ap->background,2,
				 "Median sky brightness (counts/pixel)",&status);
    (void)fits_update_key_fixflt(tptr,"SKYNOISE",ap->sigma,2,
				 "Pixel noise at sky level (counts)",&status);
    (void)fits_update_key_fixflt(tptr,"THRESHOL",ap->thresh,2,
				 "Isophotal analysis threshold (counts)",
				 &status);
    (void)fits_update_key_lng(tptr,"MINPIX",ap->ipnop,
			      "Minimum size for images (pixels)",&status);
    (void)fits_update_key_lng(tptr,"CROWDED",ap->icrowd,
			      "Crowded field analysis flag (0 none, 1 active)",
			      &status);
    (void)fits_update_key_fixflt(tptr,"RCORE",ap->rcore,2,
				 "[pixels] Core radius for default profile fit",
				 &status);
    (void)fits_update_key_fixflt(tptr,"FILTFWHM",ap->filtfwhm,2,
				 "[pixels] FWHM of smoothing kernel",
				 &status);
    (void)fits_update_key_fixflt(tptr,"SEEING",fwhm,2,"[pixels] Average FWHM",
				 &status);
    (void)fits_update_key_lng(tptr,"NBSIZE",ap->nbsize,
 			      "Smoothing box size for background",&status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to update table header: %s (status = %d)",
		errmsg,status);
	tidy();
	return(ERRCODE_FILE_IO);
    }

    /* Get out of here */

    tidy();
    return(ERRCODE_OK);
}

extern int tabclose_gen(ap_t *ap, char *errstr) {
    int status;
 
    status = 0;
    closefits(tptr);
    return(ERRCODE_OK);
}

static void tabtag(void) {
    int status,nhdu,hdutype;

    status = 0;
    (void)fits_get_hdu_num(tptr,&nhdu);
    (void)fits_movabs_hdu(tptr,1,&hdutype,&status);
    (void)fits_update_key(tptr,TINT,"CATTYPE",&cattype,"Imcore catalogue type",
			  &status);
    (void)fits_movabs_hdu(tptr,nhdu,&hdutype,&status);
}

static void tidy () {
    int i,status;

    /* Free up workspace */

    freespace(ellipt);
    freespace(pkht);
    freespace(work);
    if (areal != NULL) 
        for (i = 0; i < NAREAL; i++)
	    freespace(areal[i]);
    freespace(areal);
    closefits(iptr);
}

/*

$Log: create_table.c,v $
Revision 1.4  2012/01/04 10:17:25  jim
Added routine tagtab

Revision 1.3  2010-09-06 08:57:37  jim
Now writes NBSIZE to output headers

Revision 1.2  2010-08-02 16:03:28  jim
Fixed so that tables create from simple fits images inherit the image header
in both the primary and table extension

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.18  2010-07-19 11:05:03  jim
Fixed multiple EXTNAME problem

Revision 1.17  2010/05/05 08:23:22  jim
Traps for existence of IAMFUD in the primary so that it isn't added to the
extension headers

Revision 1.16  2010/02/05 10:04:57  jim
Added line to explicitly remove reference to header cards that start with VSA_

Revision 1.15  2010/01/19 12:01:47  jim
removed reference to colour stuff

Revision 1.14  2009/09/24 13:49:52  jim
Modified for create_table_6 issues

Revision 1.13  2009/01/29 14:02:11  jim
fix to get rid of boring compiler warning

Revision 1.12  2009/01/29 12:14:06  jim
Adds imcore_version into header

Revision 1.11  2008/04/15 18:55:25  jim
Added calls to _5 routines

Revision 1.10  2007/06/04 10:34:02  jim
Modified to add list driven routines

Revision 1.9  2006/08/21 09:06:37  jim
Modified so that MAG*, EXTINCT and BSCALE don't get copied into table
extensions

Revision 1.8  2006/07/31 13:21:09  jim
Modified imcore now allows for a smoothing kernel with variable FWHM

Revision 1.7  2005/09/08 13:15:26  jim
added code to put the image data array size into catalogue header

Revision 1.6  2005/05/09 08:43:52  jim
Fixed so that you don't accidently copy the date keyword from the input image

Revision 1.5  2004/09/07 14:18:56  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.4  2004/04/14 08:38:27  jim
Fixed bug in tabinit_gen where the wrong extension header was being copied to
the output table

Revision 1.3  2004/04/08 08:26:01  jim
Fixed memory deallocation bug in tidy routine

Revision 1.2  2004/04/07 13:34:29  jim
Fixed bug that caused a failure when a table is created for anything but
the first image extension in the input file

Revision 1.1  2004/04/02 10:54:57  jim
New version for rewrite of imcore


*/
