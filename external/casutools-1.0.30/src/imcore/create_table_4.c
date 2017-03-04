/*

$Id: create_table_4.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/

#include <stdio.h>
#include <math.h>
#include "imcore.h"
#include "errcodes.h"
#include "util.h"

static fitsfile *iptr = NULL;

extern int tabinit_4(ap_t *ap, char *infile, char *outtab, char *errstr) {
    int status,nhdus,hdunum,hdutype,nkeys,nmore,i,hdunum2,bitpix,naxis;
    int newout;
    long naxes[2];
    char errmsg[FLEN_STATUS],card[FLEN_CARD];
    unsigned char *opm;

    /* First decide whether this is a new table output file or whether we're
       just creating a new image extension. If this is a simple FITS image
       or if this is the first extension, then create a new file. Otherwise
       just tack on a new fits image. NB: this scenario assumes that there
       will be a mask generated for each image in the MEF */
  
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

    /* See if the output file exists first. If it doesn't then signal that
       we need to create a new one */

    (void)fits_open_file(&tptr,outtab,READONLY,&status);
    newout = (status != 0);
    (void)fits_close_file(tptr,&status);
    status = 0;

    /* Now create the file if necessary or just open it again */

    if (hdunum == 1 || nhdus == 1 || newout) {
        (void)fits_create_file(&tptr,outtab,&status);
	if (nhdus > 1) {
            (void)fits_get_hdu_num(iptr,&hdunum2);
	    (void)fits_movabs_hdu(iptr,1,&hdutype,&status);
	    (void)fits_get_hdrspace(iptr,&nkeys,&nmore,&status);
	    (void)fits_create_img(tptr,BYTE_IMG,0,naxes,&status);
	    for (i = 1; i <= nkeys; i++) {
		(void)fits_read_record(iptr,i,card,&status);
		if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
		    (void)fits_write_record(tptr,card,&status);
	    }
	    if (status != 0) {
		fits_get_errstatus(status,errmsg);
		sprintf(errstr,"Unable to add image PHU header cards: %s (status = %d)",                errmsg,status);
		closefits(iptr);
		return(ERRCODE_FILE_IO);
	    }
	    (void)fits_movabs_hdu(iptr,hdunum2,&hdutype,&status);
	}
    } else
        (void)fits_open_file(&tptr,outtab,READWRITE,&status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to create/open FITS container: %s (status = %d)",
		errmsg,status);
        closefits(iptr);
        return(ERRCODE_FILE_IO);
    }

    /* Get the image parameters for the input file */

    (void)fits_get_img_param(iptr,2,&bitpix,&naxis,naxes,&status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to read FITS image params: %s (status = %d)",
		errmsg,status);
        closefits(iptr);
        return(ERRCODE_FILE_IO);
    }
    
    /* Now create a similar image in the mask output file */

    bitpix = BYTE_IMG;
    (void)fits_create_img(tptr,bitpix,naxis,naxes,&status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to create new image extension: %s (status = %d)",
		errmsg,status);
        closefits(iptr);
        return(ERRCODE_FILE_IO);
    }
	
    /* Copy the relevant cards from the image extension to the mask header
       If this is a simple fits image, then this will just be the primary */

    (void)fits_movabs_hdu(iptr,hdunum,&hdutype,&status);
    (void)fits_get_hdrspace(iptr,&nkeys,&nmore,&status);
    for (i = 1; i <= nkeys; i++) {
	(void)fits_read_record(iptr,i,card,&status);
	if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
	    (void)fits_write_record(tptr,card,&status);
    }
    if (status != 0) {
	fits_get_errstatus(status,errmsg);
	sprintf(errstr,"Unable to add image EHU header cards: %s (status = %d)",
		errmsg,status);
	closefits(iptr);
	return(ERRCODE_FILE_IO);
    }

    /* Get some memory to hold the object mask */

    opm = calloc(naxes[0]*naxes[1],sizeof(*opm));
    ap->opm = opm;

    /* Get out of here */

    closefits(iptr);
    return(ERRCODE_OK);
}

extern int do_seeing_4(ap_t *ap, char *errstr) {

    /* Get out of here */

    return(ERRCODE_OK);
}
        

extern int process_results_4(ap_t *ap, char *errstr) {
    int i,j,np,iareal[NAREAL];
    long nx;
    float momresults[8];
    plstruct *plarray;
    unsigned char *opm;

    /* Do a basic moments analysis and work out the areal profiles. This
       is only done to check and make sure this object _should_ be flagged */

    moments(ap,momresults);
    
    if (momresults[0] < 0)
	return(ERRCODE_BADRES);
    areals(ap,iareal);

    /* See if this object makes the cut in terms of its size.  If not, then
       just return with good status */

    if (iareal[0] < ap->ipnop || momresults[3] < ap->xintmin)
	return(ERRCODE_OK);

    /* Loop for each object in the array */

    plarray = ap->plarray;
    np = ap->npl_pix;
    opm = ap->opm;
    nx = ap->lsiz;
    for (i = 0; i < np; i++) {
	if (plarray[i].z <= 0.0)
	    continue;
	j = nx*(plarray[i].y - 1) + plarray[i].x - 1;
	opm[j] = 1;
    }

    /* Get outta here */

    return(ERRCODE_OK);
}

extern int tabclose_4(ap_t *ap, char *errstr) {
    int status,npts,nobj,i;
    unsigned char *opm;
    short int *conf;
    char errmsg[FLEN_STATUS];

    /* Write the mask image */

    status = 0;
    opm = ap->opm;
    conf = ap->conf;
    npts = (ap->lsiz)*(ap->csiz);
    (void)fits_write_img(tptr,TBYTE,1,npts,opm,&status);
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
        sprintf(errstr,"Unable to write object mask: %s (status = %d)",
                errmsg,status);
        return(ERRCODE_FILE_IO);
    }

    /* Work out how many object pixels there were */

    nobj = 0;
    for (i = 0; i < npts; i++)
	nobj += (opm[i] == 1 && conf[i] != 0);
    (void)fits_update_key(tptr,TINT,"NOBJPIX",&nobj,"Number of object pixels",
			  &status);
    closefits(tptr);
    freespace(opm);
    return(ERRCODE_OK);
}

extern int get_ncol_4(void) {
    return(0);
}

/*

$Log: create_table_4.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.7  2008/04/30 04:47:38  jim
Modified to allow an image to be created with a different extension number
from that of the input image

Revision 1.6  2008/04/25 12:22:47  jim
Modified to take out pixels with zero confidence from the calculation
of the number of object pixels

Revision 1.5  2008/02/20 09:30:11  jim
Now adds NOBJPIX header item to output FITS file

Revision 1.4  2007/07/31 12:02:37  jim
Modified so that only pixels that match the minimum object size requirements
are flagged

Revision 1.3  2007/06/04 10:34:02  jim
Modified to add list driven routines

Revision 1.2  2004/05/05 10:54:53  jim
Fixed bug where masks weren't being written to multi-extension files
correctly

Revision 1.1  2004/04/02 10:54:58  jim
New version for rewrite of imcore


*/
