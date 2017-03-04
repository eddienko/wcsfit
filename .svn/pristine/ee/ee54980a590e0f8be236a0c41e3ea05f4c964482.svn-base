#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <tools.h>
#include <config.h>
#include <math.h>
#include <misc.h>

/*+
 *  Name:
 *      cir_nint
 *
 *  Purpose:
 *      Do a nearest integer evaluation on a value.
 *
 *  Description:
 *      Find the nearest integer to an input value
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      value = float (Given)
 *          The input value.
 *
 *  Returned values:
 *      The nearest integer to the input value.
 *
 *  Notes:
 *      None
 *              
 *  Dependencies:
 *      None
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_nint(float value) {
    int ival;

    if (value < 0.0) 
        ival = (int)(value - 0.5);
    else 
        ival = (int)(value + 0.5);
    return(ival);
}

/*+
 *  Name:
 *      cir_open_output
 *
 *  Purpose:
 *      Open an output FITS file, based on an input file with perhaps
 *      some changes to the file structure.
 *
 *  Description:
 *      Locate objects in a 2d image.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      newfile = char * (Given)
 *          This is the name of the new file to be created.  This could 
 *          either be a new file or a new image extension in an existing
 *          file.
 *      template = char * (Given)
 *          This is the template file on which you want to base your new
 *          output file.  This can either be a file or an image extension
 *          and it must exist.
 *      iptr = fitsfile ** (Returned)
 *          This is the CFITSIO pointer for the new file/extension
 *      recipe = int (Given)
 *          The name of the recipe to be used in creating the file. The
 *          allowable recipes are included in the Notes below.
 *      newbitpix = int (Given)
 *          If the recipe suggests a change in BITPIX from the template to
 *          the new file, then this is the value that the new file will have.
 *      newaxis = int (Given)
 *          If the recipe suggests that the number of axes in the data array
 *          should be different from the template, then this is the new
 *          number of axes.
 *      newaxes = long * (Given)
 *          If the recipe suggests that either the number of axes should change
 *          or the size of an axis should change, then this is the dimensions
 *          of each of the axes in the output image.
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      The parameter 'recipe' may have the following values:
 *          DIRECTCOPY
 *              The output is a direct copy of the input
 *          MAKE2D
 *              The header will be the same, but the new data array will
 *              be 2d. The dimensions of the axes will be taken from the
 *              first two dimensions of the template data array as will
 *              the value of BITPIX.
 *          NEWBITPIX
 *              The header will be the same.  The data array dimensionality
 *              will be the same as the template.  The only change will be
 *              in the data type. The parameter 'newbitpix' will specify
 *              the output data type.
 *          NEWNAXES
 *              The header will be the same.  The number of axes in the data
 *              array will be the same, but the size of each data array will
 *              now be determined by the values of 'newaxes'. Data type will
 *              remain the same.
 *          NEWDATASZ
 *              The header will remain the same.  The data type, number of
 *              axes and the size of each axis will be determined by 
 *              'newbitpix', 'newaxis' and 'newaxes' respectively.  
 *
 *      These recipes are not supposed to be an exhaustive list of the 
 *      possibilities. Rather they are just what I needed at the moment.
 *              
 *  Dependencies:
 *      cfitsio
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_open_output(char *newfile, char *template, fitsfile **iptr,
			   int recipe, int newbitpix, int newaxis,
			   long *newaxes, char *errmsg) {
    int status,isnew,retval,hdrnum;
    fitsfile *tptr;
    char rootname[BUFSIZ],temproot[BUFSIZ];

    /* See if the file exists in the first place.  If it doesn't then
       create the file.  If it does then just open it */

    status = 0;
    isnew = 0;
    *iptr = NULL;
    (void)fits_parse_rootname(newfile,rootname,&status);
    if (access(rootname,F_OK) != 0) {
        retval = cir_crfile(rootname,iptr,errmsg);
        isnew = 1;
    } else {
        (void)fits_open_file(iptr,rootname,READWRITE,&status);
    }

    /* Ok parse the extension number.  If it's zero, then assume you
       are writing to the PHU, that is you have a simple fits file */

    (void)fits_parse_extnum(newfile,&hdrnum,&status);
    if (hdrnum == -99)
        hdrnum = 1;
    hdrnum--;
    
    /* Ok, here's where it gets complicated...

       hdrnum = 0; isnew = 1 -- This output is a simple fits file
                                so just copy the header from the template
			        using the recipe specified in the
				calling parameters
       hdrnum != 0; isnew = 1 -- This is a MEF that has only just been
                                 created (hdrnum is probably 1).  Create
				 a PHU first and then create the extension
       hdrnum != 0; isnew = 0 -- This is a MEF that already exists and
                                 you're adding a new extension to it

       A lot of things don't get tested here, for example, does the extension
       already exist?  I don't have time for that now... */

    if ((hdrnum == 0 && isnew == 1) || (hdrnum != 0 && isnew == 0)) {
        (void)fits_open_file(&tptr,template,READONLY,&status);
        retval = cir_copy_hdu(*iptr,tptr,recipe,newbitpix,newaxis,newaxes,
			      errmsg);
        (void)fits_close_file(tptr,&status);
    } else if (hdrnum != 0 && isnew == 1) {
	(void)fits_parse_rootname(template,temproot,&status);
        (void)fits_open_file(&tptr,temproot,READONLY,&status);
        retval = cir_copy_hdu(*iptr,tptr,DIRECTCOPY,0,0,NULL,errmsg);
        (void)fits_close_file(tptr,&status);
        (void)fits_open_file(&tptr,template,READONLY,&status);
        retval = cir_copy_hdu(*iptr,tptr,recipe,newbitpix,newaxis,newaxes,
			      errmsg);
        (void)fits_close_file(tptr,&status);
    }	

    return(CIR_OK);
}


extern int cir_crfile(char *newfile, fitsfile **iptr, char *errmsg) {
    int status;

    /* Create an empty fits file */

    status = 0;
    (void)fits_create_file(iptr,newfile,&status);
    if (status != 0) {
        (void)sprintf(errmsg,"Couldn't open file: %s\n",newfile);
        return(CIR_FATAL);
    }

    /* Close up the file and return */

    return(CIR_OK);
}

extern int cir_copy_hdu(fitsfile *iptr, fitsfile *tptr, int recipe, 
			int newbitpix, int newnaxis, long *newnaxes, 
			char *errmsg) {
    int status,ncards,bitpix,naxis,retval;
    long naxes[3];

    /* Copy the header. Assume that the 'move' functions have been done
       before calling this routine.  First find out how many header cards
       there are in the current file template */

    status = 0;
    (void)fits_get_hdrspace(tptr,&ncards,NULL,&status);

    /* Now set this as the total space (plus a little extra for slop) in the
       new file */

    status = 0;
    (void)fits_create_hdu(iptr,&status);
    (void)fits_clear_errmsg();
    status = 0;
    (void)fits_set_hdrsize(iptr,ncards+MOREKEYS,&status);

    /* Right, now find out how big the template image is.  Create an image
       in the new file that's the same size and type, but that's all zeros */

    status = 0;
    switch (recipe) {
    case DIRECTCOPY:
        (void)fits_get_img_param(tptr,3,&bitpix,&naxis,naxes,&status);
	if (naxis == 0)
	    bitpix = BYTE_IMG;
        (void)fits_create_img(iptr,bitpix,naxis,naxes,&status);
        break;
    case MAKE2D:
	(void)fits_get_img_param(tptr,3,&bitpix,&naxis,naxes,&status);
        (void)fits_create_img(iptr,bitpix,2,naxes,&status);
        break;
    case NEWBITPIX:
	(void)fits_get_img_param(tptr,3,&bitpix,&naxis,naxes,&status);
	(void)fits_create_img(iptr,newbitpix,naxis,naxes,&status);
        break;
    case NEWNAXES:
	(void)fits_get_img_param(tptr,3,&bitpix,&naxis,naxes,&status);
	(void)fits_create_img(iptr,bitpix,naxis,newnaxes,&status);
        break;
    case NEWDATASZ:
	(void)fits_create_img(iptr,newbitpix,newnaxis,newnaxes,&status);
        break;
    }

    /* Now copy over all the header records that are not 'data specific' */

    retval = cir_copy_cards(iptr,tptr);

    /* Get out of here */

    return(CIR_OK);
}

extern int cir_copy_cards(fitsfile *iptr, fitsfile *tptr) {
    int i,status,ncards,bitpix,bzcards,keyclass,naxis,bitpix2;
    long naxes[3];
    float bscale,bzero;
    char card[81],junk[64];

    /* How many cards are there in the template? */

    status = 0;
    (void)fits_get_hdrspace(tptr,&ncards,NULL,&status);

    /* Flush just to be on the safe side */

    (void)fits_flush_file(iptr,&status);

    /* If the output file is not float or double then you want to copy
       the bscale and bzero values */

    (void)fits_get_img_param(iptr,3,&bitpix2,&naxis,naxes,&status);
    (void)fits_get_img_param(tptr,3,&bitpix,&naxis,naxes,&status);
    bzcards = (bitpix != FLOAT_IMG && bitpix != DOUBLE_IMG && 
	       bitpix2 != FLOAT_IMG && bitpix2 != DOUBLE_IMG) ? 1 : 0;

    /* Copy all cards from the template to the input HDU, except those that
       have anything to do with data */

    for (i = 1; i <= ncards; i++) {
        (void)fits_read_record(tptr,i,card,&status);
	keyclass = fits_get_keyclass(card);
        if (keyclass > TYP_CKSUM_KEY || (bzcards && keyclass == TYP_SCAL_KEY))
            (void)fits_write_record(iptr,card,&status);
    }
    if (bzcards && bitpix != BYTE_IMG && naxis != 0) {
	status = 0;
	(void)fits_read_key(tptr,TFLOAT,"BZERO",&bzero,junk,&status);
	if (status != 0) {
	    bzero = 0.0;
	    status = 0;
	}
	(void)fits_read_key(tptr,TFLOAT,"BSCALE",&bscale,junk,&status);
	if (status != 0) {
	    bscale = 1.0;
	    status = 0;
	}
	(void)fits_set_bscale(iptr,bscale,bzero,&status);
    }

    return(status);
}

/*+
 *  Name:
 *      cir_qmedsig
 *
 *  Purpose:
 *      Find a median and standard deviation in a data array.
 *
 *  Description:
 *      An estimate of the median and standard deviation of a data array our
 *      found using a quick algorithm. Solution can be iterated
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      data = float * (Given)
 *          The input data
 *      bpm = unsigned char * (Given)
 *          The input BPM.  If no BPM is to be used, this should be set to NULL
 *      npts = int (Given)
 *          The number of points in the input data arrays.
 *      thresh = float (Given)
 *          Threshold for pixel rejection during iteration.
 *      niter = int (Given)
 *          Number of iterations to be done.
 *      lowv = float (Given)
 *          The lowest allowable data value.
 *      highv = float (Given)
 *          The highest allowable data value.
 *      median = float * (Returned)
 *          The value of the median of the data array.
 *      sigma = float * (Returned)
 *          The value of the sigma of the data array.
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
 *      None.
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_qmedsig(float *data, unsigned char *bpm, int npts, 
		       float thresh, int niter, float lowv, float highv,
		       float *median,float *sigma, char *errmsg) {
    int *histo,nbins,nhist,i,ilev,iclip,nhist2,halflev,quartlev;
    int irej,jst,j;
    float mlev,qlev;

    /* Right, first thing is to histogram the data.  Cut values below
       and above the 'ceiling' values */

    nbins = cir_nint(highv - lowv + 1.0);
    histo = cir_calloc(nbins,sizeof(*histo));
    nhist = 0;
    if (bpm != NULL) {
	for (i = 0; i < npts; i++) {
	    if (bpm[i] || data[i] < lowv || data[i] > highv)
		continue;
	    ilev = cir_nint(data[i] - lowv);
	    ilev = max(0,min(nbins-1,ilev));
	    histo[ilev] += 1;
	    nhist += 1;
	}
    } else {
	for (i = 0; i < npts; i++) {
	    if (data[i] < lowv || data[i] > highv)
		continue;
	    ilev = cir_nint(data[i] - lowv);
	    ilev = max(0,min(nbins-1,ilev));
	    histo[ilev] += 1;
	    nhist += 1;
	}
    }

    /* Right, find the median value and the first quartile. */

    iclip = nbins - 1;
    nhist2 = nhist;
    for (i = 0; i <= niter; i++) {
        halflev = (nhist2 + 1)/2;
        quartlev = (nhist2 + 3)/4;
        mlev = histexam(histo,nbins,halflev);
        *median = mlev + lowv;
        qlev = histexam(histo,nbins,quartlev);
        *sigma = (mlev - qlev)*1.48;
        if (i == niter)
            break;
        irej = 0;
        jst = cir_nint(mlev + thresh*(*sigma));
        for (j = jst; j <= iclip; j++)
            irej += histo[j];
        if (irej == 0)
            break;
        iclip = jst - 1;
        nhist2 -= irej;
    }
    free(histo);
    return(CIR_OK);
}

extern float histexam(int *histo, int nhist, int level) {
    int ilev,ii;
    float value;

    ii = 0;
    ilev = -1;
    while (ii < level) 
        ii += histo[++ilev];
    value = (float)ilev - (float)(ii - level)/(float)histo[ilev] + 0.5;
    return(value);
}

/*+
 *  Name:
 *     cir_med
 *
 *  Purpose:
 *      Find the median of a float data array.  Ignore bad pixels.
 *
 *  Description:
 *      The median of a float data array is found through sorting the input
 *      array and ignoring bad pixels (if a BPM is supplied).
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      data = float * (Given)
 *          The input data
 *      bpm = unsigned char * (Given)
 *          The input BPM.  If no BPM is to be used, this should be set to NULL
 *      npts = int (Given)
 *          The number of points in the input data arrays.
 *      value = float * (Returned)
 *          A floating point value representing the median value of the
 *          input array (ignoring any bad values)
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      Care is taken not to reorder the data array.
 *              
 *  Dependencies:
 *      None.
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_med(float *data, unsigned char *bpm, int npts, float *value, 
		   char *errmsg) {
    int i,j,is_even,ilevel;
    float *buf;

    /* If there is not BPM, then just do a straight forward median */

    buf = cir_malloc(npts*sizeof(*buf));
    if (bpm == NULL) {
        is_even = !(npts & 1);
        memmove((char *)buf,(char *)data,npts*sizeof(float));
        if (is_even) {
  	    ilevel = npts/2 - 1;
	    *value = kselect(buf,npts,ilevel);
            ilevel = npts/2;
            *value = 0.5*(*value + kselect(buf,npts,ilevel));
        } else {
            ilevel = npts/2;
            *value = kselect(buf,npts,ilevel);
        }
    } else {
        j = 0;
        for (i = 0; i < npts; i++) {
            if (bpm[i] == 0)
                buf[j++] = data[i];
        }
	if (j == 0) {
	    free(buf);
	    *value = 0.0;
	    (void)sprintf(errmsg,"CIR_MED: No good pixels\n");
	    return(CIR_WARN);
	}
        is_even = !(j & 1);
        if (is_even) {
  	    ilevel = j/2 - 1;
	    *value = kselect(buf,j,ilevel);
            ilevel = j/2;
            *value = 0.5*(*value + kselect(buf,j,ilevel));
        } else {
            ilevel = j/2;
            *value = kselect(buf,j,ilevel);
        }
    }
    free(buf);

    /* Return the status value... */

    return(CIR_OK);
}

/*+
 *  Name:
 *      cir_dmed
 *
 *  Purpose:
 *      Find the median of a double data array.  Ignore bad pixels.
 *
 *  Description:
 *      The median of a double data array is found through sorting the input
 *      array and ignoring bad pixels (if a BPM is supplied).
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      data = double * (Given)
 *          The input data
 *      bpm = unsigned char * (Given)
 *          The input BPM.  If no BPM is to be used, this should be set to NULL
 *      npts = int (Given)
 *          The number of points in the input data arrays.
 *      value = double * (Returned)
 *          A double precision value representing the median value of the
 *          input array (ignoring any bad values)
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      Care is taken not to reorder the data array.
 *              
 *  Dependencies:
 *      None.
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_dmed(double *data, unsigned char *bpm, int npts, double *value,
		    char *errmsg) {
    int i,j,is_even,ilevel;
    double *buf;

    /* If there is not BPM, then just do a straight forward median */

    buf = cir_malloc(npts*sizeof(*buf));
    if (bpm == NULL) {
        is_even = !(npts & 1);
        memmove((char *)buf,(char *)data,npts*sizeof(double));
        if (is_even) {
  	    ilevel = npts/2 - 1;
	    *value = dkselect(buf,npts,ilevel);
            ilevel = npts/2;
            *value = 0.5*(*value + dkselect(buf,npts,ilevel));
        } else {
            ilevel = npts/2;
            *value = dkselect(buf,npts,ilevel);
        }
    } else {
        j = 0;
        for (i = 0; i < npts; i++) {
            if (bpm[i] == 0)
                buf[j++] = data[i];
        }
	if (j == 0) {
	    free(buf);
	    *value = 0.0;
	    (void)sprintf(errmsg,"CIR_MED: No good pixels\n");
	    return(CIR_WARN);
	}
        is_even = !(j & 1);
        if (is_even) {
  	    ilevel = j/2 - 1;
	    *value = dkselect(buf,j,ilevel);
            ilevel = j/2;
            *value = 0.5*(*value + dkselect(buf,j,ilevel));
        } else {
            ilevel = j/2;
            *value = dkselect(buf,j,ilevel);
        }
    }
    free(buf);

    /* Return the value... */

    return(CIR_OK);
}

/*+
 *  Name:
 *      cir_mean
 *
 *  Purpose:
 *      Find the mean of a float data array.  Ignore bad pixels.
 *
 *  Description:
 *      The mean of a float data array is found in the usual way, while
 *      ignoring bad pixels (if a BPM is supplied).
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      data = float * (Given)
 *          The input data
 *      bpm = unsigned char * (Given)
 *          The input BPM.  If no BPM is to be used, this should be set to NULL
 *      npts = int (Given)
 *          The number of points in the input data arrays.
 *      value = float * (Returned)
 *          A floating point value representing the mean value of the
 *          input array (ignoring any bad values)
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      A value of FLT_MAX is returned if all the pixels are bad.
 *              
 *  Dependencies:
 *      None.
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_mean(float *data, unsigned char *bpm, int npts, float *value,
		    char *errmsg) {
    int i,n;
    float sum;

    /* Separate sections depending on whether there is a BPM or not */

    sum = 0.0;
    if (bpm == NULL) {
        n = npts;
        for (i = 0; i < npts; i++) 
            sum += data[i];
    } else {
        n = 0;
        for (i = 0; i < npts; i++) {
            if (bpm[i] == 0) {
                sum += data[i];
                n++;
            }
	}
    }
    if (n > 0) {
        *value = sum/(float)n;
        return(CIR_OK);
    } else {
        *value = FLT_MAX;
        (void)sprintf(errmsg,"CIR_MEAN: All values flagged as bad\n");
        return(CIR_WARN);
    }
}

/*+
 *  Name:
 *      cir_dmean
 *
 *  Purpose:
 *      Find the mean of a double data array.  Ignore bad pixels.
 *
 *  Description:
 *      The mean of a double data array is found in the usual way, while
 *      ignoring bad pixels (if a BPM is supplied).
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      data = double * (Given)
 *          The input data
 *      bpm = unsigned char * (Given)
 *          The input BPM.  If no BPM is to be used, this should be set to NULL
 *      npts = int (Given)
 *          The number of points in the input data arrays.
 *      value = float * (Returned)
 *          A double precision value representing the mean value of the
 *          input array (ignoring any bad values)
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      A value of DBL_MAX is returned if all the pixels are bad.
 *              
 *  Dependencies:
 *      None.
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_dmean(double *data, unsigned char *bpm, int npts, double *value,
		     char *errmsg) {
    int i,n;
    double sum;

    /* Separate sections depending on whether there is a BPM or not */

    sum = 0.0;
    if (bpm == NULL) {
        n = npts;
        for (i = 0; i < npts; i++) 
            sum += data[i];
    } else {
        n = 0;
        for (i = 0; i < npts; i++) {
            if (bpm[i] == 0) {
                sum += data[i];
                n++;
            }
	}
    }
    if (n > 0) {
        *value = sum/(float)n;
        return(CIR_OK);
    } else {
        *value = DBL_MAX;
        (void)sprintf(errmsg,"CIR_DMEAN: All values flagged as bad\n");
        return(CIR_WARN);
    }
}

/*+
 *  Name:
 *     cir_medmad
 *
 *  Purpose:
 *      Find the median and the Median Absolute Deviation of an array.
 *
 *  Description:
 *      An array of values is given.  A median is found. The MAD is defined as
 *      the median absolute residual from the ensemble median for the array. 
 *      Bad pixels are ignored if a BPM is supplied
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      data = float * (Given)
 *          The input data
 *      bpm = unsigned char * (Given)
 *          The input BPM.  If no BPM is to be used, this should be set to NULL
 *      npts = int (Given)
 *          The number of points in the input data arrays.
 *      medval = float * (Returned);
 *          The output median
 *      madval = float * (Returned)
 *          A floating point value representing the median value of the
 *          input array (ignoring any bad values)
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      Care is taken not to reorder the data array.
 *              
 *  Dependencies:
 *      None.
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_medmad(float *data, unsigned char *bpm, int npts, float *medval,
		      float *madval, char *errmsg) {
    int i;
    float *wptr;

    /* First find the median value */

    (void)cir_med(data,bpm,npts,medval,errmsg);

    /* Now work out the MAD. Start by getting a bit of workspace and filling
       it with absolute residuals */

    wptr = cir_malloc(npts*sizeof(*wptr));
    for (i = 0; i < npts; i++)
        wptr[i] = (float)fabs((double)(data[i] - *medval));

    /* Now get the median value of the absolute deviations */

    (void)cir_med(wptr,bpm,npts,madval,errmsg);

    /* Tidy and exit */

    free(wptr);
    return(CIR_OK);
}

/*
 * given an array a of n elements, return the element that would be at 
 * position k, (0 <= k < n), if the array were sorted.  from Algorithms 
 * and Data Structures in C++ by Leendert Ammeraal, pg. 82.  O(n).
 *
 * NB: partially reorders data in array 
 */

/* Stolen from C. Sabbey */

extern float kselect(float *a, int n, int k) {
    while (n > 1) {
        int i = 0, j = n - 1;
        float x = a[j/2], w;

        do {
            while (a[i] < x) i++;
            while (a[j] > x) j--;
            if (i < j) { 
                w = a[i]; a[i] = a[j]; a[j] = w; 
            } else {
                if (i == j) i++;
                break;
            }
        } while (++i <= --j);

        if (k < i) 
            n = i; 
        else { 
            a += i; n -= i; k -= i; 
        }
    }

    return a[0];
}

extern double dkselect(double *a, int n, int k) {
    while (n > 1) {
        int i = 0, j = n - 1;
        double x = a[j/2], w;

        do {
            while (a[i] < x) i++;
            while (a[j] > x) j--;
            if (i < j) { 
                w = a[i]; a[i] = a[j]; a[j] = w; 
            } else {
                if (i == j) i++;
                break;
            }
        } while (++i <= --j);

        if (k < i) 
            n = i; 
        else { 
            a += i; n -= i; k -= i; 
        }
    }

    return a[0];
}

/*+
 *  Name:
 *      cir_sumbpm
 *
 *  Purpose:
 *      Find the sum of all the values in a BPM
 *
 *  Description:
 *      Find the sum of all values in a BPM. Useful for seeing of all your
 *      values are bad or not.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      bpm = unsigned char * (Given)
 *          The input BPM
 *      npts = int (Given)
 *          The number of values in the input BPM array.
 *      sumb = int * (Returned)
 *          A integer value representing the sum of the BPM
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
 *      None
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int cir_sumbpm(unsigned char *bpm, int npts, int *sumb, char *errmsg) {
    int j;

    *sumb = 0;
    for (j = 0; j < npts; j++)
        *sumb += bpm[j];
    return(CIR_OK);
}

/*+
 *  Name:
 *      fndmatch
 *
 *  Purpose:
 *      Search an ordered list of x,y coordinates to find the best match.
 *
 *  Description:
 *      A list of x and a list of y coordinates are searched to find the
 *      best match to a program value of x,y.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      x = float (Given)
 *          The X value of the object you want to match
 *      y = float (Given)
 *          The Y value of the object you want to match
 *      xlist = float * (Given)
 *          The list of possible X values
 *      ylist = float * (Given)
 *          The list of possible Y values.
 *      nlist = int (Given)
 *          The size of the xlist and ylist arrays.
 *      err = float (Given)
 *          The maximum error radius for a match.
 *
 *  Returned values:
 *      The index of the best match to x,y from the lists. A value of -1
 *      indicates that none of the items in the lists was an acceptable match.
 *
 *  Notes:
 *      The lists should be sorted so that ylist is in ascending order
 *              
 *  Dependencies:
 *      None
 *
 *  Authors:
 *      Mike Irwin (CASU, IoA)
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/  

extern int fndmatch (float x, float y, float *xlist, float *ylist, int nlist, 
    float err) {

    int isp,ifp,index,i;
    float errsq,errmin,dx,dy,poserr;
    
    /* Find lower limit index */

    isp = 0;
    ifp = nlist - 1;
    errsq = err*err;
    index = (isp + ifp)/2;
    while (ifp-isp >= 2) {
        if (ylist[index] < y - err) {
            isp = index;
            index = (index+ifp)/2;
        } else if (ylist[index] > y - err) {
            ifp = index;
            index = (index+isp)/2;
        } else {
            isp = index;
            break;
        }
    }

    /* Now find nearest one within limit */

    index = -1;
    errmin = errsq;
    for (i = isp; i < nlist; i++) {
        if (ylist[i] > y+err)
            break;
         dx = x - xlist[i];
         dy = y - ylist[i];
         poserr = dx*dx + dy*dy;
         if (poserr < errsq) {
             if (poserr <= errmin) {
                 index = i;
                 errmin = poserr;
             }
         }
    }
    return(index);
}

extern int casu_stamp(fitsfile *fptr, char *mainprog) {
    int status;

    /* Try and add the required keywords. If there is an error CFITSIO has
       inherited status and so we only need to pass the bad status back
       at the end */

    status = 0;
    (void)fits_update_key(fptr,TSTRING,"SOFTNAME",mainprog,
			  "Main program that created this file",&status);
    (void)fits_update_key(fptr,TSTRING,"SOFTVERS",PACKAGE_VERSION,
			  "Version of CASUtools",&status);
    (void)fits_update_key(fptr,TSTRING,"SOFTDATE",PACKAGE_DATE,
			  "Date this version was released",&status);
    (void)fits_update_key(fptr,TSTRING,"SOFTAUTH",PACKAGE_BUGREPORT,
			  "Contact for bug reports",&status);
    (void)fits_update_key(fptr,TSTRING,"SOFTINST",PACKAGE_URL,
			  "URL for CASU",&status);
    return(status == 0 ? CIR_OK : CIR_FATAL);
}
