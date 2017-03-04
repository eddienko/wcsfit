/*

$Id: casu_wcs.c,v 1.8 2011/06/09 12:10:33 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <fitsio.h>
#include <tools.h>

/*+
 *  Name:
 *      cir_wcsopen
 *
 *  Purpose:
 *      Opens a WCS structure.
 *
 *  Description:
 *      Reads the header of an input images and constructs a WCS structure
 *      from it.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infile = char * (Given)
 *          The name of the input FITS image. This must have a valid WCS
 *          defined in the header
 *      wcs = struct wcsprm ** (Returned)
 *          The WCSLIB structure defined from the header
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
 *      cfitsio, wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern int cir_wcsopen(char *infile, struct wcsprm **wcs, char *errmsg) {
    int status,nk,retval,nrej,nwcs,hdutype,istab;
    fitsfile *iptr = NULL;
    float val;
    char *hdr = NULL,msg[BUFSIZ];
    struct wcsprm *wwcs = NULL;

    /* Open the file and read the header out of it */

    status = 0;
    *wcs = NULL;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_get_hdu_type(iptr,&hdutype,&status);
    (void)fits_hdr2str(iptr,1,NULL,0,&hdr,&nk,&status);
    if (status != 0) {
	(void)fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"WCSOPEN: Can't read file header %s -- %s\n",
		      infile,msg);
	(void)fits_close_file(iptr,&status);
	free(hdr);
	return(2);
    }
    (void)fits_read_key(iptr,TFLOAT,"CRVAL1",&val,NULL,&status);
    if (status != 0) {
	istab = 1;
	status = 0;
    } else {
	istab = 0;
    }
    (void)fits_close_file(iptr,&status);

    /* Now get the WCS structure */

    if (hdutype == IMAGE_HDU) {
	retval = wcspih(hdr,nk,0,0,&nrej,&nwcs,&wwcs);
    } else {
	if (istab) 
	    retval = wcsbth(hdr,nk,0,0,0,NULL,&nrej,&nwcs,&wwcs);
	else
	    retval = wcspih(hdr,nk,0,0,&nrej,&nwcs,&wwcs);
    }
    if (retval != 0 || nwcs == 0) {
	(void)sprintf(errmsg,"WCSOPEN: WCSPIH failed with status: %d\n",
		      retval);
	free(hdr);
	wcsvfree(&nwcs,&wwcs);
	return(2);
    }
    free(hdr);
    if (wwcs[0].crval[0] == 0.0 && wwcs[0].crval[1] == 0.0) {
	(void)sprintf(errmsg,"WCSOPEN: crval1,2 both zero\n");
	wcsvfree(&nwcs,&wwcs);
	return(2);
    }

    /* Get the one you want (the first one in this case) */

    *wcs = (struct wcsprm *)calloc(1,sizeof(struct wcsprm));
    (*wcs)->flag = -1;
    if (wcscopy(1,wwcs,*wcs) != 0 || wcsset(*wcs) != 0) {
	(void)sprintf(errmsg,"WCSOPEN: WCSCOPY/WCSSET failed\n");
	wcsvfree(&nwcs,&wwcs);
	wcsfree(*wcs);
	*wcs = NULL;
	return(2);
    }

    /* Get out of here */

    wcsvfree(&nwcs,&wwcs);
    return(0);
}
 
/*+
 *  Name:
 *      cir_wcsclose
 *
 *  Purpose:
 *      Closes a WCS structure.
 *
 *  Description:
 *      Closes a WCS structure and frees all workspace associated with it
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcsfree.
 *
 *  Dependencies:
 *      wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/
                                                                               
extern void cir_wcsclose(struct wcsprm *wcs) {
    int retval;

    /* Free the workspace */

    retval = wcsfree(wcs);
    free(wcs);
}

/*+
 *  Name:
 *      cir_xytoradec
 *
 *  Purpose:
 *      Convert X,Y coordinates in pixels to RA and Dec in degrees.
 *
 *  Description:
 *      A WCS structure is used to convert cartesian coordinates into
 *      equatorial.  The structure must exist before entering this routine.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *      x = double (Given)
 *          The X coordinate to be converted.
 *      y = double (Given)
 *          The Y coordinate to be converted.
 *      ra = double * (Returned)
 *          The calculated RA
 *      dec = double * (Returned)
 *          The calculated Dec.
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcsp2s. The
 *      equinox of the output coordinates are defined by the wcs structure
 *      when it was created.
 *
 *  Dependencies:
 *      wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern void cir_xytoradec(struct wcsprm *wcs, double x, double y, 
			  double *ra, double *dec) {
    double xy[2],std[2],phi,theta,radec[2];
    int stat,retval;

    /* Load up the information and call the wcslib routine */

    xy[0] = x;
    xy[1] = y;
    retval = wcsp2s(wcs,1,2,xy,std,&phi,&theta,radec,&stat);

    /* Pass it back now */

    *ra = radec[0];
    *dec = radec[1];
}

/*+
 *  Name:
 *      cir_radectoxy
 *
 *  Purpose:
 *      Convert RA and Dec in degrees to X and Y in pixels.
 *
 *  Description:
 *      A WCS structure is used to convert equatorial coordiantes to
 *      cartesian  The structure must exist before entering this routine.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *      ra = double (Given)
 *          The input RA
 *      dec = double (Given)
 *          The input Dec.
 *      x = double * (Returned)
 *          The output X coordinate
 *      y = double * (Returned)
 *          The output Y coordinate
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcss2p.
 *
 *  Dependencies:
 *      wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern void cir_radectoxy(struct wcsprm *wcs, double ra, double dec,
			  double *x, double *y) {
    double xy[2],std[2],phi,theta,radec[2];
    int stat,retval;

    /* Load up the information and call the wcslib routine */

    radec[0] = ra;
    radec[1] = dec;
    retval = wcss2p(wcs,1,2,radec,&phi,&theta,std,xy,&stat);

    /* Pass it back now */

    *x = xy[0];
    *y = xy[1];
}


/*+
 *  Name:
 *      cir_radectoxieta
 *
 *  Purpose:
 *      Convert RA and Dec in degrees to standard coordinates (xi and eta)
 *      in degrees.
 *
 *  Description:
 *      A WCS structure is used to convert equatorial coordiantes to
 *      standard coordinates.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *      ra = double (Given)
 *          The input RA
 *      dec = double (Given)
 *          The input Dec.
 *      xi = double * (Returned)
 *          The output xi coordinate
 *      eta = double * (Returned)
 *          The output eta coordinate
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcss2p.
 *
 *  Dependencies:
 *      wcstools
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern void cir_radectoxieta(struct wcsprm *wcs, double ra, double dec,
			     double *xi, double *eta) {
    double xy[2],std[2],phi,theta,radec[2];
    int stat,retval;

    /* Load up the information and call the wcslib routine */

    radec[0] = ra;
    radec[1] = dec;
    retval = wcss2p(wcs,1,2,radec,&phi,&theta,std,xy,&stat);

    /* Pass it back now */

    *xi = DEGRAD*std[0];
    *eta = DEGRAD*std[1];
}

extern void cir_xytoxy_list(struct wcsprm *wcs1, struct wcsprm *wcs2, int nc,
			    double *x1, double *y1, double *ra, double *dec, 
			    double *x2, double *y2) {
    double *xy,*xieta,*phi,*theta,*radec;
    int *stat,i;

    /* Get some workspace for the wcslib arrays */

    xy = cir_malloc(nc*2*sizeof(*xy));
    xieta = cir_malloc(nc*2*sizeof(*xieta));
    radec = cir_malloc(nc*2*sizeof(*radec));
    stat = cir_malloc(nc*sizeof(*stat));
    phi = cir_malloc(nc*sizeof(*phi));
    theta = cir_malloc(nc*sizeof(*theta));

    /* Load up */

    for (i = 0; i < nc; i++) {
	xy[2*i] = x1[i];
	xy[2*i+1] = y1[i];
    }

    /* Do the transformation */

    (void)wcsp2s(wcs1,nc,2,xy,xieta,phi,theta,radec,stat);

    /* Transform the Ra Dec array to xy in the second frame */

    (void)wcss2p(wcs2,nc,2,radec,phi,theta,xieta,xy,stat);

    /* Copy the ra,dec,x,y answers over now */

    for (i = 0; i < nc; i++) {
	if (ra != NULL && dec != NULL) {
  	    ra[i] = radec[2*i];
	    dec[i] = radec[2*i+1];
	}
	x2[i] = xy[2*i];
	y2[i] = xy[2*i+1];
    }

    /* Tidy and exit */

    freespace(radec);
    freespace(xieta);
    freespace(xy);
    freespace(stat);
    freespace(phi);
    freespace(theta);
}

/*

$Log: casu_wcs.c,v $
Revision 1.8  2011/06/09 12:10:33  jim
Fixed a little bug in wcsopen to trap for situations where a WCS isn't
parsed from a header

Revision 1.7  2011-01-31 15:03:11  jim
Modified again to make sure that tables with image WCS keywords are parsed
correctly

Revision 1.6  2011-01-14 10:53:20  jim
Modified so that different wcslib routines open the header parameters
depending on whether the file is a table or an image

Revision 1.5  2011-01-12 13:17:07  jim
Added an error message into wcsopen and changed a parameter to wcsbth so that
only FITS standard WCS cards are used

Revision 1.4  2010-11-22 10:42:42  jim
Added cir_xytoxy_list


*/
