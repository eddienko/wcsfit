/*

$Id: cir_vstfix.c,v 1.1 2013/05/30 14:28:26 jim Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>

#include <vstfix.h>
#include <tools.h>
/* Columns in VST tables that have flux information */

#define VSTFIX_NFLUXCOLS 35
int fluxcols[] = {2,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,
		  35,36,37,38,39,40,41,42,43,44,45,49,50,51,52,53,54};

static fitsfile *iptr = NULL;
static fitsfile *cptr = NULL;
static fitsfile *optr = NULL;
static float *work = NULL;
static short int *confs[DEFNEXTN] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
				     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
				     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
				     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
static struct wcsprm **wcslist = NULL;
static struct wcsprm *wcsin = NULL;
static float *xpos = NULL;
static float *ypos = NULL;
static float *fluxes[VSTFIX_NFLUXCOLS];
static int ncomb = 1;

static float fraction (float x, float y, float r_out);
static void tidy(void);

extern int cir_vstfix(char *infile, char *corr_file, char *confmap, 
		      char *outcat, char *errmsg) {
    int status,nhdu,hdutype,anynul,nicomb,i,j,k,ii,jj,n,ixlim1,ixlim2;
    int iylim1,iylim2,iz,colnum,midpt,kk,imef,ix,iy;
    long naxis[2],npts,nrows;
    float *zpsmed,*bin1dx,*bin1dy,*coord,*radial,rcore,delmag,scale,confwt;
    float fluxcor,cnumb,xj,yj;
    double x1,y1,xi,eta,ra,dec,radius,x2,y2;
    char prov[FLEN_VALUE],cmapextn[FLEN_VALUE],confbase[BUFSIZ],msg[BUFSIZ];
    char inextn[BUFSIZ],*tfits;

    /* Initialise a few things */

    for (i = 0; i < VSTFIX_NFLUXCOLS; i++)
	fluxes[i] = NULL;

    /* Start by opening the input catalogue and making sure it has the
       correct number of extensions */

    status = 0;
    (void)fits_open_file(&iptr,infile,READWRITE,&status);
    if (status != 0) {
        (void)sprintf(errmsg,"VSTFIX: Can't open input catalogue %s",infile);
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_num_hdus(iptr,&nhdu,&status);
    nhdu--;
    if (nhdu != DEFNEXTN) {
        (void)sprintf(errmsg,
		      "VSTFIX: Catalogue file %s has %d extensions and must have %d",
		      infile,nhdu,DEFNEXTN);
	tidy();
	return(CIR_FATAL);
    }

    /* Is this a stack of more than one exposure? */

    (void)fits_movabs_hdu(iptr,2,&hdutype,&status);
    (void)fits_read_key(iptr,TINT,"NICOMB",&nicomb,NULL,&status);
    if (status != 0) {
	status = 0;
	nicomb = 1;
    }
    ncomb = nicomb;

    /* If this is the product of a stack then we need to read the 
       confidence map. If the user specified 'auto' then try to use the
       provenance info header to find the name of a confidence map for one
       of the constituent images. If not, then use the name that the user
       provided and hope that it's correct */

    if (nicomb > 1) {
	if (!strncmp(confmap,"auto",4)) {
	    (void)fits_read_key(iptr,TSTRING,"PROV0001",prov,NULL,&status);
	    (void)fits_open_file(&optr,prov,READONLY,&status);
	    if (status != 0) {
		(void)sprintf(errmsg,"VSTFIX: Unable to open %s using auto",
			      prov);
		tidy();
		return(CIR_FATAL);
	    }
	    (void)fits_read_key(optr,TSTRING,"CIR_CPM",cmapextn,NULL,&status);
	    (void)fits_parse_rootname(cmapextn,confbase,&status);
	    closefits(optr);
	} else {
	    (void)strcpy(confbase,confmap);
	}

	/* Ok, try and open this file now... */

	status = 0;
	(void)fits_open_file(&cptr,confbase,READWRITE,&status);
	if (status != 0) {
	    (void)sprintf(errmsg,"VSTFIX: Can't open input conf map %s",
			  confbase);
	    tidy();
	    return(CIR_FATAL);
	}
	(void)fits_get_num_hdus(cptr,&nhdu,&status);
	nhdu--;
	if (nhdu != DEFNEXTN) {
	    (void)sprintf(errmsg,
			  "VSTFIX: Conf map %s has %d extensions and must have %d",
			  confbase,nhdu,DEFNEXTN);
	    tidy();
	    return(CIR_FATAL);
	}

	/* Read each of the extensions */

	npts = DEFNX_TRIMMED*DEFNY_TRIMMED;
	for (i = 0; i < DEFNEXTN; i++) {
	    (void)fits_movabs_hdu(cptr,i+2,&hdutype,&status);
	    (void)fits_get_img_size(cptr,2,naxis,&status);
	    if (naxis[0] != DEFNX_TRIMMED || naxis[1] != DEFNY_TRIMMED) {
		(void)sprintf(errmsg,
			      "VSTFIX: Confmap %s[%d] has data size [%ld,%ld] and should be [%d,%d]",
			      confbase,i+1,naxis[0],naxis[1],DEFNX_TRIMMED,
			      DEFNY_TRIMMED);
		tidy();
		return(CIR_FATAL);
	    }
	    confs[i] = cir_malloc(npts*sizeof(short int));
	    (void)fits_read_img(cptr,TSHORT,1,npts,NULL,confs[i],&anynul,
				&status);
	}
	closefits(cptr);

	/* Allocate some memory for the wcs structures for the input 
	   tables */

	wcslist = cir_malloc(nicomb*DEFNEXTN*sizeof(struct wcsprm *));
    }

    /* Open up the correction factor FITS file */

    (void)fits_open_file(&optr,corr_file,READONLY,&status);
    if (status != 0) {
        (void)sprintf(errmsg,"VSTFIX: Can't open correction file %s",
		      corr_file);
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_num_hdus(optr,&nhdu,&status);
    nhdu--;
    if (nhdu != VSTFIX_NTABS) {
        (void)sprintf(errmsg,
		      "VSTFIX: Correction file %s has %d extensions and must have %d",
		      corr_file,nhdu,VSTFIX_NTABS);
	tidy();
	return(CIR_FATAL);
    }  
    
    /* Get some memory to hold all the info you need from the correction file */

    work = cir_malloc((DEFNEXTN+2*VSTFIX_MPBIN+2*VSTFIX_MRBIN)*sizeof(float));
    zpsmed = work;
    bin1dx = work + DEFNEXTN;
    bin1dy = work + DEFNEXTN + VSTFIX_MPBIN;
    coord = work + DEFNEXTN + 2*VSTFIX_MPBIN;
    radial = work + DEFNEXTN + 2*VSTFIX_MPBIN + VSTFIX_MRBIN;

    /* Now read the detector offsets from the first extension */

    (void)fits_movabs_hdu(optr,2,&hdutype,&status);
    (void)fits_get_num_rows(optr,&nrows,&status);
    if (nrows != DEFNEXTN) {
	sprintf(errmsg,
		"VSTFIX: First extension of corr file has %d rows and must have %d\n",
		nrows,DEFNEXTN);
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_colnum(optr,CASEINSEN,"zpsmed",&colnum,&status);
    (void)fits_read_col(optr,TFLOAT,colnum,1,1,DEFNEXTN,NULL,zpsmed,&anynul,
			&status);
    if (status != 0) {
	sprintf(errmsg,
		"VSTFIX: Unable to read column %s in corr file first extn",
		"zpsmed");
	tidy();
	return(CIR_FATAL);
    }
	
    /* Now read the bilinear coefficients from the second extension */

    (void)fits_movabs_hdu(optr,3,&hdutype,&status);
    (void)fits_get_num_rows(optr,&nrows,&status);
    if (nrows != VSTFIX_MPBIN) {
	sprintf(errmsg,
		"VSTFIX: Second extension of corr file has %d rows and must have %d\n",
		nrows,VSTFIX_MPBIN);
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_colnum(optr,CASEINSEN,"bin1dx",&colnum,&status);
    (void)fits_read_col(optr,TFLOAT,colnum,1,1,VSTFIX_MPBIN,NULL,bin1dx,&anynul,
			&status);
    if (status != 0) {
	sprintf(errmsg,
		"VSTFIX: Unable to read column %s in corr file second extn",
		"bin1dx");
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_colnum(optr,CASEINSEN,"bin1dy",&colnum,&status);
    (void)fits_read_col(optr,TFLOAT,colnum,1,1,VSTFIX_MPBIN,NULL,bin1dy,&anynul,
			&status);
    if (status != 0) {
	sprintf(errmsg,
		"VSTFIX: Unable to read column %s in corr file second extn",
		"bin2dx");
	tidy();
	return(CIR_FATAL);
    }
    
    /* Now read radial coefficients from third extension */

    (void)fits_movabs_hdu(optr,4,&hdutype,&status);
    (void)fits_get_num_rows(optr,&nrows,&status);
    if (nrows != VSTFIX_MRBIN) {
	sprintf(errmsg,
		"VSTFIX: third extension of corr file has %d rows and must have %d\n",
		nrows,VSTFIX_MRBIN);
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_colnum(optr,CASEINSEN,"coord",&colnum,&status);
    (void)fits_read_col(optr,TFLOAT,colnum,1,1,VSTFIX_MRBIN,NULL,coord,&anynul,
			&status);
    if (status != 0) {
	sprintf(errmsg,
		"VSTFIX: Unable to read column %s in corr file third extn",
		"coord");
	tidy();
	return(CIR_FATAL);
    }
    (void)fits_get_colnum(optr,CASEINSEN,"radial",&colnum,&status);
    (void)fits_read_col(optr,TFLOAT,colnum,1,1,VSTFIX_MRBIN,NULL,radial,&anynul,
			&status);
    if (status != 0) {
	sprintf(errmsg,
		"VSTFIX: Unable to read column %s in corr file third extn",
		"radial");
	tidy();
	return(CIR_FATAL);
    }
    closefits(optr);

    /* Create the output file. Start by copying the primary */

    if (access(outcat,F_OK) == 0) 
	remove(outcat);
    (void)fits_create_file(&optr,outcat,&status);
    (void)fits_movabs_hdu(iptr,1,&hdutype,&status);
    (void)fits_copy_hdu(iptr,optr,0,&status);
    i = 1;
    (void)fits_update_key(optr,TLOGICAL,"VSTFIX",&i,"VST illumcor fix done",
			  &status);

    /* If this is a stack then read the provenence keywords, translate
       them to catalogue names and read the WCS for each of them too */

    if (nicomb > 1) {
	(void)fits_movabs_hdu(iptr,2,&hdutype,&status);
	n = 0;
	for (j = 1; j <= nicomb; j++) {
	    (void)sprintf(msg,"PROV%04d",j);
	    (void)fits_read_key(iptr,TSTRING,msg,prov,NULL,&status);
	    tfits = strstr(prov,".fit");
	    for (k = 1; k <= DEFNEXTN; k++) {
		n++;
		(void)sprintf(tfits,"_cat.fits[%d]",k);
		if (cir_wcsopen(prov,wcslist+n-1,msg) != CIR_OK) {
		    sprintf(errmsg,"VSTFIX: Unable to read WCS in %s\n",
			    prov);
		    tidy();
		    return(CIR_FATAL);
		}
	    }
	}
    }

    /* Right, now back to the input catalogue. Loop for each of the 
       extensions and start out by reading the WCS info from the header */

    for (i = 1; i <= DEFNEXTN; i++) {
	(void)fits_movabs_hdu(iptr,i+1,&hdutype,&status);
	(void)sprintf(inextn,"%s[%d]",infile,i);
	if (cir_wcsopen(inextn,&wcsin,msg) != CIR_OK) {
	    sprintf(errmsg,"VSTFIX: Unable to read WCS in %s\n",inextn);
	    tidy();
	    return(CIR_FATAL);
	}
	(void)fits_read_key(iptr,TFLOAT,"RCORE",&rcore,NULL,&status);
	
	/* Copy the current input HDU to the output. Get some memory
	   for the x,y positions and read them in */

	(void)fits_copy_hdu(iptr,optr,0,&status);
	(void)fits_get_num_rows(optr,&nrows,&status);
	xpos = cir_malloc(nrows*sizeof(float));
	ypos = cir_malloc(nrows*sizeof(float));
	(void)fits_read_col(optr,TFLOAT,3,1,1,nrows,NULL,xpos,&anynul,&status);
	(void)fits_read_col(optr,TFLOAT,5,1,1,nrows,NULL,ypos,&anynul,&status);

	/* Read the flux columns now */

	for (j = 0; j < VSTFIX_NFLUXCOLS; j++) {
	    fluxes[j] = cir_malloc(nrows*sizeof(float));
	    (void)fits_read_col(optr,TFLOAT,fluxcols[j],1,1,nrows,NULL,
				fluxes[j],&anynul,&status);
	}
	    
	/* Loop for each object. If this isn't a stack, then it's we just
	   read the corrections from the tables. */

	midpt = (VSTFIX_MPBIN+1)/2;
	for (j = 0; j < nrows; j++) {
	    x1 = (double)xpos[j];
	    y1 = (double)ypos[j];
	    cir_xytoradec(wcsin,x1,y1,&ra,&dec);
	    if (nicomb == 1) {
		cir_radectoxieta(wcsin,ra,dec,&xi,&eta);
		xi /= DEGRAD;
		eta /= DEGRAD;
		radius = sqrt(xi*xi + eta*eta);
		ii = min(VSTFIX_MPBIN,max(1,cir_nint(xi*VSTFIX_MPRES+midpt)));
		jj = min(VSTFIX_MPBIN,max(1,cir_nint(eta*VSTFIX_MPRES+midpt)));
		kk = min(VSTFIX_MRBIN,max(1,cir_nint(radius*VSTFIX_MPRES)));
		delmag = bin1dx[ii-1] + bin1dy[jj-1] + radial[kk-1] + 
		    zpsmed[i-1];
		scale = pow(10.0,(double)(-0.4*delmag));
		confwt = 100.0;
		fluxcor = confwt*scale;
		cnumb = confwt;

	    /* Otherwise we need to go through all extensions of all the
	       input images to see which extensions this object came from.
	       Then we can sum up the contributions from each extension for
	       the object */

	    } else {
		fluxcor = 0.0;
		cnumb = 0.0;
		for (k = 0; k < DEFNEXTN*nicomb; k++) {
		    imef = k % DEFNEXTN;
		    cir_radectoxy(wcslist[k],ra,dec,&x2,&y2);

		    /* Only continue with this extension if it's possble
		       the object is here */

		    if (x2 < -rcore || x2 > DEFNX_TRIMMED+rcore ||
			y2 < -rcore || y2 > DEFNY_TRIMMED+rcore)
			continue;

		    /* Work out an average confidence in the aperture */

		    cir_radectoxieta(wcslist[k],ra,dec,&xi,&eta);
		    xi /= DEGRAD;
		    eta /= DEGRAD;
		    confwt = 0.0;
		    ixlim1 = max(0,(int)(x2 - rcore)-1);
		    ixlim2 = min(DEFNX_TRIMMED,(int)(x2 + rcore));
		    iylim1 = max(0,(int)(y2 - rcore)-1);
		    iylim2 = min(DEFNY_TRIMMED,(int)(y2 + rcore));
		    for (iy = iylim1; iy <= iylim2; iy++) {
			iz = iy*DEFNX_TRIMMED;
			for (ix = ixlim1; ix <= ixlim2; ix++) {
			    xj = (float)(ix+1) - (float)x2;
			    yj = (float)(iy+1) - (float)y2;
			    confwt += fraction(xj,yj,rcore)*confs[imef][iz+ix];
			}
		    }
		    confwt /= (M_PI*rcore*rcore);

		    /* Now compute the contribution of the correction from
		       this particular extension */

		    radius = sqrt(xi*xi + eta*eta);
		    ii = min(VSTFIX_MPBIN,max(1,cir_nint(xi*VSTFIX_MPRES+midpt)));
		    jj = min(VSTFIX_MPBIN,max(1,cir_nint(eta*VSTFIX_MPRES+midpt)));
		    kk = min(VSTFIX_MRBIN,max(1,cir_nint(radius*VSTFIX_MPRES)));
		    delmag = bin1dx[ii-1] + bin1dy[jj-1] + radial[kk-1] + 
			zpsmed[imef];
		    scale = pow(10.0,(double)(-0.4*delmag));
		    fluxcor += confwt*scale;
		    cnumb += confwt;
		}
	    }

	    /* Ok, if we've got some sort of contribution then we can
	       apply it to the fluxes */
		    
	    if (cnumb > 0.0) {
		scale = fluxcor/cnumb;
		for (k = 0; k < VSTFIX_NFLUXCOLS; k++) 
		    fluxes[k][j] *= scale;
	    }
	}

	/* Write out the fluxes now and release the memory for each
	   column in turn */

	for (j = 0; j < VSTFIX_NFLUXCOLS; j++) {
	    (void)fits_write_col(optr,TFLOAT,fluxcols[j],1,1,nrows,fluxes[j],
				 &status);
	    freespace(fluxes[j]);
	}
	j = 1;
	(void)fits_update_key(optr,TLOGICAL,"VSTFIX",&j,"VST illumcor fix done",
			      &status);

	/* Free up some more memory */

	freespace(xpos);
	freespace(ypos);
	freewcs(wcsin);
    } /* End of catalogue header extension loop */

    /* Free up some workspace and close up some files */

    tidy();
    return(CIR_OK);
}

static float fraction (float x, float y, float r_out) {
    float r,t,x_a,x_b,frac,tanao2,cosa,tanp2a,sqrt2o2;

    r = sqrtf(x*x + y*y);
    sqrt2o2 = 0.5*M_SQRT2;

    /* is it worth bothering? */

    if(r > r_out+sqrt2o2) 
	return(0.0);

    /* is it trivially all in? */

    if(r < r_out-sqrt2o2) 
	return(1.0);

    /* bugger - have to do some work then ... ok first ...
     * use 8-fold symmetry to convert to 0-45 degree range */

    x = fabsf(x);
    y = fabsf(y);
    if(y > x) {
	t = x;
	x = y;
	y = t;
    }

    /* If the angles are too close to cardinal points, then fudge something */

    if (x > 0.0 && y > 0.0) {
        tanao2 = 0.5*y/x;
        tanp2a = x/y;
        cosa = x/sqrt(x*x + y*y);
    } else {
        tanao2 = 0.00005;
        tanp2a = 10000.0;
        cosa = 1.0;
    }

    /* only outer radius - compute linear intersections top and bot of pixel */

    x_a = x - tanao2 + (r_out - r)/cosa;
    if(x_a < x+0.5) {

	/* intersects */

	x_b = x + tanao2 + (r_out - r)/cosa;

	/* three cases to consider */

	if(x_a < x-0.5)
	    frac = 0.5*max(0.0,x_b-(x-0.5))*max(0.0,x_b-(x-0.5))*tanp2a;
	else {
	    if(x_b > x+0.5)
		frac = 1.0 - 0.5*(x+0.5-x_a)*(x+0.5-x_a)*tanp2a;
	    else
		frac = 0.5-(x-x_a)+0.5*(x_b-x_a);
	}
    } else  /* missed entirely */
	frac = 1.0;

    return(frac);
}

static void tidy(void) {
    int i;
    closefits(iptr);
    closefits(optr);
    closefits(cptr);
    freespace(work);
    for (i = 0; i < DEFNEXTN; i++) 
	freespace(confs[i]);
    if (wcslist != NULL) {
	for (i = 0; i < ncomb*DEFNEXTN; i++)
	    freewcs(wcslist[i]);
    }
    freewcs(wcsin);
    freespace(xpos);
    freespace(ypos);
    for (i = 0; i < VSTFIX_NFLUXCOLS; i++)
	freespace(fluxes[i]);
}

/*

$Log: cir_vstfix.c,v $
Revision 1.1  2013/05/30 14:28:26  jim
New Entry


*/
