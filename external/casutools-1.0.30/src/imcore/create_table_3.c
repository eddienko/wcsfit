/*

$Id: create_table_3.c,v 1.2 2014/07/18 08:46:11 jim Exp $

*/

#include <stdio.h>
#include <math.h>
#include "imcore.h"
#include "errcodes.h"
#include "util.h"
#include "floatmath.h"

/* Column numbers for each item to be stored */

#define COL_NUMBER      1
#define COL_X           2
#define COL_Y           3
#define COL_FLUXISO     4
#define COL_PEAKHEIGHT  5
#define COL_ELLIPT      6
#define COL_SIGMA       7
#define COL_PA          8
#define COL_AREAL1      9
#define COL_AREAL2     10
#define COL_AREAL3     11
#define COL_AREAL4     12
#define COL_AREAL5     13
#define COL_AREAL6     14
#define COL_AREAL7     15
#define COL_AREAL8     16

/* Number of columns in the table */

#define NCOLS 16

/* Column definitions */

static char *ttype[NCOLS]={"No.","X_coordinate","Y_coordinate",
			   "Isophotal_flux","Peak_height","Ellipticity",
			   "Gaussian_sigma","Position_angle",
			   "Areal_1_profile","Areal_2_profile","Areal_3_profile",
			   "Areal_4_profile","Areal_5_profile","Areal_6_profile",
			   "Areal_7_profile","Areal_8_profile"};
static char *tunit[NCOLS]={" ","Pixels","Pixels","Counts","Counts"," ",
			   "Pixels","Degrees","Pixels","Pixels","Pixels",
			   "Pixels","Pixels","Pixels","Pixels","Pixels"};

static char *tform[NCOLS]={"1E","1E","1E","1E","1E","1E","1E","1E","1E","1E",
			   "1E","1E","1E","1E","1E","1E"};

static int areal_cols[NAREAL] = {COL_AREAL1,COL_AREAL2,COL_AREAL3,COL_AREAL4,
				 COL_AREAL5,COL_AREAL6,COL_AREAL7,COL_AREAL8};

extern int tabinit_3(char *infile, char *outtab, char *errstr) {
    int retval;

    /* Call the generic routine to open a new output table */

    retval = tabinit_gen(infile,outtab,NCOLS,ttype,tunit,tform,errstr);
    if (retval != ERRCODE_OK)
	return(retval);
	
    /* Get out of here */

    return(ERRCODE_OK);
}

extern int do_seeing_3(ap_t *ap, char *errstr) {
    int retval;

    /* Just call the generic seeing routine */

    retval = do_seeing_gen(ap,COL_ELLIPT,COL_PEAKHEIGHT,areal_cols,
			   errstr);
    if (retval != ERRCODE_OK)
	return(retval);

    /* Get out of here */

    return(ERRCODE_OK);
}
        

extern int process_results_3(ap_t *ap, char *errstr) {
    float momresults[8],parmall[NPAR];
    float sxx,syy,srr,sxy,ecc,temp,xx,theta,radeg,ell,nr;
    float yy,sigma,peak,areal1,iso_flux;
    float areal2,areal3,areal4,areal5,areal6,areal7,areal8,aa,bb;
    int iareal[NAREAL],i,status;
    long nrows;
    char errmsg[FLEN_STATUS];

    /* Do a basic moments analysis and work out the areal profiles*/

    moments(ap,momresults);
    if (momresults[0] < 0)
	return(ERRCODE_BADRES);
    areals(ap,iareal);

    /* See if this object makes the cut in terms of its size.  If not, then
       just return with good status */

    if (iareal[0] < ap->ipnop || momresults[3] < ap->xintmin)
	return(ERRCODE_OK);

    /* Try and deblend the images if it is requested and justified */

    parmall[0] = momresults[3];
    parmall[1] = momresults[1];
    parmall[2] = momresults[2];
    parmall[3] = ap->thresh;
    for (i = 4; i < 8; i++) 
        parmall[i] = momresults[i];
    for (i = 0; i < NAREAL; i++)
        parmall[i+8] = (float)iareal[i];

    /* Massage the results and write them to the fits table */

    radeg = 180.0/M_PI;
    sxx = parmall[4];
    sxy = parmall[5];
    syy = parmall[6];    
    if(sxy > 0.0)
      sxy = MAX(1.0e-4,MIN(sxy,sqrtf(sxx*syy)));
    else
      sxy = MIN(-1.0e-4,MAX(sxy,-sqrtf(sxx*syy)));

    srr = MAX(0.5,sxx+syy);
    ecc = sqrtf((syy-sxx)*(syy-sxx)+4.0*sxy*sxy)/srr;
    temp = MAX((1.0-ecc)/(1.0+ecc),0.0);
    ell = 1.0 - sqrtf(temp);
    ell = MIN(0.99,MAX(0.0,ell));
    xx = 0.5*(1.0+ecc)*srr-sxx;
    if(xx == 0.0)
	theta = 0.0;
    else
	theta = 90.0-radeg*atanf(sxy/xx);

    /* Create a list of values */

    status = 0;
    (void)fits_get_num_rows(tptr,&nrows,&status);
    nrows++;
    iso_flux = parmall[0];
    nr = (float)nrows;
    xx = parmall[1];
    yy = parmall[2];
    sigma = sqrt(srr);
    peak = parmall[7];
    areal1 = parmall[8];
    areal2 = parmall[9];
    areal3 = parmall[10];
    areal4 = parmall[11];
    areal5 = parmall[12];
    areal6 = parmall[13];
    areal7 = parmall[14];
    areal8 = parmall[15];

    /* Store away the results for this object */

    nr = (float)nrows;
    (void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&iso_flux,
			 &status);
    (void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&peak,
			 &status);
    (void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&ell,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&sigma,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&theta,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&areal1,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&areal2,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&areal3,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&areal4,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&areal5,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&areal6,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&areal7,&status);
    (void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&areal8,&status);
    if (status != 0) {
	fits_get_errstatus(status,errmsg);
	sprintf(errstr,"Unable to write catalogue data: %s (status = %d)",
		errmsg,status);
	return(ERRCODE_FILE_IO);
    }

    /* Write out the .ell file info... */

    if (ellfp != NULL) {
	aa = sqrtf(areal1/(M_PI*(1.0 - ell)));
	bb = aa*(1.0 - ell);
	fprintf(ellfp,"image; ellipse %9.3f %9.3f %7.2f %7.2f %6.1f\n",
		xx,yy,aa,bb,theta);
    }

    /* Get outta here */

    return(ERRCODE_OK);
}

extern int process_results_list_3(ap_t *ap, objstruct *objlist[], int nbit, 
                                  char *errstr) {
    int j,i;
    float sxx,syy,sxy,tpk,*d,parrad,tmn,radeg;
    float xsq,ysq,xy,dx,dy,t,srr,ecc,temp,ell,xxx,theta,nr;
    float iso_flux,zero;
    float sigma,areal1,areal2,areal3,areal4,areal5,areal6,areal7,areal8;
    float peak,xx,yy,xbar,ybar,skylev,skyrms,aa,bb;
    int iareal[NAREAL],x1,y1,x2,y2,ix,iy,indy,ind,nup,status,k;
    long nrows;
    char errmsg[BUFSIZ];

    /* Check to see if this can even appear on the image */

    xx = objlist[0]->x_i;
    yy = objlist[0]->y_i;
    if (xx < 1 || yy < 1 || xx > ap->lsiz || yy > ap->csiz) {
	for (i = 0; i < nbit; i++) {
	    xx = objlist[i]->x_i;
	    yy = objlist[i]->y_i;
	    zero = 0.0;
	    status = 0;
            (void)fits_get_num_rows(tptr,&nrows,&status);
            nrows++;
	    nr = (float)nrows;
	    (void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&zero,
				 &status);
	    if (status != 0) {
		fits_get_errstatus(status,errmsg);
		sprintf(errstr,"Unable to write catalogue data: %s (status = %d)",
			errmsg,status);
		return(ERRCODE_FILE_IO);
	    }
	}
    }

    /* Ok, here's a big loop where we try to refine the final parameters. */

    parrad = ap->rcore + 0.5;
    radeg = 180.0/M_PI;
    d = ap->data;
    for (k = 0; k < nbit; k++) {
        xx = objlist[k]->x_i;
        yy = objlist[k]->y_i;

        /* Start with a better estimate of the x,y positions. This is only
           done if the object is fully on the image. If this is not the case
           then zero the results */

        if (xx < ap->rcore || yy < ap->rcore || xx > (ap->lsiz - ap->rcore) ||
            yy > (ap->csiz - ap->rcore)) {
            for (j = 0; j < NAREAL; j++)
                iareal[j] = 0;
            xbar = xx;
            ybar = yy;
            sxx = 0.0;
            syy = 0.0;
            sxy = 0.0;
            tpk = 0.0;
            tmn = 0.0;
        } else {
            
            /* Get the range in the data map to examine */

            x1 = MAX(1,(int)(xx - parrad)) - 1;
            x2 = MIN(ap->lsiz,(int)(xx + parrad)) - 1;
            y1 = MAX(1,(int)(yy - parrad)) - 1;
            y2 = MIN(ap->csiz,(int)(yy + parrad)) - 1;

            /* Zero some accumulators */

            xsq = 0.0;
            ysq = 0.0;
            xy = 0.0;
            xbar = 0.0;
            ybar = 0.0;
            tmn = 0.0;
            tpk = 0.0;
            for (j = 0; j < NAREAL; j++)
                iareal[j] = 0;
            for (iy = y1; iy < y2; iy++) {
                indy = iy*(ap->lsiz);
                for (ix = x1; ix < x2; ix++) {
                    ind = indy + ix;
                    dx = (float)(ix+1) - xx;
                    dy = (float)(iy+1) - yy;
                    t = d[ind];
                    tpk = MAX(tpk,t);
                    t = fraction(dx,dy,ap->rcore)*t;
                    tmn += t;
                    xbar += dx*t;
                    ybar += dy*t;
                    xsq += (dx*dx + 0.0833333)*t;
                    ysq += (dy*dy + 0.0833333)*t;
                    xy += dx*dy*t;
                    if (t > ap->thresh) {
                        nup = MIN(NAREAL,(int)(logf(t)*(ap->fconst) - 
                                               ap->areal_offset)+1);
                        nup = MAX(1,nup);
                        for (j = 0; j < nup; j++)
                            iareal[j]++;
                    }
                }
            }
	    if (tmn != 0.0) {
		xbar /= tmn;
		ybar /= tmn;
	    } else {
		xbar = 0.0;
		ybar = 0.0;
	    }
	    if (fabs(xbar) > 5.0 || fabs(ybar) > 5.0) {
		xbar = 0.0;
		ybar = 0.0;
	    }
            sxx = MAX(0.0,xsq/MAX(1.0,tmn) - xbar*xbar);
            syy = MAX(0.0,ysq/MAX(1.0,tmn) - ybar*ybar);
            sxy = xy/MAX(1.0,tmn) - xbar*ybar;
            xbar += xx;
            ybar += yy;
        }
           
        /* Ellipse parameters */

        if (sxy > 0.0) 
            sxy = MAX(0.0001,MIN(sxy,sqrtf(sxx*syy)));
        else
            sxy = MIN(-0.0001,MAX(sxy,-sqrtf(sxx*syy)));
        srr = MAX(0.5,sxx+syy);
        ecc = sqrtf(powf((sxx-syy),2.0) + 4.0*sxy*sxy)/srr;
        temp = MAX((1.0-ecc)/(1.0+ecc),0.0);
        ell = 1.0 - sqrtf(temp);
        ell = MIN(0.99,MAX(0.0,ell));
        xxx = 0.5*(1.0+ecc)*srr-sxx;
        if(xxx == 0.0)
            theta = 0.0;
        else
            theta = 90.0 - radeg*atanf(sxy/xxx);

        /* Create a list of values */

        status = 0;
        (void)fits_get_num_rows(tptr,&nrows,&status);
        nrows++;
        nr = (float)nrows;
        iso_flux = tmn;
        peak = tpk;
        xx = xbar;
        yy = ybar;
        sigma = sqrt(srr);
        areal1 = iareal[0];
        areal2 = iareal[1];
        areal3 = iareal[2];
        areal4 = iareal[3];
        areal5 = iareal[4];
        areal6 = iareal[5];
        areal7 = iareal[6];
        areal8 = iareal[7];
        imcore_backest(ap,xx,yy,&skylev,&skyrms);

        /* Store away the results for this object */

        nr = (float)nrows;
	(void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&iso_flux,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&peak,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&ell,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&sigma,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&theta,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&areal1,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&areal2,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&areal3,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&areal4,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&areal5,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&areal6,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&areal7,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&areal8,&status);
        if (status != 0) {
            fits_get_errstatus(status,errmsg);
            sprintf(errstr,"Unable to write catalogue data: %s (status = %d)",
                    errmsg,status);
            return(ERRCODE_FILE_IO);
        }

        /* Write out the .ell file info... */

        if (ellfp != NULL) {
            /* aa = sqrtf(areal1/(M_PI*(1.0 - ell))); */
            /* bb = aa*(1.0 - ell); */
            /* fprintf(ellfp,"image; ellipse %9.3f %9.3f %7.2f %7.2f %6.1f\n", */
            /*         xx,yy,aa,bb,theta); */
            fprintf(ellfp,"image; ellipse %9.3f %9.3f %7.2f %7.2f %6.1f\n",
                    objlist[k]->x_i,objlist[k]->y_i,ap->rcore,ap->rcore,0.0);
        }
    }
    return(ERRCODE_OK);
}                   


extern int tabclose_3(ap_t *ap, char *errstr) {

    return(tabclose_gen(ap,errstr));
}

extern int readtab_list_3(char *infile, objstruct **ob, int *nr, 
			  char *errstr) {
    int status,anynul;
    long nrows,i;
    float *work;
    fitsfile *fptr;
    char msg[BUFSIZ];

    /* Start by opening the FITS table and getting the number of rows */
   
    *nr = 0;
    status = 0;
    (void)fits_open_file(&fptr,infile,READONLY,&status);
    (void)fits_get_num_rows(fptr,&nrows,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errstr,"Unable to open input table %s -- %s\n",infile,
                      msg);
        return(ERRCODE_FILE_IO);
    }

    /* Get some workspace */

    work = malloc(nrows*sizeof(*work));
    *ob = malloc(nrows*sizeof(objstruct));

    /* Read X, Y columns */

    (void)fits_read_col(fptr,TFLOAT,COL_X,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errstr,"Unable to read X from table %s -- %s\n",infile,
                      msg);
        free(ob);
        free(work);
        (void)fits_close_file(fptr,&status);
        return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
        (*ob+i)->x_t = work[i];
    (void)fits_read_col(fptr,TFLOAT,COL_Y,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errstr,"Unable to read Y from table %s -- %s\n",infile,
                      msg);
        free(ob);
        free(work);
        (void)fits_close_file(fptr,&status);
        return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
        (*ob+i)->y_t = work[i];

    /* No RA, Dec columns */

    for (i = 0; i < nrows; i++)
        (*ob+i)->ra = 0.0;
    for (i = 0; i < nrows; i++)
        (*ob+i)->dec = 0.0;

    /* No Class, Stat columns */

    for (i = 0; i < nrows; i++)
        (*ob+i)->class = 0.0;
    for (i = 0; i < nrows; i++)
        (*ob+i)->stat = 0.0;

    /* Read areal7, areal8 columns */

    (void)fits_read_col(fptr,TFLOAT,COL_AREAL7,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errstr,"Unable to read areal7 from table %s -- %s\n",
                      infile,msg);
        free(ob);
        free(work);
        (void)fits_close_file(fptr,&status);
        return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
        (*ob+i)->areal7 = work[i];
    (void)fits_read_col(fptr,TFLOAT,COL_AREAL8,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errstr,"Unable to read areal8 from table %s -- %s\n",
                      infile,msg);
        free(ob);
        free(work);
        (void)fits_close_file(fptr,&status);
        return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
        (*ob+i)->areal8 = work[i];
                      
    /* Read isoflux column */

    (void)fits_read_col(fptr,TFLOAT,COL_FLUXISO,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
        (void)sprintf(errstr,"Unable to read coreflux from table %s -- %s\n",
                      infile,msg);
        free(ob);
        free(work);
        (void)fits_close_file(fptr,&status);
        return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
        (*ob+i)->coreflux = work[i];

    /* Put some values into the remaining elements */

    for (i = 0; i < nrows; i++) {
        (*ob+i)->x_i = 0.0;
        (*ob+i)->y_i = 0.0;
    }

    /* Get out of here */

    free(work);
    (void)fits_close_file(fptr,&status);
    *nr = nrows;
    return(ERRCODE_OK);
}

extern int get_ncol_3(void) {
    return(NCOLS);
}

/*

$Log: create_table_3.c,v $
Revision 1.2  2014/07/18 08:46:11  jim
Fixed ell file so that in list mode it puts a circle where the original
coordinates were specified, rather than where it thinks it found the object

Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.10  2009/01/22 13:36:15  jim
removed unnecessary declaration

Revision 1.9  2008/04/15 18:57:49  jim
Removed code to do pixel flagging

Revision 1.8  2007/06/14 20:29:20  jim
fixed dodgy malloc

Revision 1.7  2007/06/14 09:34:59  jim
Modified so that if the give coordinates are off the image, then null results
are written to the tables.

Revision 1.6  2007/06/14 06:51:05  jim
Fixed to trap for images that 'move' too far from their initial positions

Revision 1.5  2007/06/04 10:34:02  jim
Modified to add list driven routines

Revision 1.4  2005/11/22 14:19:08  jim
iso_flux was accessing the wrong element in parmall

Revision 1.3  2005/07/08 11:14:53  jim
Modified so that the flux measurement is just the isophotal flux

Revision 1.2  2004/09/07 14:18:56  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:54:57  jim
New version for rewrite of imcore


*/
