/*

$Id: create_table_1.c,v 1.4 2014/07/30 08:12:45 jim Exp $

*/

#include <stdio.h>
#include <math.h>
#include "imcore.h"
#include "errcodes.h"
#include "util.h"
#include "floatmath.h"

/* Column numbers for each item to be stored */

#define COL_NUMBER      1
#define COL_FLUXISO     2
#define COL_FLUXTOTAL   3
#define COL_APFLUX1     4
#define COL_X           5
#define COL_Y           6
#define COL_SIGMA       7
#define COL_ELLIPT      8
#define COL_PA          9
#define COL_PEAKHEIGHT 10
#define COL_AREAL1     11
#define COL_AREAL2     12
#define COL_AREAL3     13
#define COL_AREAL4     14
#define COL_AREAL5     15
#define COL_AREAL6     16
#define COL_AREAL7     17
#define COL_AREAL8     18
#define COL_APFLUX2    19
#define COL_APFLUX3    20
#define COL_APFLUX4    21
#define COL_APFLUX5    22
#define COL_RA         23
#define COL_DEC        24
#define COL_CLASS      25
#define COL_STAT       26
#define COL_APFLUX6    27
#define COL_SKYLEV     28
#define COL_SKYRMS     29

/* Number of columns in the table */

#define NCOLS 32

/* Core radius extras */

#define NRADS 6
static float rmults[] = {0.5,1.0,M_SQRT2,2.0,2.0*M_SQRT2,4.0};
static int nrcore = 1;
static float apertures[NRADS];
static char *rcorelabs[] = {"1x","1/2x","sqrt(2)x","2x","2sqrt(2)x","4x"};
static int corecols[] = {COL_APFLUX1,COL_APFLUX2,COL_APFLUX3,COL_APFLUX4,
		       COL_APFLUX5,COL_APFLUX6};

/* Column definitions */

static char *ttype[NCOLS]={"No.","Isophotal_flux","Total_flux","Core_flux",
			   "X_coordinate","Y_coordinate","Gaussian_sigma",
			   "Ellipticity","Position_angle","Peak_height",
			   "Areal_1_profile","Areal_2_profile","Areal_3_profile",
			   "Areal_4_profile","Areal_5_profile","Areal_6_profile",
			   "Areal_7_profile","Areal_8_profile","Core1_flux",
			   "Core2_flux","Core3_flux","Core4_flux",
			   "RA","DEC","Classification","Statistic",
			   "Core5_flux","Skylev",
			   "Skyrms","Blank30","Blank31","Blank32"};
static char *tunit[NCOLS]={" ","Counts","Counts","Counts","Pixels","Pixels",
			   "Pixels"," ","Degrees","Counts","Pixels","Pixels",
			   "Pixels","Pixels","Pixels","Pixels","Pixels","Pixels",
			   "Counts","Counts","Counts","Counts",
			   "Radians","Radians","Flag","N-sigma",
			   "Counts","Counts","Counts","Blank",
			   "Blank","Blank"};
static char *tform[NCOLS]={"1E","1E","1E","1E","1E","1E","1E","1E","1E","1E",
		    "1E","1E","1E","1E","1E","1E","1E","1E","1E","1E","1E",
		    "1E","1E","1E","1E","1E","1E","1E","1E","1E","1E","1E"};

static int areal_cols[NAREAL] = {COL_AREAL1,COL_AREAL2,COL_AREAL3,COL_AREAL4,
				 COL_AREAL5,COL_AREAL6,COL_AREAL7,COL_AREAL8};

extern int tabinit_1(char *infile, char *outtab, char *errstr) {
    int retval,status,i;
    char errmsg[FLEN_STATUS],key[FLEN_KEYWORD],comment[FLEN_COMMENT];

    /* Call the generic routine to open a new output table */

    retval = tabinit_gen(infile,outtab,NCOLS,ttype,tunit,tform,errstr);
    if (retval != ERRCODE_OK)
	return(retval);
	
    /* Add in a few helpful comments */

    status = 0;
    for (i = 0; i < NRADS; i++) {
	(void)sprintf(key,"TTYPE%d",corecols[i]);
	(void)sprintf(comment,"Fitted flux within %s core radius",rcorelabs[i]);
	(void)fits_modify_comment(tptr,key,comment,&status);
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
	sprintf(errstr,"Unable to update table comments: %s (status = %d)",
		errmsg,status);
	return(ERRCODE_FILE_IO);
    }

    /* Get out of here */

    return(ERRCODE_OK);
}

extern int do_seeing_1(ap_t *ap, char *errstr) {
    int retval;

    /* Just call the generic seeing routine */

    retval = do_seeing_gen(ap,COL_ELLIPT,COL_PEAKHEIGHT,areal_cols,
			   errstr);
    if (retval != ERRCODE_OK)
	return(retval);

    /* Get out of here */

    return(ERRCODE_OK);
}
        

extern int process_results_1(ap_t *ap, char *errstr) {
    float momresults[8],ttotal,parmall[IMNUM][NPAR],ratio,cflux[NRADS*IMNUM];
    float sxx,syy,srr,sxy,ecc,temp,xx,theta,radeg,ell,nr,iso_flux,total_flux;
    float apflux1,apflux2,apflux3,apflux4,apflux5,yy,sigma,peak,areal1,apflux6;
    float areal2,areal3,areal4,areal5,areal6,areal7,areal8,aa,bb;
    float skylev,skyrms,badpix[IMNUM],avconf[IMNUM];
    int iareal[NAREAL],nbit,i,k,status,mbit,j;
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

    /* Work out the total flux */

    extend(ap,momresults[3],momresults[1],momresults[2],
	   momresults[4],momresults[5],momresults[6],
	   (float)iareal[0],momresults[7],&ttotal);
    ratio = MAX(ttotal,momresults[3])/momresults[3];

    /* Try and deblend the images if it is requested and justified */

    if (iareal[0] >= ap->mulpix && ap->icrowd)
        overlp(ap,parmall,&nbit,momresults[1],momresults[2],
	       momresults[3],iareal[0],momresults[7]);
    else
	nbit = 1;
    if (nbit == 1) {
	parmall[0][0] = momresults[3];
	parmall[0][1] = momresults[1];
	parmall[0][2] = momresults[2];
	parmall[0][3] = ap->thresh;
	for (i = 4; i < 8; i++) 
	    parmall[0][i] = momresults[i];
	for (i = 0; i < NAREAL; i++)
	    parmall[0][i+8] = (float)iareal[i];
    } else {
	mbit = 0;
	for (i = 0; i < nbit; i++) {
	    if (parmall[i][1] > 1.0 && parmall[i][1] < ap->lsiz &&
		parmall[i][2] > 1.0 && parmall[i][2] < ap->csiz) {
		for (j = 0; j < NPAR; j++) 
		    parmall[mbit][j] = parmall[i][j];
		mbit++;
	    }
	}
	nbit = mbit;
    } 

    /* Create a list of apertures */

    for (i = 0; i < NRADS; i++) 
	apertures[i] = rmults[i]*(ap->rcore);

    /* Initialise badpix and average confidence accumulators */

    for (i = 0; i < nbit; i++) {
	badpix[i] = 0.0;
	avconf[i] = 0.0;
    }

    /* Get the core fluxes in all apertures */

    phopt(ap,parmall,nbit,NRADS,apertures,cflux,badpix,nrcore,avconf);

    /* Massage the results and write them to the fits table */

    radeg = 180.0/M_PI;
    for (k = 0; k < nbit; k++) {
	sxx = parmall[k][4];
	sxy = parmall[k][5];
	syy = parmall[k][6];    
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
        nr = (float)nrows;
	iso_flux = parmall[k][0];
	total_flux = ratio*parmall[k][0];
	apflux1 = cflux[k*NRADS + 0];
	apflux2 = cflux[k*NRADS + 1];
	apflux3 = cflux[k*NRADS + 2];
	apflux4 = cflux[k*NRADS + 3];
	apflux5 = cflux[k*NRADS + 4];
	apflux6 = cflux[k*NRADS + 5];
	xx = parmall[k][1];
	yy = parmall[k][2];
	sigma = sqrt(srr);
	peak = parmall[k][7];
	areal1 = parmall[k][8];
	areal2 = parmall[k][9];
	areal3 = parmall[k][10];
	areal4 = parmall[k][11];
	areal5 = parmall[k][12];
	areal6 = parmall[k][13];
	areal7 = parmall[k][14];
	if (nbit > 1 && k == 0)
	    areal8 = 0.0;
	else 
	    areal8 = parmall[k][15];
	imcore_backest(ap,xx,yy,&skylev,&skyrms);

	/* Store away the results for this object */
    
        nr = (float)nrows;
        (void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&iso_flux,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_FLUXTOTAL,nrows,1,1,&total_flux,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_APFLUX1,nrows,1,1,&apflux2,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&sigma,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&ell,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&theta,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&peak,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&areal1,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&areal2,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&areal3,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&areal4,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&areal5,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&areal6,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&areal7,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&areal8,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_APFLUX2,nrows,1,1,&apflux1,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_APFLUX3,nrows,1,1,&apflux3,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_APFLUX4,nrows,1,1,&apflux4,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_APFLUX5,nrows,1,1,&apflux5,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_APFLUX6,nrows,1,1,&apflux6,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_SKYLEV,nrows,1,1,&skylev,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_SKYRMS,nrows,1,1,&skyrms,&status);
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
    }

    /* Get outta here */

    return(ERRCODE_OK);
}

extern int process_results_list_1(ap_t *ap, objstruct *objlist[], int nbit, 
                                  char *errstr) {
    int i,j;
    float parmall[IMNUM][NPAR],cflux[NRADS*IMNUM],badpix[IMNUM],avconf[IMNUM];
    float sxx,syy,sxy,tpk,*d,parrad,ra,dec,class,stat,tmn,radeg,total_flux;
    float xsq,ysq,xy,dx,dy,t,srr,ecc,temp,ell,xxx,theta,theta_ra,nr;
    float iso_flux,apflux1,apflux2,apflux3,apflux4,apflux5,apflux6;
    float sigma,areal1,areal2,areal3,areal4,areal5,areal6,areal7,areal8;
    float peak,xx,yy,xbar,ybar,skylev,skyrms,aa,bb,zero;
    int iareal[NAREAL],x1,y1,x2,y2,ix,iy,indy,ind,nup,status,k;
    long nrows;
    char errmsg[BUFSIZ];

    /* Check to see if this can even appear on the image */

    radeg = 180.0/M_PI;
    xx = objlist[0]->x_i;
    yy = objlist[0]->y_i;
    if (xx < 1 || yy < 1 || xx > ap->lsiz || yy > ap->csiz) {
	for (i = 0; i < nbit; i++) {
	    xx = objlist[i]->x_i;
	    yy = objlist[i]->y_i;
	    ra = (objlist[i]->ra)/radeg;
	    dec = (objlist[i]->dec)/radeg;
	    class = objlist[i]->class;
	    stat = objlist[i]->stat; 
	    zero = 0.0;
	    status = 0;
            (void)fits_get_num_rows(tptr,&nrows,&status);
            nrows++;
	    nr = (float)nrows;
	    (void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_FLUXTOTAL,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&zero,
				 &status);
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
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_SKYLEV,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_SKYRMS,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_RA,nrows,1,1,&ra,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_DEC,nrows,1,1,&dec,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_CLASS,nrows,1,1,&class,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_STAT,nrows,1,1,&stat,&status);
	    if (status != 0) {
		fits_get_errstatus(status,errmsg);
		sprintf(errstr,"Unable to write catalogue data: %s (status = %d)",
			errmsg,status);
		return(ERRCODE_FILE_IO);
	    }
	}
	return(ERRCODE_OK);
    }

    /* Create a list of apertures */

    for (i = 0; i < NRADS; i++) 
        apertures[i] = rmults[i]*(ap->rcore);

    /* Initialise the badpix and average confidence accumulators */

    for (i = 0; i < nbit; i++) {
        badpix[i] = 0.0;
	avconf[i] = 0.0;
    }

    /* Load the input parameters (these routines only use x and y) */

    for (i = 0; i < nbit; i++) {
        parmall[i][1] = objlist[i]->x_i;
        parmall[i][2] = objlist[i]->y_i;
    }

    /* Get the core fluxes in all apertures */

    phopt(ap,parmall,nbit,NRADS,apertures,cflux,badpix,nrcore,avconf);

    /* Ok, here's a big loop where we try to refine the final parameters. */

    parrad = ap->rcore + 0.5;
    d = ap->data;
    for (k = 0; k < nbit; k++) {
        xx = objlist[k]->x_i;
        yy = objlist[k]->y_i;
        ra = (objlist[k]->ra)/radeg;
        dec = (objlist[k]->dec)/radeg;
        class = objlist[k]->class;
        stat = objlist[k]->stat;

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
            for (j = 0; j < NRADS; j++)
                cflux[k*NRADS+j] = 0.0;
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
        theta_ra = theta/radeg;

        /* Create a list of values */

        status = 0;
        (void)fits_get_num_rows(tptr,&nrows,&status);
        nrows++;
        nr = (float)nrows;
        iso_flux = tmn;
	total_flux = tmn;
        apflux1 = cflux[k*NRADS + 0];
        apflux2 = cflux[k*NRADS + 1];
        apflux3 = cflux[k*NRADS + 2];
        apflux4 = cflux[k*NRADS + 3];
        apflux5 = cflux[k*NRADS + 4];
        apflux6 = cflux[k*NRADS + 5];
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
        (void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&iso_flux,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_FLUXTOTAL,nrows,1,1,&total_flux,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1,nrows,1,1,&apflux2,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&sigma,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&ell,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&theta,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&peak,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&areal1,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&areal2,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&areal3,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&areal4,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&areal5,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&areal6,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&areal7,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&areal8,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2,nrows,1,1,&apflux1,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3,nrows,1,1,&apflux3,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4,nrows,1,1,&apflux4,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5,nrows,1,1,&apflux5,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6,nrows,1,1,&apflux6,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_SKYLEV,nrows,1,1,&skylev,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_SKYRMS,nrows,1,1,&skyrms,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_RA,nrows,1,1,&ra,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_DEC,nrows,1,1,&dec,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_CLASS,nrows,1,1,&class,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_STAT,nrows,1,1,&stat,&status);
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


extern int tabclose_1(ap_t *ap, char *errstr) {

    return(tabclose_gen(ap,errstr));
}

extern int readtab_list_1(char *infile, objstruct **ob, int *nr, 
			  char *errstr) {
    int status,anynul;
    long nrows,i;
    float *work;
    double degrad;
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

    /* Read RA, Dec columns */

    degrad = 180.0/M_PI;
    (void)fits_read_col(fptr,TFLOAT,COL_RA,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
	(void)sprintf(errstr,"Unable to read RA from table %s -- %s\n",infile,
		      msg);
	free(ob);
	free(work);
	(void)fits_close_file(fptr,&status);
	return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
	(*ob+i)->ra = degrad*(double)work[i];
    (void)fits_read_col(fptr,TFLOAT,COL_DEC,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
	(void)sprintf(errstr,"Unable to read Dec from table %s -- %s\n",infile,
		      msg);
	free(ob);
	free(work);
	(void)fits_close_file(fptr,&status);
	return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
	(*ob+i)->dec = degrad*(double)work[i];

    /* Read Class, Stat columns */

    (void)fits_read_col(fptr,TFLOAT,COL_CLASS,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
	(void)sprintf(errstr,"Unable to read Class from table %s -- %s\n",
		      infile,msg);
	free(ob);
	free(work);
	(void)fits_close_file(fptr,&status);
	return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
	(*ob+i)->class = work[i];
    (void)fits_read_col(fptr,TFLOAT,COL_STAT,1,1,nrows,NULL,work,&anynul,&status);
    if (status != 0) {
        fits_get_errstatus(status,msg);
	(void)sprintf(errstr,"Unable to read stat from table %s -- %s\n",
		      infile,msg);
	free(ob);
	free(work);
	(void)fits_close_file(fptr,&status);
	return(ERRCODE_FILE_IO);
    }
    for (i = 0; i < nrows; i++)
	(*ob+i)->stat = work[i];

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
		      
    /* Read coreflux column */

    (void)fits_read_col(fptr,TFLOAT,COL_APFLUX1,1,1,nrows,NULL,work,&anynul,&status);
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

extern int get_ncol_1(void) {
    return(NCOLS);
}

/*

$Log: create_table_1.c,v $
Revision 1.4  2014/07/30 08:12:45  jim
Fixed long standing bug that means that in list mode the output equatorial
coordinates are written to the FITS tables in degrees rather than in radians
as advertised

Revision 1.3  2014/07/18 08:46:11  jim
Fixed ell file so that in list mode it puts a circle where the original
coordinates were specified, rather than where it thinks it found the object

Revision 1.2  2012/10/01 19:02:59  jim
Fixed duplicate column names

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.12  2009/12/17 11:32:38  jim
Modified to use new version of phopt

Revision 1.11  2008/04/15 18:57:49  jim
Removed code to do pixel flagging

Revision 1.10  2007/06/14 20:29:20  jim
fixed dodgy malloc

Revision 1.9  2007/06/14 09:34:59  jim
Modified so that if the give coordinates are off the image, then null results
are written to the tables.

Revision 1.8  2007/06/14 06:51:05  jim
Fixed to trap for images that 'move' too far from their initial positions

Revision 1.7  2007/06/07 14:06:10  jim
Fixed readtab routine so that RA and DEC are converted to degrees

Revision 1.6  2007/06/04 10:34:02  jim
Modified to add list driven routines

Revision 1.5  2006/12/13 12:36:21  jim
Modified to remove merged objects near the edges

Revision 1.4  2006/08/11 12:58:57  jim
Fixed so that 8th areal profile is set to zero for the start of a blend

Revision 1.3  2005/08/26 04:46:20  jim
Modified to add new radii and error estimates

Revision 1.2  2004/09/07 14:18:56  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:54:57  jim
New version for rewrite of imcore


*/
