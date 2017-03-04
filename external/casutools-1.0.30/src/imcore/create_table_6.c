/*

$Id: create_table_6.c,v 1.7 2014/07/30 08:12:45 jim Exp $

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "imcore.h"
#include "imcore_radii.h"
#include "errcodes.h"
#include "util.h"
#include "floatmath.h"

#define COL_NUMBER      1
#define COL_FLUXISO     2
#define COL_X           3
#define COL_XERR        4
#define COL_Y           5
#define COL_YERR        6
#define COL_SIGMA       7
#define COL_ELLIPT      8
#define COL_PA          9
#define COL_AREAL1     10
#define COL_AREAL2     11
#define COL_AREAL3     12
#define COL_AREAL4     13
#define COL_AREAL5     14
#define COL_AREAL6     15
#define COL_AREAL7     16
#define COL_AREAL8     17
#define COL_PEAKHEIGHT 18
#define COL_PKHTERR    19
#define COL_APFLUX1    20
#define COL_APFLUX1ERR 21
#define COL_APFLUX2    22
#define COL_APFLUX2ERR 23
#define COL_APFLUX3    24
#define COL_APFLUX3ERR 25
#define COL_APFLUX4    26
#define COL_APFLUX4ERR 27
#define COL_APFLUX5    28
#define COL_APFLUX5ERR 29
#define COL_APFLUX6    30
#define COL_APFLUX6ERR 31
#define COL_APFLUX7    32
#define COL_APFLUX7ERR 33
#define COL_APFLUX8    34
#define COL_APFLUX8ERR 35
#define COL_APFLUX9    36
#define COL_APFLUX9ERR 37
#define COL_APFLUX10    38
#define COL_APFLUX10ERR 39
#define COL_APFLUX11    40
#define COL_APFLUX11ERR 41
#define COL_APFLUX12    42
#define COL_APFLUX12ERR 43
#define COL_APFLUX13    44
#define COL_APFLUX13ERR 45
#define COL_PETRAD      46
#define COL_KRONRAD     47
#define COL_HALFRAD     48
#define COL_PETFLUX     49
#define COL_PETFLUXERR  50
#define COL_KRONFLUX    51
#define COL_KRONFLUXERR 52
#define COL_HALFFLUX    53
#define COL_HALFFLUXERR 54
#define COL_ERRFLAG     55
#define COL_SKYLEVEL    56
#define COL_SKYSIGMA    57
#define COL_AVCONF      58
#define COL_RA          59
#define COL_DEC         60
#define COL_CLASS       61
#define COL_STAT        62
 
/* Number of columns in the table */
 
#define NCOLS 80

static char *ttype[NCOLS]={"Sequence_number","Isophotal_flux",
			   "X_coordinate","X_coordinate_err",
			   "Y_coordinate","Y_coordinate_err",
			   "Gaussian_sigma","Ellipticity","Position_angle",
			   "Areal_1_profile","Areal_2_profile","Areal_3_profile",
			   "Areal_4_profile","Areal_5_profile","Areal_6_profile",
			   "Areal_7_profile","Areal_8_profile",
			   "Peak_height","Peak_height_err",
			   "Aper_flux_1","Aper_flux_1_err",
			   "Aper_flux_2","Aper_flux_2_err",
			   "Aper_flux_3","Aper_flux_3_err",
			   "Aper_flux_4","Aper_flux_4_err",
			   "Aper_flux_5","Aper_flux_5_err",
			   "Aper_flux_6","Aper_flux_6_err",
			   "Aper_flux_7","Aper_flux_7_err",
			   "Aper_flux_8","Aper_flux_8_err",
			   "Aper_flux_9","Aper_flux_9_err",
			   "Aper_flux_10","Aper_flux_10_err",
			   "Aper_flux_11","Aper_flux_11_err",
			   "Aper_flux_12","Aper_flux_12_err",
			   "Aper_flux_13","Aper_flux_13_err",
			   "Petr_radius","Kron_radius","Half_radius",
			   "Petr_flux","Petr_flux_err",
			   "Kron_flux","Kron_flux_err","Half_flux","Half_flux_err",
			   "Error_bit_flag","Sky_level","Sky_rms",
			   "Av_conf",
			   "RA","DEC","Classification","Statistic",
			   "Blank63","Blank64","Blank65","Blank66","Blank67",
			   "Blank68","Blank69","Blank70","Blank71","Blank72",
			   "Blank73","Blank74","Blank75","Blank76","Blank77",
			   "Blank78","Blank79","Blank80"};

static char *tunit[NCOLS]={"Number","ADU",
			   "Pixels","Pixels",
			   "Pixels","Pixels",
			   "Pixels","Number","Degrees",
			   "Pixels","Pixels","Pixels",
			   "Pixels","Pixels","Pixels",
			   "Pixels","Pixels",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "ADU","ADU",
			   "Pixels","Pixels","Pixels",
			   "ADU","ADU",
			   "ADU","ADU","ADU","ADU",
			   "Number","ADU","ADU","Number",
			   "Degrees","Degrees","Flag","N-sigma",
			   "Blank63","Blank64","Blank65","Blank66","Blank67",
			   "Blank68","Blank69","Blank70","Blank71","Blank72",
			   "Blank73","Blank74","Blank75","Blank76","Blank77",
			   "Blank78","Blank79","Blank80"};

static char *tform[NCOLS]={"1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E","1E",
			   "1E","1E","1E",
			   "1E","1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E",
			   "1E","1E","1E",
			   "1E","1E",
			   "1E","1E","1E","1E",
			   "1E","1E","1E","1E",
			   "1D","1D","1E","1E","1E","1E",
			   "1E","1E","1E","1E","1E","1E",
			   "1E","1E","1E","1E","1E","1E",
			   "1E","1E","1E","1E"};

#define NRADS 13
static float rmults[] = {0.5,1.0/M_SQRT2,1.0,M_SQRT2,2.0,2.0*M_SQRT2,4.0,
			 5.0,6.0,7.0,8.0,10,12.0};
static int nrcore = 2;
static int n2rcore = 4;
static float apertures[NRADS];
static char *rcorelabs[NRADS] = {"1/2*","1/sqrt(2)*","1*","sqrt(2)*","2*",
				 "2*sqrt(2)*","4*","5*","6*","7*",
				 "8*","10*","12*"};
static int corecols[NRADS] = {COL_APFLUX1,COL_APFLUX2,COL_APFLUX3,COL_APFLUX4,
			      COL_APFLUX5,COL_APFLUX6,COL_APFLUX7,COL_APFLUX8,
			      COL_APFLUX9,COL_APFLUX10,COL_APFLUX11,
			      COL_APFLUX12,COL_APFLUX13};
static int coreerrcols[NRADS] = {COL_APFLUX1ERR,COL_APFLUX2ERR,COL_APFLUX3ERR,
				 COL_APFLUX4ERR,COL_APFLUX5ERR,COL_APFLUX6ERR,
				 COL_APFLUX7ERR,COL_APFLUX8ERR,COL_APFLUX9ERR,
				 COL_APFLUX10ERR,COL_APFLUX11ERR,
				 COL_APFLUX12ERR,COL_APFLUX13ERR};
static int areal_cols[NAREAL] = {COL_AREAL1,COL_AREAL2,COL_AREAL3,COL_AREAL4,
				 COL_AREAL5,COL_AREAL6,COL_AREAL7,COL_AREAL8};

static void change_key(fitsfile *fptr, int mef, char *oldkey, char *newkey);
static int ismef(char *infile);


extern int tabinit_6(char *infile, char *outtab, char *errstr) {
    int retval,status,i,mef;
    char errmsg[FLEN_STATUS],comment[FLEN_COMMENT],key[FLEN_KEYWORD];
    char oldkey[FLEN_KEYWORD];

    /* Call the generic routine to open a new output table */

    retval = tabinit_genmef(infile,outtab,NCOLS,ttype,tunit,tform,errstr);
    if (retval != ERRCODE_OK)
	return(retval);
	
    /* Add in a few helpful comments */

    status = 0;
    for (i = 0; i < 13; i++) {
	(void)sprintf(key,"TTYPE%d",corecols[i]);
	(void)sprintf(comment,"Fitted flux within %s core radius",
		      rcorelabs[i]);
	(void)fits_modify_comment(tptr,key,comment,&status);
	(void)sprintf(key,"TTYPE%d",coreerrcols[i]);
	(void)sprintf(comment,"Error in fitted flux within %s core radius",
		      rcorelabs[i]);
	(void)fits_modify_comment(tptr,key,comment,&status);
    }
    if (status != 0) {
        fits_get_errstatus(status,errmsg);
	sprintf(errstr,"Unable to update table comments: %s (status = %d)",
		errmsg,status);
	return(ERRCODE_FILE_IO);
    }

    /* Is the input file a MEF or SEF? */

    mef = ismef(infile);

    /* Rename the image WCS keywords to the FITS standard table WCS 
       keywords */

    (void)strncpy(oldkey,"CTYPE1",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCTYP%d",COL_X);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CTYPE2",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCTYP%d",COL_Y);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CRVAL1",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCRVL%d",COL_X);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CRVAL2",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCRVL%d",COL_Y);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CRPIX1",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCRPX%d",COL_X);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CRPIX2",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCRPX%d",COL_Y);
    change_key(tptr,mef,oldkey,key);
    for (i = 1; i <= 5; i++) {
        (void)snprintf(oldkey,FLEN_KEYWORD,"PV2_%d",i);
        (void)snprintf(key,FLEN_KEYWORD,"TV%d_%d",COL_Y,i);
        change_key(tptr,mef,oldkey,key);
    }
    (void)strncpy(oldkey,"CD1_1",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TC%d_%d",COL_X,COL_X);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CD1_2",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TC%d_%d",COL_X,COL_Y);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CD2_1",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TC%d_%d",COL_Y,COL_X);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CD2_2",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TC%d_%d",COL_Y,COL_Y);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CDELT1",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCDLT%d",COL_X);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CDELT2",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCDLT%d",COL_Y);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CROTA1",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCROT%d",COL_X);
    change_key(tptr,mef,oldkey,key);
    (void)strncpy(oldkey,"CROTA2",FLEN_KEYWORD);
    (void)snprintf(key,FLEN_KEYWORD,"TCROT%d",COL_Y);
    change_key(tptr,mef,oldkey,key);

    /* Get out of here */

    return(ERRCODE_OK);
}

static void change_key(fitsfile *fptr, int mef, char *oldkey, char *newkey) {
    int status,hdutype;
    char card[81];

    /* Modify the keyword name. It's not an error if it doesn't exist */

    status = 0;
    if (! mef) {
	(void)fits_movabs_hdu(fptr,1,&hdutype,&status);
	(void)fits_read_card(fptr,oldkey,card,&status);
	(void)fits_delete_key(fptr,oldkey,&status);
	if (status != 0) {
	    status = 0;
	    (void)fits_movabs_hdu(fptr,2,&hdutype,&status);
	    return;
	}
	(void)fits_movabs_hdu(fptr,2,&hdutype,&status);
	(void)fits_update_card(fptr,oldkey,card,&status);
    }
    (void)fits_modify_name(fptr,oldkey,newkey,&status);
    return;
}

extern int do_seeing_6(ap_t *ap, char *errstr) {
    int retval;
 
    /* Just call the generic seeing routine */
 
    retval = do_seeing_gen(ap,COL_ELLIPT,COL_PEAKHEIGHT,areal_cols,
                           errstr);
    if (retval != ERRCODE_OK)
        return(retval);
 
    /* Get out of here */
 
    return(ERRCODE_OK);
}

extern int process_results_6(ap_t *ap, char *errstr) {
    float momresults[8],ttotal,parmall[IMNUM][NPAR],ratio,cflux[NRADS*IMNUM];
    float sxx,syy,srr,sxy,ecc,temp,xx,theta,radeg,ell,nr,iso_flux,total_flux;
    float apflux1,apflux2,apflux3,apflux4,apflux5,yy,sigma,peak,areal1,apflux6;
    float apflux7,apflux8,apflux9,apflux10,apflux11,apflux12,apflux13,zero;
    float areal2,areal3,areal4,areal5,areal6,areal7,areal8,aa,bb;
    float skylev,skyrms,half_rad[IMNUM],half_flux[IMNUM],kron_flux[IMNUM];
    float petr_flux[IMNUM],kron_rad[IMNUM],petr_rad[IMNUM],badpix[IMNUM];
    float theta_ra,skyvar[IMNUM],cc,dd,sigsq,xxe,yye,peake,kron_fluxe;
    float apflux1e,apflux2e,apflux3e,apflux4e,apflux5e,apflux6e,apflux7e;
    float apflux8e,apflux9e,apflux10e,apflux11e,apflux12e,apflux13e,half_fluxe;
    float petr_fluxe,r,avconf[IMNUM],rcore_area;
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

    for (i = 0; i < NRADS; i++) {
        apertures[i] = rmults[i]*(ap->rcore);
	skyvar[i] = M_PI*apertures[i]*apertures[i];
    }
    rcore_area = M_PI*pow(ap->rcore,2.0);

    /* Initialise the badpix and average confidence accumulators */

    for (i = 0; i < nbit; i++) {
	badpix[i] = 0.0;
	avconf[i] = 0.0;
    }

    /* Get the core fluxes in all apertures */

    phopt(ap,parmall,nbit,NRADS,apertures,cflux,badpix,nrcore,avconf);
    for (i = 0; i < nbit; i++) 
	avconf[i] /= rcore_area;

    /* Get half-flux radius for all images */

    for (k = 0; k < nbit; k++) {
        half_flux[k] = 0.5*(MAX(parmall[k][0],cflux[k*NRADS+n2rcore]));
        half_rad[k] = imcore_halflight(apertures,cflux+k*NRADS,half_flux[k],
				       parmall[k][7],NRADS);
    }

    /* Get Kron radius for all images and get the flux */

    for (k = 0; k < nbit; k++) {
	areal1 = parmall[k][8];
	kron_rad[k] = imcore_kronrad(areal1,apertures,cflux+k*NRADS,NRADS);
    }
    imcore_flux(ap,parmall,nbit,kron_rad,kron_flux,NRADS,apertures,cflux);

    /* Get Petrosian radius for all images and get the flux */

    for (k = 0; k < nbit; k++) {
	areal1 = parmall[k][8];
	petr_rad[k] = imcore_petrad(areal1,apertures,cflux+k*NRADS,NRADS);
    }
    imcore_flux(ap,parmall,nbit,petr_rad,petr_flux,NRADS,apertures,cflux); 

    /* Massage the results and write them to the fits table */

    sigsq = powf(ap->sigma,2.0);
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
	aa = 1.0;
	bb = aa*(1.0 - ell);
        xx = 0.5*(1.0+ecc)*srr-sxx;
        if(xx == 0.0)
            theta = 0.0;
        else
            theta = 90.0 - radeg*atanf(sxy/xx);
	theta_ra = theta/radeg;
	cc = (1.0 + ecc)*pow(cos(theta_ra),2.0) + (1.0 - ecc)*pow(sin(theta_ra),2.0);
	dd = (1.0 + ecc)*pow(sin(theta_ra),2.0) + (1.0 - ecc)*pow(cos(theta_ra),2.0);

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
        apflux7 = cflux[k*NRADS + 6];
        apflux8 = cflux[k*NRADS + 7];
        apflux9 = cflux[k*NRADS + 8];
        apflux10 = cflux[k*NRADS + 9];
        apflux11 = cflux[k*NRADS + 10];
        apflux12 = cflux[k*NRADS + 11];
        apflux13 = cflux[k*NRADS + 12];
        peak = parmall[k][7];
        xx = parmall[k][1];
        xxe = sqrtf((2.0*sigsq/(M_PI*peak*peak)) + cc/(2.0*M_PI*gain*peak) + 
                    0.0001);
        yy = parmall[k][2];
        yye = sqrtf((2.0*sigsq/(M_PI*peak*peak)) + dd/(2.0*M_PI*gain*peak) + 
                    0.0001);
        sigma = sqrt(srr);
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
        peake = sqrtf(MAX(peak,0.0)/gain + sigsq + skyrms*skyrms);
        kron_fluxe = sqrt(MAX(kron_flux[k],0.0)/gain + 
                          (sigsq + skyrms*skyrms)*M_PI*powf(kron_rad[k],2.0));
        half_fluxe = sqrt(MAX(half_flux[k],0.0)/gain + 
                         (sigsq + skyrms*skyrms)*M_PI*powf(half_rad[k],2.0));
        petr_fluxe = sqrt(MAX(petr_flux[k],0.0)/gain + 
                          (sigsq + skyrms*skyrms)*M_PI*powf(petr_rad[k],2.0));
        apflux1e = sqrt(MAX(0.0,apflux1)/gain + skyvar[0]*(sigsq + skyrms*skyrms
));
        apflux2e = sqrt(MAX(0.0,apflux2)/gain + skyvar[1]*(sigsq + skyrms*skyrms
));
	apflux3e = sqrt(MAX(0.0,apflux3)/gain + skyvar[2]*(sigsq + skyrms*skyrms));
	apflux4e = sqrt(MAX(0.0,apflux4)/gain + skyvar[3]*(sigsq + skyrms*skyrms));
	apflux5e = sqrt(MAX(0.0,apflux5)/gain + skyvar[4]*(sigsq + skyrms*skyrms));
	apflux6e = sqrt(MAX(0.0,apflux6)/gain + skyvar[5]*(sigsq + skyrms*skyrms));
	apflux7e = sqrt(MAX(0.0,apflux7)/gain + skyvar[6]*(sigsq + skyrms*skyrms));
	apflux8e = sqrt(MAX(0.0,apflux8)/gain + skyvar[7]*(sigsq + skyrms*skyrms));
	apflux9e = sqrt(MAX(0.0,apflux9)/gain + skyvar[8]*(sigsq + skyrms*skyrms));
	apflux10e = sqrt(MAX(0.0,apflux10)/gain + skyvar[9]*(sigsq + skyrms*skyrms));
	apflux11e = sqrt(MAX(0.0,apflux11)/gain + skyvar[10]*(sigsq + skyrms*skyrms));
	apflux12e = sqrt(MAX(0.0,apflux12)/gain + skyvar[11]*(sigsq + skyrms*skyrms));
	apflux13e = sqrt(MAX(0.0,apflux13)/gain + skyvar[12]*(sigsq + skyrms*skyrms));

        /* Store away the results for this object */

	zero = 0.0;
        nr = (float)nrows;
        (void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&iso_flux,
                             &status); 
        (void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_XERR,nrows,1,1,&xxe,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_YERR,nrows,1,1,&yye,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&sigma,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&ell,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&theta,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&areal1,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&areal2,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&areal3,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&areal4,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&areal5,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&areal6,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&areal7,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&areal8,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&peak,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_PKHTERR,nrows,1,1,&peake,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1,nrows,1,1,&apflux1,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1ERR,nrows,1,1,&apflux1e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2,nrows,1,1,&apflux2,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2ERR,nrows,1,1,&apflux2e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3,nrows,1,1,&apflux3,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3ERR,nrows,1,1,&apflux3e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4,nrows,1,1,&apflux4,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4ERR,nrows,1,1,&apflux4e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5,nrows,1,1,&apflux5,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5ERR,nrows,1,1,&apflux5e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6,nrows,1,1,&apflux6,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6ERR,nrows,1,1,&apflux6e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX7,nrows,1,1,&apflux7,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX7ERR,nrows,1,1,&apflux7e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX8,nrows,1,1,&apflux8,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX8ERR,nrows,1,1,&apflux8e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX9,nrows,1,1,&apflux9,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX9ERR,nrows,1,1,&apflux9e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX10,nrows,1,1,&apflux10,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX10ERR,nrows,1,1,&apflux10e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX11,nrows,1,1,&apflux11,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX11ERR,nrows,1,1,&apflux11e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX12,nrows,1,1,&apflux12,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX12ERR,nrows,1,1,&apflux12e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX13,nrows,1,1,&apflux13,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX13ERR,nrows,1,1,&apflux13e,
                             &status);
        r = 0.5*petr_rad[k];
	(void)fits_write_col(tptr,TFLOAT,COL_PETRAD,nrows,1,1,&r,&status);
	r = 0.5*kron_rad[k];
	(void)fits_write_col(tptr,TFLOAT,COL_KRONRAD,nrows,1,1,&r,&status);
	r = half_rad[k];
	(void)fits_write_col(tptr,TFLOAT,COL_HALFRAD,nrows,1,1,&r,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_PETFLUX,nrows,1,1,(petr_flux+k),
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_PETFLUXERR,nrows,1,1,&petr_fluxe,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_KRONFLUX,nrows,1,1,(kron_flux+k),
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_KRONFLUXERR,nrows,1,1,&kron_fluxe,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_HALFFLUX,nrows,1,1,(half_flux+k),
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_HALFFLUXERR,nrows,1,1,&half_fluxe,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_ERRFLAG,nrows,1,1,(badpix+k),
			     &status);
        (void)fits_write_col(tptr,TFLOAT,COL_SKYLEVEL,nrows,1,1,&skylev,
			     &status);
        (void)fits_write_col(tptr,TFLOAT,COL_SKYSIGMA,nrows,1,1,&skyrms,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_AVCONF,nrows,1,1,(avconf+k),
			     &status);
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

extern int process_results_list_6(ap_t *ap, objstruct *objlist[], int nbit, 
				  char *errstr) {
    int i,j;
    float parmall[IMNUM][NPAR],cflux[NRADS*IMNUM],skyvar[IMNUM],badpix[IMNUM];
    float peak,half_rad[IMNUM],half_flux[IMNUM],kron_rad[IMNUM],avconf[IMNUM];
    float kron_flux[IMNUM],petr_rad[IMNUM],petr_flux[IMNUM],xx,yy,xbar,ybar;
    float sxx,syy,sxy,tpk,*d,parrad,ra,dec,class,stat,tmn,radeg,rcore_area;
    float xsq,ysq,xy,dx,dy,t,srr,ecc,temp,ell,xxx,theta,theta_ra,cc,dd,nr;
    float iso_flux,apflux1,apflux2,apflux3,apflux4,apflux5,apflux6,apflux7;
    float apflux8,apflux9,apflux10,apflux11,apflux12,apflux13,sigsq,xxe,yye;
    float sigma,areal1,areal2,areal3,areal4,areal5,areal6,areal7,areal8;
    float skylev,skyrms,peake,kron_fluxe,half_fluxe,petr_fluxe,apflux1e;
    float apflux2e,apflux3e,apflux4e,apflux5e,apflux6e,apflux7e,apflux8e;
    float apflux9e,apflux10e,apflux11e,apflux12e,apflux13e,zero,r,aa,bb;
    int iareal[NAREAL],x1,y1,x2,y2,ix,iy,indy,ind,nup,status,k;
    long nrows;
    char errmsg[BUFSIZ];

    /* Check to see if this can even appear on the image */

    xx = objlist[0]->x_i;
    yy = objlist[0]->y_i;
    radeg = 180.0/M_PI;
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
	    bb = M_PI*powf(ap->rcore,2.0);
	    (void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&zero,
				 &status); 
	    (void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_XERR,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_YERR,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PKHTERR,nrows,1,1,&zero,&status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX7,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX7ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX8,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX8ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX9,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX9ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX10,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX10ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX11,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX11ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX12,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX12ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX13,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_APFLUX13ERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PETRAD,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_KRONRAD,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_HALFRAD,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PETFLUX,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_PETFLUXERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_KRONFLUX,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_KRONFLUXERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_HALFFLUX,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_HALFFLUXERR,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_ERRFLAG,nrows,1,1,&bb,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_SKYLEVEL,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_SKYSIGMA,nrows,1,1,&zero,
				 &status);
	    (void)fits_write_col(tptr,TFLOAT,COL_AVCONF,nrows,1,1,&zero,
				 &status);
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
	}
        return(ERRCODE_OK);
    }

    /* Create a list of apertures */

    for (i = 0; i < NRADS; i++) {
        apertures[i] = rmults[i]*(ap->rcore);
	skyvar[i] = M_PI*apertures[i]*apertures[i];
    }
    rcore_area = M_PI*pow(ap->rcore,2.0);

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
    for (i = 0; i < nbit; i++)
	avconf[i] /= rcore_area;

    /* Peak and 1st areal level */

    peak = expf(1.0)*(ap->thresh);
    areal1 = skyvar[NRADS-1];

    /* Get half-flux radius for all images */

    for (k = 0; k < nbit; k++) {
        half_flux[k] = 0.5*(MAX(parmall[k][0],cflux[k*NRADS+n2rcore]));
        half_rad[k] = imcore_halflight(apertures,cflux+k*NRADS,half_flux[k],
				       parmall[k][7],NRADS);
    }

    /* Get Kron radius for all images and get the flux */

    for (k = 0; k < nbit; k++) 
	kron_rad[k] = imcore_kronrad(areal1,apertures,cflux+k*NRADS,NRADS);
    imcore_flux(ap,parmall,nbit,kron_rad,kron_flux,NRADS,apertures,cflux);

    /* Get Petrosian radius for all images and get the flux */

    for (k = 0; k < nbit; k++) 
	petr_rad[k] = imcore_petrad(areal1,apertures,cflux+k*NRADS,NRADS);
    imcore_flux(ap,parmall,nbit,petr_rad,petr_flux,NRADS,apertures,cflux);

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
	cc = (1.0 + ecc)*pow(cos(theta_ra),2.0) + 
	    (1.0 - ecc)*pow(sin(theta_ra),2.0);
	dd = (1.0 + ecc)*pow(sin(theta_ra),2.0) + 
	    (1.0 - ecc)*pow(cos(theta_ra),2.0);
			
        /* Create a list of values */

        status = 0;
        (void)fits_get_num_rows(tptr,&nrows,&status);
        nrows++;
        nr = (float)nrows;
        iso_flux = tmn;
        apflux1 = cflux[k*NRADS + 0];
        apflux2 = cflux[k*NRADS + 1];
        apflux3 = cflux[k*NRADS + 2];
        apflux4 = cflux[k*NRADS + 3];
        apflux5 = cflux[k*NRADS + 4];
        apflux6 = cflux[k*NRADS + 5];
        apflux7 = cflux[k*NRADS + 6];
        apflux8 = cflux[k*NRADS + 7];
        apflux9 = cflux[k*NRADS + 8];
        apflux10 = cflux[k*NRADS + 9];
        apflux11 = cflux[k*NRADS + 10];
        apflux12 = cflux[k*NRADS + 11];
        apflux13 = cflux[k*NRADS + 12];
        peak = tpk;
        xx = xbar;
	sigsq = (ap->sigma)*(ap->sigma);
	if (peak <= 0.0) 
	    xxe = 0.0;
	else
            xxe = sqrtf((2.0*sigsq/(M_PI*peak*peak)) + 
			cc/(2.0*M_PI*gain*peak) + 0.0001);
        yy = ybar;
	if (peak <= 0.0) 
	    yye = 0.0;
	else 
            yye = sqrtf((2.0*sigsq/(M_PI*peak*peak)) + 
			dd/(2.0*M_PI*gain*peak) + 0.0001);
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
        peake = sqrtf(MAX(peak,0.0)/gain + sigsq + skyrms*skyrms);
        kron_fluxe = sqrt(MAX(kron_flux[k],0.0)/gain + 
                          (sigsq + skyrms*skyrms)*M_PI*powf(kron_rad[k],2.0));
        half_fluxe = sqrt(MAX(half_flux[k],0.0)/gain + 
                         (sigsq + skyrms*skyrms)*M_PI*powf(half_rad[k],2.0));
        petr_fluxe = sqrt(MAX(petr_flux[k],0.0)/gain + 
                          (sigsq + skyrms*skyrms)*M_PI*powf(petr_rad[k],2.0));
        apflux1e = sqrt(MAX(0.0,apflux1)/gain + skyvar[0]*(sigsq + skyrms*skyrms
));
        apflux2e = sqrt(MAX(0.0,apflux2)/gain + skyvar[1]*(sigsq + skyrms*skyrms
));
	apflux3e = sqrt(MAX(0.0,apflux3)/gain + skyvar[2]*(sigsq + skyrms*skyrms));
	apflux4e = sqrt(MAX(0.0,apflux4)/gain + skyvar[3]*(sigsq + skyrms*skyrms));
	apflux5e = sqrt(MAX(0.0,apflux5)/gain + skyvar[4]*(sigsq + skyrms*skyrms));
	apflux6e = sqrt(MAX(0.0,apflux6)/gain + skyvar[5]*(sigsq + skyrms*skyrms));
	apflux7e = sqrt(MAX(0.0,apflux7)/gain + skyvar[6]*(sigsq + skyrms*skyrms));
	apflux8e = sqrt(MAX(0.0,apflux8)/gain + skyvar[7]*(sigsq + skyrms*skyrms));
	apflux9e = sqrt(MAX(0.0,apflux9)/gain + skyvar[8]*(sigsq + skyrms*skyrms));
	apflux10e = sqrt(MAX(0.0,apflux10)/gain + skyvar[9]*(sigsq + skyrms*skyrms));
	apflux11e = sqrt(MAX(0.0,apflux11)/gain + skyvar[10]*(sigsq + skyrms*skyrms));
	apflux12e = sqrt(MAX(0.0,apflux12)/gain + skyvar[11]*(sigsq + skyrms*skyrms));
	apflux13e = sqrt(MAX(0.0,apflux13)/gain + skyvar[12]*(sigsq + skyrms*skyrms));

        /* Store away the results for this object */

	zero = 0.0;
        nr = (float)nrows;
        (void)fits_write_col(tptr,TFLOAT,COL_NUMBER,nrows,1,1,&nr,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_FLUXISO,nrows,1,1,&iso_flux,
                             &status); 
        (void)fits_write_col(tptr,TFLOAT,COL_X,nrows,1,1,&xx,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_XERR,nrows,1,1,&xxe,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_Y,nrows,1,1,&yy,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_YERR,nrows,1,1,&yye,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_SIGMA,nrows,1,1,&sigma,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_ELLIPT,nrows,1,1,&ell,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_PA,nrows,1,1,&theta,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL1,nrows,1,1,&areal1,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL2,nrows,1,1,&areal2,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL3,nrows,1,1,&areal3,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL4,nrows,1,1,&areal4,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL5,nrows,1,1,&areal5,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL6,nrows,1,1,&areal6,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL7,nrows,1,1,&areal7,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_AREAL8,nrows,1,1,&areal8,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_PEAKHEIGHT,nrows,1,1,&peak,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_PKHTERR,nrows,1,1,&peake,&status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1,nrows,1,1,&apflux1,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX1ERR,nrows,1,1,&apflux1e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2,nrows,1,1,&apflux2,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX2ERR,nrows,1,1,&apflux2e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3,nrows,1,1,&apflux3,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX3ERR,nrows,1,1,&apflux3e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4,nrows,1,1,&apflux4,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX4ERR,nrows,1,1,&apflux4e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5,nrows,1,1,&apflux5,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX5ERR,nrows,1,1,&apflux5e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6,nrows,1,1,&apflux6,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX6ERR,nrows,1,1,&apflux6e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX7,nrows,1,1,&apflux7,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX7ERR,nrows,1,1,&apflux7e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX8,nrows,1,1,&apflux8,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX8ERR,nrows,1,1,&apflux8e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX9,nrows,1,1,&apflux9,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX9ERR,nrows,1,1,&apflux9e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX10,nrows,1,1,&apflux10,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX10ERR,nrows,1,1,&apflux10e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX11,nrows,1,1,&apflux11,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX11ERR,nrows,1,1,&apflux11e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX12,nrows,1,1,&apflux12,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX12ERR,nrows,1,1,&apflux12e,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX13,nrows,1,1,&apflux13,
                             &status);
        (void)fits_write_col(tptr,TFLOAT,COL_APFLUX13ERR,nrows,1,1,&apflux13e,
                             &status);
        r = 0.5*petr_rad[k];
	(void)fits_write_col(tptr,TFLOAT,COL_PETRAD,nrows,1,1,&r,&status);
	r = 0.5*kron_rad[k];
	(void)fits_write_col(tptr,TFLOAT,COL_KRONRAD,nrows,1,1,&r,&status);
	r = half_rad[k];
	(void)fits_write_col(tptr,TFLOAT,COL_HALFRAD,nrows,1,1,&r,&status);
	(void)fits_write_col(tptr,TFLOAT,COL_PETFLUX,nrows,1,1,(petr_flux+k),
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_PETFLUXERR,nrows,1,1,&petr_fluxe,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_KRONFLUX,nrows,1,1,(kron_flux+k),
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_KRONFLUXERR,nrows,1,1,&kron_fluxe,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_HALFFLUX,nrows,1,1,(half_flux+k),
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_HALFFLUXERR,nrows,1,1,&half_fluxe,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_ERRFLAG,nrows,1,1,(badpix+k),
			     &status);
        (void)fits_write_col(tptr,TFLOAT,COL_SKYLEVEL,nrows,1,1,&skylev,
			     &status);
        (void)fits_write_col(tptr,TFLOAT,COL_SKYSIGMA,nrows,1,1,&skyrms,
			     &status);
	(void)fits_write_col(tptr,TFLOAT,COL_AVCONF,nrows,1,1,(avconf+k),
			     &status);
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

extern int tabclose_6(ap_t *ap, char *errstr) {

    return(tabclose_gen(ap,errstr));
}

extern int readtab_list_6(char *infile, objstruct **ob, int *nr, 
			  char *errstr) {
    int status,anynul;
    long nrows,i;
    float *work;
    fitsfile *fptr;
    char msg[BUFSIZ];
    double degrad;

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

    (void)fits_read_col(fptr,TFLOAT,COL_APFLUX3,1,1,nrows,NULL,work,&anynul,&status);
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

extern int get_ncol_6(void) {
    return(NCOLS);
}

static int ismef(char *infile) {
    fitsfile *fptr;
    int status,nhdu;

    /* Open the file and see how many HDUs there are. If it's 1, then
       this is a SEF otherwise it's a MEF */

    status = 0;
    (void)fits_open_file(&fptr,infile,READONLY,&status);
    (void)fits_get_num_hdus(fptr,&nhdu,&status);
    (void)fits_close_file(fptr,&status);
    return(nhdu > 1 ? 1 : 0);
}

/* 

$Log: create_table_6.c,v $
Revision 1.7  2014/07/30 08:12:45  jim
Fixed long standing bug that means that in list mode the output equatorial
coordinates are written to the FITS tables in degrees rather than in radians
as advertised

Revision 1.6  2014/07/18 08:46:11  jim
Fixed ell file so that in list mode it puts a circle where the original
coordinates were specified, rather than where it thinks it found the object

Revision 1.5  2011/02/03 10:43:54  jim
Added old style FITS WCS keywords to lookup tables

Revision 1.4  2011-02-03 10:26:38  jim
Fixed a little bug that might affect SEFs if some keywords are missing

Revision 1.3  2010-11-01 10:59:04  jim
Fixed problem in process_results_list routine where we could possibly get NULL
values for the position errors

Revision 1.2  2010-08-06 08:30:13  jim
nothing

Revision 1.1.1.1  2010-07-27 08:41:20  jim
Imported casutools

Revision 1.5  2010/02/10 11:54:03  jim
Modified calculation of half-light flux and radius

Revision 1.4  2009/12/17 11:35:19  jim
replaced hall radius/flux with half light radius/flux. Also writes out
average confidence now into column 58

Revision 1.3  2009/11/16 18:50:46  jim
Modified so that SEF input files are detected and the WCS is transfered over
to the extension header in the table

Revision 1.2  2009/10/29 15:22:05  jim
Now translates the image WCS keywords to official FITS table WCS keywords

Revision 1.1  2009/09/25 08:51:30  jim
new entry


*/
