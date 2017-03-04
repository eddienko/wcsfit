/* 

$Id: vsttab_a2f_main.c,v 1.1 2013/05/30 14:28:26 jim Exp $

*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>

#include <vstfix.h>
#include <tools.h>
static char *vstfix_tables[VSTFIX_NTABS] = {"Det_offs","Bilinear_cor",
					    "Radial_cor"};
static int vstfix_ncols[VSTFIX_NTABS] = {2,5,3};
static char *vstfix_ttype1[] = {"zpsmed","mzps"};
static char *vstfix_tform1[] = {"1E","1J"};
static char *vstfix_tunit1[] = {"",""};
static char *vstfix_ttype2[] = {"xn","bin1dx","bin1dy","nbin1dx","nbin1dy"};
static char *vstfix_tform2[] = {"1E","1E","1E","1J","1J"};
static char *vstfix_tunit2[] = {"","","","",""};
static char *vstfix_ttype3[] = {"coord","radial","nradial"};
static char *vstfix_tform3[] = {"1E","1E","1J"};
static char *vstfix_tunit3[] = {"","",""};

static char *USAGE = "Usage: vsttab_a2f inascii outfits\n";

static FILE *fd = NULL;
static fitsfile *optr = NULL;
static float *work1 = NULL;
static float *work2 = NULL;
static float *work3 = NULL;
static int *iwork1 = NULL;
static int *iwork2 = NULL;

static void tidy(void);

int main (int argc, char *argv[]) {
    int status,i,ival;
    char *inascii,*outfits,msg[BUFSIZ];

    /* Check that we have enough information */

    if (argc < 3) {
        fprintf(stderr,USAGE);
	exit(1);
    } else {
	inascii = argv[1];
	outfits = argv[2];
    }

    /* Open up the ASCII file */

    fd = fopen(inascii,"r");
    if (fd == NULL) {
	fprintf(stderr,"VSTTAB_A2F: Unable to open ASCII file %s\n",inascii);
	tidy();
	exit(1);
    }

    /* If the output file already exists, the clobber it */

    if (access((const char *)outfits,F_OK) == 0) 
	remove((const char *)outfits);

    /* Now create the output FITS file */

    status = 0;
    (void)fits_create_file(&optr,outfits,&status);
    if (status != 0) {
	fprintf(stderr,"VSTTAB_A2F: Unable to create output FITS file %s\n",
		outfits);
	tidy();
	exit(1);
    }

    /* Read the first lot of information from the file */

    work1 = cir_malloc(DEFNEXTN*sizeof(float));
    iwork1 = cir_malloc(DEFNEXTN*sizeof(int));
    (void)fscanf(fd,"%[^\n] ",msg);
    i = 0;
    while (i < DEFNEXTN && fscanf(fd,"%g %d %d",work1+i,iwork1+i,&ival) != EOF) 
	i++;
    if (i < DEFNEXTN) {
	fprintf(stderr,"VSTTAB_A2F: Only %d entries for detector zeropoints\n",
		i);
	tidy();
	exit(1);
    }
	
    /* Create the first table */

    (void)fits_create_tbl(optr,BINARY_TBL,0,vstfix_ncols[0],vstfix_ttype1,
			  vstfix_tform1,vstfix_tunit1,vstfix_tables[0],
			  &status);
    
    /* Now write the column data */

    (void)fits_write_col(optr,TFLOAT,1,1,1,DEFNEXTN,work1,&status);
    (void)fits_write_col(optr,TINT,2,1,1,DEFNEXTN,iwork1,&status);
    freespace(work1);
    freespace(iwork1);

    /* Read the second lot of information from the file */

    work1 = cir_malloc(VSTFIX_MPBIN*sizeof(float));
    work2 = cir_malloc(VSTFIX_MPBIN*sizeof(float));
    work3 = cir_malloc(VSTFIX_MPBIN*sizeof(float));
    iwork1 = cir_malloc(VSTFIX_MPBIN*sizeof(int));
    iwork2 = cir_malloc(VSTFIX_MPBIN*sizeof(int));
    (void)fscanf(fd,"%s",msg);
    i = 0;
    while (i < VSTFIX_MPBIN && fscanf(fd,"%g %g %g %d %d",
				      work1+i,work2+i,work3+i,
				      iwork1+i,iwork2+i) != EOF) 
	i++;
    if (i < VSTFIX_MPBIN) {
	fprintf(stderr,"VSTTAB_A2F: Only %d entries for 2D interpolation\n",i);
	tidy();
	exit(1);
    }
		
    /* Create the second table */

    (void)fits_create_tbl(optr,BINARY_TBL,0,vstfix_ncols[1],vstfix_ttype2,
			  vstfix_tform2,vstfix_tunit2,vstfix_tables[1],
			  &status);
    
    /* Now write the column data */

    (void)fits_write_col(optr,TFLOAT,1,1,1,VSTFIX_MPBIN,work1,&status);
    (void)fits_write_col(optr,TFLOAT,2,1,1,VSTFIX_MPBIN,work2,&status);
    (void)fits_write_col(optr,TFLOAT,3,1,1,VSTFIX_MPBIN,work3,&status);
    (void)fits_write_col(optr,TINT,4,1,1,VSTFIX_MPBIN,iwork1,&status);
    (void)fits_write_col(optr,TINT,5,1,1,VSTFIX_MPBIN,iwork2,&status);
    freespace(work1);
    freespace(work2);
    freespace(work3);
    freespace(iwork1);
    freespace(iwork2);
    
    /* Read the third lot of information from the file */

    work1 = cir_malloc(VSTFIX_MRBIN*sizeof(float));
    work2 = cir_malloc(VSTFIX_MRBIN*sizeof(float));
    iwork1 = cir_malloc(VSTFIX_MRBIN*sizeof(int));
    (void)fscanf(fd,"%s",msg);
    i = 0;
    while (i < VSTFIX_MRBIN && fscanf(fd,"%g %g %d",work1+i,work2+i,
				      iwork1+i) != EOF) 
	i++;
    if (i < VSTFIX_MRBIN) {
	fprintf(stderr,"VSTTAB_A2F: Only %d entries for radial correction\n",i);
	tidy();
	exit(1);
    }
		
    /* Create the third table */

    (void)fits_create_tbl(optr,BINARY_TBL,0,vstfix_ncols[2],vstfix_ttype3,
			  vstfix_tform3,vstfix_tunit3,vstfix_tables[2],
			  &status);
    
    /* Now write the column data */

    (void)fits_write_col(optr,TFLOAT,1,1,1,VSTFIX_MRBIN,work1,&status);
    (void)fits_write_col(optr,TFLOAT,2,1,1,VSTFIX_MRBIN,work2,&status);
    (void)fits_write_col(optr,TINT,3,1,1,VSTFIX_MRBIN,iwork1,&status);
    freespace(work1);
    freespace(work2);
    freespace(iwork1);

    /* Close things up and get out of here */

    tidy();
    return(0);
}

static void tidy(void) {
    closefits(optr);
    closefile(fd);
    freespace(work1);
    freespace(work2);
    freespace(work3);
    freespace(iwork1);
    freespace(iwork2);
}

/*

$Log: vsttab_a2f_main.c,v $
Revision 1.1  2013/05/30 14:28:26  jim
New Entry


*/
