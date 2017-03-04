/* 

$Id: seeing.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/

#include <stdio.h>
#include <math.h>
#include "imcore.h"
#include "util.h"
#include "floatmath.h"

static void sortit (float [], int);

extern void seeing(ap_t *ap, int nrows, float *ellipt, float *pkht, 
		   float **areal, float *work, float *fwhm) {
    int i,ii,iaper;
    float aper,delaper,area,logf5t,logf2,arg;

    /* Convenience variables */

    logf5t = logf(0.5/ap->thresh);
    logf2 = logf(2.0);

    /* Do the seeing calculation */

    ii = 0;
    for (i = 0; i < nrows; i++) {
        if (ellipt[i] < 0.2 && pkht[i] < 30000.0 && pkht[i] > 10.0*ap->thresh) {
	    aper = (logf5t + logf(pkht[i]))/logf2 + 1.0;
            iaper = (int)aper;
            delaper = aper - iaper;
	    if (iaper > 0 && iaper < NAREAL && areal[1][i] > 0.0) {
		area = (1.0-delaper)*areal[iaper-1][i] + delaper*areal[iaper][i];
		work[ii++] = M_2_SQRTPI*sqrtf(area);
	    }
	}
    }

    /* Sort the resulting array and choose a location that allows for
       contamination by galaxies */

    if (ii >= 3) {    
        sortit(work,ii);
        *fwhm = work[ii/3 - 1];

	/* Allow for finite pixel size */

        arg = 0.25*M_PI*powf(*fwhm,2.0) - 1;
	*fwhm = 2.0*sqrt(MAX(0.0,arg/M_PI));
    } else 
	*fwhm = 0.0;

    /* Message if in verbose mode */

    if (verbose) 
        printf("Estimated FWHM (pixels) = %8.1f\nNo. of objects used     = %8i\n",
	       *fwhm,ii);
}

static void sortit (float ia[], int n) {
    int i, j, ii, jj, ifin;
    float it;
 
    jj = 4;
    while (jj < n) 
	jj = 2 * jj;
    jj = MIN(n,(3 * jj)/4 - 1);
    while (jj > 1) {
        jj = jj/2;
        ifin = n - jj;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + jj;
            if (ia[i] <= ia[j]) 
		continue;
            it = ia[j];
            do {
                ia[j] = ia[i];
                j = i;
                i = i - jj;
                if (i < 0) 
		    break;
            } while (ia[i] > it);
            ia[j] = it;
        }
    }
    return;
}

/*

$Log: seeing.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.3  2005/04/20 07:48:27  jim
Modified estimate of seeing allowing for undersampled data

Revision 1.2  2004/09/07 14:18:58  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.1  2004/04/02 10:55:00  jim
New version for rewrite of imcore


*/
