#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "imcore.h"
#include "util.h"
#include "floatmath.h"

extern void areals(ap_t *ap, int iareal[NAREAL]) {
    int i,nup,j,np;
    float t,thresh,fconst,offset;
    plstruct *plarray;

    /* Initialise some stuff */

    np = ap->npl_pix;
    plarray = ap->plarray;
    thresh = ap->thresh;
    fconst = ap->fconst;
    offset = ap->areal_offset;

    /* Zero the areal profile array */

    (void)memset(iareal,0,NAREAL*sizeof(int));

    /* Loop through the array and form the areal profiles */

    for (i = 0; i < np; i++) {
        t = plarray[i].z;
        if (t <= thresh) 
	    continue;
	nup = MIN(NAREAL,(int)(logf(t)*fconst - offset)+1);
	nup = MAX(1,nup);
	for (j = 0; j < nup; j++)
	    iareal[j]++;
    }
}
		
