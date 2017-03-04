/*

$Id: apinit.c,v 1.2 2014/07/31 12:45:16 jim Exp $

*/
#include <stdlib.h>
#include "imcore.h"

void apinit(ap_t *ap) {
    int maxpa;

    int i;

    maxpa = ap->lsiz / 2;		/* max possible parents */
    ap->lastline = calloc(ap->lsiz + 1, sizeof(short int));
    ap->maxip = 0;
    ap->maxpa = maxpa;
    ap->pstack = malloc(maxpa*sizeof(*ap->pstack));
    ap->parent = malloc(maxpa*sizeof(*(ap->parent)));
    for(i = 0; i < maxpa; i++) {
	ap->pstack[i] = i;
	ap->parent[i].pnop = -1;	/* mark all parents inactive */
	ap->parent[i].pnbp = -1;	/* mark all parents inactive */
    }
    ap->ipstack = 1;
    ap->maxbl = MAXBL;
    ap->bstack = malloc(ap->maxbl*sizeof(*ap->bstack));
    ap->blink = malloc(ap->maxbl*sizeof(*ap->blink));
    ap->plessey = malloc(ap->maxbl*sizeof(*ap->plessey));
    for (i = 0; i < MAXBL; i++)
        ap->bstack[i] = i;
    ap->ibstack = 2;	/* block 1 will get overwritten; don't use it */
    ap->nimages = 0;

    /* set up exponential areal-profile levels: */

    ap->areal[0] = 1;
    for (i = 1; i < 8; i++)
        ap->areal[i] = ap->areal[i-1]*2;

    /* allocate some space for a processing array */

    ap->npl = ap->lsiz;
    ap->npl_pix = 0;
    ap->plarray = malloc(ap->npl*sizeof(plstruct));

    /* set these to null values as you may not need the background structure */

    ap->backmap.nby = -1;
    ap->backmap.bvals = NULL;
}

void apreinit(ap_t *ap) {
    int i;

    for (i = 0; i < ap->lsiz+1; i++)
	ap->lastline[i] = 0;
    ap->maxip = 0;
    for(i = 0; i < ap->maxpa; i++) {
	ap->pstack[i] = i;
	ap->parent[i].pnop = -1;	/* mark all parents inactive */
	ap->parent[i].pnbp = -1;	/* mark all parents inactive */
    }
    ap->ipstack = 1;
    ap->ibstack = 2;	/* block 1 will get overwritten; don't use it */
    ap->nimages = 0;
    ap->npl_pix = 0;

}

void apclose(ap_t *ap) {
    int i;
    free(ap->lastline);
    free(ap->pstack);
    free(ap->parent);
    free(ap->bstack);
    free(ap->blink);
    free(ap->plessey);
    free(ap->plarray);
    if (ap->backmap.bvals != NULL) {
        for (i = 0; i < ap->backmap.nby; i++)
            free(ap->backmap.bvals[i]);
        free(ap->backmap.bvals);
    }
}

/*

$Log: apinit.c,v $
Revision 1.2  2014/07/31 12:45:16  jim
Modified so that nbsize <= 0 gives a constant background

Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.2  2004/04/05 11:25:42  jim
Small modifications and bug fixes

Revision 1.1  2004/04/02 10:54:56  jim
New version for rewrite of imcore



*/
