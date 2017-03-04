/*

$Id: apline.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/
#include <stdio.h>
#include "imcore.h"
#include "util.h"


void apline(ap_t *ap, float dat[], float conf[], float smoothed[], 
	    float smoothedc[], int j, unsigned char *bpm) {
    int i,i1,loop,nn;
    int is;    /* parent name for image in this slice */
    int ip;    /* parent name for image on last line */
    int ib;    /* data block name */
    float i2compare,icompare;
    unsigned char *mflag;

/*     i2compare = ap->background + ap->thresh; */
    i2compare = ap->thresh;
    icompare = i2compare * ap->multiply;
    mflag = ap->mflag;

    for (i = 0; i < ap->lsiz; i++) {
/*         if (smoothedc[i] > icompare && conf[i] != 0) { */
        if (smoothedc[i] > icompare) {

            /* Pixel is above threshold, find which parent it belongs to. */

            is = ap->lastline[i];       /* Parent last pixel this line */
            ip = ap->lastline[i + 1];   /* Guess belongs to above line */
            if (ip == 0) {

                /* New parent, or, horizontal slice: */

                if (is == 0) {

                    /* Ah - new parent. */

                    ip = ap->pstack[ap->ipstack++];
                    ap->parent[ip].first = ap->bstack[ap->ibstack];
                    ap->parent[ip].pnop = 0;
                    ap->parent[ip].pnbp = 0;
                    ap->parent[ip].growing = 0;
                    if (j == 0)

                        /* It touches first line: */

                        ap->parent[ip].touch = 1;
                    else
                        ap->parent[ip].touch = 0;

                    /* For hunt thru list for terminates: */

                    if (ip > ap->maxip)
                        ap->maxip = ip;
                } else {

                    /* Slice with no vertical join: */

                    ip = is;
                }
            } else if ((ip > 0 && is > 0) && (ip != is)) {

                /* merge: Join linked lists: */

                ap->blink[ap->parent[ip].last] = ap->parent[is].first;

                /* Copy `last block': */

                ap->parent[ip].last = ap->parent[is].last;
                ap->parent[ip].pnop += ap->parent[is].pnop;
                ap->parent[ip].pnbp += ap->parent[is].pnbp;

                /* Fix `lastline' correlator array: */

                ib = ap->parent[is].first;
                loop = 1;
                while (loop) {
                    i1 = ap->plessey[ib].x;
                    if (ap->lastline[i1 + 1] == is)
                    ap->lastline[i1 + 1] = ip;
                    if (ap->parent[is].last == ib)
                        loop = 0;
                    else
                        ib = ap->blink[ib];
                }

                /* Mark parent inactive: */

                ap->parent[is].pnop = -1;
                ap->parent[is].pnbp = -1;

                /* return name to stack: */

                ap->pstack[--ap->ipstack] = is;
            }

            /* Add in pixel to linked list: */

            ib = ap->bstack[ap->ibstack++];

            /* Patch forward link into last data block: */

            if (ap->parent[ip].pnop > 0)
                ap->blink[ap->parent[ip].last] = ib;

            /* Remember last block in chain: */

            ap->parent[ip].last = ib;

            /* Store the data: */

            ap->plessey[ib].x = i;
            ap->plessey[ib].y = j;
            ap->plessey[ib].z = dat[i];
	    nn = j*ap->lsiz + i;
	    if (mflag[nn] != MF_SATURATED)
                ap->plessey[ib].zsm = MIN(ap->saturation,smoothed[i]);
	    else
		ap->plessey[ib].zsm = ap->saturation;
	    mflag[nn] = MF_POSSIBLEOBJ;

            /* increment active count: */

            ap->parent[ip].pnop++;
            if (bpm != NULL) 
                ap->parent[ip].pnbp += bpm[i];

            /* remember which parent this pixel was for next line: */

            ap->lastline[i + 1] = ip;

        } else {

            /* Pixel was below threshold, mark lastline: */

            ap->lastline[i + 1] = 0;
        }
    }

    /* Check for images touching left & right edges:
       OR the touch flag with 2 for left, 4 for right: */

    if(ap->lastline[1] > 0 )
        ap->parent[ap->lastline[1]].touch |= 2;
    if(ap->lastline[ap->lsiz] > 0)
        ap->parent[ap->lastline[ap->lsiz]].touch |= 4;
}
/*

$Log: apline.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.6  2010/02/11 21:53:35  jim
Changed declaration of confidence to float

Revision 1.5  2009/01/22 13:35:36  jim
Now flags possible object pixels

Revision 1.4  2008/04/15 18:51:27  jim
Took out check for zero confidence. This is dealt with later

Revision 1.3  2006/06/05 11:21:03  jim
Fixed apline so that it takes confidence maps into account better

Revision 1.2  2004/04/05 11:25:42  jim
Small modifications and bug fixes

Revision 1.1  2004/04/02 10:54:56  jim
New version for rewrite of imcore



*/
