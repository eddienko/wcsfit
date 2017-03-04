/*

$Id: apclust.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/
#include <stdio.h>
#include "imcore.h"
#include "util.h"

static void minmax_xy(int, plstruct *, int *, int *, int *, int *);

void apclust(ap_t *ap, int np, plstruct *plstr) {

    int i,i1,loop,ik;
    int is;    /* parent name for image in this slice */
    int ip;    /* parent name for image on last line */
    int ib;    /* data block name */
    int k,j,ix1,ix2,iy1,iy2,nwork,nx,ny,kk;
    float i2compare,icompare;
    short int *work;
 
    /* A couple of useful things */

    i2compare = ap->thresh;
    icompare = i2compare * ap->multiply;

    /* Get the min and max positions. Create a raster with the IDs of the
       pixels in the pixel list (algorithm prefers data to be in a raster) */

    minmax_xy(np,plstr,&ix1,&ix2,&iy1,&iy2);
    nx = ix2 - ix1 + 1;
    ny = iy2 - iy1 + 1;
    nwork = nx*ny;
    work = malloc(nwork*sizeof(*work));
    for (i = 0; i < nwork; i++)
        work[i] = -1;
    for (k = 0; k < np; k++) {
        i = plstr[k].x - 1;
        j = plstr[k].y - 1;
	kk = (j - iy1)*nx + i - ix1;
	work[kk] = k;
    }

    /* Now do the job */
	    
    for (j = iy1; j <= iy2; j++) {
        for (i = ix1; i <= ix2; i++) {
            kk = (j - iy1)*nx + i - ix1;
	    k = work[kk];
	    if (k < 0) {
                ap->lastline[i + 1] = 0;
            } else {
		if (plstr[k].zsm > icompare) {

		    /* Pixel is above threshold, find which parent it belongs to. */

		    is = ap->lastline[i];       /* Parent last pixel this line */
		    ip = ap->lastline[i + 1];   /* Guess belongs to above line */
		    if (ip == 0) {

			/* New parent, or, horizontal slice: */

			if (is == 0) {

			    /* Ah - new parent. */

			    if (ap->ipstack > ap->maxpa*3/4) {
 	                        for (ik = 0; ik < ap->maxpa*3/8; ik++)
		                    apfu(ap);
                            }
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
		    ap->plessey[ib].z = plstr[k].z;
		    ap->plessey[ib].zsm = plstr[k].zsm;

		    /* increment active count: */

		    ap->parent[ip].pnop++;

		    /* remember which parent this pixel was for next line: */

		    ap->lastline[i + 1] = ip;

		} else {

		    /* Pixel was below threshold, mark lastline: */

		    ap->lastline[i + 1] = 0;
		}
	    }
	}
    }

    /* Check for images touching left & right edges:
       OR the touch flag with 2 for left, 4 for right: */

    if(ap->lastline[1] > 0 )
        ap->parent[ap->lastline[1]].touch |= 2;
    if(ap->lastline[ap->lsiz] > 0)
        ap->parent[ap->lastline[ap->lsiz]].touch |= 4;
    free(work);
}

static void minmax_xy(int np, plstruct *plstr, int *ix1, int *ix2, int *iy1,
		      int *iy2) {
    int i;

    /* Get the minmax of the positions of the pixels in a plstruct.  Take
       1 away from each position so that it runs from 0 rather than 1 */

    *ix1 = plstr[0].x - 1;
    *ix2 = plstr[0].x - 1;
    *iy1 = plstr[0].y - 1;
    *iy2 = plstr[0].y - 1; 
    for (i = 1; i < np; i++) {
	*ix1 = MIN(*ix1,plstr[i].x - 1);
	*ix2 = MAX(*ix2,plstr[i].x - 1);
	*iy1 = MIN(*iy1,plstr[i].y - 1);
	*iy2 = MAX(*iy2,plstr[i].y - 1);
    }
}


/*

$Log: apclust.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.2  2005/04/15 14:41:29  jim
Added call to apfu if the stack is close to being overrun

Revision 1.1  2004/04/02 10:54:56  jim
New version for rewrite of imcore



*/
