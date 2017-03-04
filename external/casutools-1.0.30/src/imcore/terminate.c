/*

$Id: terminate.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/
#include <stdio.h>
#include "imcore.h"
#include "errcodes.h"
#include "util.h"

void restack(ap_t *ap, int ip) {
    int i,ib,nn,np;
    unsigned char *mflag;

    /* Reset the mflag */

    np = ap->parent[ip].pnop;
    ib = ap->parent[ip].first;
    mflag = ap->mflag;
    for (i = 0; i < np; i++) {
	nn = ap->plessey[ib].y*ap->lsiz + ap->plessey[ib].x;
  	mflag[nn] = MF_POSSIBLEOBJ;
        ib = ap->blink[ib];
    }

    /* Stash all blocks back in a burst: */

    ib = ap->parent[ip].first;
    for(i = ap->ibstack - ap->parent[ip].pnop; i < ap->ibstack-1;  i++) {    
        ap->bstack[i] = ib;
        ib = ap->blink[ib];
    }

    /* and the last one: */

    ap->bstack[ap->ibstack-1] = ib;
    ap->ibstack -= ap->parent[ip].pnop;

    /* Put parent name back on stack: */

    ap->pstack[--ap->ipstack] = ip;

    /* Mark that parent inactive: */

    ap->parent[ip].pnop = -1;
    ap->parent[ip].pnbp = -1;
}

void terminate(ap_t *ap, unsigned char *mpt) {
    int ip,status;
    char errstr[BUFSIZ];

    /* Search through all possible parents!  */

    for (ip = 1; ip <= ap->maxip; ip++) {
        if (ap->parent[ip].pnop != -1) {
            if (ap->parent[ip].pnop == ap->parent[ip].growing) {

                /* That's a termination: */

                if ((ap->parent[ip].pnop >= ap->ipnop &&
                    ap->parent[ip].touch == 0) &&
                    (ap->parent[ip].pnbp < (ap->parent[ip].pnop)/2)) {
		    extract_data(ap,ip);
                    
                    /* Call the processing routine */

		    status = process_results(ap,errstr);
		    if (status != ERRCODE_OK) {
			restack(ap,ip);
			continue;
		    }
                }
                restack(ap,ip);
            } else {

                /* This parent still active: */

                ap->parent[ip].growing = ap->parent[ip].pnop;
            }
        }
    }
}

void apfu(ap_t *ap) {
    int ip, big, ipbig;

    /* Search through all possible parents and just junk the biggest
       one to free space:  */

    big = 0;
    ipbig = 0;
    for (ip = 1; ip <= ap->maxip; ip++) {
        if(ap->parent[ip].pnop != -1) {
            if(ap->parent[ip].pnop > big) {
                big = ap->parent[ip].pnop;
                ipbig = ip;
            }
        }
    }
    if(big > 0) {
        restack(ap, ipbig);

        /* clearout lastline references to this parent: */

        for (ip = 0; ip <= ap->lsiz; ip++)
            if(ap->lastline[ip] == ipbig) ap->lastline[ip] = 0;
    }
}

void extract_data(ap_t *ap, int ip) {
    int ib,i,np,nn;
    unsigned char *mflag;

    /* Check the size of the workspace and see if it's big enough. If it
       isn't then increase the size until it is */

    np = ap->parent[ip].pnop;
    if (ap->npl < np) {
        ap->plarray = realloc(ap->plarray,np*sizeof(plstruct));
	ap->npl = np;
    }

    /* Pull the info out now */

    ib = ap->parent[ip].first;
    ap->npl_pix = np;
    mflag = ap->mflag;
    for (i = 0; i < np; i++) {
        ap->plarray[i].x = ap->plessey[ib].x + 1;
        ap->plarray[i].y = ap->plessey[ib].y + 1;
	ap->plarray[i].z = ap->plessey[ib].z;
	ap->plarray[i].zsm = ap->plessey[ib].zsm;
	nn = ap->plessey[ib].y*ap->lsiz + ap->plessey[ib].x;
  	mflag[nn] = MF_OBJPIX;
        ib = ap->blink[ib];
    }
}


/*

$Log: terminate.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.6  2009/01/22 14:18:07  jim
Resets the mflag for freed pixels

Revision 1.5  2008/04/15 19:10:15  jim
superficial changes

Revision 1.4  2005/03/06 19:41:21  jim
Removed unnecessary diagnostic things

Revision 1.3  2004/09/07 14:18:58  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.2  2004/04/05 11:25:42  jim
Small modifications and bug fixes

Revision 1.1  2004/04/02 10:55:01  jim
New version for rewrite of imcore


*/
