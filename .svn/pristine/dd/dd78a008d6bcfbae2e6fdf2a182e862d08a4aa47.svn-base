/*

$Id: polynm.c,v 1.1.1.1 2010/07/27 08:41:20 jim Exp $

*/
#include <stdlib.h>
#include "imcore.h"
#include "floatmath.h"
#include "util.h"

/* least-squares fit of order m polynomial to n data points */

void polynm (float xdat[], float xcor[], int n, float polycf[], int m, int ilim) {
  double a[25][25], b[25], temp;
  int i, j, k;

/*   if(n < m)  bomboutx(1, "polynm: too few data points"); */
/*   if(m > 25) bomboutx(1, "polynm: order of polynomial too large"); */

  /* clear arrays */
  for(i = 0; i < 25; i++) {
    b[i] = 0.0;
    for(j = 0; j < 25; j++) a[i][j] = 0.0;
  }

  /* cumulate sums */
  for(i = 0; i < n; i++) {
    for(k = 0; k < m; k++) {
      temp = 1.0;
      if(k+ilim != 0)temp = pow(xcor[i], (float) (k+ilim));
      b[k] += xdat[i]*temp;

      for(j = 0; j <= k; j++) {
	temp = 1.0;
	if(k+j+2*ilim != 0)temp = pow(xcor[i], (float) (k+j+2*ilim));
	a[j][k] += temp;
      }
    }
  }

  for(k = 1; k < m; k++) {
    for(j = 0; j < k; j++) a[k][j] = a[j][k];
  }

  /* solve linear equations */
  solve(a, b, m);

  for(i = 0; i < m; i++) polycf[i] = b[i];
}

/*

$Log: polynm.c,v $
Revision 1.1.1.1  2010/07/27 08:41:20  jim
Imported casutools

Revision 1.1  2004/04/02 10:55:00  jim
New version for rewrite of imcore


*/
