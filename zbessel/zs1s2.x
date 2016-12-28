#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
void zs1s2(double zrr, double zri, double *__restrict__ s1r, double *__restrict__ s1i, double *__restrict__ s2r,
           double *__restrict__ s2i, int *__restrict__ nz, double ascle, double alim, int *__restrict__ iuf) {
  /* Local variables */
  double aa, c1i, as1, as2, c1r, aln, s1di, s1dr;

  /* ***BEGIN PROLOGUE  ZS1S2 */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZAIRY and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CS1S2-A, ZS1S2-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE */
  /*     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON- */
  /*     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION. */
  /*     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF */
  /*     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER */
  /*     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE */
  /*     PRECISION ABOVE THE UNDERFLOW LIMIT. */

  /* ***SEE ALSO  ZAIRY, ZBESK */
  /* ***ROUTINES CALLED  ZABS, ZEXP, ZLOG */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZS1S2 */
  /*     COMPLEX CZERO,C1,S1,S1D,S2,ZR */
  /* ***FIRST EXECUTABLE STATEMENT  ZS1S2 */
  *nz = 0;
  as1 = zabs(*s1r, *s1i);
  as2 = zabs(*s2r, *s2i);
  if (*s1r == 0. && *s1i == 0.) {
    goto L10;
  }
  if (as1 == 0.) {
    goto L10;
  }
  aln = -zrr - zrr + std::log(as1);
  s1dr = *s1r;
  s1di = *s1i;
  *s1r = 0.;
  *s1i = 0.;
  as1 = 0.;
  if (aln < -alim) {
    goto L10;
  }
  zlog(s1dr, s1di, &c1r, &c1i);
  c1r = c1r - zrr - zrr;
  c1i = c1i - zri - zri;
  zexp(c1r, c1i, s1r, s1i);
  as1 = zabs(*s1r, *s1i);
  ++(*iuf);
L10:
  aa = std::max(as1, as2);
  if (aa > ascle) {
    return;
  }
  *s1r = 0.;
  *s1i = 0.;
  *s2r = 0.;
  *s2i = 0.;
  *nz = 1;
  *iuf = 0;
}

}  // namespace zbessel
