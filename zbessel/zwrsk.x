#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <limits>
#include <cmath>

namespace zbessel {

template <class>
void zwrsk(double zrr, double zri, double fnu, int kode, int n, double *__restrict__ yr,
           double *__restrict__ yi, int *__restrict__ nz, double *__restrict__ cwr, double *__restrict__ cwi, double tol,
           double elim, double alim) {
  /* Local variables */
  int i__, nw;
  double c1i, c2i, c1r, c2r, act, acw, cti, ctr, pti, sti, ptr, str, ract;
  double ascle, csclr, cinui, cinur;

  /* ***BEGIN PROLOGUE  ZWRSK */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CWRSK-A, ZWRSK-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY */
  /*     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZBKNU, ZRATI */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZWRSK */
  /*     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR */
  /* ***FIRST EXECUTABLE STATEMENT  ZWRSK */
  /* ----------------------------------------------------------------------- */
  /*     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS */
  /*     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE */
  /*     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU. */
  /* ----------------------------------------------------------------------- */

  /* Parameter adjustments */
  --yi;
  --yr;
  --cwr;
  --cwi;

  /* Function Body */
  *nz = 0;
  zbknu(zrr, zri, fnu, kode, 2, &cwr[1], &cwi[1], &nw, tol, elim, alim);
  if (nw != 0) {
    goto L50;
  }
  zrati(zrr, zri, fnu, n, &yr[1], &yi[1], tol);
  /* ----------------------------------------------------------------------- */
  /*     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z), */
  /*     R(FNU+J-1,Z)=Y(J),  J=1,...,N */
  /* ----------------------------------------------------------------------- */
  cinur = 1.;
  cinui = 0.;
  if (kode == 1) {
    goto L10;
  }
  cinur = std::cos(zri);
  cinui = std::sin(zri);
L10:
  /* ----------------------------------------------------------------------- */
  /*     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH */
  /*     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE */
  /*     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT */
  /*     THE RESULT IS ON SCALE. */
  /* ----------------------------------------------------------------------- */
  acw = zabs(cwr[2], cwi[2]);
  ascle = std::numeric_limits<double>::min() * 1e3 / tol;
  csclr = 1.;
  if (acw > ascle) {
    goto L20;
  }
  csclr = 1. / tol;
  goto L30;
L20:
  ascle = 1. / ascle;
  if (acw < ascle) {
    goto L30;
  }
  csclr = tol;
L30:
  c1r = cwr[1] * csclr;
  c1i = cwi[1] * csclr;
  c2r = cwr[2] * csclr;
  c2i = cwi[2] * csclr;
  str = yr[1];
  sti = yi[1];
  /* ----------------------------------------------------------------------- */
  /*     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0D0/ABS(CT) PREVENTS */
  /*     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT) */
  /* ----------------------------------------------------------------------- */
  ptr = str * c1r - sti * c1i;
  pti = str * c1i + sti * c1r;
  ptr += c2r;
  pti += c2i;
  ctr = zrr * ptr - zri * pti;
  cti = zrr * pti + zri * ptr;
  act = zabs(ctr, cti);
  ract = 1. / act;
  ctr *= ract;
  cti = -cti * ract;
  ptr = cinur * ract;
  pti = cinui * ract;
  cinur = ptr * ctr - pti * cti;
  cinui = ptr * cti + pti * ctr;
  yr[1] = cinur * csclr;
  yi[1] = cinui * csclr;
  if (n == 1) {
    return;
  }
  for (i__ = 2; i__ <= n; ++i__) {
    ptr = str * cinur - sti * cinui;
    cinui = str * cinui + sti * cinur;
    cinur = ptr;
    str = yr[i__];
    sti = yi[i__];
    yr[i__] = cinur * csclr;
    yi[i__] = cinui * csclr;
    /* L40: */
  }
  return;
L50:
  *nz = -1;
  if (nw == -2) {
    *nz = -2;
  }
}

}  // namespace zbessel
