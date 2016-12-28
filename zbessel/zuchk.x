#pragma once
#include "zbsubr.h"
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
void zuchk(double yr, double yi, int *__restrict__ nz, double ascle, double tol) {
  double wi, ss, st, wr;

  /* ***BEGIN PROLOGUE  ZUCHK */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to SERI, ZUOIK, ZUNK1, ZUNK2, ZUNI1, ZUNI2 and */
  /*            ZKSCL */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUCHK-A, ZUCHK-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN */
  /*      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE */
  /*      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW */
  /*      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED */
  /*      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE */
  /*      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE */
  /*      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED. */

  /* ***SEE ALSO  SERI, ZKSCL, ZUNI1, ZUNI2, ZUNK1, ZUNK2, ZUOIK */
  /* ***ROUTINES CALLED  (NONE) */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   ??????  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZUCHK */

  /*     COMPLEX Y */
  /* ***FIRST EXECUTABLE STATEMENT  ZUCHK */
  *nz = 0;
  wr = std::fabs(yr);
  wi = std::fabs(yi);
  st = std::min(wr, wi);
  if (st > ascle)
    return;
  ss = std::max(wr, wi);
  st /= tol;
  if (ss < st)
    *nz = 1;
}

}  // namespace zbessel
