#pragma once
#include "zbsubr.h"
#include <cmath>

namespace zbessel {

template <class>
void zshch(double zr, double zi, double *__restrict__ cshr, double *__restrict__ cshi, double *__restrict__ cchr,
           double *__restrict__ cchi) {
  /* Local variables */
  double ch, cn, sh, sn;

  /* ***BEGIN PROLOGUE  ZSHCH */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESH and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CSHCH-A, ZSHCH-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y) */
  /*     AND CCH=COSH(X+I*Y), WHERE I**2=-1. */

  /* ***SEE ALSO  ZBESH, ZBESK */
  /* ***ROUTINES CALLED  (NONE) */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZSHCH */

  /* ***FIRST EXECUTABLE STATEMENT  ZSHCH */
  sh = std::sinh(zr);
  ch = std::cosh(zr);
  sn = std::sin(zi);
  cn = std::cos(zi);
  *cshr = sh * cn;
  *cshi = ch * sn;
  *cchr = ch * cn;
  *cchi = sh * sn;
}

}  // namespace zbessel
