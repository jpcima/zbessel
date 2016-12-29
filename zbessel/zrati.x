#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
void zrati(double zr, double zi, double fnu, int n, double *__restrict__ cyr,
           double *__restrict__ cyi, double tol) {
  /* Initialized data */

  const double rt2 = 1.41421356237309514547L;  // sqrt(2)

  /* Local variables */
  int i__, k;
  double ak;
  int id, kk;
  double az, ap1, ap2, p1i, p2i, t1i, p1r, p2r, t1r, arg, rak, rho;
  int inu;
  double pti, tti, rzi, ptr, ttr, rzr, rap1, flam, dfnu, fdnu;
  int magz;
  int idnu;
  double fnup;
  double test, test1, amagz;
  int itime;
  double cdfnui, cdfnur;

  /* ***BEGIN PROLOGUE  ZRATI */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CRATI-A, ZRATI-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD */
  /*     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD */
  /*     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B, */
  /*     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973, */
  /*     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER, */
  /*     BY D. J. SOOKNE. */

  /* ***SEE ALSO  ZBESH, ZBESI, ZBESK */
  /* ***ROUTINES CALLED  ZABS, ZDIV */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZRATI */
  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZRATI */
  az = zabs(zr, zi);
  inu = (int)fnu;
  idnu = inu + n - 1;
  magz = (int)az;
  amagz = (double)(magz + 1);
  fdnu = (double)idnu;
  fnup = std::max(amagz, fdnu);
  id = idnu - magz - 1;
  itime = 1;
  k = 1;
  ptr = 1. / az;
  rzr = ptr * (zr + zr) * ptr;
  rzi = -ptr * (zi + zi) * ptr;
  t1r = rzr * fnup;
  t1i = rzi * fnup;
  p2r = -t1r;
  p2i = -t1i;
  p1r = 1.;
  p1i = 0.;
  t1r += rzr;
  t1i += rzi;
  if (id > 0) {
    id = 0;
  }
  ap2 = zabs(p2r, p2i);
  ap1 = zabs(p1r, p1i);
  /* ----------------------------------------------------------------------- */
  /*     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU */
  /*     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT */
  /*     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR */
  /*     PREMATURELY. */
  /* ----------------------------------------------------------------------- */
  arg = (ap2 + ap2) / (ap1 * tol);
  test1 = std::sqrt(arg);
  test = test1;
  rap1 = 1. / ap1;
  p1r *= rap1;
  p1i *= rap1;
  p2r *= rap1;
  p2i *= rap1;
  ap2 *= rap1;
L10:
  ++k;
  ap1 = ap2;
  ptr = p2r;
  pti = p2i;
  p2r = p1r - (t1r * ptr - t1i * pti);
  p2i = p1i - (t1r * pti + t1i * ptr);
  p1r = ptr;
  p1i = pti;
  t1r += rzr;
  t1i += rzi;
  ap2 = zabs(p2r, p2i);
  if (ap1 <= test) {
    goto L10;
  }
  if (itime == 2) {
    goto L20;
  }
  ak = zabs(t1r, t1i) * .5;
  flam = ak + std::sqrt(ak * ak - 1.);
  /* Computing MIN */
  rho = std::min(ap2 / ap1, flam);
  test = test1 * std::sqrt(rho / (rho * rho - 1.));
  itime = 2;
  goto L10;
L20:
  kk = k + 1 - id;
  ak = (double)kk;
  t1r = ak;
  t1i = 0.;
  dfnu = fnu + (n - 1);
  p1r = 1. / ap2;
  p1i = 0.;
  p2r = 0.;
  p2i = 0.;
  for (i__ = 1; i__ <= kk; ++i__) {
    ptr = p1r;
    pti = p1i;
    rap1 = dfnu + t1r;
    ttr = rzr * rap1;
    tti = rzi * rap1;
    p1r = ptr * ttr - pti * tti + p2r;
    p1i = ptr * tti + pti * ttr + p2i;
    p2r = ptr;
    p2i = pti;
    t1r -= 1.;
    /* L30: */
  }
  if (p1r != 0. || p1i != 0.) {
    goto L40;
  }
  p1r = tol;
  p1i = tol;
L40:
  zdiv(p2r, p2i, p1r, p1i, &cyr[n], &cyi[n]);
  if (n == 1) {
    return;
  }
  k = n - 1;
  ak = (double)k;
  t1r = ak;
  t1i = 0.;
  cdfnur = fnu * rzr;
  cdfnui = fnu * rzi;
  for (i__ = 2; i__ <= n; ++i__) {
    ptr = cdfnur + (t1r * rzr - t1i * rzi) + cyr[k + 1];
    pti = cdfnui + (t1r * rzi + t1i * rzr) + cyi[k + 1];
    ak = zabs(ptr, pti);
    if (ak != 0.) {
      goto L50;
    }
    ptr = tol;
    pti = tol;
    ak = tol * rt2;
  L50:
    rak = 1. / ak;
    cyr[k] = rak * ptr * rak;
    cyi[k] = -rak * pti * rak;
    t1r -= 1.;
    --k;
    /* L60: */
  }
}

}  // namespace zbessel
