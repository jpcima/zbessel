#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
void zmlri(double zr, double zi, double fnu, int kode, int n, double *__restrict__ yr,
           double *__restrict__ yi, int *__restrict__ nz, double tol) {
  /* Local variables */
  int i__, k, m;
  double ak, bk, ap, at;
  int kk, km;
  double az, p1i, p2i, p1r, p2r, ack, cki, fnf, fkk, ckr;
  int iaz;
  double rho;
  int inu;
  double pti, raz, sti, rzi, ptr, str, tst, rzr, rho2, flam, fkap, scle, tfnf;
  int ifnu;
  double sumi, sumr;
  int itime;
  double cnormi, cnormr;

  /* ***BEGIN PROLOGUE  ZMLRI */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CMLRI-A, ZMLRI-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE */
  /*     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES. */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZMLRI */
  /*     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZMLRI */
  scle = std::numeric_limits<double>::min() / tol;
  *nz = 0;
  az = zabs(zr, zi);
  iaz = (int)az;
  ifnu = (int)fnu;
  inu = ifnu + n - 1;
  at = iaz + 1.;
  raz = 1. / az;
  str = zr * raz;
  sti = -zi * raz;
  ckr = str * at * raz;
  cki = sti * at * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  p1r = 0.;
  p1i = 0.;
  p2r = 1.;
  p2i = 0.;
  ack = (at + 1.) * raz;
  rho = ack + std::sqrt(ack * ack - 1.);
  rho2 = rho * rho;
  tst = (rho2 + rho2) / ((rho2 - 1.) * (rho - 1.));
  tst /= tol;
  /* ----------------------------------------------------------------------- */
  /*     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES */
  /* ----------------------------------------------------------------------- */
  ak = at;
  for (i__ = 1; i__ <= 80; ++i__) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r - (ckr * ptr - cki * pti);
    p2i = p1i - (cki * ptr + ckr * pti);
    p1r = ptr;
    p1i = pti;
    ckr += rzr;
    cki += rzi;
    ap = zabs(p2r, p2i);
    if (ap > tst * ak * ak) {
      goto L20;
    }
    ak += 1.;
    /* L10: */
  }
  goto L110;
L20:
  ++i__;
  k = 0;
  if (inu < iaz) {
    goto L40;
  }
  /* ----------------------------------------------------------------------- */
  /*     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS */
  /* ----------------------------------------------------------------------- */
  p1r = 0.;
  p1i = 0.;
  p2r = 1.;
  p2i = 0.;
  at = inu + 1.;
  str = zr * raz;
  sti = -zi * raz;
  ckr = str * at * raz;
  cki = sti * at * raz;
  ack = at * raz;
  tst = std::sqrt(ack / tol);
  itime = 1;
  for (k = 1; k <= 80; ++k) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r - (ckr * ptr - cki * pti);
    p2i = p1i - (ckr * pti + cki * ptr);
    p1r = ptr;
    p1i = pti;
    ckr += rzr;
    cki += rzi;
    ap = zabs(p2r, p2i);
    if (ap < tst) {
      goto L30;
    }
    if (itime == 2) {
      goto L40;
    }
    ack = zabs(ckr, cki);
    flam = ack + std::sqrt(ack * ack - 1.);
    fkap = ap / zabs(p1r, p1i);
    rho = std::min(flam, fkap);
    tst *= std::sqrt(rho / (rho * rho - 1.));
    itime = 2;
  L30:;
  }
  goto L110;
L40:
  /* ----------------------------------------------------------------------- */
  /*     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION */
  /* ----------------------------------------------------------------------- */
  ++k;
  /* Computing MAX */
  kk = std::max(i__ + iaz, k + inu);
  fkk = (double)kk;
  p1r = 0.;
  p1i = 0.;
  /* ----------------------------------------------------------------------- */
  /*     SCALE P2 AND SUM BY SCLE */
  /* ----------------------------------------------------------------------- */
  p2r = scle;
  p2i = 0.;
  fnf = fnu - ifnu;
  tfnf = fnf + fnf;
  bk = std::lgamma(fkk + tfnf + 1.) - std::lgamma(fkk + 1.) - std::lgamma(tfnf + 1.);
  bk = std::exp(bk);
  sumr = 0.;
  sumi = 0.;
  km = kk - inu;
  for (i__ = 1; i__ <= km; ++i__) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
    p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
    p1r = ptr;
    p1i = pti;
    ak = 1. - tfnf / (fkk + tfnf);
    ack = bk * ak;
    sumr += (ack + bk) * p1r;
    sumi += (ack + bk) * p1i;
    bk = ack;
    fkk += -1.;
    /* L50: */
  }
  yr[n] = p2r;
  yi[n] = p2i;
  if (n == 1) {
    goto L70;
  }
  for (i__ = 2; i__ <= n; ++i__) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
    p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
    p1r = ptr;
    p1i = pti;
    ak = 1. - tfnf / (fkk + tfnf);
    ack = bk * ak;
    sumr += (ack + bk) * p1r;
    sumi += (ack + bk) * p1i;
    bk = ack;
    fkk += -1.;
    m = n - i__ + 1;
    yr[m] = p2r;
    yi[m] = p2i;
    /* L60: */
  }
L70:
  if (ifnu <= 0) {
    goto L90;
  }
  for (i__ = 1; i__ <= ifnu; ++i__) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
    p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
    p1r = ptr;
    p1i = pti;
    ak = 1. - tfnf / (fkk + tfnf);
    ack = bk * ak;
    sumr += (ack + bk) * p1r;
    sumi += (ack + bk) * p1i;
    bk = ack;
    fkk += -1.;
    /* L80: */
  }
L90:
  ptr = zr;
  pti = zi;
  if (kode == 2) {
    ptr = 0.;
  }
  zlog(rzr, rzi, &str, &sti);
  p1r = -fnf * str + ptr;
  p1i = -fnf * sti + pti;
  ap = std::lgamma(fnf + 1.);
  ptr = p1r - ap;
  pti = p1i;
  /* ----------------------------------------------------------------------- */
  /*     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW */
  /*     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES */
  /* ----------------------------------------------------------------------- */
  p2r += sumr;
  p2i += sumi;
  ap = zabs(p2r, p2i);
  p1r = 1. / ap;
  zexp(ptr, pti, &str, &sti);
  ckr = str * p1r;
  cki = sti * p1r;
  ptr = p2r * p1r;
  pti = -p2i * p1r;
  zmlt(ckr, cki, ptr, pti, &cnormr, &cnormi);
  for (i__ = 1; i__ <= n; ++i__) {
    str = yr[i__] * cnormr - yi[i__] * cnormi;
    yi[i__] = yr[i__] * cnormi + yi[i__] * cnormr;
    yr[i__] = str;
    /* L100: */
  }
  return;
L110:
  *nz = -2;
}

}  // namespace zbessel
