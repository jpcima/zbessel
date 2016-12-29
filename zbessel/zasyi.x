#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <limits>
#include <cmath>

namespace zbessel {

template <class>
void zasyi(double zr, double zi, double fnu, int kode,
           int n, double *__restrict__ yr, double *__restrict__ yi, int *__restrict__ nz,
           double rl, double tol, double elim, double alim) {
  /* Initialized data */
  const double pi = 3.14159265358979323846L;
  const double rtpi = 0.159154943091895335769L; // 1/(2*pi)

  /* Local variables */
  int i__, j, k, m;
  double s, aa, bb;
  int ib;
  double ak, bk;
  int il, jl;
  double az;
  int nn;
  double p1i, s2i, p1r, s2r, cki, dki, fdn, arg, aez, arm, ckr, dkr,
      czi, ezi, sgn;
  int inu;
  double raz, czr, ezr, sqk, sti, rzi, tzi, str, rzr, tzr, ak1i,
      ak1r, cs1i, cs2i, cs1r, cs2r, dnu2, rtr1, dfnu, atol;
  int koded;

  /* ***BEGIN PROLOGUE  ZASYI */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CASYI-A, ZASYI-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY */
  /*     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE */
  /*     REGION ABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN. */
  /*     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1. */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZEXP, ZMLT, ZSQRT */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZASYI */
  /*     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZASYI */
  *nz = 0;
  az = zabs(zr, zi);
  arm = std::numeric_limits<double>::min() * 1e3;
  rtr1 = std::sqrt(arm);
  il = std::min(2, n);
  dfnu = fnu + (n - il);
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW TEST */
  /* ----------------------------------------------------------------------- */
  raz = 1. / az;
  str = zr * raz;
  sti = -zi * raz;
  ak1r = rtpi * str * raz;
  ak1i = rtpi * sti * raz;
  zsqrt(ak1r, ak1i, &ak1r, &ak1i);
  czr = zr;
  czi = zi;
  if (kode != 2) {
    goto L10;
  }
  czr = 0.;
  czi = zi;
L10:
  if (std::fabs(czr) > elim) {
    goto L100;
  }
  dnu2 = dfnu + dfnu;
  koded = 1;
  if (std::fabs(czr) > alim && n > 2) {
    goto L20;
  }
  koded = 0;
  zexp(czr, czi, &str, &sti);
  zmlt(ak1r, ak1i, str, sti, &ak1r, &ak1i);
L20:
  fdn = 0.;
  if (dnu2 > rtr1) {
    fdn = dnu2 * dnu2;
  }
  ezr = zr * 8.;
  ezi = zi * 8.;
  /* ----------------------------------------------------------------------- */
  /*     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE */
  /*     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE */
  /*     EXPANSION FOR THE IMAGINARY PART. */
  /* ----------------------------------------------------------------------- */
  aez = az * 8.;
  s = tol / aez;
  jl = (int)(rl + rl + 2);
  p1r = 0.;
  p1i = 0.;
  if (zi == 0.) {
    goto L30;
  }
  /* ----------------------------------------------------------------------- */
  /*     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF */
  /*     SIGNIFICANCE WHEN FNU OR N IS LARGE */
  /* ----------------------------------------------------------------------- */
  inu = (int)(fnu);
  arg = (fnu - inu) * pi;
  inu = inu + n - il;
  ak = -std::sin(arg);
  bk = std::cos(arg);
  if (zi < 0.) {
    bk = -bk;
  }
  p1r = ak;
  p1i = bk;
  if (inu % 2 == 0) {
    goto L30;
  }
  p1r = -p1r;
  p1i = -p1i;
L30:
  for (k = 1; k <= il; ++k) {
    sqk = fdn - 1.;
    atol = s * std::fabs(sqk);
    sgn = 1.;
    cs1r = 1.;
    cs1i = 0.;
    cs2r = 1.;
    cs2i = 0.;
    ckr = 1.;
    cki = 0.;
    ak = 0.;
    aa = 1.;
    bb = aez;
    dkr = ezr;
    dki = ezi;
    for (j = 1; j <= jl; ++j) {
      zdiv(ckr, cki, dkr, dki, &str, &sti);
      ckr = str * sqk;
      cki = sti * sqk;
      cs2r += ckr;
      cs2i += cki;
      sgn = -sgn;
      cs1r += ckr * sgn;
      cs1i += cki * sgn;
      dkr += ezr;
      dki += ezi;
      aa = aa * std::fabs(sqk) / bb;
      bb += aez;
      ak += 8.;
      sqk -= ak;
      if (aa <= atol) {
        goto L50;
      }
      /* L40: */
    }
    goto L110;
  L50:
    s2r = cs1r;
    s2i = cs1i;
    if (zr + zr >= elim) {
      goto L60;
    }
    tzr = zr + zr;
    tzi = zi + zi;
    zexp(-tzr, -tzi, &str, &sti);
    zmlt(str, sti, p1r, p1i, &str, &sti);
    zmlt(str, sti, cs2r, cs2i, &str, &sti);
    s2r += str;
    s2i += sti;
  L60:
    fdn = fdn + dfnu * 8. + 4.;
    p1r = -p1r;
    p1i = -p1i;
    m = n - il + k;
    yr[m] = s2r * ak1r - s2i * ak1i;
    yi[m] = s2r * ak1i + s2i * ak1r;
    /* L70: */
  }
  if (n <= 2) {
    return;
  }
  nn = n;
  k = nn - 2;
  ak = (double)k;
  str = zr * raz;
  sti = -zi * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  ib = 3;
  for (i__ = ib; i__ <= nn; ++i__) {
    yr[k] = (ak + fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
    yi[k] = (ak + fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
    ak += -1.;
    --k;
    /* L80: */
  }
  if (koded == 0) {
    return;
  }
  zexp(czr, czi, &ckr, &cki);
  for (i__ = 1; i__ <= nn; ++i__) {
    str = yr[i__] * ckr - yi[i__] * cki;
    yi[i__] = yr[i__] * cki + yi[i__] * ckr;
    yr[i__] = str;
    /* L90: */
  }
  return;
L100:
  *nz = -1;
  return;
L110:
  *nz = -2;
}

}  // namespace zbessel
