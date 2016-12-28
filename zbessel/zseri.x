#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <limits>
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
void zseri(double zr, double zi, double fnu, int kode, int n, double *__restrict__ yr,
           double *__restrict__ yi, int *__restrict__ nz, double tol, double elim, double alim) {
  /* Local variables */
  int i__, k, l, m;
  double s, aa;
  int ib;
  double ak;
  int il;
  double az;
  int nn;
  double wi[2], rs, ss;
  int nw;
  double wr[2], s1i, s2i, s1r, s2r, cki, acz, arm, ckr, czi, hzi, raz,
      czr, sti, hzr, rzi, str, rzr, ak1i, ak1r, rtr1, dfnu;
  double atol;
  double fnup;
  int iflag;
  double coefi, ascle, coefr, crscr;

  /* ***BEGIN PROLOGUE  ZSERI */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CSERI-A, ZSERI-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY */
  /*     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE */
  /*     REGION ABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN. */
  /*     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO */
  /*     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE */
  /*     CONDITION ABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE */
  /*     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ). */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZDIV, ZLOG, ZMLT, ZUCHK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added ZLOG to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZSERI */
  /*     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZSERI */
  *nz = 0;
  az = zabs(zr, zi);
  if (az == 0.) {
    goto L160;
  }
  arm = std::numeric_limits<double>::min() * 1e3;
  rtr1 = std::sqrt(arm);
  crscr = 1.;
  iflag = 0;
  if (az < arm) {
    goto L150;
  }
  hzr = zr * .5;
  hzi = zi * .5;
  czr = 0.;
  czi = 0.;
  if (az <= rtr1) {
    goto L10;
  }
  zmlt(hzr, hzi, hzr, hzi, &czr, &czi);
L10:
  acz = zabs(czr, czi);
  nn = n;
  zlog(hzr, hzi, &ckr, &cki);
L20:
  dfnu = fnu + (nn - 1);
  fnup = dfnu + 1.;
  /* ----------------------------------------------------------------------- */
  /*     UNDERFLOW TEST */
  /* ----------------------------------------------------------------------- */
  ak1r = ckr * dfnu;
  ak1i = cki * dfnu;
  ak = std::lgamma(fnup);
  ak1r -= ak;
  if (kode == 2) {
    ak1r -= zr;
  }
  if (ak1r > -elim) {
    goto L40;
  }
L30:
  ++(*nz);
  yr[nn] = 0.;
  yi[nn] = 0.;
  if (acz > dfnu) {
    goto L190;
  }
  --nn;
  if (nn == 0) {
    return;
  }
  goto L20;
L40:
  if (ak1r > -alim) {
    goto L50;
  }
  iflag = 1;
  ss = 1. / tol;
  crscr = tol;
  ascle = arm * ss;
L50:
  aa = std::exp(ak1r);
  if (iflag == 1) {
    aa *= ss;
  }
  coefr = aa * std::cos(ak1i);
  coefi = aa * std::sin(ak1i);
  atol = tol * acz / fnup;
  il = std::min(2, nn);
  for (i__ = 1; i__ <= il; ++i__) {
    dfnu = fnu + (nn - i__);
    fnup = dfnu + 1.;
    s1r = 1.;
    s1i = 0.;
    if (acz < tol * fnup) {
      goto L70;
    }
    ak1r = 1.;
    ak1i = 0.;
    ak = fnup + 2.;
    s = fnup;
    aa = 2.;
  L60:
    rs = 1. / s;
    str = ak1r * czr - ak1i * czi;
    sti = ak1r * czi + ak1i * czr;
    ak1r = str * rs;
    ak1i = sti * rs;
    s1r += ak1r;
    s1i += ak1i;
    s += ak;
    ak += 2.;
    aa = aa * acz * rs;
    if (aa > atol) {
      goto L60;
    }
  L70:
    s2r = s1r * coefr - s1i * coefi;
    s2i = s1r * coefi + s1i * coefr;
    wr[i__ - 1] = s2r;
    wi[i__ - 1] = s2i;
    if (iflag == 0) {
      goto L80;
    }
    zuchk(s2r, s2i, &nw, ascle, tol);
    if (nw != 0) {
      goto L30;
    }
  L80:
    m = nn - i__ + 1;
    yr[m] = s2r * crscr;
    yi[m] = s2i * crscr;
    if (i__ == il) {
      goto L90;
    }
    zdiv(coefr, coefi, hzr, hzi, &str, &sti);
    coefr = str * dfnu;
    coefi = sti * dfnu;
  L90:;
  }
  if (nn <= 2) {
    return;
  }
  k = nn - 2;
  ak = (double)k;
  raz = 1. / az;
  str = zr * raz;
  sti = -zi * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  if (iflag == 1) {
    goto L120;
  }
  ib = 3;
L100:
  for (i__ = ib; i__ <= nn; ++i__) {
    yr[k] = (ak + fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
    yi[k] = (ak + fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
    ak += -1.;
    --k;
    /* L110: */
  }
  return;
/* ----------------------------------------------------------------------- */
/*     RECUR BACKWARD WITH SCALED VALUES */
/* ----------------------------------------------------------------------- */
L120:
  /* ----------------------------------------------------------------------- */
  /*     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE */
  /*     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3 */
  /* ----------------------------------------------------------------------- */
  s1r = wr[0];
  s1i = wi[0];
  s2r = wr[1];
  s2i = wi[1];
  for (l = 3; l <= nn; ++l) {
    ckr = s2r;
    cki = s2i;
    s2r = s1r + (ak + fnu) * (rzr * ckr - rzi * cki);
    s2i = s1i + (ak + fnu) * (rzr * cki + rzi * ckr);
    s1r = ckr;
    s1i = cki;
    ckr = s2r * crscr;
    cki = s2i * crscr;
    yr[k] = ckr;
    yi[k] = cki;
    ak += -1.;
    --k;
    if (zabs(ckr, cki) > ascle) {
      goto L140;
    }
    /* L130: */
  }
  return;
L140:
  ib = l + 1;
  if (ib > nn) {
    return;
  }
  goto L100;
L150:
  *nz = n;
  if (fnu == 0.) {
    --(*nz);
  }
L160:
  yr[1] = 0.;
  yi[1] = 0.;
  if (fnu != 0.) {
    goto L170;
  }
  yr[1] = 1.;
  yi[1] = 0.;
L170:
  if (n == 1) {
    return;
  }
  for (i__ = 2; i__ <= n; ++i__) {
    yr[i__] = 0.;
    yi[i__] = 0.;
    /* L180: */
  }
  return;
/* ----------------------------------------------------------------------- */
/*     RETURN WITH NZ.LT.0 IF ABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE */
/*     THE CALCULATION IN CBINU WITH N=N-ABS(NZ) */
/* ----------------------------------------------------------------------- */
L190:
  *nz = -(*nz);
}

}  // namespace zbessel
