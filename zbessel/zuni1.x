#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <limits>
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
void zuni1(double zr, double zi, double fnu, int kode, int n, double *__restrict__ yr,
           double *__restrict__ yi, int *__restrict__ nz, int *__restrict__ nlast, double fnul, double tol,
           double elim, double alim) {
  /* Local variables */
  int i__, k, m, nd;
  double fn;
  int nn, nw;
  double c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, cyi[2];
  int nuf;
  double bry[3], cyr[2], sti, rzi, str, rzr, aphi, cscl, phii, crsc;
  double phir;
  int init;
  double csrr[3], cssr[3], rast, sumi, sumr;
  int iflag;
  double ascle, cwrki[16];
  double cwrkr[16];
  double zeta1i, zeta2i, zeta1r, zeta2r;

  /* ***BEGIN PROLOGUE  ZUNI1 */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUNI1-A, ZUNI1-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC */
  /*     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3. */

  /*     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC */
  /*     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET. */
  /*     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER */
  /*     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL. */
  /*     Y(I)=CZERO FOR I=NLAST+1,N */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZUCHK, ZUNIK, ZUOIK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZUNI1 */
  /*     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1, */
  /*    *S2,Y,Z,ZETA1,ZETA2 */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZUNI1 */
  *nz = 0;
  nd = n;
  *nlast = 0;
  /* ----------------------------------------------------------------------- */
  /*     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG- */
  /*     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE, */
  /*     EXP(ALIM)=EXP(ELIM)*TOL */
  /* ----------------------------------------------------------------------- */
  cscl = 1. / tol;
  crsc = tol;
  cssr[0] = cscl;
  cssr[1] = 1.;
  cssr[2] = crsc;
  csrr[0] = crsc;
  csrr[1] = 1.;
  csrr[2] = cscl;
  bry[0] = std::numeric_limits<double>::min() * 1e3 / tol;
  /* ----------------------------------------------------------------------- */
  /*     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER */
  /* ----------------------------------------------------------------------- */
  fn = std::max(fnu, 1.);
  init = 0;
  zunik(zr, zi, fn, 1, 1, tol, init, &phir, &phii, &zeta1r, &zeta1i,
        &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
  if (kode == 1) {
    goto L10;
  }
  str = zr + zeta2r;
  sti = zi + zeta2i;
  rast = fn / zabs(str, sti);
  str = str * rast * rast;
  sti = -sti * rast * rast;
  s1r = -zeta1r + str;
  s1i = -zeta1i + sti;
  goto L20;
L10:
  s1r = -zeta1r + zeta2r;
  s1i = -zeta1i + zeta2i;
L20:
  rs1 = s1r;
  if (std::fabs(rs1) > elim) {
    goto L130;
  }
L30:
  nn = std::min(2, nd);
  for (i__ = 1; i__ <= nn; ++i__) {
    fn = fnu + (nd - i__);
    init = 0;
    zunik(zr, zi, fn, 1, 0, tol, init, &phir, &phii, &zeta1r,
          &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
    if (kode == 1) {
      goto L40;
    }
    str = zr + zeta2r;
    sti = zi + zeta2i;
    rast = fn / zabs(str, sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zeta1r + str;
    s1i = -zeta1i + sti + zi;
    goto L50;
  L40:
    s1r = -zeta1r + zeta2r;
    s1i = -zeta1i + zeta2i;
  L50:
    // -----------------------------------------------------------------------/
    //     TEST FOR UNDERFLOW AND OVERFLOW */
    // -----------------------------------------------------------------------/
    rs1 = s1r;
    if (std::fabs(rs1) > elim) {
      goto L110;
    }
    if (i__ == 1) {
      iflag = 2;
    }
    if (std::fabs(rs1) < alim) {
      goto L60;
    }
    // -----------------------------------------------------------------------/
    //     REFINE  TEST AND SCALE */
    // -----------------------------------------------------------------------/
    aphi = zabs(phir, phii);
    rs1 += std::log(aphi);
    if (std::fabs(rs1) > elim) {
      goto L110;
    }
    if (i__ == 1) {
      iflag = 1;
    }
    if (rs1 < 0.) {
      goto L60;
    }
    if (i__ == 1) {
      iflag = 3;
    }
  L60:
    /* -----------------------------------------------------------------------
     */
    /*     SCALE S1 IF ABS(S1).LT.ASCLE */
    /* -----------------------------------------------------------------------
     */
    s2r = phir * sumr - phii * sumi;
    s2i = phir * sumi + phii * sumr;
    str = std::exp(s1r) * cssr[iflag - 1];
    s1r = str * std::cos(s1i);
    s1i = str * std::sin(s1i);
    str = s2r * s1r - s2i * s1i;
    s2i = s2r * s1i + s2i * s1r;
    s2r = str;
    if (iflag != 1) {
      goto L70;
    }
    zuchk(s2r, s2i, &nw, bry[0], tol);
    if (nw != 0) {
      goto L110;
    }
  L70:
    cyr[i__ - 1] = s2r;
    cyi[i__ - 1] = s2i;
    m = nd - i__ + 1;
    yr[m] = s2r * csrr[iflag - 1];
    yi[m] = s2i * csrr[iflag - 1];
    /* L80: */
  }
  if (nd <= 2) {
    goto L100;
  }
  rast = 1. / zabs(zr, zi);
  str = zr * rast;
  sti = -zi * rast;
  rzr = (str + str) * rast;
  rzi = (sti + sti) * rast;
  bry[1] = 1. / bry[0];
  bry[2] = std::numeric_limits<double>::max();
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  c1r = csrr[iflag - 1];
  ascle = bry[iflag - 1];
  k = nd - 2;
  fn = (double)k;
  for (i__ = 3; i__ <= nd; ++i__) {
    c2r = s2r;
    c2i = s2i;
    s2r = s1r + (fnu + fn) * (rzr * c2r - rzi * c2i);
    s2i = s1i + (fnu + fn) * (rzr * c2i + rzi * c2r);
    s1r = c2r;
    s1i = c2i;
    c2r = s2r * c1r;
    c2i = s2i * c1r;
    yr[k] = c2r;
    yi[k] = c2i;
    --k;
    fn += -1.;
    if (iflag >= 3) {
      goto L90;
    }
    str = std::abs(c2r);
    sti = std::abs(c2i);
    c2m = std::max(str, sti);
    if (c2m <= ascle) {
      goto L90;
    }
    ++iflag;
    ascle = bry[iflag - 1];
    s1r *= c1r;
    s1i *= c1r;
    s2r = c2r;
    s2i = c2i;
    s1r *= cssr[iflag - 1];
    s1i *= cssr[iflag - 1];
    s2r *= cssr[iflag - 1];
    s2i *= cssr[iflag - 1];
    c1r = csrr[iflag - 1];
  L90:;
  }
L100:
  return;
/* ----------------------------------------------------------------------- */
/*     SET UNDERFLOW AND UPDATE PARAMETERS */
/* ----------------------------------------------------------------------- */
L110:
  if (rs1 > 0.) {
    goto L120;
  }
  yr[nd] = 0.;
  yi[nd] = 0.;
  ++(*nz);
  --nd;
  if (nd == 0) {
    goto L100;
  }
  zuoik(zr, zi, fnu, kode, 1, nd, &yr[1], &yi[1], &nuf, tol, elim, alim);
  if (nuf < 0) {
    goto L120;
  }
  nd -= nuf;
  *nz += nuf;
  if (nd == 0) {
    goto L100;
  }
  fn = fnu + (nd - 1);
  if (fn >= fnul) {
    goto L30;
  }
  *nlast = nd;
  return;
L120:
  *nz = -1;
  return;
L130:
  if (rs1 > 0.) {
    goto L120;
  }
  *nz = n;
  for (i__ = 1; i__ <= n; ++i__) {
    yr[i__] = 0.;
    yi[i__] = 0.;
    /* L140: */
  }
}

}  // namespace zbessel
