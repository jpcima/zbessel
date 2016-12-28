#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
void zkscl(double zrr, double zri, double fnu, int n, double *__restrict__ yr,
           double *__restrict__ yi, int *__restrict__ nz, double rzr, double rzi, double ascle,
           double tol, double elim) {
  /* Local variables */
  int i__, ic;
  double as, fn;
  int kk, nn, nw;
  double s1i, s2i, s1r, s2r, acs, cki, elm, csi, ckr, cyi[2], zdi, csr,
      cyr[2], zdr, str, alas;
  double helim, celmr;

  /* ***BEGIN PROLOGUE  ZKSCL */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CKSCL-A, ZKSCL-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE */
  /*     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN */
  /*     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL. */

  /* ***SEE ALSO  ZBESK */
  /* ***ROUTINES CALLED  ZABS, ZLOG, ZUCHK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added ZLOG to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZKSCL */
  /*     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZKSCL */
  *nz = 0;
  ic = 0;
  nn = std::min(2, n);
  for (i__ = 1; i__ <= nn; ++i__) {
    s1r = yr[i__];
    s1i = yi[i__];
    cyr[i__ - 1] = s1r;
    cyi[i__ - 1] = s1i;
    as = zabs(s1r, s1i);
    acs = -zrr + std::log(as);
    ++(*nz);
    yr[i__] = 0.;
    yi[i__] = 0.;
    if (acs < -elim) {
      goto L10;
    }
    zlog(s1r, s1i, &csr, &csi);
    csr -= zrr;
    csi -= zri;
    str = std::exp(csr) / tol;
    csr = str * std::cos(csi);
    csi = str * std::sin(csi);
    zuchk(csr, csi, &nw, ascle, tol);
    if (nw != 0) {
      goto L10;
    }
    yr[i__] = csr;
    yi[i__] = csi;
    ic = i__;
    --(*nz);
  L10:;
  }
  if (n == 1) {
    return;
  }
  if (ic > 1) {
    goto L20;
  }
  yr[1] = 0.;
  yi[1] = 0.;
  *nz = 2;
L20:
  if (n == 2) {
    return;
  }
  if (*nz == 0) {
    return;
  }
  fn = fnu + 1.;
  ckr = fn * rzr;
  cki = fn * rzi;
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  helim = elim * .5;
  elm = std::exp(-elim);
  celmr = elm;
  zdr = zrr;
  zdi = zri;

  /*     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF */
  /*     S2 GETS LARGER THAN EXP(ELIM/2) */

  for (i__ = 3; i__ <= n; ++i__) {
    kk = i__;
    csr = s2r;
    csi = s2i;
    s2r = ckr * csr - cki * csi + s1r;
    s2i = cki * csr + ckr * csi + s1i;
    s1r = csr;
    s1i = csi;
    ckr += rzr;
    cki += rzi;
    as = zabs(s2r, s2i);
    alas = std::log(as);
    acs = -zdr + alas;
    ++(*nz);
    yr[i__] = 0.;
    yi[i__] = 0.;
    if (acs < -elim) {
      goto L25;
    }
    zlog(s2r, s2i, &csr, &csi);
    csr -= zdr;
    csi -= zdi;
    str = std::exp(csr) / tol;
    csr = str * std::cos(csi);
    csi = str * std::sin(csi);
    zuchk(csr, csi, &nw, ascle, tol);
    if (nw != 0) {
      goto L25;
    }
    yr[i__] = csr;
    yi[i__] = csi;
    --(*nz);
    if (ic == kk - 1) {
      goto L40;
    }
    ic = kk;
    goto L30;
  L25:
    if (alas < helim) {
      goto L30;
    }
    zdr -= elim;
    s1r *= celmr;
    s1i *= celmr;
    s2r *= celmr;
    s2i *= celmr;
  L30:;
  }
  *nz = n;
  if (ic == n) {
    *nz = n - 1;
  }
  goto L45;
L40:
  *nz = kk - 2;
L45:
  for (i__ = 1; i__ <= *nz; ++i__) {
    yr[i__] = 0.;
    yi[i__] = 0.;
    /* L50: */
  }
}

}  // namespace zbessel
