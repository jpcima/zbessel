#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <algorithm>

namespace zbessel {

template <class>
void zbinu(double zr, double zi, double fnu, int kode, int n, double *__restrict__ cyr,
           double *__restrict__ cyi, int *__restrict__ nz, double rl, double fnul, double tol,
           double elim, double alim) {
  /* Local variables */
  int i__;
  double az;
  int nn, nw;
  double cwi[2], cwr[2];
  int nui, inw;
  double dfnu;
  int nlast;

  /* ***BEGIN PROLOGUE  ZBINU */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CBINU-A, ZBINU-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE */

  /* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY */
  /* ***ROUTINES CALLED  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZBINU */
  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZBINU */
  *nz = 0;
  az = zabs(zr, zi);
  nn = n;
  dfnu = fnu + (n - 1);
  if (az <= 2.) {
    goto L10;
  }
  if (az * az * .25 > dfnu + 1.) {
    goto L20;
  }
L10:
  /* ----------------------------------------------------------------------- */
  /*     POWER SERIES */
  /* ----------------------------------------------------------------------- */
  zseri(zr, zi, fnu, kode, nn, &cyr[1], &cyi[1], &nw, tol, elim, alim);
  inw = std::abs(nw);
  *nz += inw;
  nn -= inw;
  if (nn == 0) {
    return;
  }
  if (nw >= 0) {
    goto L120;
  }
  dfnu = fnu + (nn - 1);
L20:
  if (az < rl) {
    goto L40;
  }
  if (dfnu <= 1.) {
    goto L30;
  }
  if (az + az < dfnu * dfnu) {
    goto L50;
  }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z */
/* ----------------------------------------------------------------------- */
L30:
  zasyi(zr, zi, fnu, kode, nn, &cyr[1], &cyi[1], &nw, rl, tol, elim, alim);
  if (nw < 0) {
    goto L130;
  }
  goto L120;
L40:
  if (dfnu <= 1.) {
    goto L70;
  }
L50:
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM */
  /* ----------------------------------------------------------------------- */
  zuoik(zr, zi, fnu, kode, 1, nn, &cyr[1], &cyi[1], &nw, tol, elim, alim);
  if (nw < 0) {
    goto L130;
  }
  *nz += nw;
  nn -= nw;
  if (nn == 0) {
    return;
  }
  dfnu = fnu + (nn - 1);
  if (dfnu > fnul) {
    goto L110;
  }
  if (az > fnul) {
    goto L110;
  }
L60:
  if (az > rl) {
    goto L80;
  }
L70:
  /* ----------------------------------------------------------------------- */
  /*     MILLER ALGORITHM NORMALIZED BY THE SERIES */
  /* ----------------------------------------------------------------------- */
  zmlri(zr, zi, fnu, kode, nn, &cyr[1], &cyi[1], &nw, tol);
  if (nw < 0) {
    goto L130;
  }
  goto L120;
L80:
  /* ----------------------------------------------------------------------- */
  /*     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN */
  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN */
  /* ----------------------------------------------------------------------- */
  zuoik(zr, zi, fnu, kode, 2, 2, cwr, cwi, &nw, tol, elim, alim);
  if (nw >= 0) {
    goto L100;
  }
  *nz = nn;
  for (i__ = 1; i__ <= nn; ++i__) {
    cyr[i__] = 0.;
    cyi[i__] = 0.;
    /* L90: */
  }
  return;
L100:
  if (nw > 0) {
    goto L130;
  }
  zwrsk(zr, zi, fnu, kode, nn, &cyr[1], &cyi[1], &nw, cwr, cwi, tol, elim,
        alim);
  if (nw < 0) {
    goto L130;
  }
  goto L120;
L110:
  /* ----------------------------------------------------------------------- */
  /*     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD */
  /* ----------------------------------------------------------------------- */
  nui = (int)(fnul - dfnu + 1);
  nui = std::max(nui, 0);
  zbuni(zr, zi, fnu, kode, nn, &cyr[1], &cyi[1], &nw, nui, &nlast, fnul, tol,
        elim, alim);
  if (nw < 0) {
    goto L130;
  }
  *nz += nw;
  if (nlast == 0) {
    goto L120;
  }
  nn = nlast;
  goto L60;
L120:
  return;
L130:
  *nz = -1;
  if (nw == -2) {
    *nz = -2;
  }
}

}  // namespace zbessel
