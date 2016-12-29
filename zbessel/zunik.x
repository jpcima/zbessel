#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <limits>
#include <cmath>

namespace zbessel {

template <class>
void zunik(double zrr, double zri, double fnu, int ikflg, int ipmtr,
           double tol, int init, double *__restrict__ phir, double *__restrict__ phii, double *__restrict__ zeta1r,
           double *__restrict__ zeta1i, double *__restrict__ zeta2r, double *__restrict__ zeta2i, double *__restrict__ sumr,
           double *__restrict__ sumi, double *__restrict__ cwrkr, double *__restrict__ cwrki) {
  /* Initialized data */

  const double con[2] = {
    0.39894228040143267794L,  // 1/sqrt(2*pi)
    1.25331413731550025121L};  // sqrt(pi/2)
  const double c__[120] = {1., -.208333333333333333, .125, .334201388888888889, -.401041666666666667, .0703125, -1.02581259645061728, 1.84646267361111111, -.8912109375, .0732421875, 4.66958442342624743, -11.2070026162229938, 8.78912353515625, -2.3640869140625, .112152099609375, -28.2120725582002449, 84.6362176746007346, -91.8182415432400174, 42.5349987453884549, -7.3687943594796317, .227108001708984375, 212.570130039217123, -765.252468141181642, 1059.99045252799988, -699.579627376132541, 218.19051174421159, -26.4914304869515555, .572501420974731445, -1919.457662318407, 8061.72218173730938, -13586.5500064341374, 11655.3933368645332, -5305.64697861340311, 1200.90291321635246, -108.090919788394656, 1.7277275025844574, 20204.2913309661486, -96980.5983886375135, 192547.001232531532, -203400.177280415534, 122200.46498301746, -41192.6549688975513, 7109.51430248936372, -493.915304773088012, 6.07404200127348304, -242919.187900551333, 1311763.6146629772, -2998015.91853810675, 3763271.297656404, -2813563.22658653411, 1268365.27332162478, -331645.172484563578, 45218.7689813627263, -2499.83048181120962, 24.3805296995560639, 3284469.85307203782, -19706819.1184322269, 50952602.4926646422, -74105148.2115326577, 66344512.2747290267, -37567176.6607633513, 13288767.1664218183, -2785618.12808645469, 308186.404612662398, -13886.0897537170405, 110.017140269246738, -49329253.664509962, 325573074.185765749, -939462359.681578403, 1553596899.57058006, -1621080552.10833708, 1106842816.82301447, -495889784.275030309, 142062907.797533095, -24474062.7257387285, 2243768.17792244943, -84005.4336030240853, 551.335896122020586, 814789096.118312115, -5866481492.05184723, 18688207509.2958249, -34632043388.1587779, 41280185579.753974, -33026599749.8007231, 17954213731.1556001, -6563293792.61928433, 1559279864.87925751, -225105661.889415278, 17395107.5539781645, -549842.327572288687, 3038.09051092238427, -14679261247.6956167, 114498237732.02581, -399096175224.466498, 819218669548.577329, -1098375156081.22331, 1008158106865.38209, -645364869245.376503, 287900649906.150589, -87867072178.0232657, 17634730606.8349694, -2167164983.22379509, 143157876.718888981, -3871833.44257261262, 18257.7554742931747, 286464035717.679043, -2406297900028.50396, 9109341185239.89896, -20516899410934.4374, 30565125519935.3206, -31667088584785.1584, 23348364044581.8409, -12320491305598.2872, 4612725780849.13197, -1196552880196.1816, 205914503232.410016, -21822927757.5292237, 1247009293.51271032, -29188388.1222208134, 118838.426256783253};

  /* Local variables */
  int i__, j, k, l;
  double ac, si, ti, sr, tr, t2i, t2r, rfn, sri, sti, zni, srr, str, znr;
  double test, crfni, crfnr;

  /* ***BEGIN PROLOGUE  ZUNIK */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUNIK-A, ZUNIK-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC */
  /*        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2 */
  /*        RESPECTIVELY BY */

  /*        W(FNU,ZR) = PHI*EXP(ZETA)*SUM */

  /*        WHERE       ZETA=-ZETA1 + ZETA2       OR */
  /*                          ZETA1 - ZETA2 */

  /*        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE */
  /*        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG= */
  /*        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK */
  /*        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI, */
  /*        ZETA1,ZETA2. */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZDIV, ZLOG, ZSQRT */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added EXTERNAL statement with ZLOG and ZSQRT.  (RWC) */
  /* ***END PROLOGUE  ZUNIK */
  /*     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1, */
  /*    *ZETA2,ZN,ZR */
  /* Parameter adjustments */
  --cwrki;
  --cwrkr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZUNIK */
  if (init != 0) {
    goto L40;
  }
  /* ----------------------------------------------------------------------- */
  /*     INITIALIZE ALL VARIABLES */
  /* ----------------------------------------------------------------------- */
  rfn = 1. / fnu;
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW TEST (ZR/FNU TOO SMALL) */
  /* ----------------------------------------------------------------------- */
  test = std::numeric_limits<double>::min() * 1e3;
  ac = fnu * test;
  if (std::fabs(zrr) > ac || std::fabs(zri) > ac) {
    goto L15;
  }
  *zeta1r = std::fabs(std::log(test)) * 2. + fnu;
  *zeta1i = 0.;
  *zeta2r = fnu;
  *zeta2i = 0.;
  *phir = 1.;
  *phii = 0.;
  return;
L15:
  tr = zrr * rfn;
  ti = zri * rfn;
  sr = 1. + (tr * tr - ti * ti);
  si = (tr * ti + ti * tr);
  zsqrt(sr, si, &srr, &sri);
  str = 1. + srr;
  sti = sri;
  zdiv(str, sti, tr, ti, &znr, &zni);
  zlog(znr, zni, &str, &sti);
  *zeta1r = fnu * str;
  *zeta1i = fnu * sti;
  *zeta2r = fnu * srr;
  *zeta2i = fnu * sri;
  zdiv(1., 0., srr, sri, &tr, &ti);
  srr = tr * rfn;
  sri = ti * rfn;
  zsqrt(srr, sri, &cwrkr[16], &cwrki[16]);
  *phir = cwrkr[16] * con[ikflg - 1];
  *phii = cwrki[16] * con[ikflg - 1];
  if (ipmtr != 0) {
    return;
  }
  zdiv(1., 0., sr, si, &t2r, &t2i);
  cwrkr[1] = 1.;
  cwrki[1] = 0.;
  crfnr = 1.;
  crfni = 0.;
  ac = 1.;
  l = 1;
  for (k = 2; k <= 15; ++k) {
    sr = 0.;
    si = 0.;
    for (j = 1; j <= k; ++j) {
      ++l;
      str = sr * t2r - si * t2i + c__[l - 1];
      si = sr * t2i + si * t2r;
      sr = str;
      /* L10: */
    }
    str = crfnr * srr - crfni * sri;
    crfni = crfnr * sri + crfni * srr;
    crfnr = str;
    cwrkr[k] = crfnr * sr - crfni * si;
    cwrki[k] = crfnr * si + crfni * sr;
    ac *= rfn;
    test = std::fabs(cwrkr[k]) + std::fabs(cwrki[k]);
    if (ac < tol && test < tol) {
      goto L30;
    }
    /* L20: */
  }
  k = 15;
L30:
  init = k;
L40:
  if (ikflg == 2) {
    goto L60;
  }
  /* ----------------------------------------------------------------------- */
  /*     COMPUTE SUM FOR THE I FUNCTION */
  /* ----------------------------------------------------------------------- */
  sr = 0.;
  si = 0.;
  for (i__ = 1; i__ <= init; ++i__) {
    sr += cwrkr[i__];
    si += cwrki[i__];
    /* L50: */
  }
  *sumr = sr;
  *sumi = si;
  *phir = cwrkr[16] * con[0];
  *phii = cwrki[16] * con[0];
  return;
L60:
  /* ----------------------------------------------------------------------- */
  /*     COMPUTE SUM FOR THE K FUNCTION */
  /* ----------------------------------------------------------------------- */
  sr = 0.;
  si = 0.;
  tr = 1.;
  for (i__ = 1; i__ <= init; ++i__) {
    sr += tr * cwrkr[i__];
    si += tr * cwrki[i__];
    tr = -tr;
    /* L70: */
  }
  *sumr = sr;
  *sumi = si;
  *phir = cwrkr[16] * con[1];
  *phii = cwrki[16] * con[1];
}

}  // namespace zbessel
