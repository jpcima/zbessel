#pragma once
#include "zbsubr.h"
#include "zops.h"
#include <limits>
#include <algorithm>
#include <cmath>

namespace zbessel {

template <class>
int zbesh(double zr, double zi, double fnu, int kode, int m,
          int n, double *__restrict__ cyr, double *__restrict__ cyi, int *__restrict__ nz) {
  static const double r1m5 = std::log10(std::numeric_limits<double>::radix);

  /* Initialized data */

  const double hpi = 1.570796326794896558L;  // pi/2

  /* Local variables */
  int ierr;
  int i__, k, k1, k2;
  double aa, bb, fn;
  int mm;
  double az;
  int ir, nn;
  double rl;
  int mr, nw;
  double dig, arg, aln, fmm, ufl, sgn;
  int nuf, inu;
  double tol, sti, zni, zti, str, znr, alim, elim;
  double atol, rhpi;
  int inuh;
  double fnul, rtol, ascle, csgni;
  double csgnr;

  /* ***BEGIN PROLOGUE  ZBESH */
  /* ***PURPOSE  Compute a sequence of the Hankel functions H(m,a,z) */
  /*            for superscript m=1 or 2, real nonnegative orders a=b, */
  /*            b+1,... where b>0, and nonzero complex argument z.  A */
  /*            scaling option is available to help avoid overflow. */
  /* ***LIBRARY   SLATEC */
  /* ***CATEGORY  C10A4 */
  /* ***TYPE      COMPLEX (CBESH-C, ZBESH-C) */
  /* ***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, */
  /*             BESSEL FUNCTIONS OF THE THIRD KIND, H BESSEL FUNCTIONS, */
  /*             HANKEL FUNCTIONS */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*                      ***A DOUBLE PRECISION ROUTINE*** */
  /*         On KODE=1, ZBESH computes an N member sequence of complex */
  /*         Hankel (Bessel) functions CY(L)=H(M,FNU+L-1,Z) for super- */
  /*         script M=1 or 2, real nonnegative orders FNU+L-1, L=1,..., */
  /*         N, and complex nonzero Z in the cut plane -pi<arg(Z)<=pi */
  /*         where Z=ZR+i*ZI.  On KODE=2, CBESH returns the scaled */
  /*         functions */

  /*            CY(L) = H(M,FNU+L-1,Z)*exp(-(3-2*M)*Z*i),  i**2=-1 */

  /*         which removes the exponential behavior in both the upper */
  /*         and lower half planes.  Definitions and notation are found */
  /*         in the NBS Handbook of Mathematical Functions (Ref. 1). */

  /*         Input */
  /*           ZR     - DOUBLE PRECISION real part of nonzero argument Z */
  /*           ZI     - DOUBLE PRECISION imag part of nonzero argument Z */
  /*           FNU    - DOUBLE PRECISION initial order, FNU>=0 */
  /*           KODE   - A parameter to indicate the scaling option */
  /*                    KODE=1  returns */
  /*                            CY(L)=H(M,FNU+L-1,Z), L=1,...,N */
  /*                        =2  returns */
  /*                            CY(L)=H(M,FNU+L-1,Z)*exp(-(3-2M)*Z*i), */
  /*                            L=1,...,N */
  /*           M      - Superscript of Hankel function, M=1 or 2 */
  /*           N      - Number of terms in the sequence, N>=1 */

  /*         Output */
  /*           CYR    - DOUBLE PRECISION real part of result vector */
  /*           CYI    - DOUBLE PRECISION imag part of result vector */
  /*           NZ     - Number of underflows set to zero */
  /*                    NZ=0    Normal return */
  /*                    NZ>0    CY(L)=0 for NZ values of L (if M=1 and */
  /*                            Im(Z)>0 or if M=2 and Im(Z)<0, then */
  /*                            CY(L)=0 for L=1,...,NZ; in the com- */
  /*                            plementary half planes, the underflows */
  /*                            may not be in an uninterrupted sequence) */
  /*           IERR   - Error flag */
  /*                    IERR=0  Normal return     - COMPUTATION COMPLETED */
  /*                    IERR=1  Input error       - NO COMPUTATION */
  /*                    IERR=2  Overflow          - NO COMPUTATION */
  /*                            (abs(Z) too small and/or FNU+N-1 */
  /*                            too large) */
  /*                    IERR=3  Precision warning - COMPUTATION COMPLETED */
  /*                            (Result has half precision or less */
  /*                            because abs(Z) or FNU+N-1 is large) */
  /*                    IERR=4  Precision error   - NO COMPUTATION */
  /*                            (Result has no precision because */
  /*                            abs(Z) or FNU+N-1 is too large) */
  /*                    IERR=5  Algorithmic error - NO COMPUTATION */
  /*                            (Termination condition not met) */

  /* *Long Description: */

  /*         The computation is carried out by the formula */

  /*            H(m,a,z) = (1/t)*exp(-a*t)*K(a,z*exp(-t)) */
  /*                   t = (3-2*m)*i*pi/2 */

  /*         where the K Bessel function is computed as described in the */
  /*         prologue to CBESK. */

  /*         Exponential decay of H(m,a,z) occurs in the upper half z */
  /*         plane for m=1 and the lower half z plane for m=2.  Exponential */
  /*         growth occurs in the complementary half planes.  Scaling */
  /*         by exp(-(3-2*m)*z*i) removes the exponential behavior in the */
  /*         whole z plane as z goes to infinity. */

  /*         For negative orders, the formula */

  /*            H(m,-a,z) = H(m,a,z)*exp((3-2*m)*a*pi*i) */

  /*         can be used. */

  /*         In most complex variable computation, one must evaluate ele- */
  /*         mentary functions.  When the magnitude of Z or FNU+N-1 is */
  /*         large, losses of significance by argument reduction occur. */
  /*         Consequently, if either one exceeds U1=SQRT(0.5/UR), then */
  /*         losses exceeding half precision are likely and an error flag */
  /*         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double */
  /*         precision unit roundoff limited to 18 digits precision.  Also, */
  /*         if either is larger than U2=0.5/UR, then all significance is */
  /*         lost and IERR=4.  In order to use the INT function, arguments */
  /*         must be further restricted not to exceed the largest machine */
  /*         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1 */
  /*         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and */
  /*         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision */
  /*         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This */
  /*         makes U2 limiting in single precision and U3 limiting in */
  /*         double precision.  This means that one can expect to retain, */
  /*         in the worst cases on IEEE machines, no digits in single pre- */
  /*         cision and only 6 digits in double precision.  Similar con- */
  /*         siderations hold for other machines. */

  /*         The approximate relative error in the magnitude of a complex */
  /*         Bessel function can be expressed as P*10**S where P=MAX(UNIT */
  /*         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre- */
  /*         sents the increase in error due to argument reduction in the */
  /*         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))), */
  /*         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF */
  /*         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may */
  /*         have only absolute accuracy.  This is most likely to occur */
  /*         when one component (in magnitude) is larger than the other by */
  /*         several orders of magnitude.  If one component is 10**K larger */
  /*         than the other, then one can expect only MAX(ABS(LOG10(P))-K, */
  /*         0) significant digits; or, stated another way, when K exceeds */
  /*         the exponent of P, no significant digits remain in the smaller */
  /*         component.  However, the phase angle retains absolute accuracy */
  /*         because, in complex arithmetic with precision P, the smaller */
  /*         component will not (as a rule) decrease below P times the */
  /*         magnitude of the larger component.  In these extreme cases, */
  /*         the principal phase angle is on the order of +P, -P, PI/2-P, */
  /*         or -PI/2+P. */

  /* ***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe- */
  /*                 matical Functions, National Bureau of Standards */
  /*                 Applied Mathematics Series 55, U. S. Department */
  /*                 of Commerce, Tenth Printing (1972) or later. */
  /*               2. D. E. Amos, Computation of Bessel Functions of */
  /*                 Complex Argument, Report SAND83-0086, Sandia National */
  /*                 Laboratories, Albuquerque, NM, May 1983. */
  /*               3. D. E. Amos, Computation of Bessel Functions of */
  /*                 Complex Argument and Large Order, Report SAND83-0643, */
  /*                 Sandia National Laboratories, Albuquerque, NM, May */
  /*                 1983. */
  /*               4. D. E. Amos, A Subroutine Package for Bessel Functions */
  /*                 of a Complex Argument and Nonnegative Order, Report */
  /*                 SAND85-1018, Sandia National Laboratory, Albuquerque, */
  /*                 NM, May 1985. */
  /*               5. D. E. Amos, A portable package for Bessel functions */
  /*                 of a complex argument and nonnegative order, ACM */
  /*                 Transactions on Mathematical Software, 12 (September */
  /*                 1986), pp. 265-273. */

  /* ***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZACON, ZBKNU, ZBUNK, ZUOIK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   890801  REVISION DATE from Version 3.2 */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   920128  Category corrected.  (WRB) */
  /*   920811  Prologue revised.  (DWL) */
  /* ***END PROLOGUE  ZBESH */

  /*     COMPLEX CY,Z,ZN,ZT,CSGN */

  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */

  /* ***FIRST EXECUTABLE STATEMENT  ZBESH */
  ierr = 0;
  *nz = 0;
  if (zr == 0. && zi == 0.) {
    ierr = 1;
  }
  if (fnu < 0.) {
    ierr = 1;
  }
  if (m < 1 || m > 2) {
    ierr = 1;
  }
  if (kode < 1 || kode > 2) {
    ierr = 1;
  }
  if (n < 1) {
    ierr = 1;
  }
  if (ierr != 0) {
    return ierr;
  }
  nn = n;
  /* ----------------------------------------------------------------------- */
  /*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
  /*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. */
  /*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
  /*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
  /*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
  /*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
  /*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
  /*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
  /*     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU */
  /* ----------------------------------------------------------------------- */
  /* Computing MAX */
  tol = std::max(std::numeric_limits<double>::epsilon(), 1e-18);
  k1 = std::numeric_limits<double>::min_exponent;
  k2 = std::numeric_limits<double>::max_exponent;
  /* Computing MIN */
  k = std::min(std::abs(k1), std::abs(k2));
  elim = (k * r1m5 - 3.) * 2.303;
  k1 = std::numeric_limits<double>::digits - 1;
  aa = r1m5 * k1;
  dig = std::min(aa, 18.);
  aa *= 2.303;
  /* Computing MAX */
  alim = elim + std::max(-aa, -41.45);
  fnul = (dig - 3.) * 6. + 10.;
  rl = dig * 1.2 + 3.;
  fn = fnu + (nn - 1);
  mm = 3 - m - m;
  fmm = (double)mm;
  znr = fmm * zi;
  zni = -fmm * zr;
  /* ----------------------------------------------------------------------- */
  /*     TEST FOR PROPER RANGE */
  /* ----------------------------------------------------------------------- */
  az = zabs(zr, zi);
  aa = .5 / tol;
  bb = std::numeric_limits<int>::max() * .5;
  aa = std::min(aa, bb);
  if (az > aa) {
    goto L260;
  }
  if (fn > aa) {
    goto L260;
  }
  aa = std::sqrt(aa);
  if (az > aa) {
    ierr = 3;
  }
  if (fn > aa) {
    ierr = 3;
  }
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE */
  /* ----------------------------------------------------------------------- */
  ufl = std::numeric_limits<double>::min() * 1e3;
  if (az < ufl) {
    goto L230;
  }
  if (fnu > fnul) {
    goto L90;
  }
  if (fn <= 1.) {
    goto L70;
  }
  if (fn > 2.) {
    goto L60;
  }
  if (az > tol) {
    goto L70;
  }
  arg = az * .5;
  aln = -fn * std::log(arg);
  if (aln > elim) {
    goto L230;
  }
  goto L70;
L60:
  zuoik(znr, zni, fnu, kode, 2, nn, &cyr[1], &cyi[1], &nuf, tol, elim, alim);
  if (nuf < 0) {
    goto L230;
  }
  *nz += nuf;
  nn -= nuf;
  /* ----------------------------------------------------------------------- */
  /*     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK */
  /*     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I */
  /* ----------------------------------------------------------------------- */
  if (nn == 0) {
    goto L140;
  }
L70:
  if (znr < 0. || (znr == 0. && zni < 0. && m == 2)) {
    goto L80;
  }
  /* ----------------------------------------------------------------------- */
  /*     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR. */
  /*     YN.GE.0. .OR. M=1) */
  /* ----------------------------------------------------------------------- */
  zbknu(znr, zni, fnu, kode, nn, &cyr[1], &cyi[1], nz, tol, elim, alim);
  goto L110;
/* ----------------------------------------------------------------------- */
/*     LEFT HALF PLANE COMPUTATION */
/* ----------------------------------------------------------------------- */
L80:
  mr = -mm;
  zacon(znr, zni, fnu, kode, mr, nn, &cyr[1], &cyi[1], &nw, rl, fnul,
        tol, elim, alim);
  if (nw < 0) {
    goto L240;
  }
  *nz = nw;
  goto L110;
L90:
  /* ----------------------------------------------------------------------- */
  /*     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL */
  /* ----------------------------------------------------------------------- */
  mr = 0;
  if (znr >= 0. && (znr != 0. || zni >= 0. || m != 2)) {
    goto L100;
  }
  mr = -mm;
  if (znr != 0. || zni >= 0.) {
    goto L100;
  }
  znr = -znr;
  zni = -zni;
L100:
  zbunk(znr, zni, fnu, kode, mr, nn, &cyr[1], &cyi[1], &nw, tol, elim, alim);
  if (nw < 0) {
    goto L240;
  }
  *nz += nw;
L110:
  /* ----------------------------------------------------------------------- */
  /*     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT) */

  /*     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2 */
  /* ----------------------------------------------------------------------- */
  sgn = std::copysign(hpi, -fmm);
  /* ----------------------------------------------------------------------- */
  /*     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
  /*     WHEN FNU IS LARGE */
  /* ----------------------------------------------------------------------- */
  inu = (int)fnu;
  inuh = inu / 2;
  ir = inu - (inuh << 1);
  arg = (fnu - (inu - ir)) * sgn;
  rhpi = 1. / sgn;
  /*     ZNI = RHPI*COS(ARG) */
  /*     ZNR = -RHPI*SIN(ARG) */
  csgni = rhpi * std::cos(arg);
  csgnr = -rhpi * std::sin(arg);
  if (inuh % 2 == 0) {
    goto L120;
  }
  /*     ZNR = -ZNR */
  /*     ZNI = -ZNI */
  csgnr = -csgnr;
  csgni = -csgni;
L120:
  zti = -fmm;
  rtol = 1. / tol;
  ascle = ufl * rtol;
  for (i__ = 1; i__ <= nn; ++i__) {
    /*       STR = CYR(I)*ZNR - CYI(I)*ZNI */
    /*       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR */
    /*       CYR(I) = STR */
    /*       STR = -ZNI*ZTI */
    /*       ZNI = ZNR*ZTI */
    /*       ZNR = STR */
    aa = cyr[i__];
    bb = cyi[i__];
    atol = 1.;
    /* Computing MAX */
    if (std::max(std::fabs(aa), std::fabs(bb)) > ascle) {
      goto L135;
    }
    aa *= rtol;
    bb *= rtol;
    atol = tol;
  L135:
    str = aa * csgnr - bb * csgni;
    sti = aa * csgni + bb * csgnr;
    cyr[i__] = str * atol;
    cyi[i__] = sti * atol;
    str = -csgni * zti;
    csgni = csgnr * zti;
    csgnr = str;
    /* L130: */
  }
  return ierr;
L140:
  if (znr < 0.) {
    goto L230;
  }
  return ierr;
L230:
  *nz = 0;
  ierr = 2;
  return ierr;
L240:
  if (nw == -1) {
    goto L230;
  }
  *nz = 0;
  ierr = 5;
  return ierr;
L260:
  *nz = 0;
  ierr = 4;
  return ierr;
}

}  // namespace zbessel
