#pragma once
#include <limits>
#include <cmath>

namespace zbessel {

template <class R>
R zabs(R r, R i) {
  R u = std::fabs(r);
  R v = std::fabs(i);
  R s = u + v;
  // s * 1.0 makes an unnormalized underflow on CDC machines into a
  // true floating zero
  s *= R(1);
  if (s == 0)
    return 0;
  if (u > v) {
    R q = v / u;
    return u * std::sqrt(q * q + R(1));
  }
  R q = u / v;
  return v * std::sqrt(q * q + R(1));
}

template <class R>
void zmlt(R ar, R ai, R br, R bi, R *cr, R *ci) {
  *cr = ar * br - ai * bi;
  *ci = ar * bi + ai * br;
}

template <class R>
void zdiv(R ar, R ai, R br, R bi, R *cr, R *ci) {
  R bm = 1 / zabs(br, bi);
  R cc = br * bm;
  R cd = bi * bm;
  *cr = (ar * cc + ai * cd) * bm;
  *ci = (ai * cc - ar * cd) * bm;
}

template <class R>
void zsqrt(R ar, R ai, R *br, R *bi) {
  const R drt = 0.707106781186547461715L;  // 1/sqrt(2)
  const R dpi = 3.14159265358979323846L;

  R zm;
  R dtheta;

  zm = std::sqrt(zabs(ar, ai));
  if (ar == 0)
    goto L10;
  if (ai == 0)
    goto L20;
  dtheta = std::atan(ai / ar);
  if (dtheta <= 0)
    goto L40;
  if (ar < 0)
    dtheta -= dpi;
  goto L50;
L10:
  if (ai > 0)
    goto L60;
  if (ai < 0)
    goto L70;
  *br = 0;
  *bi = 0;
  return;
L20:
  if (ar > 0)
    goto L30;
  *br = 0;
  *bi = std::sqrt(std::abs(ar));
  return;
L30:
  *br = std::sqrt(ar);
  *bi = 0;
  return;
L40:
  if (ar < 0)
    dtheta += dpi;
L50:
  dtheta *= R(.5);
  *br = zm * std::cos(dtheta);
  *bi = zm * std::sin(dtheta);
  return;
L60:
  *br = zm * drt;
  *bi = zm * drt;
  return;
L70:
  *br = zm * drt;
  *bi = -zm * drt;
}

template <class R>
void zexp(R ar, R ai, R *br, R *bi) {
  R zm = std::exp(ar);
  *br = zm * std::cos(ai);
  *bi = zm * std::sin(ai);
}

template <class R>
void zlog(R ar, R ai, R *br, R *bi) {
  const R dpi = 3.14159265358979323846L;
  const R dhpi = 1.570796326794896558L;  // pi/2

  R zm;
  R dtheta;

  if (ar == 0)
    goto L10;
  if (ai == 0)
    goto L20;
  dtheta = std::atan(ai / ar);
  if (dtheta <= 0)
    goto L40;
  if (ar < 0)
    dtheta -= dpi;
  goto L50;
L10:
  if (ai == 0)
    goto L60;
  *bi = dhpi;
  *br = std::log(std::abs(ai));
  if (ai < 0)
    *bi = -(*bi);
  return;
L20:
  if (ar > 0)
    goto L30;
  *br = std::log(std::abs(ar));
  *bi = dpi;
  return;
L30:
  *br = std::log(ar);
  *bi = 0;
  return;
L40:
  if (ar < 0)
    dtheta += dpi;
L50:
  zm = zabs(ar, ai);
  *br = std::log(zm);
  *bi = dtheta;
  return;
L60:
  *br = -std::numeric_limits<R>::infinity();
  *bi = 0;
}

}  // namespace zbessel
