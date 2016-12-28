#pragma once

#if 1
#include "zops.x"

#else
#include <complex>

namespace zbessel {

template <class R>
inline R zabs(R r, R i) {
  return std::abs(std::complex<R>(r, i));
}

template <class R>
inline void zmlt(R ar, R ai, R br, R bi, R *cr, R *ci) {
  std::complex<R> a(ar, ai);
  std::complex<R> b(br, bi);
  std::complex<R> c = a * b;
  *cr = c.real();
  *ci = c.imag();
}

template <class R>
inline void zdiv(R ar, R ai, R br, R bi, R *cr, R *ci) {
  std::complex<R> a(ar, ai);
  std::complex<R> b(br, bi);
  std::complex<R> c = a / b;
  *cr = c.real();
  *ci = c.imag();
}

template <class R>
inline void zsqrt(R ar, R ai, R *br, R *bi) {
  std::complex<R> b = std::sqrt(std::complex<R>(ar, ai));
  *br = b.real();
  *bi = b.imag();
}

template <class R>
inline void zexp(R ar, R ai, R *br, R *bi) {
  std::complex<R> b = std::exp(std::complex<R>(ar, ai));
  *br = b.real();
  *bi = b.imag();
}

template <class R>
inline void zlog(R ar, R ai, R *br, R *bi) {
  std::complex<R> b = std::log(std::complex<R>(ar, ai));
  *br = b.real();
  *bi = b.imag();
}

}  // namespace zbessel

#endif
