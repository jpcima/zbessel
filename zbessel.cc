#include "zbessel.h"
#include "zbessel.hh"

[[gnu::visibility("default")]]
int zbesh(double zr, double zi, double fnu, int kode, int m,
          int n, double *cyr, double *cyi, int *nz) {
  return zbessel::zbesh(zr, zi, fnu, kode, m, n, cyr, cyi, nz);
}

[[gnu::visibility("default")]]
int zbesi(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz) {
  return zbessel::zbesi(zr, zi, fnu, kode, n, cyr, cyi, nz);
}

[[gnu::visibility("default")]]
int zbesj(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz) {
  return zbessel::zbesj(zr, zi, fnu, kode, n, cyr, cyi, nz);
}

[[gnu::visibility("default")]]
int zbesk(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz) {
  return zbessel::zbesk(zr, zi, fnu, kode, n, cyr, cyi, nz);
}

[[gnu::visibility("default")]]
int zbesy(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz, double *cwrkr, double *cwrki) {
  return zbessel::zbesy(zr, zi, fnu, kode, n, cyr, cyi, nz, cwrkr, cwrki);
}

[[gnu::visibility("default")]]
int zairy(double zr, double zi, int id, int kode, double *air, double *aii,
          int *nz) {
  return zbessel::zairy(zr, zi, id, kode, air, aii, nz);
}

[[gnu::visibility("default")]]
int zbiry(double zr, double zi, int id, int kode, double *bir, double *bii) {
  return zbessel::zbiry(zr, zi, id, kode, bir, bii);
}
