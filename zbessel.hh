#pragma once

namespace zbessel {

template <class = void>
int zbesh(double zr, double zi, double fnu, int kode, int m,
          int n, double *cyr, double *cyi, int *nz);

template <class = void>
int zbesi(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz);

template <class = void>
int zbesj(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz);

template <class = void>
int zbesk(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz);

template <class = void>
int zbesy(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz, double *cwrkr, double *cwrki);

template <class = void>
int zairy(double zr, double zi, int id, int kode, double *air, double *aii,
          int *nz);

template <class = void>
int zbiry(double zr, double zi, int id, int kode, double *bir, double *bii);

}  // namespace zbessel

#include "zbessel/zbesh.x"
#include "zbessel/zbesi.x"
#include "zbessel/zbesj.x"
#include "zbessel/zbesk.x"
#include "zbessel/zbesy.x"
#include "zbessel/zairy.x"
#include "zbessel/zbiry.x"
