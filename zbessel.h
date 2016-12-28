#pragma once

#if defined(__cplusplus)
extern "C" {
#endif

int zbesh(double zr, double zi, double fnu, int kode, int m,
          int n, double *cyr, double *cyi, int *nz);

int zbesi(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz);

int zbesj(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz);

int zbesk(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz);

int zbesy(double zr, double zi, double fnu, int kode, int n, double *cyr,
          double *cyi, int *nz, double *cwrkr, double *cwrki);

int zairy(double zr, double zi, int id, int kode, double *air, double *aii,
          int *nz);

int zbiry(double zr, double zi, int id, int kode, double *bir, double *bii);

#if defined(__cplusplus)
}  // extern "C"
#endif
