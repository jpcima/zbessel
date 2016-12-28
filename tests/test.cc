#include <zbessel.hh>
#include <random>
#include <complex>
#include <cstdio>
#include <ctime>
#include <cmath>

namespace slatec {
extern "C" {
  void zbesh_(double *zr, double *zi, double *fnu, int *kode, int *m,
              int *n, double *cyr, double *cyi, int *nz, int *ierr);
  void zbesi_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr,
              double *cyi, int *nz, int *ierr);
  void zbesj_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr,
              double *cyi, int *nz, int *ierr);
  void zbesk_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr,
              double *cyi, int *nz, int *ierr);
  void zbesy_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr,
              double *cyi, int *nz, double *cwrkr, double *cwrki, int *ierr);
  void zairy_(double *zr, double *zi, int *id, int *kode, double *air, double *aii,
              int *nz, int *ierr);
  void zbiry_(double *zr, double *zi, int *id, int *kode, double *bir, double *bii,
              int *ierr);
}
}  // namespace slatec

static std::default_random_engine rng;

static unsigned rndmod(unsigned n) {
  return (rng() - rng.min()) % n;
}

static double rndfloat() {
  double n = rng() - rng.min();
  double d = rng.max() - rng.min();
  return n / d;
}

double relerr(double val, double ref) {
  if (ref == 0) return 0;  // ignore
  return std::abs(val - ref) / std::abs(ref);
}

double abserr(double val, double ref) {
  return std::abs(val) - std::abs(ref);
}

std::complex<double> randz() {
  double range = 1000;
  return std::complex<double>(range * (rndfloat() * 2.0 - 1.0),
                              range * (rndfloat() * 2.0 - 1.0));
}

static int errorcheckgood = 0;
static int errorchecktotal = 0;
static int computationgood = 0;
static int computationtotal = 0;

void check_conformity(int ierr1, int ierr2, double err) {
  if (ierr1 != 0 || ierr2 != 0) {
    ++errorchecktotal;
    if (ierr1 == ierr2)
      ++errorcheckgood;
  }
  const double epsilon = 1e-10;
  ++computationtotal;
  if (err < epsilon)
    ++computationgood;
}

void test_zbesh() {
  constexpr int nmax = 16;

  double err = 0;

  std::complex<double> z = randz();
  double zr = z.real(), zi = z.imag();
  double fnu = rndfloat() * 1.0;
  int kode = 1 + rndmod(2);
  int m = 1 + rndmod(2);
  int n = 1 + rndmod(nmax);

  printf("ZBESH <- z (%g,%g) fnu %g kode %d m %d n %d\n", zr, zi, fnu, kode, m, n);

  double cyr1[nmax] {}, cyr2[nmax] {};
  double cyi1[nmax] {}, cyi2[nmax] {};
  int nz1 {}, nz2 {};
  int ierr1 {}, ierr2 {};

  ierr1 = zbessel::zbesh(zr, zi, fnu, kode, m, n, cyr1, cyi1, &nz1);
  slatec::zbesh_(&zr, &zi, &fnu, &kode, &m, &n, cyr2, cyi2, &nz2, &ierr2);

  for (int i = 0; i < n; ++i) {
    printf("ZBESH -> cy[%u] (%g,%g) (%g,%g)\n", i, cyr1[i], cyi1[i], cyr2[i], cyi2[i]);
    if (ierr1 == ierr2 && (ierr1 == 0 || ierr1 == 3)) {
      err = std::max(err, relerr(cyr1[i], cyr2[i]));
      err = std::max(err, relerr(cyi1[i], cyi2[i]));
    }
  }
  printf("ZBESH -> nz %d %d\n", nz1, nz2);
  printf("ZBESH -> ierr %d %d\n", ierr1, ierr2);

  printf("ZBESH relerr %g\n", err);

  check_conformity(ierr1, ierr2, err);
}

void test_zbesi() {
  constexpr int nmax = 16;

  double err = 0;

  std::complex<double> z = randz();
  double zr = z.real(), zi = z.imag();
  double fnu = rndfloat() * 1.0;
  int kode = 1 + rndmod(2);
  int n = 1 + rndmod(nmax);

  printf("ZBESI <- z (%g,%g) fnu %g kode %d n %d\n", zr, zi, fnu, kode, n);

  double cyr1[nmax] {}, cyr2[nmax] {};
  double cyi1[nmax] {}, cyi2[nmax] {};
  int nz1 {}, nz2 {};
  int ierr1 {}, ierr2 {};

  ierr1 = zbessel::zbesi(zr, zi, fnu, kode, n, cyr1, cyi1, &nz1);
  slatec::zbesi_(&zr, &zi, &fnu, &kode, &n, cyr2, cyi2, &nz2, &ierr2);

  for (int i = 0; i < n; ++i) {
    printf("ZBESI -> cy[%u] (%g,%g) (%g,%g)\n", i, cyr1[i], cyi1[i], cyr2[i], cyi2[i]);
    if (ierr1 == ierr2 && (ierr1 == 0 || ierr1 == 3)) {
      err = std::max(err, relerr(cyr1[i], cyr2[i]));
      err = std::max(err, relerr(cyi1[i], cyi2[i]));
    }
  }
  printf("ZBESI -> nz %d %d\n", nz1, nz2);
  printf("ZBESI -> ierr %d %d\n", ierr1, ierr2);

  printf("ZBESI relerr %g\n", err);

  check_conformity(ierr1, ierr2, err);
}

void test_zbesj() {
  constexpr int nmax = 16;

  double err = 0;

  std::complex<double> z = randz();
  double zr = z.real(), zi = z.imag();
  double fnu = rndfloat() * 1.0;
  int kode = 1 + rndmod(2);
  int n = 1 + rndmod(nmax);

  printf("ZBESJ <- z (%g,%g) fnu %g kode %d n %d\n", zr, zi, fnu, kode, n);

  double cyr1[nmax] {}, cyr2[nmax] {};
  double cyi1[nmax] {}, cyi2[nmax] {};
  int nz1 {}, nz2 {};
  int ierr1 {}, ierr2 {};

  ierr1 = zbessel::zbesj(zr, zi, fnu, kode, n, cyr1, cyi1, &nz1);
  slatec::zbesj_(&zr, &zi, &fnu, &kode, &n, cyr2, cyi2, &nz2, &ierr2);

  for (int i = 0; i < n; ++i) {
    printf("ZBESJ -> cy[%u] (%g,%g) (%g,%g)\n", i, cyr1[i], cyi1[i], cyr2[i], cyi2[i]);
    if (ierr1 == ierr2 && (ierr1 == 0 || ierr1 == 3)) {
      err = std::max(err, relerr(cyr1[i], cyr2[i]));
      err = std::max(err, relerr(cyi1[i], cyi2[i]));
    }
  }
  printf("ZBESJ -> nz %d %d\n", nz1, nz2);
  printf("ZBESJ -> ierr %d %d\n", ierr1, ierr2);

  printf("ZBESJ relerr %g\n", err);

  check_conformity(ierr1, ierr2, err);
}

void test_zbesk() {
  constexpr int nmax = 16;

  double err = 0;

  std::complex<double> z = randz();
  double zr = z.real(), zi = z.imag();
  double fnu = rndfloat() * 1.0;
  int kode = 1 + rndmod(2);
  int n = 1 + rndmod(nmax);

  printf("ZBESK <- z (%g,%g) fnu %g kode %d n %d\n", zr, zi, fnu, kode, n);

  double cyr1[nmax] {}, cyr2[nmax] {};
  double cyi1[nmax] {}, cyi2[nmax] {};
  int nz1 {}, nz2 {};
  int ierr1 {}, ierr2 {};

  ierr1 = zbessel::zbesk(zr, zi, fnu, kode, n, cyr1, cyi1, &nz1);
  slatec::zbesk_(&zr, &zi, &fnu, &kode, &n, cyr2, cyi2, &nz2, &ierr2);

  for (int i = 0; i < n; ++i) {
    printf("ZBESK -> cy[%u] (%g,%g) (%g,%g)\n", i, cyr1[i], cyi1[i], cyr2[i], cyi2[i]);
    if (ierr1 == ierr2 && (ierr1 == 0 || ierr1 == 3)) {
      err = std::max(err, relerr(cyr1[i], cyr2[i]));
      err = std::max(err, relerr(cyi1[i], cyi2[i]));
    }
  }
  printf("ZBESK -> nz %d %d\n", nz1, nz2);
  printf("ZBESK -> ierr %d %d\n", ierr1, ierr2);

  printf("ZBESK relerr %g\n", err);

  check_conformity(ierr1, ierr2, err);
}

void test_zbesy() {
  constexpr int nmax = 16;

  double err = 0;

  std::complex<double> z = randz();
  double zr = z.real(), zi = z.imag();
  double fnu = rndfloat() * 1.0;
  int kode = 1 + rndmod(2);
  int n = 1 + rndmod(nmax);

  printf("ZBESY <- z (%g,%g) fnu %g kode %d n %d\n", zr, zi, fnu, kode, n);

  double cyr1[nmax] {}, cyr2[nmax] {};
  double cyi1[nmax] {}, cyi2[nmax] {};
  int nz1 {}, nz2 {};
  int ierr1 {}, ierr2 {};

  double cwrkr[nmax], cwrki[nmax];

  ierr1 = zbessel::zbesy(zr, zi, fnu, kode, n, cyr1, cyi1, &nz1, cwrkr, cwrki);
  slatec::zbesy_(&zr, &zi, &fnu, &kode, &n, cyr2, cyi2, &nz2, cwrkr, cwrki, &ierr2);

  for (int i = 0; i < n; ++i) {
    printf("ZBESY -> cy[%u] (%g,%g) (%g,%g)\n", i, cyr1[i], cyi1[i], cyr2[i], cyi2[i]);
    if (ierr1 == ierr2 && (ierr1 == 0 || ierr1 == 3)) {
      err = std::max(err, relerr(cyr1[i], cyr2[i]));
      err = std::max(err, relerr(cyi1[i], cyi2[i]));
    }
  }
  printf("ZBESY -> nz %d %d\n", nz1, nz2);
  printf("ZBESY -> ierr %d %d\n", ierr1, ierr2);

  printf("ZBESY relerr %g\n", err);

  check_conformity(ierr1, ierr2, err);
}

void test_zairy() {
  double err = 0;

  std::complex<double> z = randz();
  double zr = z.real(), zi = z.imag();
  int id = rndmod(2);
  int kode = 1 + rndmod(2);

  printf("ZAIRY <- z (%g,%g) id %d kode %d\n", zr, zi, id, kode);

  double air1 {}, air2 {};
  double aii1 {}, aii2 {};
  int nz1 {}, nz2 {};
  int ierr1 {}, ierr2 {};

  ierr1 = zbessel::zairy(zr, zi, id, kode, &air1, &aii1, &nz1);
  slatec::zairy_(&zr, &zi, &id, &kode, &air2, &aii2, &nz2, &ierr2);

  printf("ZAIRY -> ai (%g,%g) (%g,%g)\n", air1, aii1, air2, aii2);
  if (ierr1 == ierr2 && (ierr1 == 0 || ierr1 == 3)) {
    err = std::max(err, relerr(air1, air2));
    err = std::max(err, relerr(aii1, aii2));
  }

  printf("ZAIRY -> nz %d %d\n", nz1, nz2);
  printf("ZAIRY -> ierr %d %d\n", ierr1, ierr2);

  printf("ZAIRY relerr %g\n", err);

  check_conformity(ierr1, ierr2, err);
}

void test_zbiry() {
  double err = 0;

  std::complex<double> z = randz();
  double zr = z.real(), zi = z.imag();
  int id = rndmod(2);
  int kode = 1 + rndmod(2);

  printf("ZBIRY <- z (%g,%g) id %d kode %d\n", zr, zi, id, kode);

  double air1 {}, air2 {};
  double aii1 {}, aii2 {};
  int ierr1 {}, ierr2 {};

  ierr1 = zbessel::zbiry(zr, zi, id, kode, &air1, &aii1);
  slatec::zbiry_(&zr, &zi, &id, &kode, &air2, &aii2, &ierr2);

  printf("ZBIRY -> ai (%g,%g) (%g,%g)\n", air1, aii1, air2, aii2);
  if (ierr1 == ierr2 && (ierr1 == 0 || ierr1 == 3)) {
    err = std::max(err, relerr(air1, air2));
    err = std::max(err, relerr(aii1, aii2));
  }

  printf("ZBIRY -> ierr %d %d\n", ierr1, ierr2);

  printf("ZBIRY relerr %g\n", err);

  check_conformity(ierr1, ierr2, err);
}

int main() {
  rng.seed(time(nullptr));

  const unsigned runcount = 1000;
  for (unsigned i = 0; i < runcount; ++i) {
    test_zbesh();
    printf("--------------------------------\n");
    test_zbesi();
    printf("--------------------------------\n");
    test_zbesj();
    printf("--------------------------------\n");
    test_zbesk();
    printf("--------------------------------\n");
    test_zbesy();
    printf("--------------------------------\n");
    test_zairy();
    printf("--------------------------------\n");
    test_zbiry();
    printf("--------------------------------\n");
  }

  printf("Run count %u\n", runcount);
  printf("Error checks passed %u/%u\n", errorcheckgood, errorchecktotal);
  printf("Computations passed %u/%u\n", computationgood, computationtotal);

  return 0;
}
