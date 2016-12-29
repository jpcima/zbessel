# zbessel

This package provides C and C++ versions of complex Bessel functions, which were originally designed and implemented by D.E. Amos.
The original Fortran 77 implementation is available as part of the SLATEC mathematical library.

This is a port based on the f2c translation of the original source, improved for modern programming.

The characteristics of zbessel are:
- API for C and C++, provided by `zbessel.h`
- header-only API for C++, provided by `zbessel.hh`
- no dependency on f2c or gfortran runtimes
- modern calling conventions for C and C++
- thread-safety
- easy inclusion into projects

A test program provided as well; it computes the functions provided by this library with random inputs, and compares them with the original Fortran implementations in the SLATEC library.
