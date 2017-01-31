# Numo::LAPACK : LAPACK binding with Numo::NArray

## Introduction

This is a binding of LAPACK with Numo::NArray .

## About 

* [GitHub](https://github.com/metanest/numo-lapack)

## About column-major order of LAPACK

Numo::LAPACK doesn't change of matrix row or column-major order.
NArray uses row-major, LAPACK uses column-major, user script must
convert the argument matrix and returned matrix.

(except for norm)

Or see lib/numo/linalg-compat.rb , this is a tiny sample compatible layer.

## Implemented Methods

* geev, heev, heevd - eigenvalues (for complex Hermitian matrix / using divide and conquer algorithm)
* gels - finding a least squares solution / minimum norm solution
* gemm - matrix multiplication
* gesv - solves linear systems
* gesvd, gesdd - singular value decomposition (using divide and conquer algorithm)
* geqrf - QR Factorization
* getrf, potrf - linear equations
* matrix_rank - matrix rank
* norm - various norm of a vector or matrix, with optional axes

## Installation

* Install [Numo::NArray](https://github.com/ruby-numo/narray)
* Install [LAPACK](http://www.netlib.org/lapack/) or [OpenBLAS](http://www.openblas.net/)
  * Yum intall:
  ```shell
$ yum install openblas-devel
```
  * or, install OpenBLAS from source and enjoy multithread acceleration.

* Install Numo::LAPACK
  ```shell
$ git clone git://github.com/metanest/numo-lapack.git
$ cd numo-lapack
$ rake build
$ gem install pkg/numo-lapack-*.gem -- --with-openblas --with-opt-dir=/opt/OpenBLAS
```

## ToDo

* Documentation complete
