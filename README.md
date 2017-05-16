# Numo::Linalg : Linear Algebra library with BLAS/LAPACK binding to Numo::NArray

[GitHub](https://github.com/ruby-numo/linalg)

## Introduction

This is a binding of LAPACK with Numo::NArray .

### Numo::Linalg

Matrix and vector products
* dot, matmul
Decomposition
* cholesky, qr, svd,
Matrix eigenvalues
* eig
* eigh
* eigvals
* eigvalsh
Norms and other numbers
* norm
* det
* matrix_rank
* slogdet
Solving equations and inverting matrices
* solve
* lstsq
* inv

### Low-level modules

* Numo::Linalg.Blas - BLAS functions
* Numo::Linalg.Lapack - LAPACK functions

## Installation

* Install [Numo::NArray](https://github.com/ruby-numo/narray)

* Install [LAPACK](http://www.netlib.org/lapack/) or compatible packages.
  * Numo::Linalg requires C-interface
    [CBLAS](http://www.netlib.org/blas/#_cblas) and
    [Lapacke](http://www.netlib.org/lapack/lapacke.html) interface.
    These are included in LAPACK package.
  * It is recommended to use one of following faster libraries:
    * [ATLAS](https://sourceforge.net/projects/math-atlas/)
    * [OpenBLAS](http://www.openblas.net/)
    * [Intel MKL](https://software.intel.com/intel-mkl)

* Install Numo::Linalg

```shell
$ git clone https://github.com/ruby-numo/linalg.git
$ cd linalg
$ rake build
$ gem install pkg/numo-linalg-*.gem
```

## Usage

* Loading backend library

  * Numo::Linalg loads BLAS,LAPACK libraries dynamically at runtime.
    This design allows you to change backend libraries without re-compiling
    Numo::Linalg.

  * Example of dynamical loading:

```ruby
require "numo/linalg"
require "fiddle"
Fiddle.dlopen("libblas.so")
Fiddle.dlopen("liblapack.so")
Numo::Linalg::Blas.dlopen("libcblas.so")
Numo::Linalg::Lapack.dlopen("liblapacke.so")
```

  * This default behavior is defined in "numo/linalg/use/lapack.rb"

```ruby
require "numo/linalg/use/lapack"
```

  * Load Atlas:

```ruby
require "numo/linalg/use/satlas"
```

  * Load OpenBLAS:

```ruby
require "numo/linalg/use/openblas"
```

## Authors

* Masahiro TANAKA
* Makoto KISHIMOTO
* This work is partially supported by 2016 Ruby Association Grant.

## ToDo

* wrap more functions
* write test
* Documentation
