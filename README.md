# Numo::Linalg : Linear Algebra library with BLAS/LAPACK binding to Numo::NArray

[GitHub](https://github.com/ruby-numo/linalg)

Under development!

## Introduction

This is a binding of LAPACK with Numo::NArray .

### [Numo::Linalg API](http://ruby-numo.github.io/linalg/yard/Numo/Linalg.html)

* Matrix and vector products
    * dot, matmul
* Decomposition
    * cholesky, qr, svd
* Matrix eigenvalues
    * eig, eigh, eigvals, eigvalsh
* Norms and other numbers
    * norm, det, matrix_rank, slogdet
* Solving equations and inverting matrices
    * solve, lstsq, inv

More functions to come

### Low-level modules

* [Numo::Linalg::Blas](http://ruby-numo.github.io/linalg/yard/Numo/Linalg/Blas.html) - Low-level BLAS functions
* [Numo::Linalg::Lapack](http://ruby-numo.github.io/linalg/yard/Numo/Linalg/Lapack.html) - Low-level LAPACK functions

## Installation

* Install [Numo::NArray](https://github.com/ruby-numo/narray)

* Install [LAPACK](http://www.netlib.org/lapack/) or compatible packages.

    * Numo::Linalg requires C-interface
      [CBLAS](http://www.netlib.org/blas/#_cblas) and
      [LAPACKE](http://www.netlib.org/lapack/lapacke.html) interface.
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

## Loading backend library

* Numo::Linalg opens dynamic libraries of BLAS and LAPACK at runtime.
  This design allows you to change backend libraries without
  re-compiling Numo::Linalg.

* Note that the performance depends on the backend libraries as shown in
  [benchmark](https://github.com/ruby-numo/linalg/tree/master/bench).

* Dynamic loading example:

```ruby
require "numo/linalg"
require "fiddle"
Fiddle.dlopen("libblas.so")  # referenced from CBLAS
Fiddle.dlopen("liblapack.so")  # referenced from LAPACKE
Numo::Linalg::Blas.dlopen("libcblas.so")  # contains cblas_* function
Numo::Linalg::Lapack.dlopen("liblapacke.so")  # containes LAPACKE_* function
```

This default behavior is defined in "numo/linalg/use/lapack.rb"

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
* This work is partly supported by 2016 Ruby Association Grant.

## ToDo

* More functions
* write test
* Documentation
