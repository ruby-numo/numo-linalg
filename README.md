# Numo::Linalg : Linear Algebra library with LAPACK

## Implemented Methods

    x = Numo::Linalg.matmul(a,b)   (_gemm) Matrix multiply
    x = Numo::Linalg.solve(a,b)    (_gesv) Solve Linear equation with LU factorization
    x,y = Numo::Linalg.eigen(a,b)  (_geev) Eigen value and Eigen vector

* [GitHub](https://github.com/ruby-numo/linalg)

## Installation

* Install [Numo::NArray](https://github.com/ruby-numo/narray)
* Install [LAPACK](http://www.netlib.org/lapack/) or [OpenBLAS](http://www.openblas.net/)
  * Yum intall:
  ```shell
$ yum install openblas-devel
```

* Install Numo::Linalg
  ```shell
$ git clone git clone git://github.com/ruby-numo/linalg.git
$ cd linalg
$ rake build
$ gem install pkg/numo-linalg-*.gem
```

## ToDo

* rank
* norm
* det
* lstsq
* inv
* pinv
* cholesky
* qr
* svd
* eigh
* eigvals
* eigvalsh
