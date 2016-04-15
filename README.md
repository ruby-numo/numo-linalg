# NArray::Linalg : Linear Algebra librart with LAPACK

## Implemented Methods

    x = NArray::Linalg.matmul(a,b)   (_gemm) Matrix multiply
    x = NArray::Linalg.solve(a,b)    (_gesv) Solve Linear equation usin LU
    x,y = NArray::Linalg.eigen(a,b)  (_geev) Eigen value and Eigen vector

* [GitHub](https://github.com/masa16/numo-linalg)

## Installation

* Install [Numo::NArray](https://github.com/masa16/numo-narray)
* Install [LAPACK](http://www.netlib.org/lapack/) or [OpenBlas](http://www.openblas.net/)
  * Yum intall:
  ```shell
$ yum install openblas-devel
```

* Install Numo::Linalg
  ```shell
$ git clone git clone git://github.com/masa16/numo-linalg.git
$ cd numo-linalg
$ rake build
$ gem install pkg/numo-linalg-*.gem
```
