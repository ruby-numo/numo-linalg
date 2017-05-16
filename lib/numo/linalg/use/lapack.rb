require 'fiddle'
require 'numo/linalg'

Fiddle.dlopen("libblas.so")
Fiddle.dlopen("liblapack.so")
Numo::Linalg::Blas.dlopen("libcblas.so")
Numo::Linalg::Lapack.dlopen("liblapacke.so")
