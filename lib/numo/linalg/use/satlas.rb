require 'fiddle'
require 'numo/linalg'

Fiddle.dlopen("libsatlas.so")
Numo::Linalg::Blas.dlopen("libsatlas.so")
Numo::Linalg::Lapack.dlopen("liblapacke.so")
