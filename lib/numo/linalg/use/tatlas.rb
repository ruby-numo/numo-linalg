require 'fiddle'
require 'numo/linalg'

Fiddle.dlopen("libtatlas.so")
Numo::Linalg::Blas.dlopen("libtatlas.so")
Numo::Linalg::Lapack.dlopen("liblapacke.so")
