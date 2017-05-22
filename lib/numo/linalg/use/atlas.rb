require 'fiddle'
require 'numo/linalg'

Fiddle.dlopen("libatlas.so")
Numo::Linalg::Blas.dlopen("libatlas.so")
Numo::Linalg::Lapack.dlopen("liblapacke.so")
