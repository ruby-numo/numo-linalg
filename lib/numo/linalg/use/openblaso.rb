require 'numo/linalg'

Numo::Linalg::Blas.dlopen("libopenblaso.so")
Numo::Linalg::Lapack.dlopen("libopenblaso.so")
