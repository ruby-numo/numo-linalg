require 'numo/linalg'

Numo::Linalg::Blas.dlopen("libopenblasp.so")
Numo::Linalg::Lapack.dlopen("libopenblasp.so")
