require 'numo/linalg'

Numo::Linalg::Blas.dlopen("libopenblas.so")
Numo::Linalg::Lapack.dlopen("libopenblas.so")
