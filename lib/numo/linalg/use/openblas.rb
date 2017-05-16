require 'numo/linalg'

Numo::Linalg::Blae.dlopen("libopenblas.so")
Numo::Linalg::Lapack.dlopen("libopenblas.so")
