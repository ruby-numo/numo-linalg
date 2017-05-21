require 'numo/linalg'
require 'fiddle'
require 'benchmark'

def load_libs(use)
  case use
  when "blas-pkg"
    Fiddle.dlopen("/usr/lib64/libblas.so")
    Fiddle.dlopen("/usr/lib64/liblapack.so")
    Numo::Linalg::Blas.dlopen("/usr/lib64/libcblas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "blas-build"
    Fiddle.dlopen("/opt/LAPACK/lib64/libblas.so")
    Fiddle.dlopen("/opt/LAPACK/lib64/liblapack.so")
    Numo::Linalg::Blas.dlopen("/opt/LAPACK/lib64/libcblas.so")
    Numo::Linalg::Lapack.dlopen("/opt/LAPACK/lib64/liblapacke.so")

  when "satlas-pkg"
    Fiddle.dlopen("/usr/lib64/atlas/libsatlas.so")
    Numo::Linalg::Blas.dlopen("/usr/lib64/atlas/libsatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "satlas-build"
    Fiddle.dlopen("/opt/ATLAS/lib/libsatlas.so")
    Numo::Linalg::Blas.dlopen("/opt/ATLAS/lib/libsatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "tatlas-pkg"
    Fiddle.dlopen("/usr/lib64/atlas/libtatlas.so")
    Numo::Linalg::Blas.dlopen("/usr/lib64/atlas/libtatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "tatlas-build"
    Fiddle.dlopen("/opt/ATLAS/lib/libtatlas.so")
    Numo::Linalg::Blas.dlopen("/opt/ATLAS/lib/libtatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")

  when "openblas-pkg"
    Numo::Linalg::Blas.dlopen("/usr/lib64/libopenblas.so")
    Numo::Linalg::Lapack.dlopen("/usr/lib64/libopenblas.so")
  when "openblas-build"
    Numo::Linalg::Blas.dlopen("/opt/OpenBLAS/lib/libopenblas.so")
    Numo::Linalg::Lapack.dlopen("/opt/OpenBLAS/lib/libopenblas.so")
  end
end
