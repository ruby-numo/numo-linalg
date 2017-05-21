require 'numo/linalg'
require 'fiddle'
require 'benchmark'

def load_libs(use)
  lapack_dir="/opt/LAPACK/lib64/"
  atlas_dir="/opt/ATLAS/lib/"
  openblas_dir="/opt/OpenBLAS/lib/"

  case use
  when "blas-pkg","lapack-pkg"
    Fiddle.dlopen("/usr/lib64/libblas.so")
    Fiddle.dlopen("/usr/lib64/liblapack.so")
    Numo::Linalg::Blas.dlopen("/usr/lib64/libcblas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "blas-build","lapack-build"
    Fiddle.dlopen(lapack_dir+"libblas.so")
    Fiddle.dlopen(lapack_dir+"liblapack.so")
    Numo::Linalg::Blas.dlopen(lapack_dir+"libcblas.so")
    Numo::Linalg::Lapack.dlopen(lapack_dir+"liblapacke.so")

  when "satlas-pkg"
    Fiddle.dlopen("/usr/lib64/atlas/libsatlas.so")
    Numo::Linalg::Blas.dlopen("/usr/lib64/atlas/libsatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "satlas-build"
    Fiddle.dlopen(atlas_dir+"libsatlas.so")
    Numo::Linalg::Blas.dlopen(atlas_dir+"libsatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "tatlas-pkg"
    Fiddle.dlopen("/usr/lib64/atlas/libtatlas.so")
    Numo::Linalg::Blas.dlopen("/usr/lib64/atlas/libtatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")
  when "tatlas-build"
    Fiddle.dlopen(atlas_dir+"libtatlas.so")
    Numo::Linalg::Blas.dlopen(atlas_dir+"libtatlas.so")
    Numo::Linalg::Lapack.dlopen("liblapacke.so")

  when "openblas-pkg"
    Numo::Linalg::Blas.dlopen("/usr/lib64/libopenblas.so")
    Numo::Linalg::Lapack.dlopen("/usr/lib64/libopenblas.so")
  when "openblas-build"
    Numo::Linalg::Blas.dlopen(openblas_dir+"libopenblas.so")
    Numo::Linalg::Lapack.dlopen(openblas_dir+"libopenblas.so")
  end
end
