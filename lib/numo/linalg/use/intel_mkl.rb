require 'fiddle'
require 'numo/linalg'

root = "/opt/intel/composerxe/"
Fiddle.dlopen(root+"lib/intel64/libiomp5.so")

dir = root+"mkl/lib/intel64/"
Fiddle.dlopen(dir+"libmkl_core.so")
Fiddle.dlopen(dir+"libmkl_intel_thread.so")
#Fiddle.dlopen(dir+"libmkl_gnu_thread.so")
#Fiddle.dlopen(dir+"libmkl_sequential.so")

# LP64: 32bit integer
Numo::Linalg::Blas.dlopen(dir+"libmkl_intel_lp64.so")
Numo::Linalg::Lapack.dlopen(dir+"libmkl_intel_lp64.so")
