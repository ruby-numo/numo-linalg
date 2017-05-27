# Selecting Backend Library

* Numo::Linalg accesses BLAS/LAPACK library throuh dynamic linking
  loader at runtime.  This design allows you to change backend
  libraries without re-compiling Numo::Linalg.

* Note that the performance depends on the backend libraries as shown in
  [benchmark](https://github.com/ruby-numo/linalg/tree/master/bench).

## Install time option to specify default library path

    gem install pkg/numo-linalg-0.0.1.gem -- \
      --with-backend=mkl \
      --with-mkl-lib=/opt/intel/composerxe/mkl/lib/intel64 \
      --with-openblas-dir=/opt/OpenBLAS \
      --with-atlas-dir=/opt/ATLAS \
      --with-lapack-lib=/opt/LAPACK/lib64

* Note that these options are for setting default path to load library
  at runtime, not required for building Numo::Linalg library.

## Use Numo::Linalg

### Use default backend

```ruby
require "numo/linalg"
```

* if --with-backend option is specified, use it.
* if --with-backend option is not specified, try to load mkl,openblas,atlas,lapack.

### Specify backend

* 1st level (specify backend library)

```ruby
require "numo/linalg/use/openblas"
```

* 2nd level (specify library path)

```ruby
require "numo/linalg/linalg"
Numo::Linalg::Loader.load_openblas "/opt/OpenBLAS/lib"
```

* 3rd level (specify library files to load)

```ruby
require "numo/linalg/linalg"
Numo::Linalg::Blas.dlopen("/opt/OpenBLAS/lib/libopenblas.so")
Numo::Linalg::Lapack.dlopen("/opt/OpenBLAS/lib/libopenblas.so")
```

* Numo::Linalg requires C-interface CBLAS and LAPACKE interface.
  If they are provided as separate libraries,
  you need to open all of them.

```ruby
require "numo/linalg/linalg"
require "fiddle"
Fiddle.dlopen("libblas.so")  # referenced from CBLAS
Fiddle.dlopen("liblapack.so")  # referenced from LAPACKE
Numo::Linalg::Blas.dlopen("libcblas.so")  # contains cblas_* function
Numo::Linalg::Lapack.dlopen("liblapacke.so")  # containes LAPACKE_* function
```

* Reopen (change backend at runtime) is not supported.

## Show libraries

```ruby
p Numo::Linalg::Loader.libs
```
