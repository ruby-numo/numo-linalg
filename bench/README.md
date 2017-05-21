# Numo::Linalg benchmark

* Environment:
    * VirtualBox 5.1.22 r115126
    * CentOS Linux release 7.3.1611 (Core)
    * gcc 4.8.5 20150623 (Red Hat 4.8.5-11) (GCC)
    * Intel(R) Core(TM) i7-3612QM CPU @ 2.10GHz

* -pkg : yum-installed package
* -build : build & installed from source code

### Result

* Check and compare *real* time.

```
$ sh bench.sh

== DGEMM benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
blas-build         9.380000   0.020000   9.400000 (  9.397385)
satlas-pkg         2.120000   0.010000   2.130000 (  2.128070)
satlas-build       1.250000   0.010000   1.260000 (  1.256298)
openblas-pkg       1.060000   0.010000   1.070000 (  1.070094)
openblas-build     1.060000   0.010000   1.070000 (  1.070562)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         2.360000   0.080000   2.440000 (  0.644536)
tatlas-build       1.510000   0.070000   1.580000 (  0.422579)
openblas-pkg       1.050000   0.020000   1.070000 (  1.063671)
openblas-build     1.220000   0.130000   1.350000 (  0.337856)

== DGELS benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
lapack-build      10.690000   0.010000  10.700000 ( 10.707058)
satlas-pkg         1.930000   0.030000   1.960000 (  1.960371)
satlas-build       1.390000   0.020000   1.410000 (  1.401828)
openblas-pkg       1.310000   0.020000   1.330000 (  1.338973)
openblas-build     1.340000   0.020000   1.360000 (  1.363780)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         2.130000   0.430000   2.560000 (  1.011448)
tatlas-build       1.610000   0.260000   1.870000 (  0.889226)
openblas-pkg       1.350000   0.020000   1.370000 (  1.370657)
openblas-build     1.930000   1.430000   3.360000 (  0.842967)
```
