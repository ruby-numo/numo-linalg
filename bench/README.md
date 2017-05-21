# Numo::Linalg benchmark

* Environment:
    * VirtualBox 5.1.22 r115126
    * CentOS Linux release 7.3.1611 (Core)
    * gcc 4.8.5 20150623 (Red Hat 4.8.5-11) (GCC)
    * Intel(R) Core(TM) i7-3612QM CPU @ 2.10GHz

* -pkg : yum-installed package
* -build : build & installed from source code

### Result

```
$ sh bench.sh

== DGEMM benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
blas-build         9.530000   0.010000   9.540000 (  9.543595)
                       user     system      total        real
satlas-pkg         2.160000   0.020000   2.180000 (  2.172579)
                       user     system      total        real
satlas-build       1.250000   0.010000   1.260000 (  1.261787)
                       user     system      total        real
openblas-pkg       1.050000   0.040000   1.090000 (  1.092658)
                       user     system      total        real
openblas-build     1.060000   0.020000   1.080000 (  1.088277)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         2.310000   0.060000   2.370000 (  0.643968)
                       user     system      total        real
tatlas-build       1.430000   0.040000   1.470000 (  0.405767)
                       user     system      total        real
openblas-pkg       1.060000   0.030000   1.090000 (  1.083399)
                       user     system      total        real
openblas-build     1.210000   0.190000   1.400000 (  0.359571)

== DGELS benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
blas-build        14.340000   0.040000  14.380000 ( 14.390546)
                       user     system      total        real
satlas-pkg         2.660000   0.060000   2.720000 (  2.714880)
                       user     system      total        real
satlas-build       1.920000   0.040000   1.960000 (  1.962759)
                       user     system      total        real
openblas-pkg       1.840000   0.040000   1.880000 (  1.867554)
                       user     system      total        real
openblas-build     1.850000   0.050000   1.900000 (  1.903175)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         3.140000   0.520000   3.660000 (  1.476800)
                       user     system      total        real
tatlas-build       2.300000   0.340000   2.640000 (  1.386921)
                       user     system      total        real
openblas-pkg       1.840000   0.040000   1.880000 (  1.886264)
                       user     system      total        real
openblas-build     2.620000   2.230000   4.850000 (  1.274298)
```
