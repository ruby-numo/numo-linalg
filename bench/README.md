# Numo::Linalg benchmark

* Environment:
    * VirtualBox 5.1.22 r115126
    * CentOS Linux release 7.3.1611 (Core)
    * gcc 4.8.5 20150623 (Red Hat 4.8.5-11) (GCC)
    * Intel(R) Core(TM) i7-3612QM CPU @ 2.10GHz

* -pkg : yum-installed package
* -build : built & installed from source code

### Result

* Check and compare *real* time.

```
$ sh bench.sh

== DGEMM benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
blas-build         9.450000   0.010000   9.460000 (  9.467362)
satlas-pkg         2.150000   0.010000   2.160000 (  2.168709)
satlas-build       1.240000   0.010000   1.250000 (  1.252929)
openblas-pkg       1.080000   0.010000   1.090000 (  1.089870)
openblas-build     1.070000   0.000000   1.070000 (  1.073924)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         2.350000   0.090000   2.440000 (  0.656182)
tatlas-build       1.420000   0.040000   1.460000 (  0.395914)
openblaso-pkg      1.270000   0.090000   1.360000 (  0.358807)
openblasp-pkg      1.200000   0.130000   1.330000 (  0.334057)
openblas-build     1.180000   0.150000   1.330000 (  0.338615)

== DGELS benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
lapack-build      10.740000   0.020000  10.760000 ( 10.757387)
satlas-pkg         1.940000   0.020000   1.960000 (  1.948135)
satlas-build       1.430000   0.020000   1.450000 (  1.444768)
openblas-pkg       1.320000   0.020000   1.340000 (  1.345208)
openblas-build     1.390000   0.020000   1.410000 (  1.401886)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         2.250000   0.360000   2.610000 (  1.051549)
tatlas-build       1.620000   0.240000   1.860000 (  0.888229)
openblaso-pkg      2.280000   0.130000   2.410000 (  0.778036)
openblasp-pkg      1.770000   1.210000   2.980000 (  0.750317)
openblas-build     1.850000   1.460000   3.310000 (  0.867976)
```

## Conclusion

Recommendation: OpenBLAS (pkg, built) or ATLAS (built).
