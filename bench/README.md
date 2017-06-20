# Numo::Linalg benchmark

* Environment:
    * VirtualBox 5.1.22 r115126
    * CentOS Linux release 7.3.1611 (Core)
    * gcc 4.8.5 20150623 (Red Hat 4.8.5-11) (GCC)
    * Intel(R) Core(TM) i7-3612QM CPU @ 2.10GHz

* -pkg : yum-installed package
* -build : built & installed from source code

    * atlas-pkg : atlas-3.10.1-10.el7.x86_64
    * atlas-build : atlas3.10.3
    * openblas-pkg : openblas-0.2.19-4.el7.x86_64
    * openblas-build : openblas-0.2.19
    * intel\_mkl : 2017.3.196

### Result

* Check and compare *real* time.

```
$ sh bench.sh

== DGEMM benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
blas-build        31.690000   0.020000  31.710000 ( 31.717582)
satlas-pkg         7.060000   0.010000   7.070000 (  7.087082)
satlas-build       4.010000   0.020000   4.030000 (  4.036639)
openblas-pkg       3.500000   0.030000   3.530000 (  3.525150)
openblas-build     3.500000   0.020000   3.520000 (  3.523451)
intel_mkl-seq      3.620000   0.010000   3.630000 (  3.638581)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         7.970000   0.260000   8.230000 (  2.136432)
tatlas-pkg         7.950000   0.200000   8.150000 (  2.111306)
tatlas-build       4.670000   0.110000   4.780000 (  1.257292)
tatlas-build       4.600000   0.110000   4.710000 (  1.232838)
openblaso-pkg      4.110000   0.240000   4.350000 (  1.117955)
openblaso-pkg      4.090000   0.140000   4.230000 (  1.076477)
openblasp-pkg      4.300000   0.310000   4.610000 (  1.161065)
openblasp-pkg      4.510000   0.420000   4.930000 (  1.245599)
openblas-build     4.090000   0.350000   4.440000 (  1.121971)
openblas-build     4.130000   0.360000   4.490000 (  1.129311)
intel_mkl-thread   4.390000   0.100000   4.490000 (  1.194790)
intel_mkl-thread   4.250000   0.080000   4.330000 (  1.093744)
intel_mkl-thread   4.200000   0.070000   4.270000 (  1.073778)

== DGELS benchmark ==

OMP_NUM_THREADS=1
                       user     system      total        real
lapack-build      36.310000   0.040000  36.350000 ( 36.344356)
satlas-pkg         6.530000   0.050000   6.580000 (  6.574038)
satlas-build       4.600000   0.040000   4.640000 (  4.640768)
openblas-pkg       4.310000   0.020000   4.330000 (  4.342447)
openblas-build     4.420000   0.050000   4.470000 (  4.471163)
intel_mkl-seq      3.310000   0.030000   3.340000 (  3.333907)

OMP_NUM_THREADS=4
                       user     system      total        real
tatlas-pkg         7.360000   1.220000   8.580000 (  3.373866)
tatlas-pkg         7.360000   1.180000   8.540000 (  3.334231)
tatlas-build       5.540000   0.760000   6.300000 (  2.969774)
tatlas-build       5.570000   0.780000   6.350000 (  3.067050)
openblaso-pkg      7.260000   0.350000   7.610000 (  2.460789)
openblaso-pkg      7.480000   0.400000   7.880000 (  2.513503)
openblasp-pkg      5.950000   3.790000   9.740000 (  2.445235)
openblasp-pkg      6.000000   3.770000   9.770000 (  2.457592)
openblas-build     6.190000   4.360000  10.550000 (  2.646929)
openblas-build     6.120000   4.010000  10.130000 (  2.543951)
intel_mkl-thread   5.230000   0.150000   5.380000 (  1.374205)
intel_mkl-thread   4.610000   0.100000   4.710000 (  1.184957)
intel_mkl-thread   4.520000   0.090000   4.610000 (  1.161126)
```

## Conclusion

Recommendation: OpenBLAS (pkg, built) or ATLAS (built).
