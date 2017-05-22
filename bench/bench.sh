echo "
== DGEMM benchmark =="

n=1; echo "
OMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgemm.rb blas-build cap
OMP_NUM_THREADS=$n ruby dgemm.rb satlas-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb satlas-build
OMP_NUM_THREADS=$n ruby dgemm.rb openblas-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb openblas-build

n=4; echo "
OMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgemm.rb tatlas-pkg cap
OMP_NUM_THREADS=$n ruby dgemm.rb tatlas-build
OMP_NUM_THREADS=$n ruby dgemm.rb openblaso-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb openblasp-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb openblas-build

echo "
== DGELS benchmark =="

n=1; echo "
OMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgels.rb lapack-build cap
OMP_NUM_THREADS=$n ruby dgels.rb satlas-pkg
OMP_NUM_THREADS=$n ruby dgels.rb satlas-build
OMP_NUM_THREADS=$n ruby dgels.rb openblas-pkg
OMP_NUM_THREADS=$n ruby dgels.rb openblas-build

n=4; echo "
OMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgels.rb tatlas-pkg cap
OMP_NUM_THREADS=$n ruby dgels.rb tatlas-build
OMP_NUM_THREADS=$n ruby dgels.rb openblaso-pkg
OMP_NUM_THREADS=$n ruby dgels.rb openblasp-pkg
OMP_NUM_THREADS=$n ruby dgels.rb openblas-build
