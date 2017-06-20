echo -e "\n== DGEMM benchmark =="

n=1; echo -e "\nOMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgemm.rb 1 blas-build cap
OMP_NUM_THREADS=$n ruby dgemm.rb 1 satlas-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb 1 satlas-build
OMP_NUM_THREADS=$n ruby dgemm.rb 1 openblas-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb 1 openblas-build
OMP_NUM_THREADS=$n ruby dgemm.rb 1 intel_mkl-seq

n=4; echo -e "\nOMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgemm.rb 2 tatlas-pkg cap
OMP_NUM_THREADS=$n ruby dgemm.rb 2 tatlas-build
OMP_NUM_THREADS=$n ruby dgemm.rb 2 openblaso-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb 2 openblasp-pkg
OMP_NUM_THREADS=$n ruby dgemm.rb 2 openblas-build
OMP_NUM_THREADS=$n ruby dgemm.rb 3 intel_mkl-thread

echo -e "\n== DGELS benchmark =="

n=1; echo -e "\nOMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgels.rb 1 lapack-build cap
OMP_NUM_THREADS=$n ruby dgels.rb 1 satlas-pkg
OMP_NUM_THREADS=$n ruby dgels.rb 1 satlas-build
OMP_NUM_THREADS=$n ruby dgels.rb 1 openblas-pkg
OMP_NUM_THREADS=$n ruby dgels.rb 1 openblas-build
OMP_NUM_THREADS=$n ruby dgels.rb 1 intel_mkl-seq

n=4; echo -e "\nOMP_NUM_THREADS=$n"

OMP_NUM_THREADS=$n ruby dgels.rb 2 tatlas-pkg cap
OMP_NUM_THREADS=$n ruby dgels.rb 2 tatlas-build
OMP_NUM_THREADS=$n ruby dgels.rb 2 openblaso-pkg
OMP_NUM_THREADS=$n ruby dgels.rb 2 openblasp-pkg
OMP_NUM_THREADS=$n ruby dgels.rb 2 openblas-build
OMP_NUM_THREADS=$n ruby dgels.rb 3 intel_mkl-thread
