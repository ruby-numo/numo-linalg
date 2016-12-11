module Numo::Linalg
    class << self

        # module methods

        ## Matrix and vector products

        # not linalg
        def dot a, b
            a.dot b
        end

        # not linalg
        def vdot a, b
            raise NotImplementedError.new
        end

        # not linalg
        def inner a, b
            a.mulsum b
        end

        # not linalg
        def outer a, b
            a[false, :new] * b
        end

        def matmul a, b
            result = Numo::LAPACK.gemm a.transpose, b.transpose
            result.transpose
        end

        # not linalg
        def tensordot a, b #, axes
            raise NotImplementedError.new
        end

        # not linalg
        def einsum *a
            raise NotImplementedError.new
        end

        # not linalg
        def matrix_power m, n
            raise NotImplementedError.new
        end

        # not linalg
        def kron a, b
            raise NotImplementedError.new
        end


        ## Decompositions

        def cholesky a, upper:false
            uflag = upper
            Numo::LAPACK.potrf(a.transpose, upper:uflag).transpose
        end

        def qr a #, mode
            result = Numo::LAPACK.geqrf a.transpose
            q, r = result
            result[0] = q.transpose
            result[1] = r.transpose
            result
        end

        def svd a, turbo:false
            a = a.transpose
            result =
                if turbo then
                    Numo::LAPACK.gesdd a, job: :FULL
                else
                    Numo::LAPACK.gesvd a, job: :FULL
                end
            u, _, v = result
            result[0] = u.transpose
            result[2] = v.transpose
            result
        end

        def svdvals a, turbo:false
            a = a.transpose
            if turbo then
                Numo::LAPACK.gesdd a, job: :VALS_ONLY
            else
                Numo::LAPACK.gesvd a, job: :VALS_ONLY
            end
        end


        ## Matrix eigenvalues

        def eig a
            result = Numo::LAPACK.geev a.transpose, vals_only:false
            w, v = result
            result[1] = v.transpose
            result
        end

        def eigh a, upper:false, turbo:false
            uflag = upper
            a = a.transpose
            result =
                if turbo then
                    Numo::LAPACK.heevd a, vals_only:false, upper:uflag
                else
                    Numo::LAPACK.heev a, vals_only:false, upper:uflag
                end
            w, v = result
            result[1] = v.transpose
            result
        end

        def eigvals a
            Numo::LAPACK.geev a.transpose, vals_only:true
        end

        def eigvalsh a, upper:false, turbo:false
            uflag = upper
            a = a.transpose
            if turbo then
                Numo::LAPACK.heevd a, vals_only:true, upper:uflag
            else
                Numo::LAPACK.heev a, vals_only:true, upper:uflag
            end
        end


        ## Norms and other numbers

        def norm *x
            raise NotImplementedError.new
        end

        def cond *x
            raise NotImplementedError.new
        end

        def det *x
            raise NotImplementedError.new
        end

        def matrix_rank m, tol=nil, turbo:false
            tflag = turbo
            m = m.transpose
            if tol then
                Numo::LAPACK.matrix_rank m, tol, turbo:tflag
            else
                Numo::LAPACK.matrix_rank m, turbo:tflag
            end
        end

        def slogdet *a
            raise NotImplementedError.new
        end

        # not linalg
        def trace *a
            raise NotImplementedError.new
        end


        ## Solving equations and inverting matrices

        def solve a, b
            Numo::LAPACK.gesv(a.transpose, b).transpose
        end

        def tensorsolve a, b, *_
            raise NotImplementedError.new
        end

        def lstsq a, b
            Numo::LAPACK.gels a.transpose, b
        end

        def inv a
            sh_a = a.shape
            raise "ndim(a) != 2" if sh_a.size != 2
            x, y = sh_a
            raise "argument matrix is not square" if x != y
            b = a.class.new(y, x).eye
            Numo::LAPACK.gesv(a.transpose, b).transpose
        end

        def pinv a, turbo:false
            tflag = turbo
            a = a.conj.transpose
            args = a, {job: :THIN}
            if tflag then
                result = Numo::LAPACK.gesdd *args
            else
                result = Numo::LAPACK.gesvd *args
            end
            u, s, vt = result
            n = s.shape[0]
            s_mat = s.class.new(n, n).fill 0.0
            (0 ... n).each {|i|
                s_mat[i, i] = 1.0 / s[i]
            }
            vt.dot(s_mat).dot(u)
        end

        def tensorinv a, *_
            raise NotImplementedError.new
        end
    end
end
