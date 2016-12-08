module Numo::Linalg
    class << self

        # module methods

        ## Matrix and vector products

        def dot a, b
            a.dot b
        end

        def vdot a, b
            raise NotImplementedError.new
        end

        def inner a, b
            a.mulsum b
        end

        def outer a, b
            a[false, :new] * b
        end

        def matmul a, b
            result = Numo::LAPACK.gemm a.transpose, b.transpose
            result.transpose
        end

        def tensordot a, b #, axes
            raise NotImplementedError.new
        end

        def einsum *a
            raise NotImplementedError.new
        end

        def matrix_power m, n
            raise NotImplementedError.new
        end

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
                    Numo::LAPACK.gesdd a, vals_only:false
                else
                    Numo::LAPACK.gesvd a, vals_only:false
                end
            u, _, v = result
            result[0] = u.transpose
            result[2] = v.transpose
            result
        end

        def svdvals a, turbo:false
            a = a.transpose
            if turbo then
                Numo::LAPACK.gesdd a, vals_only:true
            else
                Numo::LAPACK.gesvd a, vals_only:true
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

        def matrix_rank m, tol=nil
            m = m.transpose
            if tol then
                Numo::LAPACK.matrix_rank m, tol
            else
                Numo::LAPACK.matrix_rank m
            end
        end

        def slogdet *a
            raise NotImplementedError.new
        end

        def trace *a
            raise NotImplementedError.new
        end


        ## Solving equations and inverting matrices

        def solve a, b
            Numo::LAPACK.gesv a.transpose, b
        end

        def tensorsolve a, b, *_
            raise NotImplementedError.new
        end

        def lstsq a, b
            Numo::LAPACK.gels a.transpose, b
        end

        def inv a
            raise NotImplementedError.new
        end

        def pinv a, *_
            raise NotImplementedError.new
        end

        def tensorinv a, *_
            raise NotImplementedError.new
        end
    end
end
