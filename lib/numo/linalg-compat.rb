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

        def lu a
            result = Numo::LAPACK.getrf(a.transpose)
            lu, _ = result
            result[0] = lu.transpose
            result
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

        def det a
            lu, piv = Numo::LAPACK.getrf a
            piv_sgn = 1
            (0 ... piv.shape[0]).each {|i|
                if piv[i] != i then
                    piv_sgn = (-piv_sgn)
                end
            }
            acc = 1.0
            (0 ... lu.shape.min).each {|i|
                acc *= lu[i, i]
            }
            if acc.zero? then
                if acc.real? then
                    0.0
                else
                    (0.0+0.0i)
                end
            else
                piv_sgn*acc
            end
        end

        def slogdet a
            lu, piv = Numo::LAPACK.getrf a
            if Numo::DFloat === a or Numo::SFloat === a then
                sgn = 1
                (0 ... piv.shape[0]).each {|i|
                    if piv[i] != i then
                        sgn = (-sgn)
                    end
                }
                sum_arr = []
                (0 ... lu.shape.min).each {|i|
                    x = lu[i, i]
                    if x.zero? then
                        return 0, (-Float::INFINITY)
                    else
                        if x.negative? then
                            sgn = (-sgn)
                            x = -x
                        end
                        sum_arr << Math.log(x)
                    end
                }
                acc = enum_sum sum_arr
                if acc == (-Float::INFINITY) then
                    sgn = 0
                end
                return sgn, acc
            elsif Numo::DComplex === a or Numo::SComplex === a then
                sgn = 1.0
                (0 ... piv.shape[0]).each {|i|
                    if piv[i] != i then
                        sgn = (-sgn)
                    end
                }
                sum_arr = []
                (0 ... lu.shape.min).each {|i|
                    x = lu[i, i]
                    x_abs = x.abs
                    if x_abs.zero? then
                        return (0.0+0.0i), (-Float::INFINITY)
                    else
                        sgn *= x/x_abs
                        sum_arr << Math.log(x_abs)
                    end
                }
                acc = enum_sum sum_arr
                if acc == (-Float::INFINITY) then
                    sgn = (0.0+0.0i)
                end
                return sgn, acc
            else
                raise
            end
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

        private

        def enum_sum e
            if e.respond_to? :sum then
                return e.sum
            end
            e.reduce :+
        end
    end
end
