require "numo/narray"
require "numo/linalg/blas"
require "numo/linalg/lapack"

module Numo; module Linalg

  module Blas
    FIXNAME =
    {
     znrm2: :dznrm2, cnrm2: :csnrm2,
    }
    # Call BLAS function prefixed with BLAS char ([sdcz])
    # defined from data-types of arguments.
    # @param [Symbol] func  function name without BLAS char.
    # @param args  arguments passed to Blas function.
    # @example
    #    c = Numo::Linalg::Blas.call(:gemm, a, b)
    def self.call(func,*args)
      fn = (Linalg.blas_char(*args) + func.to_s).to_sym
      fn = FIXNAME[fn] || fn
      send(fn,*args)
    end
  end

  module Lapack
    # Call LAPACK function prefixed with BLAS char ([sdcz])
    # defined from data-types of arguments.
    # @param [Symbol,String] func  function name without BLAS char.
    # @param args  arguments passed to Lapack function.
    # @example
    #    s = Numo::Linalg::Lapack.call(:gesv, a)
    def self.call(func,*args)
      fn = (Linalg.blas_char(*args) + func.to_s).to_sym
      send(fn,*args)
    end
  end

  BLAS_CHAR =
  {
   SFloat => "s",
   DFloat => "d",
   SComplex => "c",
   DComplex => "z",
  }

  module_function

  def blas_char(*args)
    t = Float
    args.each do |a|
      k =
        case a
        when NArray
          a.class
        when Array
          NArray.array_type(a)
        end
      if k && k < NArray
        t = k::UPCAST[t]
      end
    end
    BLAS_CHAR[t] || raise(TypeError,"invalid data type for BLAS/LAPACK")
  end

  # module methods

  ## Matrix and vector products

  # not linalg
  def dot(a, b)
    if a.ndim >= 2
      if b.ndim >= 2
        return Blas.call(:gemm, a, b)
      else
        return Blas.call(:gemv, a, b)
      end
    else
      if b.ndim >= 2
        return Blas.call(:gemv, b, a, trans:'t')
      else
        return Blas.call(:dot, a, b)
      end
    end
  end

  # not linalg
  def vdot(a, b)
    raise NotImplementedError
  end

  # not linalg
  def inner(a, b)
    Blas.call(:gemv, a, b)
  end

  # not linalg
  def outer(a, b)
    a[false, :new] * b
  end

  def matmul(a, b)
    Blas.call(:gemm, a, b)
  end

  # not linalg
  def matrix_power(m, n)
    raise NotImplementedError
  end

  # not linalg
  def kron(a, b)
    raise NotImplementedError
  end


  ## Decompositions

  def cholesky(a, uplo:false)
    Lapack.call(:potrf, a, uplo:false)[0]
  end

  def qr(a) #, mode)
    Lapack.call(:geqrf, a)[0]
  end

  def svd(a, turbo:false)
    (turbo ?
     Lapack.call(:gesdd, a, jobz:'A') :
     Lapack.call(:gesvd, a, jobu:'A', jobvt:'A')
    )[0]
  end

  def svdvals(a, turbo:false)
    (turbo ?
     Lapack.call(:gesdd, a, jobz:'N') :
     Lapack.call(:gesvd, a, jobu:'N', jobvt:'N')
    )[0]
  end

  def lu(a)
    Lapack.call(:getrf, a)[0]
  end


  ## Matrix eigenvalues

  def eig(a, left:false, right:true)
    jobvl, jobvr = left, right
    case blas_char(a)
    when /cz/
      w, vl, vr, info = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
    else
      wr, wi, vl, vr, info = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
      w  = wr + wi * Complex::I
      vl = _make_complex_eigvecs(w,vl) if left
      vr = _make_complex_eigvecs(w,vr) if right
    end
    [w,vl,vr].compact
  end

  def eigh(a, vals_only:false, uplo:false, turbo:false)
    jobz = vals_only ? 'N' : 'V' # jobz: Compute eigenvalues and eigenvectors.
    case blas_char(a)
    when /cz/
      func = turbo ? :hegv : :heev
    else
      func = turbo ? :sygv : :syev
    end
    w, v, = Lapack.call(func, a, uplo:uplo, jobz:jobz)
    [w,v].compact
  end

  def eigvals(a)
    jobvl, jobvr = 'N','N'
    case blas_char(a)
    when /cz/
      w, = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
    else
      wr, wi, = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
      w  = wr + wi * Complex::I
    end
    w
  end

  def eigvalsh(a, uplo:false, turbo:false)
    jobz = 'N' # jobz: Compute eigenvalues and eigenvectors.
    case blas_char(a)
    when /cz/
      func = turbo ? :hegv : :heev
    else
      func = turbo ? :sygv : :syev
    end
    Lapack.call(func, a, uplo:uplo, jobz:jobz)[0]
  end


  ## Norms and other numbers

  def norm(a, ord=nil, axis:nil, keepdims:false)
    a = Numo::NArray.asarray(a)

    # check axis
    if axis
      case axis
      when Integer
        axis = [axis]
      when Array
        if axis.size < 1 || axis.size > 2
          raise ArgmentError, "axis option should be 1- or 2-element array"
      end
      else
        raise ArgumentError, "invalid option for axis: #{axis}"
      end
      # swap axes
      if a.ndim > 1
        idx = (0...a.ndim).to_a
        tmp = []
        (axis.size-1).downto(0) do |i|
          tmp.unshift( idx.delete_at(axis[i]) )
        end
        idx.concat(tmp)
        a = a.transpose(*idx)
      end
    else
      case a.ndim
      when 0
        raise ArgumentError, "zero-dimensional array"
      when 1
        axis = [-1]
      else
        axis = [-2,-1]
      end
    end

    # calculate norm
    case axis.size

    when 1  # vector
      k = keepdims
      ord ||= 2  # default
      case ord.to_s
      when "0"
        r = a.class.cast(a.ne(0)).sum(axis:-1, keepdims:k)
      when "1"
        r = a.abs.sum(axis:-1, keepdims:k)
      when "2"
        r = Blas.call(:nrm2, a, keepdims:k)
      when /^-?\d+$/
        o = ord.to_i
        r = (a.abs**o).sum(axis:-1, keepdims:k)**(1.0/o)
      when /^inf(inity)?$/i
        r = a.abs.max(axis:-1, keepdims:k)
      when /^-inf(inity)?$/i
        r = a.abs.min(axis:-1, keepdims:k)
      else
        raise ArgumentError, "ord (#{ord}) is invalid for vector norm"
      end

    when 2  # matrix
      if keepdims
        fixdims = [true] * a.ndim
        axis.each do |i|
          if i < -a.ndim || i >= a.ndim
            raise ArgmentError, "axis (%d) is out of range", i
          end
          fixdims[i] = :new
        end
      end
      ord ||= "fro"  # default
      case ord.to_s
      when "1"
        r, = Lapack.call(:lange, a, '1')
      when "-1"
        r = a.abs.sum(axis:-2).min(axis:-1)
      when "2"
        svd, = Lapack.call(:gesvd, a, jobu:'N', jobvt:'N')
        r = svd.max(axis:-1)
      when "-2"
        svd, = Lapack.call(:gesvd, a, jobu:'N', jobvt:'N')
        r = svd.min(axis:-1)
      when /^f(ro)?$/i
        r, = Lapack.call(:lange, a, 'F')
      when /^inf(inity)?$/i
        r, = Lapack.call(:lange, a, 'I')
      when /^-inf(inity)?$/i
        r = a.abs.sum(axis:-1).min(axis:-1)
      else
        raise ArgumentError, "ord (#{ord}) is invalid for matrix norm"
      end
      if keepdims
        if NArray===r
          r = r[*fixdims]
        else
          r = a.class.new(1,1).store(r)
        end
      end
    end
    return r
  end

  def cond(*x)
    raise NotImplementedError
  end

  def det(a)
    lu, piv, = Lapack.call(:getrf, a)
    idx = piv.new_empty.store(1..piv.shape[-1])
    m = piv.eq(idx).count_false(axis:-1) % 2
    sign = m * -2 + 1
    lu.diagonal.prod(axis:-1) * sign
  end

  def slogdet(a)
    lu, piv, = Lapack.call(:getrf, a)
    idx = piv.new_empty.store(1..piv.shape[-1])
    m = piv.eq(idx).count_false(axis:-1) % 2
    sign = m.inplace * -2 + 1

    lud = lu.diagonal
    if (lud.eq 0).any?
      return 0, (-Float::INFINITY)
    end
    lud_abs = lud.abs
    sign *= (lud/lud_abs).prod
    NMath.log(lud_abs).sum(axis:-1) * sign
  end

  def matrix_rank(m, tol:nil, turbo:true)
    m = Numo::NArray.asarray(m)
    if m.ndim < 2
      m.ne(0).any? ? 1 : 0
    else
      s = (turbo ?
           Lapack.call(:gesdd, m, jobz:'N') :
           Lapack.call(:gesvd, m, jobu:'N', jobvt:'N')
          )[0]
      tol ||= s.max(axis:-1, keepdims:true) *
              (m.shape[-2..-1].max * s.class::EPSILON)
      (s > tol).count(axis:-1)
    end
  end

  # not linalg
  def trace(*a)
    raise NotImplementedError
  end


  ## Solving equations and inverting matrices

  def solve(a, b)
    Lapack.call(:gesv, a, b)[0]
  end

  def lstsq(a, b)
    Lapack.call(:gels, a, b)[0]
  end

  def inv(a)
    b = a.new_zeros.eye
    Lapack.call(:gesv, a, b)[0]
  end

  # @!visibility private
  def _make_complex_eigvecs(w, vin) # :nodoc:
    v = w.class.cast(vin)
    # broadcast to vin.shape
    m = (w.imag > 0 | Bit.zeros(*vin.shape)).where
    v[m].imag = vin[m+1]
    v[m+1] = v[m].conj
    v
  end

=begin
  # numpy.linalg.pinv and scipy.linalg.pinv are different
  def pinv(a, turbo:false)
    tflag = turbo
    a = a.conj
    args = [a, job: :THIN]
    if tflag
      result = Lapack.gesdd(*args)
    else
      result = Lapack.gesvd(*args)
    end
    u, s, vt = result
    n = s.shape[0]
    s_mat = s.class.new(n, n).fill(0.0)
    (0 ... n).each {|i|
      s_mat[i, i] = 1.0 / s[i]
    }
    vt.dot(s_mat).dot(u)
  end
=end

end
end
