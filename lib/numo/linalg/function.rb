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

  # dot product
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

  #
  def vdot(a, b)
    raise NotImplementedError
  end

  #
  def inner(a, b)
    Blas.call(:gemv, a, b)
  end

  #
  def outer(a, b)
    a[false, :new] * b
  end

  # Matrix product
  def matmul(a, b)
    Blas.call(:gemm, a, b)
  end

  #
  def matrix_power(m, n)
    raise NotImplementedError
  end

  #
  def kron(a, b)
    raise NotImplementedError
  end


  ## Decompositions

  # Cholesky decomposition
  def cholesky(a, uplo:false)
    Lapack.call(:potrf, a, uplo:false)[0]
  end

  # QR factorization
  def qr(a) #, mode)
    Lapack.call(:geqrf, a)[0]
  end

  # Singular Value Decomposition
  def svd(a, driver:'svd')
    case driver.to_s
    when /^(ge)?sdd$/i, "turbo"
      Lapack.call(:gesdd, a, jobz:'A')[0]
    when /^(ge)?svd$/i
      Lapack.call(:gesvd, a, jobu:'A', jobvt:'A')[0]
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
  end

  # Singular values
  def svdvals(a, driver:'svd')
    case driver.to_s
    when /^(ge)?sdd$/i, "turbo"
      Lapack.call(:gesdd, a, jobz:'N')[0]
    when /^(ge)?svd$/i
      Lapack.call(:gesvd, a, jobu:'N', jobvt:'N')[0]
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
  end

  # LU decomposition
  def lu(a)
    Lapack.call(:getrf, a)[0]
  end


  ## Matrix eigenvalues

  # Computes the eigenvalues and, optionally, the left and/or right eigenvectors
  # for a square nonsymmetric matrix A.
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

  # Computes the eigenvalues and, optionally, the left and/or right eigenvectors
  # for a square symmetric/hermitian matrix A.
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

  # Computes the eigenvalues
  # for a square nonsymmetric matrix A.
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

  # Computes the eigenvalues
  # for a square symmetric/hermitian matrix A.
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
  #
  #     |  ord  |  matrix norm            | vector norm                 |
  #     | ----- | ----------------------- | --------------------------- |
  #     |  nil  | Frobenius norm          | 2-norm                      |
  #     | 'fro' | Frobenius norm          |  -                          |
  #     | 'inf' | x.abs.sum(axis:-1).max  | x.abs.max                   |
  #     |    0  |  -                      | (x.ne 0).sum                |
  #     |    1  | x.abs.sum(axis:-2).max  | same as below               |
  #     |    2  | 2-norm (max sing-value) | same as below               |
  #     | other |  -                      | (x.abs**ord).sum**(1.0/ord) |

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

  # Determinant of a matrix
  # @param a [Numo::NArray] matrix (>= 2-dimensional NArray)
  # @return [Float or Complex or Numo::NArray]

  def det(a)
    lu, piv, = Lapack.call(:getrf, a)
    idx = piv.new_narray.store(piv.class.new(piv.shape[-1]).seq(1))
    m = piv.eq(idx).count_false(axis:-1) % 2
    sign = m * -2 + 1
    lu.diagonal.prod(axis:-1) * sign
  end

  # Natural logarithm of the determinant of a matrix
  # @param a [Numo::NArray] matrix (>= 2-dimensional NArray)
  # @return [[sign,logdet]]
  #  sign: A number representing the sign of the determinant
  #  logdet: The natural log of the absolute value of the determinant.

  def slogdet(a)
    lu, piv, = Lapack.call(:getrf, a)
    idx = piv.new_narray.store(piv.class.new(piv.shape[-1]).seq(1))
    m = piv.eq(idx).count_false(axis:-1) % 2
    sign = m * -2 + 1

    lud = lu.diagonal
    if (lud.eq 0).any?
      return 0, (-Float::INFINITY)
    end
    lud_abs = lud.abs
    sign *= (lud/lud_abs).prod
    [sign, NMath.log(lud_abs).sum(axis:-1)]
  end

  # matrix rank using SVD method
  def matrix_rank(m, tol:nil, driver:'svd')
    m = Numo::NArray.asarray(m)
    if m.ndim < 2
      m.ne(0).any? ? 1 : 0
    else
      case driver.to_s
      when /^(ge)?sdd$/, "turbo"
        s = Lapack.call(:gesdd, m, jobz:'N')[0]
      when /^(ge)?svd$/
        s = Lapack.call(:gesvd, m, jobu:'N', jobvt:'N')[0]
      else
        raise ArgumentError, "invalid driver: #{driver}"
      end
      tol ||= s.max(axis:-1, keepdims:true) *
        (m.shape[-2..-1].max * s.class::EPSILON)
      (s > tol).count(axis:-1)
    end
  end

  #
  def trace(*a)
    raise NotImplementedError
  end


  ## Solving equations and inverting matrices


  # Solves linear equation ```a * x = b``` for ```x```
  # from square matrix ```a```
  # @param a [Numo::NArray] n-by-n square matrix
  # @param b [Numo::NArray] n-by-nrhs right-hand-side matrix
  # @param driver [Numo::NArray] optional, default='gen'. one of 'gen','sym','her','pos'
  # @param uplo [String or Symbol] optional, default='U'. Access upper or ('U') lower ('L') triangle.

  def solve(a, b, driver:"gen", uplo:'U')
    case driver.to_s
    when /^gen?(sv)?$/i
      # returns lu, x, ipiv, info
      Lapack.call(:gesv, a, b)[1]
    when /^(sym?|her?|pos?)(sv)?$/i
      func = driver[0..2].downcase+"sv"
      Lapack.call(func, a, b, uplo:uplo)[1]
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
  end


  # Inverse matrix from square matrix ```a```
  # @param a [Numo::NArray] n-by-n square matrix
  # @param b [Numo::NArray] n-by-nrhs right-hand-side matrix
  # @param driver [Numo::NArray] optional, default='gen'. one of 'gen','sym','her','pos'
  # @param uplo [String or Symbol] optional, default='U'. Access upper or ('U') lower ('L') triangle.

  def inv(a, driver:"gen", uplo:'U')
    b = a.new_zeros.eye
    solve(a, b, driver:driver, uplo:uplo)
  end


  # Computes the minimum-norm solution to a linear least squares
  # problem:
  #         minimize 2-norm(| b - A*x |)
  # using the singular value decomposition (SVD) of A.
  # A is an M-by-N matrix which may be rank-deficient.

  # @param a [Numo::NArray] m-by-n matrix
  # @param b [Numo::NArray] m-by-nrhs right-hand-side matrix
  # @param driver [String or Symbol] (optional, default='lsd') one of 'lsd','lss','lsy'
  # @param rcond [Float] (optional, default=-1) RCOND is used to determine the effective rank of A. Singular values S(i) <= RCOND*S(1) are treated as zero. If RCOND < 0, machine precision is used instead.

  def lstsq(a, b, driver:'lsd', rcond:-1)
    a = NArray.asarray(a)
    b = NArray.asarray(b)
    b_orig = nil
    if b.shape.size==1
      b_orig = b
      b = b_orig[true,:new]
    end
    m = a.shape[-2]
    n = a.shape[-1]
    #nrhs = b.shape[-1]
    if m != b.shape[-2]
      raise NArray::ShapeError, "size mismatch: A-row and B-row"
    end
    if m < n   # need to extend b matrix
      shp = b.shape
      shp[-2] = n
      b2 = b.class.zeros(*shp)
      b2[false,0...m,true] = b
      b = b2
    end
    case driver.to_s
    when /^(ge)?lsd$/i
      # x, s, rank, info
      x, s, rank, = Lapack.call(:gelsd, a, b, rcond:rcond)
    when /^(ge)?lss$/i
      # v, x, s, rank, info
      _, x, s, rank, = Lapack.call(:gelss, a, b, rcond:rcond)
    when /^(ge)?lsy$/i
      jpvt = Int32.zeros(*a[false,0,true].shape)
      # v, x, jpvt, rank, info
      _, x, _, rank, = Lapack.call(:gelsy, a, b, jpvt, rcond:rcond)
      s = nil
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
    resids = nil
    if m > n
      if /ls(d|s)$/i =~ driver
        case rank
        when n
          resids = (x[n..-1,true].abs**2).sum(axis:0)
        when NArray
          if true
            resids = (x[false,n..-1,true].abs**2).sum(axis:-2)
          else
            resids = x[false,0,true].new_zeros
            mask = rank.eq(n)
            # NArray does not suppurt this yet.
            resids[mask,true] = (x[mask,n..-1,true].abs**2).sum(axis:-2)
          end
        end
      end
      x = x[false,0...n,true]
    end
    if b_orig && b_orig.shape.size==1
      x = x[true,0]
      resids &&= resids[false,0]
    end
    [x, resids, rank, s]
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
