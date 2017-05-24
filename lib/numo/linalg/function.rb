module Numo; module Linalg

  module Blas

    FIXNAME =
    {
     cnrm2: :csnrm2,
     znrm2: :dznrm2,
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

    FIXNAME =
    {
     corgqr: :cungqr,
     zorgqr: :zungqr,
    }

    # Call LAPACK function prefixed with BLAS char ([sdcz])
    # defined from data-types of arguments.
    # @param [Symbol,String] func  function name without BLAS char.
    # @param args  arguments passed to Lapack function.
    # @example
    #    s = Numo::Linalg::Lapack.call(:gesv, a)
    def self.call(func,*args)
      fn = (Linalg.blas_char(*args) + func.to_s).to_sym
      fn = FIXNAME[fn] || fn
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

  # Dot product.
  # @param a [Numo::NArray] matrix or vector (>= 1-dimensinal NArray)
  # @param b [Numo::NArray] matrix or vector (>= 1-dimensinal NArray)
  # @return [Numo::NArray] result of dot product
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

  # Matrix product.
  # @param a [Numo::NArray] matrix (>= 2-dimensinal NArray)
  # @param b [Numo::NArray] matrix (>= 2-dimensinal NArray)
  # @return [Numo::NArray] result of matrix product
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


  ## factorization

  # Upper triangular matrix.
  def triu(a,k=0)
    if a.ndim < 2
      raise NArray::ShapeError, "ivalid shape"
    end
    *shp,m,n = a.shape
    b = a.reshape(*shp,m*n)
    x = Numo::Int64.new(m,1).seq + k
    y = Numo::Int64.new(1,n).seq
    b[false,(x>y).where] = 0
    b.reshape(*shp,m,n)
  end

  # Computes a QR factorization of a complex M-by-N matrix A: A = Q \* R.
  #
  # @param a [Numo::NArray] m-by-n matrix A (>= 2-dimensinal NArray)
  # @param mode [String]
  #   - "reduce"  -- returns both Q and R,
  #   - "r"       -- returns only R,
  #   - "economy" -- returns both Q and R but computed in economy-size,
  #   - "raw"     -- returns QR and TAU used in LAPACK.
  # @return [r]        if mode:"r"
  # @return [[q,r]]    if mode:"reduce" or "economic"
  # @return [[qr,tau]] if mode:"raw" (LAPACK geqrf result)

  def qr(a, mode:"reduce")
    qr,tau, = Lapack.call(:geqrf, a)
    *shp,m,n = qr.shape
    r = (m >= n && %w[economic raw].include?(mode)) ?
      triu(qr[false, 0...n, true]) : triu(qr)
    mode = mode.to_s.downcase
    case mode
    when "r"
      return r
    when "raw"
      return [qr,tau]
    when "reduce","economic"
      # skip
    else
      raise ArgumentError, "invalid mode:#{mode}"
    end
    if m < n
      q, = Lapack.call(:orgqr, qr[false, 0...m], tau)
    elsif mode == "economic"
      q, = Lapack.call(:orgqr, qr, tau)
    else
      qqr = qr.class.zeros(*(shp+[m,m]))
      qqr[false,0...n] = qr
      q, = Lapack.call(:orgqr, qqr, tau)
    end
    return [q,r]
  end

  # Computes the Singular Value Decomposition (SVD) of a M-by-N matrix A,
  # and the left and/or right singular vectors.  The SVD is written
  #
  #     A = U * SIGMA * transpose(V)
  #
  # where SIGMA is an M-by-N matrix which is zero except for its
  # min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  # V is an N-by-N orthogonal matrix. The diagonal elements of SIGMA
  # are the singular values of A; they are real and non-negative, and
  # are returned in descending order. The first min(m,n) columns of U
  # and V are the left and right singular vectors of A. Note that the
  # routine returns V**T, not V.
  #
  # @param a [Numo::NArray] m-by-n matrix A (>= 2-dimensinal NArray)
  # @param driver [String or Symbol] choose LAPACK solver from 'svd',
  #  'sdd'. (optional, default='svd')
  # @return [[sigma,u,vt]] SVD result. Array<Numo::NArray>

  def svd(a, driver:'svd')
    case driver.to_s
    when /^(ge)?sdd$/i, "turbo"
      Lapack.call(:gesdd, a, jobz:'A')[0..2]
    when /^(ge)?svd$/i
      Lapack.call(:gesvd, a, jobu:'A', jobvt:'A')[0..2]
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
  end

  # Computes the Singular Values of a M-by-N matrix A.
  # The SVD is written
  #
  #     A = U * SIGMA * transpose(V)
  #
  # where SIGMA is an M-by-N matrix which is zero except for its
  # min(m,n) diagonal elements. The diagonal elements of SIGMA
  # are the singular values of A; they are real and non-negative, and
  # are returned in descending order.
  #
  # @param a [Numo::NArray] m-by-n matrix A (>= 2-dimensinal NArray)
  # @param driver [String or Symbol] choose LAPACK solver from 'svd',
  #  'sdd'. (optional, default='svd')
  # @return [Numo::NArray] returns SIGMA (singular values).

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

  # Computes an LU factorization of a M-by-N matrix A
  # using partial pivoting with row interchanges.
  #
  # The factorization has the form
  #
  #     A = P * L * U
  #
  # where P is a permutation matrix, L is lower triangular with unit
  # diagonal elements (lower trapezoidal if m > n), and U is upper
  # triangular (upper trapezoidal if m < n).
  #
  # @param a [Numo::NArray] m-by-n matrix A (>= 2-dimensinal NArray)
  # @param driver [String or Symbol] choose LAPACK diriver from
  #  'gen','sym','her'. (optional, default='gen')
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle. (omitted when driver:"gen")
  # @return [[lu, ipiv]]
  #   - **lu** [Numo::NArray] -- The factors L and U from the factorization
  #     `A = P*L*U`; the unit diagonal elements of L are not stored.
  #   - **ipiv** [Numo::NArray] -- The pivot indices; for 1 <= i <= min(M,N),
  #      row i of the matrix was interchanged with row IPIV(i).

  def lu_fact(a, driver:"gen", uplo:"U")
    case driver.to_s
    when /^gen?(trf)?$/i
      Lapack.call(:getrf, a)[0..1]
    when /^(sym?|her?)(trf)?$/i
      func = driver[0..2].downcase+"trf"
      Lapack.call(func, a, uplo:uplo)[0..1]
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
  end

  # LU decomposition
  # Computes the inverse of a matrix using the LU factorization
  # computed by Numo::Linalg.lu_fact.
  #
  # This method inverts U and then computes inv(A) by solving the system
  #
  #     inv(A)*L = inv(U)
  #
  # for inv(A).
  #
  # @param lu [Numo::NArray] matrix containing the factors L and U
  #  from the factorization `A = P*L*U` as computed by
  #  Numo::Linalg.lu_fact.
  # @param ipiv [Numo::NArray] The pivot indices from
  #  Numo::Linalg.lu_fact; for 1<=i<=N, row i of the matrix was
  #  interchanged with row IPIV(i).
  # @param driver [String or Symbol] choose LAPACK diriver from
  #  'gen','sym','her'. (optional, default='gen')
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle. (omitted when driver:"gen")
  # @return [Numo::NArray]  the inverse of the original matrix A.

  def lu_inv(lu, ipiv, driver:"gen", uplo:"U")
    case driver.to_s
    when /^gen?(tri)?$/i
      Lapack.call(:getri, lu, ipiv)[0]
    when /^(sym?|her?)(tri)?$/i
      func = driver[0..2].downcase+"tri"
      Lapack.call(func, lu, ipiv, uplo:uplo)[0]
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
  end

  # Solves a system of linear equations
  #
  #     A * X = B  or  A**T * X = B
  #
  # with a N-by-N matrix A using the LU factorization computed by
  # Numo::Linalg.lu_fact
  #
  # @param lu [Numo::NArray] matrix containing the factors L and U
  #  from the factorization `A = P*L*U` as computed by
  #  Numo::Linalg.lu_fact.
  # @param ipiv [Numo::NArray] The pivot indices from
  #  Numo::Linalg.lu_fact; for 1<=i<=N, row i of the matrix was
  #  interchanged with row IPIV(i).
  # @param b [Numo::NArray] the right hand side matrix B.
  # @param driver [String or Symbol] choose LAPACK diriver from
  #  'gen','sym','her'. (optional, default='gen')
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle. (omitted when driver:"gen")
  # @param trans [String or Symbol]
  #   Specifies the form of the system of equations
  #   (omitted if not driver:"gen"):
  #
  #     - If 'N': `A * X = B` (No transpose).
  #     - If 'T': `A*\*T* X = B` (Transpose).
  #     - If 'C': `A*\*T* X = B` (Conjugate transpose = Transpose).
  # @return [Numo::NArray]  the solution matrix X.

  def lu_solve(lu, ipiv, b, driver:"gen", uplo:"U", trans:"N")
    case driver.to_s
    when /^gen?(trs)?$/i
      Lapack.call(:getrs, lu, ipiv, b, trans:trans)[0]
    when /^(sym?|her?)(trs)?$/i
      func = driver[0..2].downcase+"trs"
      Lapack.call(func, lu, ipiv, b, uplo:uplo)[0]
    else
      raise ArgumentError, "invalid driver: #{driver}"
    end
  end


  # Computes the Cholesky factorization of a symmetric/Hermitian
  # positive definite matrix A. The factorization has the form
  #
  #     A = U**H * U,  if UPLO = 'U', or
  #     A = L  * L**H,  if UPLO = 'L',
  #
  # where U is an upper triangular matrix and L is lower triangular
  # @param a [Numo::NArray] n-by-n symmetric matrix A (>= 2-dimensinal NArray)
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle.
  # @return [Numo::NArray] the factor U or L.

  def cho_fact(a, uplo:'U')
    Lapack.call(:potrf, a, uplo:uplo)[0]
  end
  #alias cholesky cho_fact

  # Computes the inverse of a symmetric/Hermitian
  # positive definite matrix A using the Cholesky factorization
  # `A = U**T*U` or `A = L*L**T` computed by Linalg.cho_fact.
  #
  # @param a [Numo::NArray] the triangular factor U or L from the
  #  Cholesky factorization, as computed by Linalg.cho_fact.
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle.
  # @return [Numo::NArray] the upper or lower triangle of the
  #  (symmetric) inverse of A.

  def cho_inv(a, uplo:'U')
    Lapack.call(:potri, a, uplo:uplo)[0]
  end

  # Solves a system of linear equations
  #     A*X = B
  # with a symmetric/Hermitian positive definite matrix A
  # using the Cholesky factorization
  # `A = U**T*U` or `A = L*L**T` computed by Linalg.cho_fact.
  # @param a [Numo::NArray] the triangular factor U or L from the
  #  Cholesky factorization, as computed by Linalg.cho_fact.
  # @param b [Numo::NArray] the right hand side matrix B.
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle.
  # @return [Numo::NArray] the solution matrix X.

  def cho_solve(a, b, uplo:'U')
    Lapack.call(:potrs, a, b, uplo:uplo)[0]
  end


  ## Matrix eigenvalues

  # Computes the eigenvalues and, optionally, the left and/or right
  # eigenvectors for a square nonsymmetric matrix A.
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

  # Computes the eigenvalues and, optionally, the left and/or right
  # eigenvectors for a square symmetric/hermitian matrix A.
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
  #   - **sign** -- A number representing the sign of the determinant.
  #   - **logdet** -- The natural log of the absolute value of the determinant.

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


  # Solves linear equation `a * x = b` for `x`
  # from square matrix `a`
  # @param a [Numo::NArray] n-by-n square matrix  (>= 2-dimensinal NArray)
  # @param b [Numo::NArray] n-by-nrhs right-hand-side matrix (>=
  #  1-dimensinal NArray)
  # @param driver [String or Symbol] choose LAPACK diriver from
  #  'gen','sym','her' or 'pos'. (optional, default='gen')
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle. (omitted when driver:"gen")

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


  # Inverse matrix from square matrix `a`
  # @param a [Numo::NArray] n-by-n square matrix  (>= 2-dimensinal NArray)
  # @param driver [String or Symbol] choose LAPACK diriver
  #  'gen','sym','her' or 'pos'. (optional, default='gen')
  # @param uplo [String or Symbol] optional, default='U'. Access upper
  #  or ('U') lower ('L') triangle. (omitted when driver:"gen")

  def inv(a, driver:"gen", uplo:'U')
    b = a.new_zeros.eye
    solve(a, b, driver:driver, uplo:uplo)
  end


  # Computes the minimum-norm solution to a linear least squares
  # problem:
  #
  #         minimize 2-norm(| b - A*x |)
  #
  # using the singular value decomposition (SVD) of A.
  # A is an M-by-N matrix which may be rank-deficient.
  # @param a [Numo::NArray] m-by-n matrix A (>= 2-dimensinal NArray)
  # @param b [Numo::NArray] m-by-nrhs right-hand-side matrix b
  #   (>= 1-dimensinal NArray)
  # @param driver [String or Symbol] choose LAPACK driver from
  #   'lsd','lss','lsy' (optional, default='lsd')
  # @param rcond [Float] (optional, default=-1)
  #   RCOND is used to determine the effective rank of A.
  #   Singular values `S(i) <= RCOND*S(1)` are treated as zero.
  #   If RCOND < 0, machine precision is used instead.

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
