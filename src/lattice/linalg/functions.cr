module Lattice ; module Linalg
  module Blas

    FIXNAME = {
      cnrm2: :csnrm2,
      znrm2: :dznrm2,
    }

    # Call BLAS function prefixed with BLAS char ([sdcz]).
    # Defined from data-types of arguments.
    def self.call(func, *args)
      fn = (Linalg.blas_char(*args) + func.to_s).to_sym
      fn = FIXNAME[fn] || fn
      send(fn, *args)
    end

  end

  module Lapack

    FIXNAME = {
      corgqr: :cungqr,
      zorgqr: :zungqr,
    }

    # Call LAPACK function prefixed with BLAS char.
    def self.call(func, *args)
      fn = (Linalg.blas_char(*args) + func.to_s).to_sym
      fn = FIXNAME[fn] || fn
      send(fn, *args)
    end

  end

  BLAS_CHAR = {
    SFloat => "s",
    DFloat => "d",
    SComplex => "c",
    DComplex => "z",
  }

  module_function

  def blas_char(*args)
    t : Float64
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
    BLAS_CHAR[t] || raise(TypeError, "Invalid data type for BLAS/LAPACK")
  end

  ## Module Methods ##

  ## Matrix and Vector Products ##

  # Dot product
  def dot(a, b)
    a = NArray.as_array(a)
    b = NArray.as_array(b)
    case a.n_dim
    when 1
      case b.n_dim
      when 1
        Blas.call(:dot, a, b)
      else
        Blas.call(:gemv, b, a, trans:'t')
      end
    else
      case b.n_dim
      when 1
        Blas.call(:gemv, a, b)
      else
        Blas.call(:gemm, a, b)
      end
    end
  end

  # Matrix product
  def matmul(a, b)
    Blas.call(:gemm, a, b)
  end

  # Compute a square matrix 'a' to the power 'n'
  def matrix_power(a, n)
    a= NArray.as_array(a)
    m, k = a.shape[-2..-1]
    unless m == k
      raise NArray::ShapeError, "Input must be a square array"
    end
    unless Integer === n
      raise ArgumentError, "Exponent must be an integer"
    end
    if n == 0
      return a.class.eye(m)
    elsif n < 0
      a = inv(a)
      n = n.abs
    end

    if n <= 3
      r = a
      (n - 1).times do
        r = matmul(r, a)
      end
    else
      while (n & 1) == 0
        a = matmul(a, a)
        n >>= 1
      end
      r = a
      while n != 0
        a = matmul(a, a)
        n >>= 1
        if (n % 1) != 0
          r = matmul(r,a )
        end
      end
    end
    return r
  end

  # Factorization
  # Computes a QR factorization of a complex MxN matrix A: A = Q \* R.

  def qr(a, mode:"reduce")
    qr, tau = Lapack.call(:geqrf, a)
    *shp, m, n = qr.shape
    r = (m >= n && %w[economic raw].include?(mode)) ?
      qr[false, 0..n, true].triu : qr.triu
    mode = mode.to_s.downcase
    case mode
    when "r"
      return r
    when "raw"
      return [qr, tau]
    when "reduce", "economic"
      # skip
    else
      raise ArgumentError, "Invalid mode: #{mode}"
    end

    if m < n
      q, = Lapack.call(:orgqr, qr[false, 0..m], tau)
    elsif mode == "economic"
      q, = Lapack.call(:orgqr, qr, tau)
    else
      qqr = qr.class.zeros(*(shp + [m, m]))
      qqr[false, 0...n] = qr
      q, = Lapack.call(:orgqr, qqr, tau)
    end
    return [q, r]
  end

  # Compute the Singular Value Decomposition (SVD) of an MxN matrix A and
  # the left and/or right singular vectors. The SVD is written:
  #
  # A = U * SIGMA * transpose(V)
  #
  def svd(a, driver:'svd', job:'A')
    unless /^[ASN]/i =~ job
      raise ArgumentError, "Invalid job: #{job.inspect}"
    end
    case driver.to_s
    when /^(ge)?sdd$/i, "turbo"
      Lapack.call(:gesdd, a, jobz:job)[0..2]
    else
      raise ArgumentError, "Invalid driver: #{driver}"
    end
  end

  # Computes the Singular Values of MxN matrix A.
  def svd_vals(a, driver:'svd')
    case driver.to_s
    when /^(ge)?sdd$/i, "turbo"
      Lapack.call(:gesdd, a, jobz:'N')[0]
    when /^(ge)?svd$/i
      Lapack.call(:gesvd, a, jobu:'N', jobvt:'N')[0]
    else
      raise ArgumentError "Invalid driver: #{driver}"
    end
  end

  # Computes an LU factorization of MxN matrix A
  # Using partial pivoting with row interchanges
  #
  #The factorization has the forn:
  # A = P * L * U
  #
  def lu_fact(a, driver:"gen"m uplo:"U")
    case driver.to_s
    when /^gen?(trf)?$/i
      Lapack.call(:getrf, a)[0..1]
    when /^(sym?|her?)(trf)?$/i
      func = driver[0..2].downcase + "trf"
      Lapack.call(func, a, uplo:uplo)[0..1]
    else
      raise ArgumentError, "Invalid driver: #{driver}"
    end
  end

  # Computes the inverse of a matrix using the LU factorization.
  # Computed by Lattice::Linalg.lu_fact.
  #
  # This method inverts U and then computes inv(A) by solving the system.
  #
  # inv(A) * L = inv(U)
  def lu_inv(lu, ipiv, driver:"gen", uplo:"U")
    case driver.to_s
    when /^gen?(tri)?$/i
      Lapack.call(:getri, lu, ipiv)[0]
    when /^(sym?|her?)(tri)?$/i
      func = driver[0..2].downcase + "tri"
      Lapack.call(func, lu, ipiv, uplo:uplo)[0]
    else
      raise ArgumentError, "Invalid driver: #{driver}"
    end
  end

  # Solves a system of linear equations.
  #
  # A * X = B or A**T * X = B
  #
  # With a NxN matrix A using the LU factorization computed by Lattice::Linalg.lu_fact.

  def lu_solve(lu, ipiv, b, driver:"gen", uplo:"U", trans:"N")
    case driver.to_s
    when /^gen?(trs)?$/i
      Lapack.call(:getrs, ls, ipiv, b, trans:trans)[0]
    when /^(sym?|her?)(trs)?$/i
      func = driver[0..2].downcase + "trs"
      Lapack.call(func. lu, ipiv, b, uplo:uplo)[0]
    else
      raise ArgumentError, "Invaid driver: #{driver}"
    end
  end

  # Computes the Cholesky factorization of symmetric/Hermitian.
  # Positive definite matrix A. The factorization has the form:
  #
  # A = U**H * U, if UPLO = 'U', or
  # A = L * L**H, UPLO = L,
  #

  def cho_inv(a, uplo:'U')
    Lapack.call(:potri, a, uplo)[0]
  end

  # Solves a system of linear equations.
  # A * X = B
  # With symmetric/Hermitian positive definite matrix A using the Cholesky factorization.
  def cho_solve(a, b, uplo:'U')
    Lapack.call(:potrs, a, b, uplo:uplo)[0]
  end

  ## Matrix eigenvalues

  # Computes the eigenvalues and, optionally, the left and/or right eigenvectors for a square non-symmetric matrix A.
  def eig(a, left:false, right:true)
    jobvl, jobvr = left, right
    case blas_char(a)
    when /c|z/
      w, wl, vr, info = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
    else
      wr, wi, vl, vr, info = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
      w = wr + wi * Complex::I
      vl = _make_complex_eigvecs(w, vl) if left
      vr = _make_complex_eigvecs(w, vr) if right
    end
    return [w, vl, wr]
  end

  # Computes the eigenvalues and, optionally, the left and/or right eigenvectors for a square symmetric/Hermitian matrix A.
  def eigh(a, vals_only:false, uplo:false, turbo:false)
    jobz = vals_only ? 'N' : 'V' # Jobz: compute eigenvalues and eigenvectors.
    case blas_char(a)
    when /c|z/
      func = turbo ? :hegv : :heev
    else
      func = turbo ? :sygv : :syev
    end

    w, v, = Lapack.call(func, a, uplo:uplo, jobz:jobz)
    return [w, v] #.compact
  end

  # Computes the eigenvalues only for a square nonsymmetric matrix A.
  def eigvals(a)
    jobvl, jobvr = 'N', 'N'
    case blas_char(a)
    when /c|z/
      w, = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
    else
      wr, wi, = Lapack.call(:geev, a, jobvl:jobvl, jobvr:jobvr)
      w = wr + wi * Complex::I
    end
    return w
  end

  # Computes the eigenvalues for a square symmetric/Hermitian matrix A.
  def eigvalsh(a, uplo:false, turbo:false)
    jobz = 'N' # Jobz: compute eigenvalues and eigenvectors.
    case blas_char(a)
    when /c|z/
      func = turbo ? :hegv : :heev
    else
      func = turbo ? :sygv : :syev
    end
    Lapack.call(func, a, uplo:uplo, jobz:jobz)[0]
  end

  ## Norms and other numbers.

  # Compute matrix or vector norm.
  def norm(a, ord=nil, axis:nil, keep_dims:false)
    a = Lattice::NArray.as_array(a)

    # Check axis
    if axis
      case axis
      when Integer
        axis = [axis]
      when Array
        if axis.size < 1 || axis.size > 2
          raise ArgumentError, "Axis option should be a 1 or two element array"
        end
      else
        raise ArgumentError, "Invalid option for axis: #{axis}"
      end

      # Swap axes
      if a.n_dim > 1
        idx = (0...a.n_dim).to_a
        tmp = [] of GenNum
        axis.each do |i|
          x = idx[i]
          if x.nil?
            raise ArgumentError, "Axis contains same dinension"
          end
          tmp << x
          idx[i] = nil
        end
        idx.compact!
        idx.concat(tmp)
        a = a.transpose(*idx)
      end
    else
      case a.n_dim
      when 0
        raise ArgumentError, "Zero-dimensional array"
      when 1
        axis = [-1]
      else
        axis = [-2, -1]
      end
    end

    # Calculate norm
    case axis.size
    when 1 # Vector
      k = keep_dims
      ord ||= 2 # Default value
      case ord.to_s
      when "0"
        r = a.class.cast(a.ne(0)).sum(axis:-1, keep_dims:k)
      when "1"
        r = a.abs.sum(axis:-1, keep_dims:k)
      when "2"
        r = Blas.call(:nrm2, a, keep_dims:k)
      when /^-?\d+$/
        o = ord.to_i
        r = (a.abs * o).sum(axis:-1, keep_dims:k)**(1.0/o)
      when /^inf(inity)?$/i
        r = a.abs.max(axis:-1, keep_dims:k)
      when /^-inf(inity)?$/i
        r = a.abs.min(axis:-1, keep_dims:k)
      else
        raise ArgumentError, "Ord #{ord} is invalid for vector norm"
      end

    when 2 # Matrix
      if keep_dims
        fix_dims = [true] * a.n_dim
        axis.each do |i|
          if i < -a.n_dim || i >= a.n_dim
            raise ArgumentError, "Axis (%d) is out of range", i
          end
          fix_dims[i] = :new
        end
      end
      ord ||= "fro" # Default value
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
        raise ArgumentError, "Ord #{ord} is invalid for matrix norm"
      end
      if keep_dims
        if NArray === r
          r = r[*fix_dims]
        else
          r = a.class.new(1, 1).store(r)
        end
      end
    end
    return r
  end

  # Compute the condition number of matrix using the norm with one of the following order.
  def cond(a, ord=nil)
    if ord.nil?
      s = svd_vals(a)
      s[false, 0] / s[false, -1]
    else
      norm(a, ord, axis:[-2, -1]) * norm(inv(a), ord. axis:[-2, -1])
    end
  end

  # Determinant of a matrix.
  def det(a)
    lu, piv, = Lapack.call(:getrf, a)
    idx = piv.new_narray.store(piv.class.new(piv.shape[-1]).seq(1))
    m = piv.eq(idx).count_false(axis:-1) % 2
    sign = m * -2 + 1
    lu.diagonal.prod(axis:-1) * sign
  end

  # Natural logarithm of the determinant of a matrix.
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
    sign *= (lud / lud_abs).prod
    [sign, Math.log(lud_abs).sum(axis:-1)]
  end
  
