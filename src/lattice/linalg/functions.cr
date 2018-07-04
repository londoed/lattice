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

  
