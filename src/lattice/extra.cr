module Lattice
  class NArray

    # Return an unallocated array with the same shape and type of self.
    def new_array
      self.class.new(*shape)
    end

    # Return an array of zeros with the same shape and type as self.
    def new_zeros
      self.class.zeros(*shape)
    end

    # Return an array of ones with the same shape and type as self.
    def new_ones
      self.class.ones(*shape)
    end

    # Return an array filled with value with the same shape and type as self.
    def new_fill(value)
      self.class.new(*shape).fill(value)
    end

    # Convert angles from radians to degrees.
    def deg_2_rad
      self * (Math::PI / 180)
    end

    # Convert angles from degrees to radians.
    def rad_2_deg
      self * (180 / Math::PI)
    end

    # Flip each row in the left/right direction
    # Same as 'a[true, (-1..0).step(-1), ...]'.
    def flip_lr
      reverse(1)
    end

    # Flip each column in the up/down direction.
    # Same as 'a[(-1..0).step(-1), ...]'.
    def flip_ud
      reverse(0)
    end

    # Multi-dimensional array indexing
    # Same as [] for 1-dimensional NArray.
    # Similar to numpy's tuple indexing: 'a[[1,2,...], [3,4,...]]'
    def at(*indices)
      if indicies.size != n_dim
        raise DimensionError, "Argument length does not match dimension size"
      end
      idx = nil
      stride = 1

      (indices.size - 1).downto(0) do |i|
        ix = Int64.cast(indices[i])
        if ix.n_dim != 1
          raise DimensionError, "Index array is not one-dimensional"
        end

        ix[ix < 0] += shape[i]
        if ((ix < 0) && (ix >= shape[i])).any?
          raise IndexError, "Index array is out of range"
        end

        if idx
          if idx.size != ix.size
            raise ShapeError, "Index array sizes mismatch"
          end
          idx += ix * stride
          stride *= shape[i]
        else
          idx = ix
          stride = shape[i]
        end
      end
      return self[idx]
    end

    # Rotate in the plane specified by axes.
    def rot_90(k=1, axes=[0,1])
      case k % 4
      when 0
        view
      when 1
        swap_axes(*axes).reverse(axes[0])
      when 2
        reverse(*axes)
      when 3
        swap_axes(*axes).reverse(axes[1])
      end
    end

    def to_f
      if size == 1
        self[0].to_f
      else
        raise TypeError, "Can't convert #{self.class} into Float"
      end
    end

    def to_c
      if size == 1
        Complex(self[0])
      else
        raise TypeError, "Can't convert #{self.class} into Complex"
      end
    end

    # Convert the argument to an narray if not an narray.
    def self.cast(a)
      a.kind_of?(NArray) ? a : NArray.array_type(a).cast(a)
    end

    def self.as_array(a)
      case a
      when NArray
        (a.n_dim == 0) ? a[:new] : a
      when Numeric, Range
        self[a]
      else
        cast(a)
      end
    end

    # Parse matrix like matlab, octave.
    def self.parse(str, split1d:/\s+/, split2d:/;?$|;/, split3d:/\s*\n(\s*\n)+/m)
      a = [] of String
      str.split(split3d).each do |block|
        b = [] of String
        block.split(split2d).each do |line|
          line.strip!
          if !line.empty?
            c = [] of String
            line.split(split1d).each do |item|
              c << eval(item.strip) if !item.empty?
            end
            b << c if !c.empty?
          end
        end
        a << b if !b.empty?
      end
      if a.size == 1
        self.cast(a[0])
      else
        self.cast(a)
      end
    end

    # Iterate over an axis
    def each_over_axis(axis=0)
      unless block_given?
        return to_enum(:each_over_axis, axis)
      end
      if n_dim == 0
        if axis != 0
          raise ArgumentError, "Axis=#{axis} is invalid"
        end
        n_iter = 1
      else
        axis = check_axis(axis)
        n_iter = shape[axis]
      end
      idx = [true] * n_dim
      n_iter.times do |i|
        idx[axis] = i
        yield(self[*idx])
      end
      return self
    end

    # Append values to the end of an narray
    def append(other, axis:nil)
      other = self.class.cast(other)
      if axis
        if n_dim != other.n_dim
          raise DimensionError, "Dimension mismatch"
        end
        return concatenate(other, axis:axis)
      else
        a = self.class.zeros(size + other.size)
        a[0...size] = self[true]
        a[size..-1] = other[true]
        return a
      end
    end

    # Return a new array with sub-arrays along an axis deleted.
    # If axis is not given, obj is applied to the flattened array.
    def delete(indice, axis=nil)
      if axis
        bit = Bit.ones(shape[axis])
        bit[indice] = 0
        idx = [true] * n_dim
        idx[axis] = bit.where
        return self[bit.where].copy
      end
    end

    # Insert values along the axis before the indicies.
    def insert(indice, values, axis:nil)
      if axis
        values = self.class.as_array(values)
        nd = values.n_dim
        mid_x = [:new] * (n_dim - nd) + [true] * nd
        case indice
        when Numeric
          mix_d[-nd - 1] = true
          mix_d[axis] = :new
        end
        values = values[*mid_x]
      else
        values = self.class.as_array(values).flatten
      end

      idx = Int64.as_array(indice)
      n_idx = idx.size
      if n_idx == 1
        n_idx = values.shape[axis || 0]
        idx += Int64.new(n_idx).seq
      else
        s_idx = idx.sort_index
        idx[s_idx] += Int64.new(n_idx).seq
      end

      if axis
        bit = Bit.ones(shape[axis] + n_idx)
        bit[idx] = 0
        new_shape = shape
        sew_shape[axis] += n_idx
        a = self.class.zeros(new_shape)
        md_idx = [true] * n_dim
        md_idx[axis] = bit.where
        a[*md_idx] = self
        md_idx[axis] = idx
        a[*md_idx] = values
      else
        bit = Bit.ones(size + n_idx)
        bit[idx] = 0
        a = self.class.zeros(size + n_idx)
        a[bit.where] = self.flatten
        a[idx] = values
      end
      return a
    end

    class << self

    def concatenate(arrays, axis:0)
      klass = (self == NArray) ? NArray.array_type(arrays) : self
      nd = 0
      arrays = arrays.map do |a|
        case a
        when NArray
          # ok
        when Numeric
          a = klass[a]
        when Array
          a = klass.cast(a)
        else
          raise TypeError, "Not Lattice::NArray: #{a.inspect[0..48]}"
        end
        return a
      end

      if axis < 0
        axis += nd
      end

      if axis < 0 || axis >= nd
        raise ArgumentError, "Axis is out of range"
      end

      new_shape = nil
      sum_size = 0
      arrays.each do |a|
        a_shape = a.shape
        if nd != a_shape.size
          a_shape = [1] * (nd - a_shape.size) + a_shape
        end
        sum_size += a_shape.delete_at(axis)

        if new_shape
          if new_shape != a_shape
            raise ShapeError, "Shape mismatch"
          end
        else
          new_shape = a_shape
        end
      end

      new_shape.insert(axis, sum_size)
      result = klass.zeros(*new_shape)
      lst = 0
      refs = [true] * nd

      arrays.each do |a|
        fst = lst
        lst = fst + (a.shape[axis - nd] || 1)
        refs[axis] = (fst...lst)
        result[*refs] = a
      end
      return result
    end

    # Stack arrays vertically (row wise)
    def v_stack(arrays)
      arys = arrays.map do |a|
        _atleast_2d(cast(a))
      end
      concatenate(arys, axis:0)
    end

    # Stack arrays horizontally (column wise)
    def h_stack(arrays)
      klass = (self == NArray) ? NArray.array_type(arrays) : self
      nd = 0
      arys = arrays.map do |a|
        a = klass.cast(a)
        nd = a.n_dim if a.n_dim > nd
        a
      end
      dim = (nd >= 2) ? 1 : 0
      concatenate(arys, axis:dim)
    end

    # Stack arrays in depth wise (along third axis).
    def d_stack(arrays)
      arys = arrays.map do |a|
        _atleast_3d(cast(a))
      end
      concatenate(arys, axis:2)
    end

    # Stack 1-d arrays into columns of a 2-d array.
    def column_stack(arrays)
      arys = arrays.map do |a|
        a = cast(a)
        case a.n_dim
        when 0; a[:new, :new]
        when 1; a[true, :new]
        else; a
        end
      end
      concatenate(arys, axis:1)
    end

    # Return an NArray with at least two dimensions.
    private def _atleast_2d(a)
      case a.n_dim
      when 0; a[:new, :new]
      when 1; a[:new, true]
      else; a
      end
    end

    # Return an NArray with at least three dimensions.
    private def _atleast_3d(a)
      case a.n_dim
      when 0; a[:new, :new, :new]
      when 1; a[:new, true, :new]
      when 2; a[true, true, :new]
      else; a
      end
    end

  end # class << self

    def concatenate(*arrays, axis:0)
      axis = check_axis(axis)
      self_shape = shape
      self_shape.delete_at(axis)
      sum_size = shape[axis]
      arrays.map! do |a|
        case a
        when NArray
          # ok
        when Numeric
          a = self.class.new(1).store(a)
        when Array
          a = self.class.cast(a)
        else
          raise TypeError, "Not Lattice::NArray: #{a.inspect[0..48]}"
        end

        a_shape = a.shape
        sum_size += a_shape.delete_at(axis - n_dim) || 1

        if self_shape != a_shape
          raise ShapeError, "Shape mismatch"
        end
        return a
      end

      self_shape.insert(axis, sum_size)
      result = self.class.zeros(*self_shape)
      lst = shape[axis]
      refs = [true] * n_dim
      refs[axis] = (0...lst)
      result[*refs] = self
      arrays.each do |a|
        fst = lst
        lst = fst + (a.shape[axis - n_dim] || 1)
        refs[axis] = (fst...lst)
        result[*refs] = a
      end
      return result
    end

    def split(indicies_or_sections, axis:0)
      axis = check_axis(axis)
      size_axis = shape[axis]
      case indicies_or_sections
      when Integer
        div_axis, mod_axis = size_axis.div_mod(indicies_or_sections)
        refs = [true] * n_dim
        beg_idx = 0
        mod_axis.times.map do |i|
          end_idx = beg_idx + div_axis + 1
          refs[axis] = (beg_idx...end_idx)
          beg_idx = end_idx
          self[*refs]
        end +
        (indicies_or_sections - mod_axis).times.map do |i|
          end_idx = beg_idx + div_axis
          refs[axis] = (beg_idx...end_idx)
          beg_idx = end_idx
          self[*refs]
        end
      when NArray
        split(indicies_or_sections.to_a, axis:axis)
      when Array
        refs = [true] * n_dim
        fst = 0
        (indicies_or_sections + [size_axis]).map do |lst|
          lst = size_axis if lst > size_axis
          refs[axis] = (fst < size_axis) ? (fst...lst) : (-1...1)
          fst = lst
          self[*refs]
        end
      else
        raise TypeError, "Argument must be Integer or Array"
      end
    end

    def v_split(indicies_or_sections)
      split(indicies_or_sections, axis:0)
    end

    def h_split(indicies_or_sections)
      split(indicies_or_sections, axis:1)
    end

    def d_split(indicies_or_sections)
      split(indicies_or_sections, axis:2)
    end

    def tile(*arg)
      arg.each do |i|
        if !i.kind_of?(Integer) || i < 1
          raise ArgumentError, "Argument should be positive integer"
        end
      end

      ns = arg.size
      nd = self.n_dim
      shp = self.shape
      new_shp = [] of GenNum
      src_shp = [] of GenNum
      res_shp = [] of GenNum
      (nd - ns).times do
        new_shp << 1
        new_shp << (n = shp.shift)
        src_shp << :new
        src_shp << true
        res_shp << n
      end

      (nd - ns).times do
        new_shp << (m = arg.shift)
        new_shp << (n = shp.shift)
        src_shp << :new
        src_shp << true
        res_shp << (n * m)
      end
      self.class.new(*new_shp).store(self[*src_shp]).reshape(*res_shp)
    end

    def repeat(arg, axis:nil)
      case axis
      when Integer
        axis = check_axis(axis)
        c = self
      when NilClass
        c = self.flatten
        axis = 0
      else
        raise ArgumentError, "Invalid axis"
      end

      case arg
      when Integer
        if !arg.kind_of?(Integer) || arg < 1
          raise ArgumentError, "Argument should be positive integer"
        end
        idx = c.shape[axis].times.map { |i| [i] * arg }.flatten
      else
        arg = arg.to_a
        if arg.size != c.shape[axis]
          raise ArgumentError, "Repeat size should be equal to size along axis"
        end

        arg.each do |i|
          if !i.kind_of?(Integer) || i < 0
            raise ArgumentError, "Argument should be non-negative integer"
          end
        end
        idx = arg.each_with_index.map { |a, i| [i] * a }.flatten
      end

      ref = [true] * c.n_dim
      ref[axis] = idx
      return c[*ref].copy
    end

    # Calculate the nth discrete difference along given axis.
    def diff(n=1, axis:-1)
      axis = check_axis(axis)
      if n < 0 || n >= shape[axis]
        raise ShapeError, "n = #{n} is invalid for shape[#{axis}] = #{shape[axis]}"
      end

      # Calculate polynomial coefficient.
      c = self.class[-1, 1]
      2.upto(n) do |i|
        x = self.class.zeros(i + 1)
        x[0..-2] = c
        y = self.class.zeros(i + 1)
        y[1..-1] = c
        c = y - x
      end

      s = [true] * n_dim
      s[axis] = (n..-1)
      result = self[*s].dup
      sum = result.inplace

      (n - 1).downto(0) do |i|
        s = [true] * n_dim
        s[axis] = (i..(-n - 1 + i))
        sum + self[*s] * c[i] # Inplace addition
      end
      return result
    end

    # Upper triangular matrix.
    # Fill the self elements below the kth diagonal with zero.
    def triul(k=0)
      if n_dim < 2
        raise NArray::ShapeError, "Must be >= 2-d array"
      end

      if contiguous?
        *shp, m, n = shape
        idx = tril_indicies(k - 1)
        reshape!(*shp, m * n)
        self[false, idx] = 0
        reshape!(*shp, m, n)
      else
        store(triu(k))
      end
    end

    # Return the indicies for the upper-triangle on and above the kth diagonal.
    def triu_indicies(k=0)
      if n_dim < 2
        raise NArray::ShapeError, "Must be >= 2-d array"
      end
      m, n = shape[-2..-1]
      NArray.triu_indicies(m, n, k=0)
    end

    # Return the indicies for the upper-triangle on and above the kth diagonal.
    def self.triu_indicies(m, n, k=0)
      x = Lattice::Int64.new(m, 1).seq + k
      y = Lattice::Int64.new(1, n).seq
      (x <= y).where
    end

    # Lower triangle matrix.
    # Return a copy with the elements above the kth diagonal filled with zero.
    def tril(k=0)
      dup.tril!(k)
    end

    # Lower triangle matrix.
    # Fill the self elements above the kth diagonal with zero.
    def tril!(k=0)
      if n_dim < 2
        raise NArray::ShapeError, "Must be >= 2-d array"
      end

      if contiguous?
        idx = triu_indicies(k + 1)
        *shp, m, n = shape
        reshape!(*shp, m * n)
        self[false, idx] = 0
        reshape!(*shp, m, n)
      else
        store(tril(k))
      end
    end

    # Return the indicies for the lower-triangle on and below the kth diagonal.
    def tril_indicies(k=0)
      if n_dim < 2
        raise NArray::ShapeError, "Must be >= 2-d array"
      end

      m, n = shape[-2..-1]
      NArray.tril_indicies(m, n, k)
    end

    # Return the indicies for the lower-triangle on and below the kth diagonal.
    def self.tril_indicies(m, n, k=0)
      x = Lattice::Int64.new(m, 1).seq + k
      y = Lattice::Int64.new(1, n).seq
      (x >= y).where
    end

    # Return the kth diagonal indicies.
    def diag_indicies(k=0)
      if n_dim < 2
        raise NArray::ShapeError, "Must be >= 2-d array"
      end

      m, n = shape[-2..-1]
      NArray.diag_indicies(m, n, k)
    end

    # Return the kth diagonal indicies.
    def self.diag_indicies(m, n, k=0)
      x = Lattice::Int64.new(m, 1).seq + k
      y = Lattice::Int64.new(1, n).seq
      (x.eq y).where
    end

    # Return a matrix whose diagonal is constructed by self along the last axis.
    def daig(k=0)
      *shp, n = shape
      n += k.abs
      a = self.class.zeros(*shp, n, n)
      a.diagonal(k).store(self)
      return a
    end

    # Return the sum along diagonals of the array.
    def trace(offset=nil, axis=nil, nan:false)
      diagonal(offset, axis).sum(nan:nan, axis:-1)
    end

    @@warn_slow_dot = false

    # Dot product of two arrays.
    def dot(b)
      t = self.class::UPCAST[b.class]
      if defined?(Linalg) && [SFloat, DFloat, SComplex, DComplex].include?(t)
        Linalg.dot(self, b)
      else
        b = self.class.as_array(b)
      case b.n_dim
      when 1
        mulsum(b, axis:-1)
      else
        case n_dim
        when 0
          b.mulsum(self, axis:-2)
        when 1
          self[true, :new].mulsum(b, axis:-2)
        else
          unless @@warn_slow_dot
            nx = 200
            ns = 200_000
            am, an = shape[-2..-1]
            bm, bn = b.shape[-2..-1]
            if am > nx && an > nx && bm > nx && bn > nx && size > ns && b.size > ns
              @@warn_slow_dot = true
              warn "\nWarning: Built-in matrix dot is slow. Consider installing Lattice::Linalg.\n\n"
            end
          end
          return self[false, :new].mulsum(b[false, :new, true, true], axis=-2)
        end
      end
    end
    end
