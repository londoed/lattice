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

    
