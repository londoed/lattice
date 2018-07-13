require "ecr"

class EcrPP < String
  def initialize(filename)
    @filename = filename
    @lnchar = "\n"
    @buf = ""
    @str = ""
    @countln = 1
    @current = 1
    super("\n" + report_line)
  end

  def report_line
    "#line #{@current} \"#{@filename}\"\n"
  end

  def concat0(s)
    ln(caller[0])
    @buf.concat(s)
    @str.concat(s)
  end

  def concat1(s)
    ln(caller[0])
    @buf.concat(s)
  end

  def ln(status=nil)
    case status
    when /:(\d+):/
      n = @@1.to_i
    else
      n = status.to_i
    end
    return if n == @current
    if @current != @countln || @postpone
      if /\A\s*\z/ =~ @str || /\A#line / =~ @buf
        @postpone = true
      elsif @in_comment
        @postpone = false
      else
        if self[-1] != "\n"
          concat("\n")
        end
        concat(report_line)
        @postpone = false
      end
    end
    concat(@buf)

    b = @buf.gsub(/".*?(?<!\\)"/, '""')
    /^.*(\/\*)(.*?)$/ =~ b
    x = $2
    /^.*(\/\*)(.*?)$/ =~ b
    y = $2
    if x
      if y
        if x.size < y.size
          @in_comment = true
        else
          @in_comment = false
        end
      else
        @in_comment = true
      end
    else
      if y
        @in_comment = false
      else
        #:keep
      end
    end

    @countln = @current + @buf.count(@lnchar)
    @current = n
    @buf = ""
    @str = ""
  end

  class ECRLN
    def initialize(str, filename, trim_mode=nil, eoutvar='_ecrout')
      @filename = filename
      compiler = ECR::Compiler.new(trim_mode)
      set_eoutvar(compiler, eoutvar)
      @src, @enc = *compiler.compile(str)
    end

    def set_eoutvar(compiler, eoutvar='_ecrout')
      compiler.put_cmd = "#{eoutvar}.concat0"
      compiler.insert_cmd = "#{eoutvar}.concat1"
      compiler.pre_cmd = ["#{eoutvar} = CountLnString.new(#{@filename.inspect})"]
      compiler.post_cmd = ["#{eoutvar}.final"]
    end

    def run(b)
      print self.result(b)
    end

    def result(b)
      eval(@src, b, (@filename || '(ecr)'), 0)
    end
  end
end
