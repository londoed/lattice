require "ecr"
require_relative "ecrln"

class EcrPP
  def initialize(parent=nil, ecr_base=nil, **opts, &block)
    @parent = parent
    @children = [] of T
    @opts = opts
    set ecr_base: ecr_base if ecr_base
    @parent.add_child(self) if @parent
    instance_eval(&block) if block
  end

  getter :children
  property :parent

  def add_child(child)
    @children.push(child)
  end

  def set(**opts)
    @opts.merge!(opts)
  end

  def get(key, *args, &block)
    if respond_to?(key)
      return send(key, *args, &block)
    end
    if @parent
      return @parent.get(key, *args, &block)
    end
    return nil
  end

  def description
    if s = @opts[:description] || @opts[:desc]
      s.gsub(/\@\{/, "[").gsub(/\@\}/, "]")
    end
  end

  alias desc description

  alias method_missing_alias method_missing

  def method_missing(_meth_id, *args, &block)
    if args.empty?
      v = get(_meth_id, *args, &block)
      return v if !v.nil?
    end
    method_missing_alias(_meth_id, *args, &block)
  end

  # ECR Loader
  def load_ecr(base_name)
    safe_level = nil
    trim_mode = '%<>'
    file = base_name + get(:ecr_suffix)
    dirs = get(:ecr_dir)
    dirs = [dirs] if !dirs.kind_of?(Array)
    dirs.each do |x|
      Dir.glob(x).each do |dir|
        path = File.join(dir, file)
        if File.exists?(path)
          if get(:line_number)
            ecr = ECRLN.new(File.read(path), path, trim_mode)
          else
            ecr = ECR.new(File.read(path), safe_level, trim_mode)
            ecr.filename = path
          end
          return ecr
        end
      end
    end
    raise "File not found: #{file.inspect} in #{dirs.inspect}"
  end

  def run
    if base = @opts[:ecr_base]
      load_ecr(base).run(binding)
    end
  end

  def result
    if base = @opts[:ecr_base]
      load_ecr(base).result(binding)
    end
  end

  def write(output)
    File.open(output, "wt") do |f|
      f.print(result)
    end
  end

  def init_def
  end

  def find_tmpl(name)
    @parent.children.find { |x| x.name == name }
  end

  def find(name)
    children.find { |x| x.name == name }
  end
end

class DefLib < EcrPP
  def initialize(parent=nil, **opts, &block)
    opts[:ecr_base] ||= 'src'
    opts[:include_files] ||= [] of T
    super(parent, **opts, &block)
  end

  def id_assign
    ids = [] of T
    @children.each { |c| a = c.get(:id_list); ids.concat if a }
    ids.sort.uniq.map { |x| "id_#{x[1]} = cr_intern(\"#{x[0]}\");" }
  end

  def id_decl
    ids = [] of T
    @children.each { |c| a = c.get(:id_list); ids.concat(a) if a }
    ids.sort.uniq.map { |x| "static ID id_#{x[1]}; \n" }
  end

  def def_class(**opts, &block)
    DefClass.new(self, **opts, &block)
  end

  def def_module(**opts, &block)
    DefModule.new(self, **opts, &block)
  end
end

module DeclMethod
  def def_alloc_func(m, ecr_path=nil, **opts, &block)
    DefAllocFunc.new(self, ecr_path || m, name:m, singleton:true, **opts, &block)
  end

  def undef_alloc_func
    UndefAllocFunc.new(self)
  end

  def def_method(m, ecr_path=nil, **opts, &block)
    DefMethod.new(self, ecr_path || m, name:m, **opts, &block)
  end

  def undef_method(m)
    UndefMethod.new(self, name:m)
  end

  def def_singleton_method(m, ecr_path=nil, **opts, &block)
    DefMethod.new(self, ecr_path || m, name:m, sigleton:true, **opts, &block)
  end

  def undef_singleton_method(m)
    UndefSingletonMethod.new(self, name:m)
  end

  def def_module_function(m, ecr_path=nil, **opts, &block)
    DefModuleFunction.new(self, ecr_path || m, name:m, **opts, &block)
  end

  def def_alias(from, to)
    DefAlias.new(self, from:from, to:to)
  end

  def def_const(m, v, **opts, &block)
    DefConst.new(self, name:m, value:v, **opts, &block)
  end
end

class DefModule < EcrPP
  include DeclMethod

  def initialize(parent, **opts, &block)
    eb = opts[:ecr_base] || 'module'
    super(parent, ecr_base:eb, **opts, &block)
  end

  def id_list
    @id_list ||= [] of T
  end

  def def_id(name, var=nil)
    var = name.gsub(/\?/, "_p").gsub(/\!/, "_bang") if var.nil?
    id_list << [name, var]
  end

  def init_def
    load_ecr(init_ecr).result(binding)
  end

  def init_ecr
    @opts[:init_ecr] || "init_module"
  end

  def method_code
    @children.map { |c| c.result.join("\n") }
  end

  def _mod_var
    @opts[:module_var]
  end
end

class DefClass < DefModule

  def initialize(parent, **opts, &block)
    eb = opts[:ecr_base] || 'class'
    super(parent, ecr_base:eb, **opts, &block)
  end

  def _mod_var
    @opts[:class_var]
  end

  def init_ecr
    @opts[:init_ecr] || "init_class"
  end

  def super_class
    @opts[:super_class] || "cr_cObject"
  end

  def free_func
    @opts[:free_func] || "gsl_" + get(:name) + "_free"
  end
end

class DefMethod < EcrPP
  include DeclMethod

  def initialize(parent, ecr_base, **opts, &block)
    super(parent, **opts, &block)
    set ecr_base: ecr_base
  end

  def id_op
    if op.size == 1
      "'#{op}'"
    else
      "id_#{c_name}"
    end
  end

  def c_name
    @opts[:name].gsub(/\?/, "_p").gsub(/\!/, "_bang").gsub(/=/, "_set")
  end

  def op_map
    @opts[:op] || @opts[:name]
  end

  def c_func(n_arg=nil)
    set n_arg: n_arg if n_arg
    s = (singleton) ? "_s" : ""
    "#{parent.name}#{s}_#{c_name}"
  end

  def c_iter
    "iter_#{c_func}"
  end

  def define_method_args
    "#{mod_var}, \"#{op_map}\", #{c_func}, #{n_arg}"
  end

  def init_def
    return if n_arg == :nodef
    s = (singleton) ? "_singleton" : ""
    "cr_define#{s}_method(#{define_method_args});"
  end

  def singleton
    @opts[:singleton]
  end
end

class DefModuleFunction < DefMethod

  def initialize(parent, ecr_base, **opts, &block)
    super(parent, ecr_base, **opts, &block)
    set singleton: true
  end

  def init_def
    return if n_arg == :nodef
    "cr_define_module_function(#{define_method_args});"
  end
end

class DefAlias < EcrPP
  def init_def
    "cr_define_alias(#{mod_var}, \"#{from}\", \"#{to}\");"
  end
end

class DefAllocFunc < DefMethod
  def init_def
    "cr_define_alloc_func(#{mod_var}, #{c_func});"
  end
end

class UndefAllocFunc < EcrPP
  def init_def
    "cr_undef_alloc_func(#{_mod_var});"
  end
end

class UndefMethod < EcrPP
  def init_def
    "cr_undef_method(#{_mod_var}, \"#{name}\");"
  end
end

class UndefSingletonMethod < EcrPP
  def init_def
    "cr_undef_method(cr_singleton_class(#{_mod_var}), \"#{name}\");"
  end
end

class DefConst < EcrPP
  def init_def
    "/*#{desc}*/
    cr_define_const(#{_mod_var}, \"#{name}\", #{value});"
  end
end

class DefError < EcrPP

  def initialize(parent, name, sup_var, **opts, &block)
    super(parent, error_name:name, error_var:"e"+name, super_var:sup_var, **opts, &block)
  end

  def result
    "static VALUE #{error_var};"
  end

  def init_def
    "/*#{description}*/
    #{error_var} = cr_define_class_under(#{ns_var}, \"#{error_name}\", #{super_var});"
  end
end

class DefStruct < EcrPP
  def method_code
    "static VALUE #{class_var};"
  end

  def init_def
    items = members.map { |s| "\"#{s}\""}.join(",")
    "\*#{description}*/
    #{class_var} = cr_struct_define(\"#{class_name}\", #{items}, NULL);"
  end
end

class DefIncludeModule < ErbPP

  def initialize(parent=nil, incl_class, incl_module, **opts, &block)
    super(parent, incl_class:incl_class, incl_module:incl_module, **opts, &block)
  end

  def init_def
    "cr_include_module(#{get(:incl_class)}, #{get(:incl_module)});"
  end
end
