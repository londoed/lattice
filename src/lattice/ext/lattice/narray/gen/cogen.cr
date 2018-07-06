this_dir = File.dirname(__FILE__)
lib_path = File.absolute_path(File.dirname(__FILE__)) + "/../../../../lib"
@@LOAD_PATH.unshift(lib_path)

require_relative "narray_def"

@@line_number = false

while true
  if ARGV[0] == "-l"
    @@line_number = true
    ARGV.shift
  elsif ARGV[0] == "-o"
    ARGV.shift
    @@output = ARGV.shift
    require "fileutils"
    FileUtils.rm_f(@@output)
  else
    break
  end
end

if ARGV.size != 1
  puts "Usage:\n crystal #{@@0} [-1] ecr_base [type_file]"
  exit 1
end

type_file = ARGV
type_name = File.basename(type_file, ".cr")

ecr_dir = ["tmpl"]
ecr_dir.unshift("tmpl_bit") if (type_name == "bit")
ecr_dir.map! { |d| File.join(this_dir, d) }

code = DefLib.new do
  set line_number: @@line_number
  set ecr_dir: ecr_dir
  set ecr_suffix: ".c"
  set ns_var: "mLattice"

  set file_name: @@output || ""
  set include_files: ["lattice/types/#{type_name}.h"]
  set lib_name: "lattice_" + type_name

  def_class do
    extend LatticeMethod
    extend LatticeType
    eval File.read(type_file), binding, type_file
    eval File.read(File.join(this_dir, "spec.cr")), binding, "spec.cr"
  end
end.result

if @@output
  open(@@output, "w").write(code)
else
  @@stdout.write(code)
end
