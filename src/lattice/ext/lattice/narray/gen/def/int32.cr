set name: "int32"
set type_name: "int32"
set full_class_name: "Lattice::Int32"
set class_name: "Int32"
set class_var: "cT"
set ctype: "int32_t"

set has_math: false
set is_bit: false
set is_int: true
set is_unsigned: false
set is_float: false
set is_complex: false
set is_object: false
set is_real: true
set is_comparable: true
set is_double_precison: false
set need_align: true

upcast_cr "Integer"
upcast_cr "Float", "DFloat"
upcast_cr "Complex", "DComplex"

upcast "RObject", "RObject"
upcast "DComplex", "DComplex"
upcast "SComplex", "SComplex"
upcast "DFloat", "DFloat"
upcast "SFloat", "SFloat"
upcast "Int64", "Int64"
upcast "Int32"
upcast "Int16"
upcast "Int8"
upcast "UInt64", "Int64"
upcast "UInt32"
upcast "UInt16"
upcast "UInt8"
