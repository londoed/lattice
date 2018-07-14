set name: "bit"
set type_name: "bit"
set full_class_name: "Lattice::Bit"
set class_name: "Bit"
set class_alias: nil
set class_var: "cT"
set ctype: "BIG_DIGIT"

set has_math: false
set is_bit: true
set is_int: false
set is_unsigned: false
set is_float: false
set is_complex: false
set is_object: false
set is_real: false
set is_comparable: false
set is_double_precison: false
set need_align: false

upcast_cr "Integer"
upcast_cr "Float", "DComplex"
upcast_cr "Complex", "DComplex"

upcast "RObject", "RObject"
upcast "DComplex", "DComplex"
upcast "SComplex", "SComplex"
upcast "DFloat", "DFloat"
upcast "SFloat", "SFloat"
upcast "Int64", "Int64"
upcast "Int32", "Int32"
upcast "Int16", "Int16"
upcast "Int8", "Int8"
upcast "UInt64", "UInt64"
upcast "UInt32", "UInt32"
upcast "UInt16", "UInt16"
upcast "UInt8", "UInt8"
