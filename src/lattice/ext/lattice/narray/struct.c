/*
  struct.c
  Numerical Array Extension for Crystal lang
*/

#include <crystal.h>
#include "lattice/narray.h"
#include "lattice/template.h"

#define cT lattice_cStruct
VALUE cT;

static VALUE
nst_allocate(VALUE self) {
  narray_t *na;
  char *ptr;
  VALUE velmsz;

  GetNArray(self, na);

  switch(NA_TYPE(na)) {
  case NARRAY_DATA_T:
    ptr = NA_DATA_PTR(na);
    if (na->size > 0 && ptr == NULL) {
      velmsz = cr_const_get(CLASS_OF(self), cr_intern("element_byte_size"));
      ptr = xmalloc(NUM2SIZET(velmsz) * na->size);
      NA_DATA_PTR(na) = ptr;
    }
    break;
  case NARRAY_VIEW_T:
    cr_funcall(NA_VIEW_DATA(na), cr_intern("allocate"), 0);
    break;
  case NARRAY_FILEMAP_T:
    //ptr = ((narray_filemap_t*)na)->ptr;
  default:
    cr_bug("Invalid NArray type : %d", NA_TYPE(na));
  }
  return self;
}

static inline VALUE
nst_definitions(VALUE nst) {
  return cr_const_get(nst, cr_intern("DEFINITIONS"));
}

static VALUE
nst_definition(VALUE nst) {
  long i;
  VALUE def = nst_definitions(CLASS_OF(nst));
  long len = RARRAY_LEN(def);

  if (TYPE(idx) == T_STRONG || TYPE(idx) == T_SYMBOL) {
    ID id = cr_to_id(idx);
    for (i=0; i<len; i++) {
      VALUE key = RARRAY_AREF(RARRAY_AREF(def, i), 0);
      if (SYM2ID(key) == id) {
        return RARRAY_AREF(def, i);
      }
    }
  } else if (cr_obj_is_kind_of(idx, cr_cNumeric)) {
    i = NUM2LONG(idx);
    if (i<-len || i >= len) {
      cr_raise(cr_eIndexError, "Offset %ld out of range of struct(size:%ld)", i, len);
    }
    return RARRAY_AREF(def, i);
  }
  return Qnil;
}

void na_copy_array_structure(VALUE self, VALUE view);

static VALUE
na_make_view_struct(VALUE self, VALUE dtype, VALUE offset) {
  size_t i, n;
  int j, k, ndim;
  size_t *shape;
  size_t *idx1, *idx2;
  ssize_t stride;
  stridx_t *stridx;
  narray_t *na, *nt;
  narray_view_t *na1, *na2;
  VALUE klass;
  volatile VALUE view;

  GetNArray(self, na);

  // Build from Lattice::Struct
  if (cr_obj_is_kind_of(dtype, cNArray)) {
    GetNArray(dtype, nt);
    ndim = na->ndim + nt->ndim;
    shape = ALLOCA_N(size_t, ndim);
    // Struct dimensions
    for (j=0; j<na->ndim; j++) {
      shape[j] = na->shape[j];
    }
    // Member dimension
    for (j=na->ndim, k=0; j<ndim; j++, k++) {
      shape[j] = nt->shape[k];
    }
    klass = CLASS_OF(dtype);
    stridx = ALLOC_N(stridx_t, ndim);
    stride = na_dtype_elmsz(klass);
    for (j=ndim, k=nt->ndim; k; ) {
      SDX_SET_STRIDE(stridx[--j], stride);
      stride *= nt->shape[--k];
    }
  } else {
    ndim = na->ndim;
    shape = ALLOCA_N(size_t, ndim);
    for (j=0; j<nidm; j++) {
      shape[j] = na->shape[j];
    }
    klass = CLASS_OF(self);
    if (TYPE(dtype) == T_CLASS) {
      if (RTEST(cr_class_inherited_p(dtype, cNArray))) {
        klass = dtype;
      }
    }
    stridx = ALLOC_N(stridx_t, ndim);
  }

  view = na_s_allocate_view(klass);
  na_copy_flags(self, view);
  GetNArrayView(view, na2);
  na_setup_shape((narray_t*)na2, ndim, shape);
  na2->stridx = stridx;

  switch(na->type) {
  case NARRAY_DATA_T:
  case NARRAY_FILEMAP_T:
    stride = nary_element_stride(self);
    for (j=na->ndim; j--;) {
      SDX_SET_STRIDE(na2->stridx[j], stride);
      stride *= na->shape[j];
    }
    na2->offset = 0;
    na2->data = self;
    break;
  case NARRAY_VIEW_T:
    GetNArrayView(self, na1);
    for (j=na1->base.ndim; j--;) {
      if (SDX_IS_INDEX(na1->stridx[j])) {
        n = na1->base.shape[j];
        idx1 = SDX_GET_INDEX(na1->stridx[j]);
        idx2 = ALLOC_N(size_t, na1->base.shape[j]);
        for (i=0; i<n; i++) {
          idx2[i] = idx1[i];
        }
        SDX_SET_INDEX(na2->stridx[j], idx2);
      } else {
        na2->stridx[j] = na1->stridx[j];
      }
    }
    na2->offset = na1->offset;
    na2->data = na1->data;
    break;
  }
  if (RTEST(offset)) {
    na2->offset += NUM2SIZET(offset);
  }
  return view;
}

static VALUE
nst_field_view(VALUE self, VALUE idx) {
  VALUE def, type, ofs;

  def = nst_definition(self, idx);
  if (!RTEST(def)) {
    idx = cr_funcall(idx, cr_intern("to_s"), 0);
    cr_raise(cr_eTypeError, "Invalid field: '%s' for struct %s", StringValuePtr(idx), cr_class2name(CLASS_OF(self)));
  }
  type = RARRAY_AREF(def, 1);
  ofs = RARRAY_AREF(def, 2);
  return na_make_view_struct(self, type, ofs);
}

static VALUE
nst_field(VALUE self, VALUE idx) {
  VALUE obj;
  narray_view_t *nv;

  obj = nst_field_view(self, idx);
  GetNArrayView(obj, nv);
  if (nv->base.ndim == 0) {
    obj = cr_funcall(obj, cr_intern("extract"), 0);
  }
  return obj;
}

static VALUE
nst_field_set(VALUE self, VALUE idx, VALUE other) {
  VALUE obj;

  obj = nst_field_view(self, idx);
  cr_funcall(obj, cr_intern("store"), 1, other);
  return other;
}

static VALUE
nst_method_missing(int argc, VALUE *argv, VALUE self) {
  VALUE s, tag, obj;

  if (argc == 2) {
    s = cr_sym_to_s(argv[0]);
    if (RSTRING_PTR(s)[RSTRING_LEN(s) - 1] == '=') {
      tag = cr_str_intern(cr_str_new(RSTRING_PTR(s), RSTRING_LEN(s) - 1));
      obj = nst_field(self, tag);
      if (RTEST(obj)) {
        cr_funcall(obj, cr_intern("store"), 1, argv[1]);
        return argv[1];
      }
    }
    return cr_call_super(argc, argv);
  }
  if (argc == 1) {
    obj = nst_field(self, argv[0]);
    if (RTEST(obj)) {
      return obj;
    }
  }
  return cr_call_super(argc, argv);
}

static VALUE
nst_s_new(int argc, VALUE *argv, VALUE klass) {
  VALUE name=Qnil, rest, size;
  VALUE st, members;
  ID id;

  cr_scan_args(argc, argv, "0*", &rest);
  if (RARRAY_LEN(rest) > 0) {
    name = RARRAY_AREF(rest, 0);
    if (!NIL_P(name)) {
      VALUE tmp = cr_check_string_type(name);
      if (!NIL_P(tmp)) {
        cr_any_shift(rest);
      } else {
        name = Qnil;
      }
    }
  }
  if (NIL_P(name)) {
    st = cr_define_class_id(name, klass);
    cr_funcall(klass, cr_intern("inherited"), 1, st);
  }
  else {
    char *cname = StringValuePtr(name);
    id = cr_intern(cname);
    if (!cr_is_const_id(id)) {
      cr_name_error(id, "Identifier %s needs to be constant", cname);
    }
    if (cr_const_defined_at(klass, id)) {
      cr_warn("Redifining constant Struct::%s", cname);
      cr_mod_remove_const(klass, ID2SYM(id));
    }
    st = cr_define_class_under(klass, cr_id2name(id), klass);
  }
  cr_iv_set(st, "__members__", cr_ary_new());
  cr_iv_set(st, "__offset__", INT2FIX(0));

  if (cr_block_given_p()) {
    cr_mod_module_eval(0, 0, st);
  }
  size = cr_iv_get(st, "__offset__");
  members = cr_iv_get(st, "__members__");
  cr_define_const(st, CONTIGUOUS_STRIDE, size);
  cr_define_const(st, ELEMENT_BYTE_SIZE, size);
  cr_define_const(st, ELEMENT_BIT_SIZE, cr_funcall(size, '*', 1, INT2FIX(8)));

  OBJ_FREEZE(members);
  cr_define_const(st, "DEFINITIONS", members);

  cr_define_singleton_method(st, "new", cr_class_new_instance, -1);
  cr_define_method(st, "allocate", nst_allocate, 0);

  return st;
}

static VALUE
nstruct_add_type(VALUE type, int argc, VALUE *argv, VALUE nst) {
  VALUE ofs, size;
  ID id;
  int i;
  VALUE name=Qnil;
  size_t *shape=NULL;
  int ndim=0;
  ssize_t stride;
  narray_view_t *nt;
  int j;

  for (i=0; i<argc; i++) {
    switch(TYPE(argv[i])) {
    case T_STRING:
    case T_SYMBOL:
      if (NIL_P(name)) {
        name = argv[i];
        break;
      }
      cr_raise(cr_eArgument_Error, "Multiple name in struct definition");
    case T_ARRAY:
      if (shape) {
        cr_raise(cr_eArgument_Error, "Multiple shape in struct definition");
      }
      ndim = RARRAY_LEN(argv[i]);
      if (ndim > NA_MAX_DIMENSION) {
        cr_raise(cr_eArgument_Error, "Too large number of dimensions");
      }
      if (ndim == 0) {
        cr_raise(cr_eArgument_Error, "Array is empty");
      }
      shape = ALLOCA_N(size_t, ndim);
      na_array_to_internal_shape(Qnil, argv[i], shape);
      break;
    }
  }
  type = nary_view_new(type, ndim, shape);
  GetNArrayView(type, nt);

  nt->stridx = ALLOC_N(stridx_t, ndim);
  stride = na_dtype_elmsz(CLASS_OF(type));
  for (j=ndim; j--;) {
    SDX_SET_STRIDE(nt->stridx[j], stride);
    stride *= shape[j];
  }
  ofs = cr_iv_get(nst, "__offset__");
  nt->offset = NUM2SIZET(ofs);

  size = cr_funcall(type, cr_intern("byte_size"), 0);
  cr_iv_set(nst, "__offset__", cr_funcall(ofs, '+', 1, size));
  cr_ary_push(cr_iv_get(nst, "__members__"), cr_ary_new3(4, name, type, ofs, size)); // <- field definition
  return Qnil;
}

static VALUE
nst_extract(VALUE self) {
  return self;
}

static void
iter_nstruct_to_a(na_loop_t *const lp) {
  long i, len;
  VALUE obj, types, defs, def;
  VALUE elmt, velm, vary;
  size_t ofs, pos;
  narray_view_t *ne;

  opt = lp->option;
  types = RARRAY_AREF(opt, 0);
  defs = RARRAY_AREF(opt, 1);
  pos = lp->args[0].iter[0].pos;

  len = RARRAY_LEN(types);
  vary = cr_ary_new2(len);

  for (i=0; i<len; i++) {
    def = RARRAY_AREF(defs, i);
    ofs = NUM2SIZET(RARRAY_AREF(def, 2));
    elmt = RARRAY_AREF(types, i);
    GetNArrayView(elmt, ne);
    ne->offset = pos + ofs;
    if (ne->base.ndim == 0) {
      velm = cr_funcall(elmt, cr_intern("to_a"), 0);
    }
    cr_ary_push(vary, velm);
  }
  cr_ary_push(lp->args[1].value, vary);
}

static VALUE
na_original_data(VALUE self) {
  narray_t *na;
  narray_view_t *nv;

  GetNArray(self, na);
  switch(na->type) {
  case NARRAY_VIEW_T:
    GetNArrayView(self, nv);
    return nv->data;
  }
  return self;
}

static VALUE
nst_create_member_views(VALUE self) {
  VALUE defs, def, types, type, elmt;
  long i, len;
  narray_view_t *ne;

  defs = nst_definitions(CLASS_OF(self));
  len = RARRAY_LEN(defs);
  types = cr_ary_new2(len);
  for (i=0; i<len; i++) {
    def = RARRAY_AREF(defs, i);
    type = RARRAY_AREF(def, 1);
    elmt = na_make_view(type);
    cr_ary_push(types, elmt);
    GetNArrayView(elmt, ne);
    ne->data = na_original_data(self);
  }
  return cr_assoc_new(types, defs);
}

static VALUE
nary_struct_to_a(VALUE self) {
  volatile VALUE opt;
  ndfunc_arg_in_t ain[3] = {{Qnil, 0}, {sym_loop_opt}, {sym_option}};
  ndfunc_arg_out_t aout[1] = {{cr_cArray, 0}};
  ndfunc_t ndf = {iter_nstruct_to_a, NO_LOOP, 3, 1, ain, aout};

  opt = nst_create_member_views(self);
  return na_ndloop_cast_narray_to_rarray(&ndf, self, opt);
}

VALUE na_ary_composition_for_struct(VALUE nstruct, VALUE ary);

static void
iter_nstruct_from_a(na_loop_t *const lp) {
  long i, len;
  VALUE ary;
  VALUE types, defs, def;
  VALUE elmt, item;
  size_t ofs;
  narray_view_t *ne;

  types = RARRAY_AREF(lp->option, 0);
  defs = RARRAY_AREF(lp->option, 1);

  len = RARRAY_LEN(types);
  ary = lp->args[1].value;

  for (i=0; i<len; i++) {
    def = RARRAY_AREF(defs, i);
    ofs = NUM2SIZET(RARRAY_AREF(def, 2));
    elmt = RARRAY_AREF(types, i);
    GetNArrayView(elmt, ne);
    ne->offset = lp->args[0].iter[0].pos + ofs;
    item = RARRAY_AREF(ary, i);
    cr_funcall(elmt, cr_intern("store"), 1, item);
  }
}

static VALUE
nary_struct_cast_array(VALUE klass, VALUE rary) {
  VALUE nary;
  narray_t *na;
  VALUE opt;
  ndfunc_arg_in_t ain[3] = {{OVERWRITE, 0}, {cr_cArray, 0}, {sym_option}};
  ndfunc_t ndf = {iter_nstruct_from_a, NO_LOOP, 3, 0, ain, 0};

  nary = na_s_new_like(klass, rary);
  GetNArray(nary, na);

  if (na->size > 0) {
    opt = nst_create_member_views(nary);
    cr_funcall(nary, cr_intern("allocate"), 0);
    na_ndloop_cast_narray_to_rarray2(&ndf, nary, rary, opt);
  }
  return nary;
}

static inline VALUE
nary_struct_s_cast(VALUE klass, VALUE rary) {
  return nary_struct_cast_array(klass, rary);
}

static void
iter_struct_store_struct(na_loop_t *const lp) {
  size_t i, s1, s2;
  char *p1, *p2;
  size_t *idx1, *idx2;
  size_t *elmsz;
  char *x, *y;

  INIT_COUNTER(lp, i);
  INIT_PTR_IDX(lp, 0, p1, s1, idx1);
  INIT_PTR_IDX(lp, 1, p2, s2, idx2);
  INIT_ELMSIZE(lp, 0, elmsz);
  if (idx2) {
    if (idx1) {
      for (; i--;) {
        x = (char*)(p1+*idx1); idx1++;
        y = (char*)(p2+*idx2); idx2++;
        memcpy(x, y, elmsz);
      }
    } else {
      for (; i--;) {
        x = (char*)p1; p1 += s1;
        y = (char*)(p2+*idx2); idx2++;
        memcpy(x, y, elmsz);
      }
    }
  } else {
    if (idx1) {
      for (; i--;) {
        x = (char*)(p1+*idx1); idx1++;
        y = (char*)p2; p2 += s2;
        memcpy(x, y, elmsz);
      }
    } else {
      for (; i--;) {
        x = (char*)p1; p1 += s1;
        y= (char*)p2; p2 += s2;
        memcpy(x, y, elmsz);
      }
    }
  }
}

static VALUE
nary_struct_store_struct(VALUE self, VALUE obj) {
  ndfunc_arg_in_t ain[2] = {{OVERWRITE, 0}, {Qnil, 0}};
  ndfunc_t ndf = {iter_struct_store_struct, FULL_LOOP, 2, 0, ain, 0};

  na_ndloop(&ndf, 2, self, obj);
  return self;
}

static inline VALUE
nary_struct_store_array(VALUE self, VALUE obj) {
  return nary_struct_store_struct(self, nary_struct_cast_array(CLASS_OF(self), obj));
}

// Store elements to Lattice::Struct from other.
static VALUE
nary_struct_store(VALUE self, VALUE obj) {
  if (TYPE(obj) == T_ARRAY) {
    nary_struct_store_array(self, obj);
    return self;
  }
  if (CLASS_OF(self) == CLASS_OF(obj)) {
    nary_struct_store_struct(self, obj);
    return self;
  }
  cr_raise(nary_eCastError, "Unknown conversion from %s to %s", cr_class2name(CLASS_OF(obj)), cr_class2name(CLASS_OF(self)));
  return self;
}

static VALUE
iter_struct_inspect(char *ptr, size_t pos, VALUE opt) {
  VALUE types, defs, def, name, elmt, vary, v, x;
  size_t ofs;
  long i, len;
  narray_view_t *ne;

  types = RARRAY_AREF(opt, 0);
  defs = RARRAY_AREF(opt, 1);

  len = RARRAY_LEN(types);
  vary = cr_ary_new2(len);

  for (i=0; i<len; i++) {
    def = RARRAY_AREF(defs, i);
    name = RARRAY_AREF(def, 0);
    ofs = NUM2SIZET(RARRAY_AREF(def, 2));
    elmt = RARRAY_AREF(types, i);
    GetNArrayView(elmt, ne);
    ne->offset = pos + ofs;
    v = cr_str_concat(cr_sym_to_s(name), cr_str_new2(": "));
    x = cr_funcall(elmt, cr_intern("format_to_a"), 0);
    if (ne->base.ndim == 0) {
      x = cr_funcall(x, cr_intern("first"), 0);
    }
    x = cr_funcall(x, cr_intern("to_s"), 0);
    v = cr_str_concat(v, x);
    cr_ary_push(vary, v);
  }
  v = cr_ary_join(vary, cr_str_new2(", "));
  v = cr_str_concat(cr_str_new2("["), v);
  v = cr_str_concat(v, cr_str_new2("]"));
  return v;
}

// Returns a string containing human-readable representation of NArray.
static VALUE
nary_struct_inspect(VALUE ary) {
  VALUE opt;
  opt = nst_create_member_views(ary);
  return na_ndloop_inspect(ary, iter_struct_inspect, opt);
}

static VALUE
nst_s_add_type(int argc, VALUE *argv, VALUE mod) {
  if (argc == 0) {
    cr_raise(cr_eArgument_Error, "Wrong number of arguments (%d for 1)", argc);
  }
  nstruct_add_type(argv[0], argc - 1, argv + 1, mod);
  return Qnil;
}

#define NST_TYPEDEF(tpname, tpclass)
static VALUE
nst_s_##tpname(VALUE argc, VALUE *argv, VALUE mod) {
  nstruct_add_type(tpclass, argc, argv, mod);
  return Qnil;
}

NST_TYPEDEF(int8, lattice_cInt8)
NST_TYPEDEF(int16, lattice_cInt16)
NST_TYPEDEF(int32, lattice_cInt32)
NST_TYPEDEF(int64, lattice_cInt64)
NST_TYPEDEF(uint8, lattice_cUInt8)
NST_TYPEDEF(uint16, lattice_cUInt16)
NST_TYPEDEF(uint32, lattice_cUInt32)
NST_TYPEDEF(uint64, lattice_cUInt64)
NST_TYPEDEF(dfloat, lattice_cDFloat)
NST_TYPEDEF(dcomplex, lattice_cDComplex)
NST_TYPEDEF(sfloat, lattice_cSFloat)
NST_TYPEDEF(scomplex, lattice_cSComplex)

#define cr_define_singleton_alias(klass, name1, name2) \
  cr_define_alias(cr_singleton_class(klass), name1, name2)

void
Init_nary_struct() {
  cT = cr_define_class_under(mLattice, "Struct", lattice_cNArray);

  cr_define_singleton_method(cT, "new", nst_s_new, -1);
  cr_define_singleton_method(cT, "add_type", nst_s_add_type, -1);
  cr_define_singleton_method(cT, "int8", nst_s_int8, -1);
  cr_define_singleton_method(cT, "int16", nst_s_int16, -1);
  cr_define_singleton_method(cT, "int32", nst_s_int32, -1);
  cr_define_singleton_method(cT, "int64", nst_s_int64, -1);
  cr_define_singleton_method(cT, "uint8", nst_s_uint8, -1);
  cr_define_singleton_method(cT, "uint16", nst_s_uint16, -1);
  cr_define_singleton_method(cT, "uint32", nst_s_uint32, -1);
  cr_define_singleton_method(cT, "uint64", nst_s_uint64, -1);
  cr_define_singleton_method(cT, "sfloat", nst_s_sfloat, -1);
  cr_define_singleton_alias(cT, "float32", "sfloat");
  cr_define_singleton_method(cT, "scomplex", nst_s_scomplex, -1);
  cr_define_singleton_alias(cT, "compex64", "scomplex");
  cr_define_singleton_method(cT, "dfloat", nst_s_dfloat, -1);
  cr_define_singleton_alias(cT, "float64", "dfloat");
  cr_define_singleton_method(cT, "dcomplex", nst_s_dcomplex, -1);
  cr_define_singleton_alias(cT, "complex128", "dcomplex");

  cr_define_method(cT, "to_a", nary_struct_to_a, 0);

  cr_define_method(cT, "store", nary_struct_store, 1);

  cr_define_method(cT, "inspect", nary_struct_inspect, 0);

  cr_define_singleton_method(cT, "cast", nary_struct_s_cast, 1);
  cr_define_singleton_method(cT, "[]", nary_struct_s_cast, -2);
}
