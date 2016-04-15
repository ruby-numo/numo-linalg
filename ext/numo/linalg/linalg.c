#include "ruby.h"
#include "numo/narray.h"

extern void Init_nary_dfloat_linalg();
extern void Init_nary_sfloat_linalg();
extern void Init_nary_dcomplex_linalg();
extern void Init_nary_scomplex_linalg();

/*
  Dispatches method to Numo::Linalg module of upcasted type,
  eg, Numo::DFloat::Linalg.
  @overload method_missing(name,x,...)
  @param [Symbol] name  method name.
  @param [NArray,Numeric] x  input array.
  @return [NArray] result.
*/
VALUE nary_linalg_method_missing(int argc, VALUE *argv, VALUE mlinalg)
{
    VALUE type, mod, hash;
    if (argc>1) {
	type = nary_mathcast(argc-1,argv+1);
	hash = rb_const_get(mlinalg, rb_intern("DISPATCH"));
	mod = rb_hash_aref(hash, type);
	if (NIL_P(mod)) {
	    rb_raise(rb_eTypeError,"%s is unknown for Numo::Math",
		     rb_class2name(type));
	}
	return rb_funcall2(mod,rb_intern("send"),argc,argv);
    }
    rb_raise(rb_eArgError,"argument or method missing");
    return Qnil;
}

void
Init_linalg()
{
    VALUE mNumo, hCast, mLinalg, c;

    mNumo = rb_const_get(rb_cObject, rb_intern("Numo"));
    mLinalg = rb_define_module_under(mNumo, "Linalg");
    rb_define_singleton_method(mLinalg, "method_missing",
                               nary_linalg_method_missing, -1);

    hCast = rb_hash_new();
    rb_define_const(mLinalg, "DISPATCH", hCast);

    Init_nary_dfloat_linalg();
    Init_nary_sfloat_linalg();
    Init_nary_dcomplex_linalg();
    Init_nary_scomplex_linalg();

    c = numo_cDFloat;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("Linalg")));
    c = numo_cSFloat;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("Linalg")));
    c = numo_cDComplex;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("Linalg")));
    c = numo_cSComplex;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("Linalg")));
}
