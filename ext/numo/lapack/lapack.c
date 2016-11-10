#include "ruby.h"
#include "numo/narray.h"

extern void Init_nary_dfloat_lapack();
extern void Init_nary_sfloat_lapack();
extern void Init_nary_dcomplex_lapack();
extern void Init_nary_scomplex_lapack();

/*
  Dispatches method to Numo::LAPACK module of upcasted type,
  eg, Numo::DFloat::LAPACK.
  @overload method_missing(name,x,...)
  @param [Symbol] name  method name.
  @param [NArray,Numeric] x  input array.
  @return [NArray] result.
*/
VALUE nary_lapack_method_missing(int argc, VALUE *argv, VALUE mlapack)
{
    VALUE type, mod, hash;
    if (argc>1) {
	type = nary_mathcast(argc-1,argv+1);
	hash = rb_const_get(mlapack, rb_intern("DISPATCH"));
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
Init_lapack()
{
    VALUE mNumo, hCast, mLAPACK, c;

    mNumo = rb_const_get(rb_cObject, rb_intern("Numo"));
    mLAPACK = rb_define_module_under(mNumo, "LAPACK");
    rb_define_singleton_method(mLAPACK, "method_missing",
                               nary_lapack_method_missing, -1);

    hCast = rb_hash_new();
    rb_define_const(mLAPACK, "DISPATCH", hCast);

    Init_nary_dfloat_lapack();
    Init_nary_sfloat_lapack();
    Init_nary_dcomplex_lapack();
    Init_nary_scomplex_lapack();

    c = numo_cDFloat;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("LAPACK")));
    c = numo_cSFloat;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("LAPACK")));
    c = numo_cDComplex;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("LAPACK")));
    c = numo_cSComplex;
    rb_hash_aset(hCast, c, rb_const_get(c,rb_intern("LAPACK")));
}
