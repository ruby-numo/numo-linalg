void
Init_dummy_for_yard(void)
{
    VALUE const mNumo   = rb_define_module("Numo");
    VALUE const cNArray = rb_define_class_under(mNumo, "NArray");
    VALUE const cDFloat   = rb_define_class_under(mNumo, "DFloat",   cNArray);
    VALUE const cDComplex = rb_define_class_under(mNumo, "DComplex", cNArray);
    VALUE const cSFloat   = rb_define_class_under(mNumo, "SFloat",   cNArray);
    VALUE const cSComplex = rb_define_class_under(mNumo, "SComplex", cNArray);
}
