static void
<%=c_iter%>(na_loop_t * const lp)
{
    size_t n = lp->args[1].shape[0];
    dtype *a  = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    double *o = (double *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    {
        size_t i;
        for (i = 0; i < n; ++i) {
            dtype a_i = a[i];
            <% if is_complex %>
            o[i] = c_abs_square(a_i);
            <% else %>
            o[i] = a_i*a_i;
            <% end %>
        }
    }
}

static VALUE
norm_sub(ndfunc_t * const ndf, VALUE const a)
{
    ID const id_kahan_sum = rb_intern("kahan_sum");

    VALUE sum;
    VALUE const tmp = na_ndloop3(ndf, 0, 1, a);

    if (rb_respond_to(tmp, id_kahan_sum)) {
        sum = rb_funcall(tmp, id_kahan_sum, 0);
    } else {
        sum = rb_funcall(tmp, rb_intern("sum"), 0);
    }

    return rb_funcall(rb_mMath, rb_intern("sqrt"), 1, sum);
}

static VALUE
norm1d(VALUE const a, narray_t * const na)
{
    size_t shape[1];

    ndfunc_arg_in_t ain[1] = {{cT, 1}};
    ndfunc_arg_out_t aout[1] = {
      {numo_cDFloat, COUNT_OF_(shape), shape}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
      COUNT_OF_(ain), COUNT_OF_(aout),
      ain, aout};
    shape[0] = na->shape[0];

    return norm_sub(&ndf, a);
}

static VALUE
norm2d(VALUE const a, narray_t * const na)
{
    size_t shape[1];

    ndfunc_arg_in_t ain[1] = {{cT, 2}};
    ndfunc_arg_out_t aout[1] = {
      {numo_cDFloat, COUNT_OF_(shape), shape}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
      COUNT_OF_(ain), COUNT_OF_(aout),
      ain, aout};
    shape[0] = na->shape[1] * na->shape[0];

    return norm_sub(&ndf, a);
}

static VALUE
<%=c_func%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    VALUE a;
    narray_t *na;
#if 0
    VALUE ord = Qnil;
    {
        VALUE const h = rb_check_hash_type(argv[argc-1]);
        if ( ! NIL_P(h)) {
            --argc;
        }
        rb_check_arity(argc, 1, 1);

        if ( ! NIL_P(h)) {
            ID table[1];
            VALUE values[1];
            table[0] = rb_intern("ord");
            rb_get_kwargs(h, table, 0, 1, values);
            if (values[0] != Qundef) {
                ord = values[0];
            }
        }
    }
#endif
    rb_check_arity(argc, 1, 1);
    a = argv[0];
    GetNArray(a, na);

    switch (na->ndim) {
    case 1:
        return norm1d(a, na);
    case 2:
        return norm2d(a, na);
    default:
        rb_raise(nary_eDimensionError, "Improper number of dimensions to norm.");
    }
}
