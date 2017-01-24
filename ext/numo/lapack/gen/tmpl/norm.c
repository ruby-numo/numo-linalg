/*
 * norm
 */

typedef struct {
    size_t outer_n;
    size_t inner_n;
    unsigned axis;
} norm_opt_1d_t;

static void
<%=c_iter%>2(na_loop_t * const lp)
{
    size_t n = lp->args[1].shape[0];
    dtype  *a = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    double *o = (double *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    {
        size_t i;
        for (i = 0; i < n; ++i) {
            dtype a_i = a[i];
            <% if is_complex %>
            o[i] = c_abs_square(a_i);
            <% else %>
            double da_i = (double)a_i;
            o[i] = da_i*da_i;
            <% end %>
        }
    }
}

static inline void
<%=c_iter%>1d_sub(na_loop_t * const lp, size_t const osize
                , double (* const map_func)(dtype), void (* const reduce_func)(void *, VALUE))
{
    norm_opt_1d_t * const opt = (norm_opt_1d_t *)lp->opt_ptr;
    unsigned const ax = opt->axis;
    size_t const ax_n = lp->args[0].shape[ax];
    unsigned const outer_n = opt->outer_n;
    unsigned const inner_n = opt->inner_n;
    size_t i, j, k;
    dtype *a = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    char  *o =            lp->args[1].ptr + lp->args[1].iter[0].pos ;
    {
        size_t na_shape[1] = { ax_n };
        VALUE const na_ptr = rb_narray_new(numo_cDFloat, 1, na_shape);
        double * const ptr = (double *)na_get_pointer_for_write(na_ptr);
        dtype *pa = a, *pa2, *pa3;
        char  *po = o;
        for (i = 0;
             i < outer_n;
             ++i, pa+=ax_n*inner_n) {
            pa2 = pa;
            for (j = 0; j < inner_n; ++j, ++pa2) {
                pa3 = pa2;
                for (k = 0;
                     k < ax_n;
                     ++k, pa3+=inner_n) {

                    ptr[k] = (*map_func)(*pa3);

                }

                (*reduce_func)(po, na_ptr);

                po += osize;
            }
        }
    }
}

static inline double
map_abs_squ(dtype const x)
{
    <% if is_complex %>
    double dxr = REAL(x);
    double dxi = IMAG(x);
    return dxr*dxr + dxi*dxi;
    <% else %>
    double dx = x;
    return dx*dx;
    <% end %>
}
static inline void
reduce_sum_sqrt(void * const o, VALUE const v)
{
    ID const id_kahan_sum = rb_intern("kahan_sum");
    double * const oo = o;
    VALUE sum;
    if (rb_respond_to(v, id_kahan_sum)) {
        sum = rb_funcall(v, id_kahan_sum, 0);
    } else {
        sum = rb_funcall(v, rb_intern("sum"), 0);
    }
    *oo = NUM2DBL(rb_funcall(rb_mMath, rb_intern("sqrt"), 1, sum));
}
static void
<%=c_iter%>1d_2norm(na_loop_t * const lp)
{
    size_t const osize = 8;
    <%=c_iter%>1d_sub(lp, osize, &map_abs_squ, &reduce_sum_sqrt);
}

static inline double
map_abs(dtype const x)
{
    <% if is_complex %>
    double dxr = REAL(x);
    double dxi = IMAG(x);
    return hypot(dxr, dxi);
    <% else %>
    double dx = x;
    return fabs(dx);
    <% end %>
}
static inline void
reduce_max(void * const o, VALUE const v)
{
    double * const oo = o;
    *oo = NUM2DBL(rb_funcall(v, rb_intern("max"), 0));
}
static void
<%=c_iter%>1d_abs_max(na_loop_t * const lp)
{
    size_t const osize = 8;
    <%=c_iter%>1d_sub(lp, osize, &map_abs, &reduce_max);
}

static inline void
reduce_min(void * const o, VALUE const v)
{
    double * const oo = o;
    *oo = NUM2DBL(rb_funcall(v, rb_intern("min"), 0));
}
static void
<%=c_iter%>1d_abs_min(na_loop_t * const lp)
{
    size_t const osize = 8;
    <%=c_iter%>1d_sub(lp, osize, &map_abs, &reduce_min);
}

#if SIZEOF_INT == 2
# define UINTSUM_TYPE numo_cUInt16
#elif SIZEOF_INT == 4
# define UINTSUM_TYPE numo_cUInt32
#elif SIZEOF_INT == 8
# define UINTSUM_TYPE numo_cUInt64
#else
# error
#endif

static inline double
map_notzero(dtype const x)
{
    <% if is_complex %>
    return (double)((REAL(x) != 0.0) || (IMAG(x) != 0.0));
    <% else %>
    return (double)(x != 0.0);
    <% end %>
}
static inline void
reduce_uintsum(void * const o, VALUE const v)
{
    unsigned * const oo = o;
    *oo = NUM2UINT(rb_funcall(v, rb_intern("sum"), 0));
}
static void
<%=c_iter%>1d_notzero_uintsum(na_loop_t * const lp)
{
    size_t const osize = SIZEOF_INT;
    <%=c_iter%>1d_sub(lp, osize, &map_notzero, &reduce_uintsum);
}

static inline void
reduce_sum(void * const o, VALUE const v)
{
    ID const id_kahan_sum = rb_intern("kahan_sum");
    double * const oo = o;
    VALUE sum;
    if (rb_respond_to(v, id_kahan_sum)) {
        sum = rb_funcall(v, id_kahan_sum, 0);
    } else {
        sum = rb_funcall(v, rb_intern("sum"), 0);
    }
    *oo = NUM2DBL(sum);
}
static void
<%=c_iter%>1d_abs_sum(na_loop_t * const lp)
{
    size_t const osize = 8;
    <%=c_iter%>1d_sub(lp, osize, &map_abs, &reduce_sum);
}

static VALUE
norm2_sub(ndfunc_t * const ndf, VALUE const a)
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
norm2(VALUE const a, narray_t * const na)
{
    size_t shape[1];
    unsigned const d = na->ndim;
    ndfunc_arg_in_t ain[1] = {{cT, d}};
    ndfunc_arg_out_t aout[1] = {
      {numo_cDFloat, COUNT_OF_(shape), shape}};
    ndfunc_t ndf = {&<%=c_iter%>2, NO_LOOP,
      COUNT_OF_(ain), COUNT_OF_(aout),
      ain, aout};
    {
        unsigned i;
        size_t n = na->shape[0];
        for (i = 1; i < d; ++i) {
            n *= na->shape[i];
        }
        shape[0] = n;
    }
    return norm2_sub(&ndf, a);
}

static VALUE
norm1d(VALUE const a, narray_t * const na, VALUE const ord, unsigned ax)
{
    VALUE otype = Qundef;
    void (*iter_func)(na_loop_t *);
    VALUE inf;
    // inf, -inf, 0, 1, 2
    if (NIL_P(ord) || (ord == INT2FIX(2))) {
        otype = numo_cDFloat;
        iter_func = &<%=c_iter%>1d_2norm;
    } else if (rb_obj_is_kind_of(ord, rb_cNumeric)) {
        if ( ! NIL_P(inf = rb_funcall(ord, rb_intern("infinite?"), 0))) {
            if (FIX2INT(inf) == 1) {
                otype = numo_cDFloat;
                iter_func = &<%=c_iter%>1d_abs_max;
            } else if (FIX2INT(inf) == -1) {
                otype = numo_cDFloat;
                iter_func = &<%=c_iter%>1d_abs_min;
            } else {
                assert(0);
            }
        } else {
            if (ord == INT2FIX(0)) {
                otype = UINTSUM_TYPE;
                iter_func = &<%=c_iter%>1d_notzero_uintsum;
            } else if (ord == INT2FIX(1)) {
                otype = numo_cDFloat;
                iter_func = &<%=c_iter%>1d_abs_sum;
            } else {
                rb_notimplement();
            }
        }
    } else {
        rb_notimplement();
    }
    {
        unsigned const d = na->ndim;
        size_t * const shape = ALLOCA_N(size_t, d-1);
        ndfunc_arg_in_t ain[1] = {{cT, d}};
        ndfunc_arg_out_t aout[1] = {
          {otype, d-1, shape}};
        ndfunc_t ndf = {iter_func, NO_LOOP,
          COUNT_OF_(ain), COUNT_OF_(aout),
          ain, aout};
        {
            size_t outer_n = 1, inner_n = 1;
            unsigned i, j;
            for (i = 0, j = 0; i < ax; ++i, ++j) {
                size_t const n = na->shape[i];
                shape[j] = n;
                outer_n *= n;
            }
            ++i;
            for (; i < d; ++i, ++j) {
                size_t const n = na->shape[i];
                shape[j] = n;
                inner_n *= n;
            }
            {
                norm_opt_1d_t opt = {
                    .outer_n = outer_n
                  , .inner_n = inner_n
                  , .axis    = ax
                };
                return na_ndloop3(&ndf, &opt, 1, a);
            }
        }
    }
}

#undef UINTSUM_TYPE

static VALUE
norm2d(VALUE const a, narray_t * const na, VALUE const ord, unsigned const ax0, unsigned const ax1)
{
    return Qnil;
}

// inf, -inf, 0, 1, 2

static VALUE
<%=c_func%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    VALUE a;
    narray_t *na;
    unsigned ndim;
    VALUE ord = Qnil, axis = Qnil;
    {
        VALUE const h = rb_check_hash_type(argv[argc-1]);
        if ( ! NIL_P(h)) {
            --argc;
        }
        rb_check_arity(argc, 1, 1);

        if ( ! NIL_P(h)) {
            ID table[2];
            VALUE values[2];
            table[0] = rb_intern("ord");
            table[1] = rb_intern("axis");
            rb_get_kwargs(h, table, 0, 2, values);
            if (values[0] != Qundef) {
                ord = values[0];
            }
            if (values[1] != Qundef) {
                axis = values[1];
            }
        }
    }
    rb_check_arity(argc, 1, 1);

    a = argv[0];
    GetNArray(a, na);
    ndim = na->ndim;

    if (NIL_P(axis)) {
        if ( NIL_P(ord)
          || ((ord == ID2SYM(rb_intern("FRO"))) && (ndim == 2))
          || ((ord == INT2FIX(2)) && (ndim == 1)) ) {
            return norm2(a, na);
        }
        axis = rb_ary_new_capa(ndim);
        {
            unsigned i;
            for (i = 0; i < ndim; ++i) {
                rb_ary_store(axis, i, UINT2NUM(i));
            }
        }
    } else if (RB_INTEGER_TYPE_P(axis)) {
        axis = rb_ary_new_from_args(1, axis);
        // Does need checking for Integer-like object ?
    } else if ( ! RB_TYPE_P(axis, T_ARRAY)) {
        /*
        VALUE tmp = rb_check_array_type(axis);
        if (NIL_P(tmp)) {
            tmp = rb_check_convert_type(axis, T_ARRAY, "Array", "to_a");
            if (NIL_P(tmp)) {
                rb_raise(rb_eTypeError, "axis: must be an Array-like of Integers or an Integer or nil");
            }
        }
        axis = tmp;
        */
        rb_raise(rb_eTypeError, "axis: must be an Array of Integers or an Integer or nil");
    }

    assert(RB_TYPE_P(axis, T_ARRAY));

    if ((RARRAY_LEN(axis) < 1) || (RARRAY_LEN(axis) > 2)) {
        rb_raise(rb_eArgError, "Improper number of dimensions to norm.");
    }
    assert((RARRAY_LEN(axis) >= 1) && (RARRAY_LEN(axis) <= 2));

    if ( ! RB_INTEGER_TYPE_P(rb_ary_entry(axis, 0))) {
        rb_raise(rb_eTypeError, "axis: must be an Array of Integers or an Integer or nil");
    }
    switch (RARRAY_LEN(axis)) {
    case 1:
        {
            int ax = NUM2INT(rb_ary_entry(axis, 0));
            unsigned uax;
            if (ax < 0) {
                ax += ndim;
            }
            uax = ax;
            if (uax >= ndim) {
                rb_raise(rb_eArgError, "Invalid axis");
            }
            return norm1d(a, na, ord, uax);
        }
    case 2:
        if ( ! RB_INTEGER_TYPE_P(rb_ary_entry(axis, 1))) {
            rb_raise(rb_eTypeError, "axis: must be an Array of Integers or an Integer or nil");
        }
        {
            int ax0 = NUM2INT(rb_ary_entry(axis, 0));
            int ax1 = NUM2INT(rb_ary_entry(axis, 1));
            unsigned uax0, uax1;
            if (ax0 < 0) {
                ax0 += ndim;
            }
            if (ax1 < 0) {
                ax1 += ndim;
            }
            uax0 = ax0; uax1 = ax1;
            if ((uax0 >= ndim) || (uax1 >= ndim)) {
                rb_raise(rb_eArgError, "Invalid axis");
            }
            if (uax0 == uax1) {
                rb_raise(rb_eArgError, "Duplicate axes given.");
            }
            return norm2d(a, na, ord, uax0, uax1);
        }
    default:
        assert(0);
    }
    assert(0);
}
