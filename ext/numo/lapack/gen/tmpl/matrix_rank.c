#if !defined DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131E-16
#endif
#if !defined FLT_EPSILON
#define FLT_EPSILON 1.19209290E-07F
#endif

typedef struct {
    double *th_tmp;
    double *th;
} matrix_rank_opt_t;

static void
<%=c_iter%>_2(na_loop_t * const lp)
{
    size_t  m =           lp->args[0].shape[1];
    size_t  n =           lp->args[0].shape[0];
    rtype *a0 = (rtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    VALUE  *b = (VALUE *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    matrix_rank_opt_t *opt = lp->opt_ptr;
    {
        size_t i, j;
        for (j=0; j<n; ++j) {
            double th;
            size_t count;
            rtype *a = a0 + m*j;
            if (opt->th) {
                th = *opt->th;
            } else {
                double mx = a[0];
                for (i=1; i<m; ++i) {
                    if (a[i] > mx) {
                       mx = a[i];
                    }
                }
                th = *opt->th_tmp * mx;
            }
            count = 0;
            for (i=0; i<m; ++i) {
                if (a[i] > th) {
                    ++count;
                }
            }
            b[j] = INT2FIX(count);
        }
    }
}

<% if is_complex %>

static void
<%=c_iter%>_1d_c(na_loop_t * const lp)
{
    size_t   n =             lp->args[0].shape[0];
    dtype   *a =   (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    uint8_t *b = (uint8_t *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    {
        size_t i;
        for (i=0; i<n; ++i) {
            b[i] = ((REAL(a[i]) != 0.0) | (IMAG(a[i]) != 0.0));
        }
    }
}

<% end %>

static void
<%=c_iter%>_1d_r(na_loop_t * const lp)
{
    double th;
    size_t   n =             lp->args[0].shape[0];
    rtype   *a =   (rtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    uint8_t *b = (uint8_t *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    th         = *(double *) lp->opt_ptr;
    {
        size_t i;
        for (i=0; i<n; ++i) {
            b[i] = (a[i] > th);
        }
    }
}

<% if is_complex %>

static VALUE
<%=c_func%>_sub_1d_c(VALUE const a)
{
    volatile VALUE ans;  // !!! DONT REMOVE: volatile qualification !!!
    VALUE mN, u8;
    mN = rb_const_get(rb_cObject, rb_intern("Numo"));
    u8 = rb_const_get(mN, rb_intern("UInt8"));
    {
        narray_t *na;
        size_t shape[1];
        ndfunc_arg_in_t ain[1] = {{cT, 1}};
        ndfunc_arg_out_t aout[1] = {{u8, COUNT_OF_(shape), shape}};
        ndfunc_t ndf = {<%=c_iter%>_1d_c, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};
        GetNArray(a, na);
        shape[0] = na->shape[0];
        ans = na_ndloop3(&ndf, 0, 1, a);
        return INT2FIX( (FIX2LONG(rb_funcall(ans, rb_intern("max"), 0))) ? 1 : 0 );
    }
}

<% end %>

static VALUE
<%=c_func%>_sub_1d_r(VALUE const a, double th, int flg)
{
    volatile VALUE ans;  // !!! DONT REMOVE: volatile qualification !!!
    VALUE mN, u8;
    mN = rb_const_get(rb_cObject, rb_intern("Numo"));
    u8 = rb_const_get(mN, rb_intern("UInt8"));
    {
        narray_t *na;
        size_t shape[1];
        ndfunc_arg_in_t ain[1] = {{cRT, 1}};
        ndfunc_arg_out_t aout[1] = {{u8, COUNT_OF_(shape), shape}};
        ndfunc_t ndf = {<%=c_iter%>_1d_r, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};
        GetNArray(a, na);
        shape[0] = na->shape[0];
        ans = na_ndloop3(&ndf, &th, 1, a);
        if (flg) {
            return INT2FIX( (FIX2LONG(rb_funcall(ans, rb_intern("max"), 0))) ? 1 : 0 );
        } else {
            size_t n;
            GetNArray(ans, na);
            n = na->shape[0];
            if (n < 256) {
                return rb_funcall(ans, rb_intern("sum"), 0);
            } else {
                VALUE id_index;
                size_t i;
                long sum=0;
                id_index = rb_intern("[]");
                for (i=0; i<n; ++i) {
                    sum += FIX2LONG(rb_funcall(ans, id_index, 1, INT2FIX(i)));
                }
                return INT2FIX(sum);
            }
        }
    }
}

/*
  @overload matrix_ran(a)
  @param ***TBD
  @return ***TBD***
  @raise

  ***TBD***
*/
static VALUE
<%=c_func%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    volatile VALUE ans;  // !!! DONT REMOVE: volatile qualification !!!
    VALUE tol = Qnil, a;
    int flg_turbo = 0;
    narray_t *na;

    if (argc > 0) {
        int const last = argc-1;
        VALUE h;
        if ( ! NIL_P(h = rb_check_hash_type(argv[last]))) {
            ID tbl;
            VALUE v;
            tbl = rb_intern("turbo");
            rb_get_kwargs(h, &tbl, 0, 1, &v);
            if (v != Qundef) {
                flg_turbo = RTEST(v);
            }
            --argc;
        }
    }
    rb_check_arity(argc, 1, 2);
    if (argc == 2) {
        tol = argv[1];
    }

    a = argv[0];
    GetNArray(a, na);
    CHECK_DIM_GE(na, 0);

    if (NA_NDIM(na) < 2) {
        if (NA_NDIM(na) == 0) {
            a = rb_funcall(a, rb_intern("reshape"), 1, INT2FIX(1));
        }
        <% if is_complex %>
        ans = <%=c_func%>_sub_1d_c(a);
        <% else %>
        ans = <%=c_func%>_sub_1d_r(a, 0.0, 1);
        <% end %>
    } else {
        double dmax;
        <% if real_ctype == "double" %>
        double const eps = DBL_EPSILON;
        <% elsif real_ctype == "float" %>
        double const eps = FLT_EPSILON;
        <% else %>
        double const eps = NAN;
        <% end %>
        {
            size_t const d = NA_NDIM(na);
            size_t const a = na->shape[d-1];
            size_t const b = na->shape[d-2];
            dmax = (double)max_(a, b);
        }
        {
            VALUE (* const f)(VALUE, svd_job) = (
                flg_turbo
                ? numo_<%=type_name%>_s_gesdd_sub
                : numo_<%=type_name%>_s_gesvd_sub ) ;
            a = (*f)(a, SVD_VALS_ONLY);
        }
        GetNArray(a, na);
        if (NA_NDIM(na) < 2) {
            double th;
            if (tol == Qnil) {
                double smax = NUM2DBL(rb_funcall(a, rb_intern("max"), 0));
                th = smax * dmax * eps;
            } else {
                th = NUM2DBL(tol);
            }
            ans = <%=c_func%>_sub_1d_r(a, th, 0);
        } else {
            VALUE mN, rO;
            mN = rb_const_get(rb_cObject, rb_intern("Numo"));
            rO = rb_const_get(mN, rb_intern("RObject"));
            {
                matrix_rank_opt_t opt;
                double th_tmp;
                size_t shape[1];
                ndfunc_arg_in_t ain[1] = {{cRT, 2}};
                ndfunc_arg_out_t aout[1] = {{rO, COUNT_OF_(shape), shape}};
                ndfunc_t ndf = {<%=c_iter%>_2, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};
                shape[0] = na->shape[NA_NDIM(na)-2];
                if (tol == Qnil) {
                    th_tmp = dmax * eps;
                    opt.th_tmp = &th_tmp;
                    opt.th = 0;
                } else {
                    th_tmp = NUM2DBL(tol);
                    opt.th_tmp = 0;
                    opt.th = &th_tmp;
                }
                ans = na_ndloop3(&ndf, &opt, 1, a);
            }
        }
    }
    return ans;
}
