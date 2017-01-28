/*
 * norm
 */

typedef struct {
    size_t outer_n, inner_n;
    unsigned axis;
    int extra;  /* extra field for ord > 2 */
} norm_opt_1d_t;

typedef struct {
    size_t outer_n, mid_n, inner_n;
    unsigned axis0, axis1;
} norm_opt_2d_t;

static ID id_kahan_sum_or_sum = (ID)0xDEADBEAF;

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
<%=c_iter%>1dg_sub(na_loop_t * const lp
                 , double (* const map_func)(dtype, void *)
                 , void (* const reduce_func)(void *, VALUE, void *))
{
    norm_opt_1d_t * const opt = (norm_opt_1d_t *)lp->opt_ptr;
    unsigned const ax = opt->axis;
    size_t const ax_n = lp->args[0].shape[ax];
    unsigned const outer_n = opt->outer_n;
    unsigned const inner_n = opt->inner_n;
    dtype * const a = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    char  * const o =            lp->args[1].ptr + lp->args[1].iter[0].pos ;
    size_t const osize =         lp->args[1].elmsz;

    size_t na_shape[1] = { ax_n };
    VALUE const na_ptr = rb_narray_new(numo_cDFloat, COUNT_OF_(na_shape), na_shape);
    double * const ptr = (double *)na_get_pointer_for_write(na_ptr);
    {
        size_t i;
        dtype *pa = a;
        char  *po = o;
        for (i = 0; i < outer_n; ++i, pa+=ax_n*inner_n) {
            size_t j;
            dtype *pa2 = pa;
            for (j = 0; j < inner_n; ++j, ++pa2) {
                size_t k;
                dtype *pa3 = pa2;
                double *tmp = ptr;
                for (k = 0; k < ax_n; ++k, pa3+=inner_n) {

                    *tmp = (*map_func)(*pa3, lp->opt_ptr);
                    ++tmp;

                }
                (*reduce_func)(po, na_ptr, lp->opt_ptr);
                po += osize;
            }
        }
    }
}

static inline void
<%=c_iter%>1ds_sub(na_loop_t * const lp
                 , double (* const map_func)(dtype, void *)
                 , void (* const reduce_func)(void *, VALUE, void *))
{
    size_t const n = lp->args[0].shape[0];
    dtype *a = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    char  *o =            lp->args[1].ptr + lp->args[1].iter[0].pos ;

    size_t na_shape[1] = { n };
    VALUE const na_ptr = rb_narray_new(numo_cDFloat, 1, na_shape);
    double * const ptr = (double *)na_get_pointer_for_write(na_ptr);
    {
        size_t i;
        for (i = 0; i < n; ++i) {
            ptr[i] = (*map_func)(a[i], lp->opt_ptr);
        }
    }
    (*reduce_func)(o, na_ptr, lp->opt_ptr);
}

static inline void
<%=c_iter%>2dg_sub(na_loop_t * const lp
                 , double (* const map_func)(dtype, void *)
                 , double (* const mid_func)(VALUE)
                 , void (* const reduce_func)(void *, VALUE, void *)
                 , int reverse)
{
    norm_opt_2d_t * const opt = (norm_opt_2d_t *)lp->opt_ptr;
    unsigned const ax0 = opt->axis0;
    unsigned const ax1 = opt->axis1;
    size_t const ax0_n = lp->args[0].shape[ax0];
    size_t const ax1_n = lp->args[0].shape[ax1];
    unsigned const outer_n = opt->outer_n;
    unsigned const mid_n   = opt->mid_n;
    unsigned const inner_n = opt->inner_n;
    dtype * const a = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    char  * const o =            lp->args[1].ptr + lp->args[1].iter[0].pos ;
    size_t const osize =         lp->args[1].elmsz;

    size_t na_shape0[1] = { ax0_n };
    size_t na_shape1[1] = { ax1_n };
    VALUE const na_ptr0 = rb_narray_new(numo_cDFloat, COUNT_OF_(na_shape0), na_shape0);
    VALUE const na_ptr1 = rb_narray_new(numo_cDFloat, COUNT_OF_(na_shape1), na_shape1);
    double * const ptr0 = (double *)na_get_pointer_for_write(na_ptr0);
    double * const ptr1 = (double *)na_get_pointer_for_write(na_ptr1);
    {
        size_t i;
        dtype *pa = a;
        char  *po = o;
        size_t step_1, step_2, step_4, step_5;
        if ( ! reverse) {
            step_5 = inner_n;
            step_2 = ax1_n*step_5;
            step_4 = mid_n*step_2;
            step_1 = ax0_n*step_4;
        } else {
            step_4 = inner_n;
            step_2 = ax0_n*step_4;
            step_5 = mid_n*step_2;
            step_1 = ax1_n*step_5;
        }
        for (i = 0; i < outer_n; ++i, pa+=step_1) {
            size_t j;
            dtype *pa2 = pa;
            for (j = 0; j < mid_n; ++j, pa2+=step_2) {
                size_t k;
                dtype *pa3 = pa2;
                for (k = 0; k < inner_n; ++k, ++pa3) {
                    size_t l;
                    dtype *pa4 = pa3;
                    double *tmp0 = ptr0;
                    for (l = 0; l < ax0_n; ++l, pa4+=step_4) {
                        size_t m;
                        dtype *pa5 = pa4;
                        double *tmp1 = ptr1;
                        for (m = 0; m < ax1_n; ++m, pa5+=step_5) {

                            *tmp1 = (*map_func)(*pa5, 0);
                            ++tmp1;

                        }
                        *tmp0 = (*mid_func)(na_ptr1);
                        ++tmp0;
                    }
                    (*reduce_func)(po, na_ptr0, 0);
                    po += osize;
                }
            }
        }
    }
}

static inline void
<%=c_iter%>2ds_sub(na_loop_t * const lp
                 , double (* const map_func)(dtype, void *)
                 , double (* const mid_func)(VALUE)
                 , void (* const reduce_func)(void *, VALUE, void *)
                 , int reverse)
{
    unsigned const ax0 = reverse;
    unsigned const ax1 = !reverse;
    size_t const ax0_n = lp->args[0].shape[ax0];
    size_t const ax1_n = lp->args[0].shape[ax1];
    dtype * const a = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    char  * const o =            lp->args[1].ptr + lp->args[1].iter[0].pos ;

    size_t na_shape0[1] = { ax0_n };
    size_t na_shape1[1] = { ax1_n };
    VALUE const na_ptr0 = rb_narray_new(numo_cDFloat, COUNT_OF_(na_shape0), na_shape0);
    VALUE const na_ptr1 = rb_narray_new(numo_cDFloat, COUNT_OF_(na_shape1), na_shape1);
    double * const ptr0 = (double *)na_get_pointer_for_write(na_ptr0);
    double * const ptr1 = (double *)na_get_pointer_for_write(na_ptr1);
    {
        size_t i;
        dtype *pa = a;
        size_t step_1, step_2;
        if ( ! reverse) {
            step_1 = ax1_n;
            step_2 = 1;
        } else {
            step_1 = 1;
            step_2 = ax0_n;
        }
        double *tmp0 = ptr0;
        for (i = 0; i < ax0_n; ++i, pa+=step_1) {
            size_t j;
            dtype *pa2 = pa;
            double *tmp1 = ptr1;
            for (j = 0; j < ax1_n; ++j, pa2+=step_2) {

                *tmp1 = (*map_func)(*pa2, 0);
                ++tmp1;

            }
            *tmp0 = (*mid_func)(na_ptr1);
            ++tmp0;
        }
        (*reduce_func)(o, na_ptr0, 0);
    }
}

static inline void
<%=c_iter%>2dg_svd(na_loop_t * const lp
                 , void (* const reduce_func)(void *, VALUE, void *)
                 , int reverse)
{
    norm_opt_2d_t * const opt = (norm_opt_2d_t *)lp->opt_ptr;
    unsigned const ax0 = opt->axis0;
    unsigned const ax1 = opt->axis1;
    size_t const ax0_n = lp->args[0].shape[ax0];
    size_t const ax1_n = lp->args[0].shape[ax1];
    unsigned const outer_n = opt->outer_n;
    unsigned const mid_n   = opt->mid_n;
    unsigned const inner_n = opt->inner_n;
    dtype * const a = (dtype  *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    char  * const o =            lp->args[1].ptr + lp->args[1].iter[0].pos ;
    size_t const osize =         lp->args[1].elmsz;

    size_t na_shape[2] = { ax0_n, ax1_n };
    <% if is_complex %>
    VALUE const na_ptr = rb_narray_new(numo_cDComplex, COUNT_OF_(na_shape), na_shape);
    dcomplex * const ptr = (dcomplex *)na_get_pointer_for_write(na_ptr);
    <% else %>
    VALUE const na_ptr = rb_narray_new(numo_cDFloat, COUNT_OF_(na_shape), na_shape);
    double * const ptr = (double *)na_get_pointer_for_write(na_ptr);
    <% end %>
    {
        size_t i;
        dtype *pa = a;
        char  *po = o;
        size_t step_1, step_2, step_4, step_5;
        if ( ! reverse) {
            step_5 = inner_n;
            step_2 = ax1_n*step_5;
            step_4 = mid_n*step_2;
            step_1 = ax0_n*step_4;
        } else {
            step_4 = inner_n;
            step_2 = ax0_n*step_4;
            step_5 = mid_n*step_2;
            step_1 = ax1_n*step_5;
        }
        for (i = 0; i < outer_n; ++i, pa+=step_1) {
            size_t j;
            dtype *pa2 = pa;
            for (j = 0; j < mid_n; ++j, pa2+=step_2) {
                size_t k;
                dtype *pa3 = pa2;
                for (k = 0; k < inner_n; ++k, ++pa3) {
                    size_t l;
                    dtype *pa4 = pa3;
                    <% if is_complex %>
                    dcomplex *tmp = ptr;
                    <% else %>
                    double *tmp = ptr;
                    <% end %>
                    for (l = 0; l < ax0_n; ++l, pa4+=step_4) {
                        size_t m;
                        dtype *pa5 = pa4;
                        for (m = 0; m < ax1_n; ++m, pa5+=step_5) {
                            <% if is_complex %>
                            dtype x = *pa5;
                            REAL(*tmp) = REAL(x);
                            IMAG(*tmp) = IMAG(x);
                            <% else %>
                            *tmp = *pa5;
                            <% end %>
                            ++tmp;
                        }
                    }
                    {
                        VALUE const kwargs = rb_hash_new();
                        VALUE const mLAPACK = rb_const_get(mNumo, rb_intern("LAPACK"));
                        rb_hash_aset(kwargs, ID2SYM(rb_intern("job")), ID2SYM(rb_intern("VALS_ONLY")));
                        (*reduce_func)(po
                                     , rb_funcall(mLAPACK, rb_intern("gesvd"), 2, na_ptr, kwargs)
                                     , 0);
                    }
                    po += osize;
                }
            }
        }
    }
}

static void
<%=c_iter%>2ds_svd_special(na_loop_t *lp
                         , void (* const reduce_func)(void *, VALUE, void *))
{
    VALUE const kwargs = rb_hash_new();
    VALUE const argv[] = { lp->option, kwargs } ;
    rb_hash_aset(kwargs, ID2SYM(rb_intern("job")), ID2SYM(rb_intern("VALS_ONLY")));
    (*reduce_func)(&(lp->option)
                 , numo_<%=type_name%>_s_gesvd(COUNT_OF_(argv), argv, Qnil)
                 , 0);
}

static inline double
map_abs_squ(dtype const x, void *p)
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
reduce_sum_sqrt(void * const o, VALUE const v, void *p)
{
    double * const oo = o;
    VALUE const sum = rb_funcall(v, id_kahan_sum_or_sum, 0);
    *oo = NUM2DBL(rb_funcall(rb_mMath, rb_intern("sqrt"), 1, sum));
}
static inline double
mid_sum(VALUE const v)
{
    return NUM2DBL(rb_funcall(v, id_kahan_sum_or_sum, 0));
}
static void
<%=c_iter%>1dg_2norm(na_loop_t * const lp)
{
    <%=c_iter%>1dg_sub(lp, &map_abs_squ, &reduce_sum_sqrt);
}
static void
<%=c_iter%>1ds_2norm(na_loop_t * const lp)
{
    <%=c_iter%>1ds_sub(lp, &map_abs_squ, &reduce_sum_sqrt);
}
static void
<%=c_iter%>2dg_fro(na_loop_t * const lp)
{
    <%=c_iter%>2dg_sub(lp, &map_abs_squ, &mid_sum, &reduce_sum_sqrt, 0);
}
static void
<%=c_iter%>2dg_fro_x(na_loop_t * const lp)
{
    <%=c_iter%>2dg_sub(lp, &map_abs_squ, &mid_sum, &reduce_sum_sqrt, 1);
}
static void
<%=c_iter%>2ds_fro(na_loop_t * const lp)
{
    <%=c_iter%>2ds_sub(lp, &map_abs_squ, &mid_sum, &reduce_sum_sqrt, 0);
}
static void
<%=c_iter%>2ds_fro_x(na_loop_t * const lp)
{
    <%=c_iter%>2ds_sub(lp, &map_abs_squ, &mid_sum, &reduce_sum_sqrt, 1);
}

static inline double
map_abs(dtype const x, void *p)
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
reduce_max_dtype(void * const o, VALUE const v, void *p)
{
    <% if is_complex %>
    double * const oo = o;
    <% else %>
    dtype * const oo = o;
    <% end %>
    *oo = NUM2DBL(rb_funcall(v, rb_intern("max"), 0));
}
static void
<%=c_iter%>1dg_abs_max(na_loop_t * const lp)
{
    <%=c_iter%>1dg_sub(lp, &map_abs, &reduce_max_dtype);
}
static void
<%=c_iter%>1ds_abs_max(na_loop_t * const lp)
{
    <%=c_iter%>1ds_sub(lp, &map_abs, &reduce_max_dtype);
}
static void
<%=c_iter%>2dg_abs_sum_max(na_loop_t * const lp)
{
    <%=c_iter%>2dg_sub(lp, &map_abs, &mid_sum, &reduce_max_dtype, 0);
}
static void
<%=c_iter%>2dg_abs_sum_max_x(na_loop_t * const lp)
{
    <%=c_iter%>2dg_sub(lp, &map_abs, &mid_sum, &reduce_max_dtype, 1);
}
static void
<%=c_iter%>2ds_abs_sum_max(na_loop_t * const lp)
{
    <%=c_iter%>2ds_sub(lp, &map_abs, &mid_sum, &reduce_max_dtype, 0);
}
static void
<%=c_iter%>2ds_abs_sum_max_x(na_loop_t * const lp)
{
    <%=c_iter%>2ds_sub(lp, &map_abs, &mid_sum, &reduce_max_dtype, 1);
}

static inline void
reduce_min_dtype(void * const o, VALUE const v, void *p)
{
    <% if is_complex %>
    double * const oo = o;
    <% else %>
    dtype * const oo = o;
    <% end %>
    *oo = NUM2DBL(rb_funcall(v, rb_intern("min"), 0));
}
static void
<%=c_iter%>1dg_abs_min(na_loop_t * const lp)
{
    <%=c_iter%>1dg_sub(lp, &map_abs, &reduce_min_dtype);
}
static void
<%=c_iter%>1ds_abs_min(na_loop_t * const lp)
{
    <%=c_iter%>1ds_sub(lp, &map_abs, &reduce_min_dtype);
}
static void
<%=c_iter%>2dg_abs_sum_min(na_loop_t * const lp)
{
    <%=c_iter%>2dg_sub(lp, &map_abs, &mid_sum, &reduce_min_dtype, 0);
}
static void
<%=c_iter%>2dg_abs_sum_min_x(na_loop_t * const lp)
{
    <%=c_iter%>2dg_sub(lp, &map_abs, &mid_sum, &reduce_min_dtype, 1);
}
static void
<%=c_iter%>2ds_abs_sum_min(na_loop_t * const lp)
{
    <%=c_iter%>2ds_sub(lp, &map_abs, &mid_sum, &reduce_min_dtype, 0);
}
static void
<%=c_iter%>2ds_abs_sum_min_x(na_loop_t * const lp)
{
    <%=c_iter%>2ds_sub(lp, &map_abs, &mid_sum, &reduce_min_dtype, 1);
}

static inline double
map_notzero(dtype const x, void *p)
{
    <% if is_complex %>
    return (double)((REAL(x) != 0.0) || (IMAG(x) != 0.0));
    <% else %>
    return (double)(x != 0.0);
    <% end %>
}
static inline void
reduce_uintsum(void * const o, VALUE const v, void *p)
{
    unsigned * const oo = o;
    *oo = NUM2UINT(rb_funcall(v, id_kahan_sum_or_sum, 0));
}
static void
<%=c_iter%>1dg_notzero_uintsum(na_loop_t * const lp)
{
    <%=c_iter%>1dg_sub(lp, &map_notzero, &reduce_uintsum);
}
static void
<%=c_iter%>1ds_notzero_uintsum(na_loop_t * const lp)
{
    <%=c_iter%>1ds_sub(lp, &map_notzero, &reduce_uintsum);
}

static inline void
reduce_sum_dtype(void * const o, VALUE const v, void *p)
{
    <% if is_complex %>
    double * const oo = o;
    <% else %>
    dtype * const oo = o;
    <% end %>
    VALUE const sum = rb_funcall(v, id_kahan_sum_or_sum, 0);
    *oo = NUM2DBL(sum);
}
static void
<%=c_iter%>1dg_abs_sum(na_loop_t * const lp)
{
    <%=c_iter%>1dg_sub(lp, &map_abs, &reduce_sum_dtype);
}
static void
<%=c_iter%>1ds_abs_sum(na_loop_t * const lp)
{
    <%=c_iter%>1ds_sub(lp, &map_abs, &reduce_sum_dtype);
}

static inline void
reduce_sum_double(void * const o, VALUE const v, void *p)
{
    double * const oo = o;
    *oo = NUM2DBL(rb_funcall(v, id_kahan_sum_or_sum, 0));
}

static inline void
reduce_sum_value(void * const o, VALUE const v, void *p)
{
    VALUE * const oo = o;
    *oo = rb_funcall(v, id_kahan_sum_or_sum, 0);
}

static inline void
reduce_max_double(void * const o, VALUE const v, void *p)
{
    double * const oo = o;
    *oo = NUM2DBL(rb_funcall(v, rb_intern("max"), 0));
}

static inline void
reduce_max_value(void * const o, VALUE const v, void *p)
{
    VALUE * const oo = o;
    *oo = rb_funcall(v, rb_intern("max"), 0);
}

static inline void
reduce_min_double(void * const o, VALUE const v, void *p)
{
    double * const oo = o;
    *oo = NUM2DBL(rb_funcall(v, rb_intern("min"), 0));
}

static inline void
reduce_min_value(void * const o, VALUE const v, void *p)
{
    VALUE * const oo = o;
    *oo = rb_funcall(v, rb_intern("min"), 0);
}

static void
<%=c_iter%>2dg_nuc(na_loop_t * const lp)
{
    <%=c_iter%>2dg_svd(lp, &reduce_sum_double, 0);
}
static void
<%=c_iter%>2dg_nuc_x(na_loop_t * const lp)
{
    <%=c_iter%>2dg_svd(lp, &reduce_sum_double, 1);
}
static void
<%=c_iter%>2ds_nuc(na_loop_t * const lp)
{
    <%=c_iter%>2ds_svd_special(lp, &reduce_sum_value);
}

static void
<%=c_iter%>2dg_svd_max(na_loop_t * const lp)
{
    <%=c_iter%>2dg_svd(lp, &reduce_max_double, 0);
}
static void
<%=c_iter%>2dg_svd_max_x(na_loop_t * const lp)
{
    <%=c_iter%>2dg_svd(lp, &reduce_max_double, 1);
}
static void
<%=c_iter%>2ds_svd_max(na_loop_t * const lp)
{
    <%=c_iter%>2ds_svd_special(lp, &reduce_max_value);
}

static void
<%=c_iter%>2dg_svd_min(na_loop_t * const lp)
{
    <%=c_iter%>2dg_svd(lp, &reduce_min_double, 0);
}
static void
<%=c_iter%>2dg_svd_min_x(na_loop_t * const lp)
{
    <%=c_iter%>2dg_svd(lp, &reduce_min_double, 1);
}
static void
<%=c_iter%>2ds_svd_min(na_loop_t * const lp)
{
    <%=c_iter%>2ds_svd_special(lp, &reduce_min_value);
}

static inline double
map_extra(dtype const x, void * const p)
{
    norm_opt_1d_t * const opt = p;
    <% if is_complex %>
    double dxr = REAL(x);
    double dxi = IMAG(x);
    double tmp = hypot(dxr, dxi);
    <% else %>
    double dx = x;
    double tmp = fabs(dx);
    <% end %>
    return pow(tmp, opt->extra);
}
static inline void
reduce_extra(void * const o, VALUE const v, void * const p)
{
    double * const oo = o;
    norm_opt_1d_t * const opt = p;
    VALUE const tmp = rb_funcall(v, id_kahan_sum_or_sum, 0);
    *oo = NUM2DBL(rb_funcall(tmp, rb_intern("**"), 1, DBL2NUM(1.0/(double)opt->extra)));
}
static void
<%=c_iter%>1dg_extra(na_loop_t * const lp)
{
    <%=c_iter%>1dg_sub(lp, &map_extra, &reduce_extra);
}
static void
<%=c_iter%>1ds_extra(na_loop_t * const lp)
{
    <%=c_iter%>1ds_sub(lp, &map_extra, &reduce_extra);
}

static VALUE
norm2dg_svd_optimized(VALUE const a, narray_t * const na, void (*iter_func)(na_loop_t *))
{
    unsigned const d_2 = na->ndim - 2;
    VALUE const na_ptr = rb_narray_new(cT, d_2, na->shape);
    rtype * const o = (rtype *)na_get_pointer_for_write(na_ptr);
    void (*reduce_func)(void *, VALUE, void *);
    VALUE * const indices = ALLOCA_N(VALUE, d_2 + 1);
    int dim;  // this must be signed

    if ((iter_func == <%=c_iter%>2dg_nuc) || (iter_func == <%=c_iter%>2dg_nuc_x)) {
        reduce_func = &reduce_sum_dtype;
    } else if ((iter_func == <%=c_iter%>2dg_svd_max) || (iter_func == <%=c_iter%>2dg_svd_max_x))  {
        reduce_func = &reduce_max_dtype;
    } else if ((iter_func == <%=c_iter%>2dg_svd_min) || (iter_func == <%=c_iter%>2dg_svd_min_x))  {
        reduce_func = &reduce_min_dtype;
    } else {
        assert(0);
    }

    for (dim = 0; dim < (int)d_2; ++dim) {
        indices[dim] = INT2FIX(0);
    }
    indices[d_2] = Qfalse;

    rtype *po = o;
    do {
        VALUE const view = rb_funcallv(a, rb_intern("[]"), d_2 + 1, indices);
        VALUE const kwargs = rb_hash_new();
        VALUE const argv[] = { view, kwargs } ;
        rb_hash_aset(kwargs, ID2SYM(rb_intern("job")), ID2SYM(rb_intern("VALS_ONLY")));
        (*reduce_func)(po
                     , numo_<%=type_name%>_s_gesvd(COUNT_OF_(argv), argv, Qnil)
                     , 0);
        ++po;
        for (dim = d_2 - 1; dim >= 0; --dim) {
            size_t const tmp = FIX2LONG(indices[dim]) + 1;
            if (tmp < na->shape[dim]) {
                indices[dim] = LONG2FIX(tmp);
                break;
            }
            indices[dim] = INT2FIX(0);
        }
    } while (dim >= 0);

    return na_ptr;
}

static VALUE
norm2_sub(ndfunc_t * const ndf, VALUE const a)
{
    VALUE const tmp = na_ndloop3(ndf, 0, 1, a);
    VALUE const sum = rb_funcall(tmp, id_kahan_sum_or_sum, 0);
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

static int
norm1d_ord_idx(VALUE const ord)
{
    if (NIL_P(ord)) {
        return 2;
    } else if (rb_obj_is_kind_of(ord, rb_cNumeric)) {
        VALUE inf;
        if ( ! NIL_P(inf = rb_funcall(ord, rb_intern("infinite?"), 0))) {
            if (FIX2LONG(inf) == 1) {
                return -1;
            } else if (FIX2LONG(inf) == -1) {
                return -2;
            }
            assert(0);
        } else if (FIXNUM_P(ord)) {
            int idx;
            if ((idx = FIX2LONG(ord)) >= 0) {
                return idx;
            } else {
                rb_raise(rb_eArgError, "Invalid norm order for vectors.");
            }
            assert(0);
        } else {
            rb_raise(rb_eArgError, "Invalid norm order for vectors.");
        }
        assert(0);
    } else {
        rb_raise(rb_eArgError, "Invalid norm order for vectors.");
    }
    assert(0);
}

typedef struct {
    VALUE otype;
    void (*iter_func)(na_loop_t *);
    int extra;
} norm1d_info;

#if SIZEOF_INT == 2
# define UINTSUM_TYPE numo_cUInt16
#elif SIZEOF_INT == 4
# define UINTSUM_TYPE numo_cUInt32
#elif SIZEOF_INT == 8
# define UINTSUM_TYPE numo_cUInt64
#else
# error
#endif

static void
norm1d_setup_info(norm1d_info * const info, VALUE const ord, int const simple)
{
    <% if is_complex %>
    VALUE const otype_tbl0[] = { numo_cDFloat, numo_cDFloat };
    VALUE const otype_tbl1[] = { UINTSUM_TYPE, numo_cDFloat, numo_cDFloat };
    <% else %>
    VALUE const otype_tbl0[] = { cT, cT };
    VALUE const otype_tbl1[] = { UINTSUM_TYPE, cT, numo_cDFloat };
    <% end %>

#define F(SUFFIX) { &<%=c_iter%>1dg_##SUFFIX, &<%=c_iter%>1ds_##SUFFIX }

    void (* const func_tbl0[][2])(na_loop_t *) = { F(abs_max), F(abs_min) };
    void (* const func_tbl1[][2])(na_loop_t *) = { F(notzero_uintsum), F(abs_sum), F(2norm) };

#undef F

    int const idx = norm1d_ord_idx(ord);

    if (idx < 0) {
        int const i = -(1+idx);
        info->otype     = otype_tbl0[i];
        info->iter_func = func_tbl0[i][simple];
    } else if (idx < 3) {
        info->otype     = otype_tbl1[idx];
        info->iter_func = func_tbl1[idx][simple];
    } else {
        info->extra = idx;
    }
}

#undef UINTSUM_TYPE

static VALUE
norm1d_general(VALUE const a, narray_t * const na, VALUE const ord, unsigned ax)
{
    norm1d_info info = { .otype=Qundef, .iter_func=0, .extra=0 };

    norm1d_setup_info(&info, ord, 0);

    if (info.extra) {
        info.otype = numo_cDFloat;
        info.iter_func = &<%=c_iter%>1dg_extra;
    }
    {
        unsigned const d = na->ndim;
        unsigned const d_1 = d-1;
        size_t * const shape = ALLOCA_N(size_t, d_1);
        ndfunc_arg_in_t ain[1] = {{cT, d}};
        ndfunc_arg_out_t aout[1] = {
          {info.otype, d_1, shape}};
        ndfunc_t ndf = {info.iter_func, NO_LOOP,
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
                norm_opt_1d_t opt =
                  { .outer_n = outer_n
                  , .inner_n = inner_n
                  , .axis    = ax
                  , .extra   = info.extra } ;
                return na_ndloop3(&ndf, &opt, 1, a);
            }
        }
    }
}

static VALUE
norm1d_simple(VALUE const a, narray_t * const na, VALUE const ord)
{
    norm1d_info info = { .otype=Qundef, .iter_func=0, .extra=0 };

    norm1d_setup_info(&info, ord, 1);

    {
        norm_opt_1d_t opt;
        norm_opt_1d_t *opt_p = 0;
        if (info.extra) {
            opt_p = &opt;
            opt.extra = info.extra;

            info.otype = numo_cDFloat;
            info.iter_func = &<%=c_iter%>1ds_extra;
        }
        {
            ndfunc_arg_in_t ain[1] = {{cT, 1}};
            ndfunc_arg_out_t aout[1] = {
              {info.otype, 0, 0}};
            ndfunc_t ndf = {info.iter_func, NO_LOOP,
              COUNT_OF_(ain), COUNT_OF_(aout),
              ain, aout};
            VALUE const tmp = na_ndloop3(&ndf, opt_p, 1, a);
            return rb_funcall(tmp, rb_intern("[]"), 0);
        }
    }
}

static VALUE
norm2d_get_otype(int const idx)
{
    VALUE const tbl[8] = {
      <% if is_complex %>
        numo_cDFloat
      , numo_cDFloat
      , numo_cDFloat
      , numo_cDFloat
      , numo_cDFloat
      , numo_cDFloat
      , numo_cDFloat
      , numo_cDFloat
      <% else %>
        numo_cDFloat
      , numo_cDFloat
      , cT
      , cT
      , cT
      , cT
      , numo_cDFloat
      , numo_cDFloat
      <% end %>
    };
    return tbl[idx];
}

static void (*
norm2d_get_func(int const idx, int const simple, int const reverse)
)(na_loop_t *)
{

#define F(SUFFIX) \
  { {&<%=c_iter%>2dg_##SUFFIX, &<%=c_iter%>2dg_##SUFFIX##_x} \
  , {&<%=c_iter%>2ds_##SUFFIX, &<%=c_iter%>2ds_##SUFFIX##_x} \
  }

#define G(SUFFIX) \
  { {&<%=c_iter%>2dg_##SUFFIX, &<%=c_iter%>2dg_##SUFFIX##_x} \
  , {&<%=c_iter%>2ds_##SUFFIX, &<%=c_iter%>2ds_##SUFFIX} \
  } // special case for svd

    void (* const tbl[8][2][2])(na_loop_t *) = {
        F(fro)
      , G(nuc)
      , F(abs_sum_max)
      , F(abs_sum_min)
      , F(abs_sum_max)
      , F(abs_sum_min)
      , G(svd_max)
      , G(svd_min)
    };

#undef F
#undef G

    return tbl[idx][simple][reverse];
}

static int
norm2d_ord_idx(VALUE const ord)
{
    static VALUE idx_hash = Qnil;

    if (NIL_P(idx_hash)) {
        VALUE const h = rb_hash_new();

        rb_hash_aset(h,                     Qnil, INT2FIX(0));//nil=FrobeniusNorm
        rb_hash_aset(h, ID2SYM(rb_intern("FRO")), INT2FIX(0));
        rb_hash_aset(h, ID2SYM(rb_intern("NUC")), INT2FIX(1));
        rb_hash_aset(h,       DBL2NUM( INFINITY), INT2FIX(2));
        rb_hash_aset(h,       DBL2NUM(-INFINITY), INT2FIX(3));
        rb_hash_aset(h,              INT2FIX( 1), INT2FIX(4));
        rb_hash_aset(h,              INT2FIX(-1), INT2FIX(5));
        rb_hash_aset(h,              INT2FIX( 2), INT2FIX(6));
        rb_hash_aset(h,              INT2FIX(-2), INT2FIX(7));

        rb_hash_freeze(h);

        idx_hash = h;
    }
    {
        VALUE const v = rb_hash_aref(idx_hash, ord);
        if ( ! NIL_P(v)) {
            return FIX2LONG(v);
        }
    }
    rb_raise(rb_eArgError, "Invalid norm order for matrices.");
    RB_GC_GUARD(idx_hash);
}

typedef struct {
    VALUE otype;
    void (*iter_func)(na_loop_t *);
    int swapaxis, svd;
} norm2d_info;

static void
norm2d_setup_info(norm2d_info * const info, VALUE const ord, int const simple, unsigned axis0, unsigned axis1)
{
    int const idx = norm2d_ord_idx(ord);
    if ((idx == 4) || (idx == 5)) {
        unsigned tmp = axis0;
        axis0 = axis1;
        axis1 = tmp;
        info->swapaxis = 1;
    }
    if ((idx == 1) || (idx == 6) || (idx == 7)) {
        info->svd = 1;
    }
    info->otype = norm2d_get_otype(idx);
    info->iter_func = norm2d_get_func(idx, simple, axis0 > axis1);
}


static VALUE
norm2d_general(VALUE const a, narray_t * const na, VALUE const ord, unsigned axis0, unsigned axis1)
{
    unsigned const d = na->ndim;
    norm2d_info info = { .otype=Qundef, .iter_func=0, .swapaxis=0, .svd=0 };

    norm2d_setup_info(&info, ord, 0, axis0, axis1);
    if ((info.svd) && (axis0 >= d-2) && (axis1 >= d-2)) {
        return norm2dg_svd_optimized(a, na, info.iter_func);
    }
    if (info.swapaxis) {
        unsigned tmp = axis0;
        axis0 = axis1;
        axis1 = tmp;
    }
    {
        unsigned const d_2 = d-2;
        size_t * const shape = ALLOCA_N(size_t, d_2);
        ndfunc_arg_in_t ain[1] = {{cT, d}};
        ndfunc_arg_out_t aout[1] = {
          {info.otype, d_2, shape}};
        ndfunc_t ndf = {info.iter_func, NO_LOOP,
          COUNT_OF_(ain), COUNT_OF_(aout),
          ain, aout};
        {
            unsigned a0, a1;
            unsigned i, j;
            size_t outer_n = 1, mid_n = 1, inner_n = 1;
            if (axis0 < axis1) {
                a0 = axis0; a1 = axis1;
            } else {
                a0 = axis1; a1 = axis0;
            }
            for (i = 0, j = 0; i < a0; ++i, ++j) {
                size_t const n = na->shape[i];
                shape[j] = n;
                outer_n *= n;
            }
            ++i;
            for (; i < a1; ++i, ++j) {
                size_t const n = na->shape[i];
                shape[j] = n;
                mid_n *= n;
            }
            ++i;
            for (; i < d; ++i, ++j) {
                size_t const n = na->shape[i];
                shape[j] = n;
                inner_n *= n;
            }
            {
                norm_opt_2d_t opt =
                  { .outer_n = outer_n
                  , .mid_n   = mid_n
                  , .inner_n = inner_n
                  , .axis0   = axis0
                  , .axis1   = axis1 } ;
                return na_ndloop3(&ndf, &opt, 1, a);
            }
        }
    }
}

static VALUE
norm2d_simple(VALUE const a, narray_t * const na_, VALUE const ord, unsigned axis0, unsigned axis1)
{
    norm2d_info info = { .otype=Qundef, .iter_func=0, .swapaxis=0, .svd=0 };

    norm2d_setup_info(&info, ord, 1, axis0, axis1);
    if (info.svd) {
        na_loop_t lp;
        lp.option = a;
        (*info.iter_func)(&lp);
        return lp.option;
    }
    {
        ndfunc_arg_in_t ain[1] = {{cT, 2}};
        ndfunc_arg_out_t aout[1] = {
          {info.otype, 0, 0}};
        ndfunc_t ndf = {info.iter_func, NO_LOOP,
          COUNT_OF_(ain), COUNT_OF_(aout),
          ain, aout};
        VALUE const tmp = na_ndloop3(&ndf, 0, 1, a);
        return rb_funcall(tmp, rb_intern("[]"), 0);
    }
}

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
            unsigned uax;
            {
                int ax = NUM2INT(rb_ary_entry(axis, 0));
                if (ax < 0) {
                    ax += ndim;
                }
                if ((ax == 0) && (ndim == 1)) {
                    return norm1d_simple(a, na, ord);
                }
                uax = ax;
            }
            if (uax >= ndim) {
                rb_raise(rb_eArgError, "Invalid axis");
            }
            return norm1d_general(a, na, ord, uax);
        }
    case 2:
        if ( ! RB_INTEGER_TYPE_P(rb_ary_entry(axis, 1))) {
            rb_raise(rb_eTypeError, "axis: must be an Array of Integers or an Integer or nil");
        }
        {
            unsigned uax0, uax1;
            {
                int ax0 = NUM2INT(rb_ary_entry(axis, 0));
                int ax1 = NUM2INT(rb_ary_entry(axis, 1));
                if (ax0 < 0) {
                    ax0 += ndim;
                }
                if (ax1 < 0) {
                    ax1 += ndim;
                }
                uax0 = ax0; uax1 = ax1;
                if (ndim == 2) {
                    if (uax0 == uax1) {
                        rb_raise(rb_eArgError, "Duplicate axes given.");
                    } else {
                        return norm2d_simple(a, na, ord, uax0, uax1);
                    }
                }
            }
            if ((uax0 >= ndim) || (uax1 >= ndim)) {
                rb_raise(rb_eArgError, "Invalid axis");
            }
            if (uax0 == uax1) {
                rb_raise(rb_eArgError, "Duplicate axes given.");
            } else {
                return norm2d_general(a, na, ord, uax0, uax1);
            }
            assert(0);
        }
    default:
        assert(0);
    }
    assert(0);
}
