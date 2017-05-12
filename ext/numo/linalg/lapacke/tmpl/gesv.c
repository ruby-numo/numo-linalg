/*
<%
 uplo = (/^?ge/=~name) ? nil : "g->uplo,"
 ipiv = (/^?po/=~name) ? nil : "ipiv,"

 tp = "Numo::"+class_name
 iary = "Numo::Int"
 iscal = "Integer"
 if uplo
   a = "a, b [, order:'r', uplo:'u']"
 else
   a = "a, b [, order:'r']"
 end
 if ipiv
   n = "lu, x, piv, info"
   t = [tp,tp,iary,iscal]
 else
   n = "lu, x, info"
   t = [tp,tp,iscal]
 end
 return_type = t.join(", ")
 return_name = n
 params = a
%>
*/
#define UPLO <%=(/^?ge/!~name) ? "1":"0"%>
#define IPIV <%=(/^?po/!~name) ? "1":"0"%>
#define args_t <%=func_name%>_args_t
#define func_p <%=func_name%>_p

typedef struct {
    int order;
    char uplo;
} args_t;

static <%=func_name%>_t <%=func_name%>_p = 0;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a, *b;
    int   *info;
    int    n, nrhs;
    int    lda, ldb;
    args_t *g;
#if IPIV
    int *ipiv;
    ipiv = (int*)NDL_PTR(lp,2);
    info = (int*)NDL_PTR(lp,3);
#else
    info = (int*)NDL_PTR(lp,2);
#endif
    a = (dtype*)NDL_PTR(lp,0);
    b = (dtype*)NDL_PTR(lp,1);
    g = (args_t*)(lp->opt_ptr);

    n = lp->args[0].shape[0];
    lda = lp->args[0].iter[0].step / sizeof(dtype);
    if (lp->args[1].ndim == 1) {
        nrhs = 1;
        ldb = (g->order==LAPACK_COL_MAJOR) ? n : 1;
    } else {
        nrhs = lp->args[1].shape[1];
        ldb = lp->args[1].iter[0].step / sizeof(dtype);
    }
    //printf("order=%d n=%d nrhs=%d lda=%d ldb=%d b.ndim=%d\n",
    //       g->order,n,nrhs,lda,ldb,lp->args[1].ndim);
    *info = (*func_p)(g->order, <%=uplo%>
                      n, nrhs, a, lda, <%=ipiv%> b, ldb);
    CHECK_ERROR(*info);
}

/*
  @overload <%=name%>(<%=params%>)
  @param [<%=tp%>] a  >=2-dimentional NArray.
  @param [<%=tp%>] b  >=1-dimentional NArray.
  @return [[<%=return_type%>]] array of [<%=return_name%>]

  <%=description%>
*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE const mod)
{
    VALUE a, b, ans;
    narray_t *na1, *na2;
    size_t n11, n12, n21, n22;
    ndfunc_arg_in_t ain[2] = {{OVERWRITE,2},{OVERWRITE,2}};
    size_t shape[2];
    ndfunc_arg_out_t aout[2] = {{cInt,1,shape},{cInt,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 2, 2, ain, aout};
    args_t g;
    VALUE kw_hash = Qnil;
    ID kw_table[2] = {id_order,id_uplo};
    VALUE opts[2] = {Qundef,Qundef};

    CHECK_FUNC(func_p,"<%=func_name%>");

    rb_scan_args(argc, argv, "2:", &a, &b, &kw_hash);
    rb_get_kwargs(kw_hash, kw_table, 0, 2, opts);
    g.order = option_order(opts[0]);
    g.uplo = option_uplo(opts[1]);

    a = rb_funcall(cT,rb_intern("cast"),1,a);
    if (!TEST_INPLACE(a)) {a = na_copy(a);}
    b = rb_funcall(cT,rb_intern("cast"),1,b);
    if (!TEST_INPLACE(b)) {b = na_copy(b);}

    GetNArray(a, na1);
    GetNArray(b, na2);
    CHECK_DIM_GE(na1, 2);
    CHECK_DIM_GE(na2, 1);
    n11 = na1->shape[na1->ndim-2]; // n
    n12 = na1->shape[na1->ndim-1]; // n
    if (NA_NDIM(na2) == 1) {
        ain[1].dim = 1;
        n21 = na2->shape[na2->ndim-1]; // n
        n22 = 1;                       // nrhs
    } else {
        n21 = na2->shape[na2->ndim-2]; // n
        n22 = na2->shape[na2->ndim-1]; // nrhs
    }
    if ((n11 != n12) || (n11 != n21)) {
        rb_raise(nary_eShapeError, "matrix dimension mismatch: "
                 "n11=%"SZF"u n12=%"SZF"u n21=%"SZF"u", n11, n12, n21);
    }
    shape[0] = n11; // n
    shape[1] = n22; // nrhs
#if !IPIV
    aout[0] = aout[1];
    ndf.nout--;
#endif

    ans = na_ndloop3(&ndf, &g, 2, a, b);
#if IPIV
    return rb_ary_concat(rb_assoc_new(a,b),ans);
#else
    return rb_ary_new3(3,a,b,ans);
#endif
}

#undef func_p
#undef args_t
#undef UPLO
#undef IPIV
