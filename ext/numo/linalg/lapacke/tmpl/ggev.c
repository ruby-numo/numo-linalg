/*<%
 aout = [         "{cT,1,shape},{cT,1,shape}",
   !is_complex && "{cT,1,shape}",
                  "{cT,2,shape},{cT,2,shape},{cInt,0}"
 ].select{|x| x}.join(",")

 func_args = [
   "g->order, g->jobvl, g->jobvr, n, a, lda, b, ldb",
   is_complex ? "alpha" : "alphar, alphai",
   "beta, vl, ldvl, vr, ldvr"
 ].join(",")

 tp = "Numo::"+class_name
 return_type = ([tp]*(is_complex ? 4 : 5)+["Integer"]).join(", ")
 return_name = (is_complex ? "alpha,":"alphar, alphai,")+" beta, vl, vr, info"
%>*/
#define args_t <%=func_name%>_args_t
#define func_p <%=func_name%>_p

typedef struct {
    int order;
    char jobvl, jobvr;
} args_t;

static <%=func_name%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a, *b, *beta, *vl, *vr;
#if IS_COMPLEX
    dtype *alpha;
#else
    dtype *alphar, *alphai;
#endif
    int   *info;
    int    n, lda, ldb, ldvl, ldvr;
    args_t *g;

    a = (dtype*)NDL_PTR(lp,0);
    b = (dtype*)NDL_PTR(lp,1);
#if IS_COMPLEX
#define N 0
    alpha = (dtype*)NDL_PTR(lp,2);
#else
#define N 1
    alphar = (dtype*)NDL_PTR(lp,2);
    alphai = (dtype*)NDL_PTR(lp,3);
#endif
    beta = (dtype*)NDL_PTR(lp,N+3);
    vl = (dtype*)NDL_PTR(lp,N+4);
    vr = (dtype*)NDL_PTR(lp,N+5);
    info = (int*)NDL_PTR(lp,N+6);
    g = (args_t*)(lp->opt_ptr);

    n = lp->args[0].shape[1];
    lda = lp->args[0].iter[0].step / sizeof(dtype);
    ldb = lp->args[1].iter[0].step / sizeof(dtype);
    ldvl = lp->args[N+4].iter[0].step / sizeof(dtype);
    if (ldvl == 0) { ldvl = n; } // jobvl == 'N'
    ldvr = lp->args[N+5].iter[0].step / sizeof(dtype);
    if (ldvr == 0) { ldvr = n; } // jobvr == 'N'

    //printf("order=%d jobvl=%c jobvr=%c n=%d lda=%d ldb=%d ldvl=%d ldvr=%d\n",g->order,g->jobvl, g->jobvr, n, lda,ldb,ldvl,ldvr);

    *info = (*func_p)( <%=func_args%> );
    CHECK_ERROR(*info);
}

/*
  @overload <%=name%>(a, b [,jobvl:'v', jobvr:'v', order:'r'] )
  @param [Numo::<%=class_name%>] a >=2-dimentional NArray.
  @param [Numo::<%=class_name%>] b >=2-dimentional NArray.
  @param [String,Symbol] jobvl
    jobvl='N':  do not compute the left generalized eigenvectors;
    jobvl='V':  compute the left generalized eigenvectors.
  @param [String,Symbol] jobvr
    jobvr='N':  do not compute the left generalized eigenvectors;
    jobvr='V':  compute the left generalized eigenvectors.
  @return [[<%=return_type%>]] array of [<%=return_name%>]

 <%= description %>

*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    VALUE a, b, ans;
    int   n, nb;
    narray_t *na1, *na2;
    size_t shape[2];
    ndfunc_arg_in_t ain[2] = {{OVERWRITE,2},{OVERWRITE,2}};
#if IS_COMPLEX
    ndfunc_arg_out_t aout[5] = {{cT,1,shape},{cT,1,shape},
                                {cT,2,shape},{cT,2,shape},{cInt,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 2, 5, ain, aout};
#else
    ndfunc_arg_out_t aout[6] = {{cT,1,shape},{cT,1,shape},{cT,1,shape},
                                {cT,2,shape},{cT,2,shape},{cInt,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 2, 6, ain, aout};
#endif

    args_t g;
    VALUE opts[3] = {Qundef,Qundef,Qundef};
    VALUE kw_hash = Qnil;
    ID kw_table[3] = {id_order,id_jobvl,id_jobvr};

    CHECK_FUNC(func_p,"<%=func_name%>");

    rb_scan_args(argc, argv, "2:", &a, &b, &kw_hash);
    rb_get_kwargs(kw_hash, kw_table, 0, 3, opts);
    g.order = option_order(opts[0]);
    g.jobvl = option_job(opts[1],'V','N');
    g.jobvr = option_job(opts[2],'V','N');

    COPY_OR_CAST_TO(a,cT);
    GetNArray(a, na1);
    CHECK_DIM_GE(na1, 2);

    COPY_OR_CAST_TO(b,cT);
    GetNArray(b, na2);
    CHECK_DIM_GE(na2, 2);

    CHECK_SQUARE("matrix a",na1);
    n  = COL_SIZE(na1);
    CHECK_SQUARE("matrix b",na2);
    nb = COL_SIZE(na2);
    if (n != nb) {
        rb_raise(nary_eShapeError,"matrix a and b must have same size");
    }
    shape[0] = shape[1] = n;
    if (g.jobvl=='N') { aout[N+2].dim = 0; }
    if (g.jobvr=='N') { aout[N+3].dim = 0; }

    ans = na_ndloop3(&ndf, &g, 2, a, b);

    if (aout[N+3].dim == 0) { RARRAY_ASET(ans,N+3,Qnil); }
    if (aout[N+2].dim == 0) { RARRAY_ASET(ans,N+2,Qnil); }
    return ans;
}

#undef N
#undef args_t
#undef func_p
