/*<%
 has_uplo = (/^.(sy|he|po)/ =~ name)
 has_piv  = ((/tr(f|i)/ =~ name && /potr/ !~ name) || /geqp3/ =~ name)
 has_tau  = (/q/ =~ name || /tzrzf/ =~ name)
 symmetric = (has_uplo || /getri/=~name)
%>*/
<% %>
#define UPLO <%= has_uplo ? "1":"0" %>
#define PIV  <%= has_piv ? "1":"0" %>
#define TAU  <%= has_tau ? "1":"0" %>
#define QP3  <%= /geqp3/ =~ name ? "1":"0" %>
#define SYMMETRIC <%= symmetric ? "1":"0" %>
#define args_t <%=func_name%>_args_t
#define func_p <%=func_name%>_p

typedef struct {
    int order;
    char uplo;
} args_t;

static <%=func_name%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a;
#if PIV
    int   *ipiv;
#endif
#if TAU
    dtype *tau;
#endif
    int   *info;
    int    m, n, lda, tmp;
    args_t *g;

    a = (dtype*)NDL_PTR(lp,0);
#if PIV
    ipiv = (int*)NDL_PTR(lp,1);
#endif
#if TAU
    tau = (dtype*)NDL_PTR(lp,1+PIV);
#endif
    info = (int*)NDL_PTR(lp,1+PIV+TAU);
    g = (args_t*)(lp->opt_ptr);

    m = NDL_SHAPE(lp,0)[0];
    n = NDL_SHAPE(lp,0)[1];
    SWAP_IFCOL(g->order,m,n);
#if SYMMETRIC
    n = min_(m,n);
#endif
    lda = NDL_STEP(lp,0) / sizeof(dtype);

    //printf("order=%d m=%d n=%d lda=%d \n",g->order,m,n,lda);

    /*<%
    func_args = [   "g->order",
      has_uplo   && "g->uplo",
      !symmetric && "m",
                    "n, a, lda",
      has_piv    && "ipiv",
      has_tau    && "tau",
    ].select{|x| x}.join(", ")
    %>*/
    *info = (*func_p)(<%=func_args%>);
    CHECK_ERROR(*info);
}

/*<%
 args_opt = [
   has_uplo && "uplo:'U'",
   "order:'R'",
 ].select{|x| x}.join(", ")

 params = [
   param("a",2),
   has_uplo && param("uplo"),
   param("order"),
 ].select{|x| x}.join("\n  ")

 return_type = [
   class_name,
   has_piv && "Numo::Int",
   has_tau && class_name,
   "Integer"
 ].select{|x| x}.join(", ")

 return_name = [
   "a",
   has_piv && "piv",
   has_tau && "tau",
   "info"
 ].select{|x| x}.join(", ")
%>
  @overload <%=name%>(a, [<%=args_opt%>])
  <%=params%>
  @return [[<%=return_type%>]] array of [<%=return_name%>]

 <%=description%>

*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    VALUE a, ans;
    int   m, n, tmp;
    narray_t *na1;
    /*<%
    aout = [
      has_piv && "{cInt,1,shape_piv}",
      has_tau && "{cT,1,shape_tau}",
                 "{cInt,0}",
    ].select{|x| x}.join(",")
    %>*/
#if PIV
    size_t shape_piv[1];
#endif
#if TAU
    size_t shape_tau[1];
#endif
    ndfunc_arg_in_t ain[1] = {{OVERWRITE,2}};
    ndfunc_arg_out_t aout[PIV+TAU+1] = {<%=aout%>};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 1,PIV+TAU+1,ain,aout};

    args_t g = {0,1};
    VALUE opts[2] = {Qundef,Qundef};
    VALUE kw_hash = Qnil;
    ID kw_table[2] = {id_order,id_uplo};

    CHECK_FUNC(func_p,"<%=func_name%>");

    rb_scan_args(argc, argv, "1:", &a, &kw_hash);
    rb_get_kwargs(kw_hash, kw_table, 0, 1+UPLO, opts);
    g.order = option_order(opts[0]);
#if UPLO
    g.uplo = option_uplo(opts[1]);
#endif

    COPY_OR_CAST_TO(a,cT);
    GetNArray(a, na1);
    CHECK_DIM_GE(na1, 2);
    m = ROW_SIZE(na1);
    n = COL_SIZE(na1);
    SWAP_IFCOL(g.order,m,n);
#if PIV
#if QP3
    shape_piv[0] = n;
#else
    shape_piv[0] = min_(m,n);
#endif
#endif
#if TAU
    shape_tau[0] = min_(m,n);
#endif

    ans = na_ndloop3(&ndf, &g, 1, a);

    rb_ary_unshift(ans, a);
    return ans;
}

#undef args_t
#undef func_p
#undef UPLO
#undef PIV
#undef TAU
#undef QP3
#undef SYMMETRIC
