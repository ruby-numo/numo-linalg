#define args_t <%=func_name%>_args_t
#define func_p <%=func_name%>_p

typedef struct {
    int   order;
    char  norm;
} args_t;

static <%=func_name%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a;
    rtype *norm;
    int    m, n, lda, tmp;
    args_t *g;

    a = (dtype*)NDL_PTR(lp,0);
    norm = (rtype*)NDL_PTR(lp,1);
    g = (args_t*)(lp->opt_ptr);

    m = NDL_SHAPE(lp,0)[0];
    n = NDL_SHAPE(lp,0)[1];
    SWAP_IFCOL(g->order,m,n);
    lda = NDL_STEP(lp,0) / sizeof(dtype);

    printf("order=%d m=%d n=%d lda=%d\n",g->order,m,n,lda);

    *norm = (*func_p)(g->order, g->norm, m, n, a, lda);
}

/*
<% tp = "Numo::"+class_name %>
  @overload <%=name%>(a, norm, [axis:nil, keepdims:false, order:'r'])
  @param [<%=tp%>]  a  >=2-dimentional NArray.
  @return [String or Symbol] norm option
  @return [<%=tp%>] <%=name%> =
                        ( max(abs(A(i,j))), NORM = 'M' or 'm'
                        ( norm1(A),         NORM = '1', 'O' or 'o'
                        ( normI(A),         NORM = 'I' or 'i'
                        ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

  <%=description%>

*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    VALUE a, norm, reduce, ans;
    int   m, n, tmp;
    narray_t *na1;
    ndfunc_arg_in_t ain[2] = {{cT,2},{sym_reduce,0}};
    ndfunc_arg_out_t aout[1] = {{cRT,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 2, 1, ain, aout};

    args_t g;
    VALUE opts[3] = {Qundef,Qundef,Qundef};
    ID kw_table[3] = {id_order,id_axis,id_keepdims};
    VALUE kw_hash = Qnil;

    CHECK_FUNC(func_p,"<%=func_name%>");

    rb_scan_args(argc, argv, "2:", &a, &norm, &kw_hash);
    rb_get_kwargs(kw_hash, kw_table, 0, 3, opts);
    g.order = option_order(opts[0]);
    g.norm  = option_job(norm,'F','F');
    reduce = nary_reduce_options(Qnil, &opts[1], 1, &a, &ndf);
    //A is DOUBLE PRECISION array, dimension (LDA,N)
    //On entry, the M-by-N matrix A.
    //COPY_OR_CAST_TO(a,cT); // not overwrite
    GetNArray(a, na1);
    CHECK_DIM_GE(na1, 2);

    ans = na_ndloop3(&ndf, &g, 2, a, reduce);
    return ans;
}

#undef func_p
#undef args_t
