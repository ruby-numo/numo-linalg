#define func_p <%=func_name%>_p

static <%=func_name%>_t func_p = 0;

#undef result_dtype
#define result_dtype <%=result_dtype%>

static void
<%=c_iter%>(na_loop_t *const lp)
{
    char *p1, *p2;
    size_t n;
    ssize_t s1;

    INIT_PTR(lp,0,p1,s1);
    p2 = NDL_PTR(lp,1);
    n  = NDL_SHAPE(lp,0)[0];

    *(result_dtype*)p2 = (*func_p)(n, (dtype*)p1, s1/sizeof(dtype));
}

/*
 *  @overload <%=name%>( x )
 *  @param [Numo::NArray] x  >= 1-dimentional NArray.
 *  @return [Numo::NArray]
 *  @raise

<%=description%>

*/
static VALUE
<%=c_func(1)%>(VALUE mod, VALUE x)
{
    VALUE     ans;
    narray_t *na1;
    ndfunc_arg_in_t ain[1] = {{cT,1}};
    ndfunc_arg_out_t aout[1] = {{cT,0}};
    ndfunc_t ndf = {<%=c_iter%>, NDF_EXTRACT, 1,1, ain,aout};

    CHECK_FUNC(func_p,"<%=func_name%>");

    GetNArray(x,na1);
    CHECK_DIM_GE(na1,1);
    CHECK_NON_EMPTY(na1);

    ans = na_ndloop(&ndf, 1, x);

    return ans;
}

#undef func_p
