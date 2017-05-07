#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t *const lp)
{
    char   *p1, *p2;
    size_t  n;
    ssize_t s1, s2;
    dtype  *g;

    INIT_COUNTER(lp,n);
    INIT_PTR(lp,0,p1,s1);
    INIT_PTR(lp,1,p2,s2);
    g = (dtype*)(lp->opt_ptr);

    (*func_p)(n, DP(*g), (dtype*)p1, s1/sizeof(dtype),
                         (dtype*)p2, s2/sizeof(dtype));
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
*
*       .. Scalar Arguments ..
*       REAL SA
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       REAL SX(*),SY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    SAXPY constant times a vector plus a vector.
*>    uses unrolled loops for increments equal to one.
*> \endverbatim
*/
/*
 * @overload <%=name%>( alpha, x, y )
 * @param [Numo::DFloat] x  >=1-dimentional NArray.
 * @param [Numo::DFloat] y  >=1-dimentional NArray.
 * @return [Numo::DFloat]  return y
 * @raise
 *
 *    constant times a vector plus a vector.
 *    uses unrolled loops for increments equal to one.
*/
static VALUE
<%=c_func(3)%>(VALUE mod, VALUE alpha, VALUE x, VALUE y)
{
    dtype g[1] = {m_one};
    narray_t *na1, *na2;
    ndfunc_arg_in_t ain[2] = {{cT,0},{OVERWRITE,0}};
    ndfunc_t ndf = {<%=c_iter%>, STRIDE_LOOP, 2, 0, ain, 0};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    if (RTEST(alpha)) {g[0] = m_num_to_data(alpha);}

    CHECK_NARRAY_TYPE(y,cT);
    GetNArray(x,na1);
    GetNArray(y,na2);
    CHECK_DIM_GE(na1,1);
    CHECK_DIM_GE(na2,1);
    CHECK_NON_EMPTY(na1);
    CHECK_NON_EMPTY(na2);
    CHECK_SAME_SHAPE(na1,na2);

    na_ndloop3(&ndf, g, 2, x, y);
    return y;
}

#undef func_p
