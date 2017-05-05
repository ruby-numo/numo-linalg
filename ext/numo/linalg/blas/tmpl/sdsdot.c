#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t *const lp)
{
    char *p1, *p2, *p3;
    size_t n;
    ssize_t s1, s2;
    dtype  *g;

    INIT_PTR(lp,0,p1,s1);
    INIT_PTR(lp,1,p2,s2);
    p3 = NDL_PTR(lp,2);
    n  = lp->args[0].shape[0];
    g  = (dtype*)(lp->opt_ptr);

    *(dtype*)p3 = (*func_p)(n, *g, (dtype*)p1, s1/sizeof(dtype),
                                   (dtype*)p2, s2/sizeof(dtype));
}

/*
*  Definition:
*  ===========
*
*       REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
*
*       .. Scalar Arguments ..
*       REAL SB
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       REAL SX(*),SY(*)
*       ..
*
*    PURPOSE
*    =======
*
*    Compute the inner product of two vectors with extended
*    precision accumulation.
*
*    Returns S.P. result with dot product accumulated in D.P.
*    SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
*    where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
*    defined in a similar way using INCY.
*/
/*
 *  @overload <%=name%>( sb, sx, xy )
 *  @param [Float] sb
 *  @param [Numo::NArray] sx  >= 1-dimentional NArray.
 *  @param [Numo::NArray] sy  >= 1-dimentional NArray.
 *  @return [Numo::NArray]
 *  @raise
 *
 *  <%=name%> forms the dot product of two vectors.
 *  uses unrolled loops for increments equal to one.
 */
static VALUE
<%=c_func(2)%>(VALUE mod, VALUE sb, VALUE x, VALUE y)
{
    VALUE     ans;
    dtype     g[1] = {m_zero};
    narray_t *na1, *na2;
    size_t    nx, ny, shape[1]={1};
    ndfunc_arg_in_t ain[2] = {{cT,1},{cT,1}};
    ndfunc_arg_out_t aout[1] = {{<%=result_class%>,0,shape}};
    ndfunc_t ndf = {<%=c_iter%>, NDF_EXTRACT, 2,1, ain,aout};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    if (RTEST(sb)) {g[0] = m_num_to_data(sb);}

    GetNArray(x,na1);
    GetNArray(y,na2);
    CHECK_DIM_GE(na1,1);
    CHECK_DIM_GE(na2,1);
    CHECK_NON_EMPTY(na1);
    CHECK_NON_EMPTY(na2);
    nx = na1->shape[na1->ndim-1];
    ny = na2->shape[na2->ndim-1];
    CHECK_SIZE_EQ(nx,ny);

    ans = na_ndloop3(&ndf, g, 2, x, y);

    return ans;
}
#undef func_p
