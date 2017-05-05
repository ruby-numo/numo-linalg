#define func_p <%=name%>_p

static <%=cblas_func%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t *const lp)
{
    char *p1, *p2;
    size_t n;
    ssize_t s1, s2;

    INIT_COUNTER(lp,n);
    INIT_PTR(lp,0,p1,s1);
    INIT_PTR(lp,1,p2,s2);

    (*func_p)(n, (dtype*)p1, s1/sizeof(dtype),
                 (dtype*)p2, s2/sizeof(dtype));
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
*
*       .. Scalar Arguments ..
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
*>    SCOPY copies a vector, x, to a vector, y.
*>    uses unrolled loops for increments equal to 1.
*> \endverbatim
*/
/*
 *  @overload <%=name%>( x, y )
 *  @param [Numo::NArray] x  1-dimentional NArray.
 *  @param [Numo::NArray] y  1-dimentional NArray.
 *  @return [nil]
 *  @raise
 *
 *    copies a vector, x, to a vector, y.
 *    uses unrolled loops for increments equal to 1.
 */
static VALUE
<%=c_func(2)%>(VALUE mod, VALUE x, VALUE y)
{
    narray_t *na1, *na2;
    ndfunc_arg_in_t ain[2] = {{cT,0},{OVERWRITE,0}};
    ndfunc_t ndf = {<%=c_iter%>, STRIDE_LOOP, 2,0, ain,0};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    CHECK_NARRAY_TYPE(x,cT);
    CHECK_NARRAY_TYPE(y,cT);
    GetNArray(x,na1);
    GetNArray(y,na2);
    CHECK_DIM_GE(na1,1);
    CHECK_DIM_GE(na2,1);
    CHECK_NON_EMPTY(na1);
    CHECK_NON_EMPTY(na2);
    CHECK_SAME_SHAPE(na1,na2);

    na_ndloop(&ndf, 2, x, y);

    return Qnil;
}

#undef func_p
