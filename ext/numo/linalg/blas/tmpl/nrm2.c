#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

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
    n  = lp->args[0].shape[0];

    *(result_dtype*)p2 = (*func_p)(n, (dtype*)p1, s1/sizeof(dtype));
}

/*
*  Definition:
*  ===========
*
*       REAL FUNCTION SNRM2(N,X,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       REAL X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    SNRM2 := sqrt( x'*x ).
*> \endverbatim
*/
/*
*  Definition:
*  ===========
*
*       REAL FUNCTION SCNRM2(N,X,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SCNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    SCNRM2 := sqrt( x**H*x )
*> \endverbatim
*/
/*
 *  @overload <%=name%>( x )
 *  @param [Numo::NArray] x  >= 1-dimentional NArray.
 *  @return [Numo::NArray]
 *  @raise
 *
 *  <%=name%> returns the euclidean norm of a vector via the function
 *  name, so that <%=name%> := sqrt( x'*x ).
 */
static VALUE
<%=c_func(1)%>(VALUE mod, VALUE x)
{
    VALUE     ans;
    narray_t *na1;
    size_t    shape[1]={1};
    ndfunc_arg_in_t ain[1] = {{cT,1}};
    ndfunc_arg_out_t aout[1] = {{<%=result_class%>,0,shape}};
    ndfunc_t ndf = {<%=c_iter%>, NDF_EXTRACT, 1,1, ain,aout};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    GetNArray(x,na1);
    CHECK_DIM_GE(na1,1);
    CHECK_NON_EMPTY(na1);

    ans = na_ndloop(&ndf, 1, x);

    return ans;
}

#undef func_p
