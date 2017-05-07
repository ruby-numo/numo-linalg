#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

<% if /^(cs|zd)scal/ =~ name %>
#define scal_t rtype
<% else %>
#define scal_t dtype
<% end %>

static void
<%=c_iter%>(na_loop_t *const lp)
{
    char *p1;
    size_t n;
    ssize_t s1;
    scal_t *g;

    INIT_COUNTER(lp,n);
    INIT_PTR(lp,0,p1,s1);
    g = (scal_t*)(lp->opt_ptr);

  <% if /^[cz]scal/ =~ name %>
    (*func_p)(n, g, (dtype*)p1, s1/sizeof(dtype));
  <% else %>
    (*func_p)(n, *g, (dtype*)p1, s1/sizeof(dtype));
  <% end %>
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE SSCAL(N,SA,SX,INCX)
*
*       .. Scalar Arguments ..
*       REAL SA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       REAL SX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    scales a vector by a constant.
*>    uses unrolled loops for increment equal to 1.
*> \endverbatim
*/
/*
 *  @overload <%=name%>( a, x )
 *  @param [Float]        a  scale factor
 *  @param [Numo::NArray] x  1-dimentional NArray. [in/out]
 *  @return x
 *  @raise
 *
 *    applies a plane rotation.
 */
static VALUE
<%=c_func(2)%>(VALUE mod, VALUE a, VALUE x)
{
    scal_t g[1];
    narray_t *na1;
    ndfunc_arg_in_t ain[1] = {{OVERWRITE,0}};
    ndfunc_t ndf = {<%=c_iter%>, STRIDE_LOOP, 1,0, ain,0};

    check_func((void*)(&func_p),"<%=cblas_func%>");

  <% if /^(cs|zd)scal/ =~ name %>
    if (RTEST(a)) {g[0] = NUM2DBL(a);} else {g[0]=1;}
  <% else %>
    if (RTEST(a)) {g[0] = m_num_to_data(a);} else {g[0]=m_one;}
  <% end %>
    CHECK_NARRAY_TYPE(x,cT);
    GetNArray(x,na1);
    CHECK_DIM_GE(na1,1);
    CHECK_NON_EMPTY(na1);

    na_ndloop3(&ndf, g, 1, x);

    return x;
}

#undef func_p
#undef scal_t
