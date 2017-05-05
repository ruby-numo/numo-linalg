#define args_t <%=name%>_args_t

typedef struct {
    enum CBLAS_ORDER order;
    enum CBLAS_UPLO uplo;
    dtype alpha;
    blasint n, lda;
} args_t;

#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t *const lp)
{
    dtype *a;
    char *p1, *p2;
    ssize_t s1, s2;
    args_t *g;

    INIT_PTR(lp,0,p1,s1);
    INIT_PTR(lp,1,p2,s2);
    a = (dtype*)NDL_PTR(lp,2);
    g = (args_t*)(lp->opt_ptr);

  <% if /[cz]/ =~ blas_char %>
    (*func_p)(g->order, g->uplo, g->n,
              &g->alpha, (dtype*)p1, s1/sizeof(dtype),
              (dtype*)p2, s2/sizeof(dtype), a, g->lda);
  <% else %>
    (*func_p)(g->order, g->uplo, g->n,
              g->alpha, (dtype*)p1, s1/sizeof(dtype),
              (dtype*)p2, s2/sizeof(dtype), a, g->lda);
  <% end %>
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER INCX,INCY,LDA,N
*       CHARACTER UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYR2  performs the symmetric rank 2 operation
*>
*>    A := alpha*x*y**T + alpha*y*x**T + A,
*>
*> where alpha is a scalar, x and y are n element vectors and A is an n
*> by n symmetric matrix.
*> \endverbatim
*
*/
/*
 * @overload <%=name%>( alpha, x, y [,a ,uplo, order] )
 * @param [Numeric]      alpha
 * @param [Numo::DFloat] x  vector (>=1-dimentional NArray)
 * @param [Numo::DFloat] y  vector (>=1-dimentional NArray)
 * @param [Numo::DFloat] a  n-by-n symmetric matrix [in/out] (>=2-dimentional NArray)
 * @param [option] uplo  (default='upper')
 * @param [option] order (default='rowmajor')
 * @return [Numo::DFloat]
 * @raise
 *
 * performs one of the symmetric rank 2k operations
 *
 *      C := alpha*A*B**T + alpha*B*A**T + beta*C,
 *    or
 *      C := alpha*A**T*B + alpha*B**T*A + beta*C,
 *
 * where alpha and beta are scalars, C is an n by n symmetric matrix
 * and A and B are n by k matrices in the first case and k by n
 * matrices in the second case.
 */
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE mod)
{
    int       i;
    VALUE     ans, arg_a = Qnil;
    VALUE     x, y, a, alpha, uplo, order;
    narray_t *na1, *na2, *na3;
    blasint   nx, ny, na;
    size_t    shape[2];
    ndfunc_arg_in_t ain[4] = {{cT,1},{cT,1},{OVERWRITE,2},{sym_init,0}};
    ndfunc_arg_out_t aout[1] = {{cT,2,shape}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, 3, 0, ain, aout};
    args_t g = {CblasRowMajor, CblasUpper, m_one};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    i = rb_scan_args(argc, argv, "33",
                     &alpha, &x, &y, &a, &uplo, &order);
    switch (i) {
    case 6: g.order = option_order(order);
    case 5: g.uplo  = option_uplo(uplo);
    case 4: arg_a = a;
    }
    if (RTEST(alpha)) {g.alpha = m_num_to_data(alpha);}

    GetNArray(x,na1);
    GetNArray(y,na2);
    CHECK_DIM_GE(na1,1);
    CHECK_DIM_GE(na2,1);
    nx = na1->shape[na1->ndim-1]; // n
    ny = na2->shape[na2->ndim-1]; // n
    CHECK_INT_EQ("nx",nx,"ny",ny);
    g.n = nx;

    if (arg_a == Qnil) { // c is not given.
        ndf.nout = 1;
        ain[2] = ain[3];
        a = INT2FIX(0);
        shape[0] = shape[1] = g.lda = nx;
    } else {
        a = rb_funcall(cT,rb_intern("cast"),1,arg_a);
        if (!TEST_INPLACE(a)) {
            a = na_copy(a);
        }
        GetNArray(a,na3);
        CHECK_DIM_GE(na3,2);
        CHECK_SQUARE("a",na3);
        g.lda = na = na3->shape[na3->ndim-1]; // n (lda)
        CHECK_SIZE_EQ(na,nx);
    }

    ans = na_ndloop3(&ndf, &g, 3, x, y, a);

    if (arg_a == Qnil) { // a is not given.
        return ans;
    } else {
        return a;
    }
}

#undef func_p
#undef args_t
