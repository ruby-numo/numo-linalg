#define args_t <%=name%>_args_t

typedef struct {
    enum CBLAS_ORDER order;
    enum CBLAS_UPLO uplo;
    rtype alpha;
    blasint n, lda;
} args_t;

#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t *const lp)
{
    dtype *a;
    char *p1;
    ssize_t s1;
    args_t *g;

    INIT_PTR(lp,0,p1,s1);
    a = (dtype*)NDL_PTR(lp,1);
    g = (args_t*)(lp->opt_ptr);

    (*func_p)(g->order, g->uplo, g->n, g->alpha,
              (dtype*)p1, s1/sizeof(dtype), a, g->lda);
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER INCX,LDA,N
*       CHARACTER UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYR   performs the symmetric rank 1 operation
*>
*>    A := alpha*x*x**T + A,
*>
*> where alpha is a real scalar, x is an n element vector and A is an
*> n by n symmetric matrix.
*> \endverbatim

*  Definition:
*  ===========
*
*       SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER INCX,LDA,N
*       CHARACTER UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZHER   performs the hermitian rank 1 operation
*>
*>    A := alpha*x*x**H + A,
*>
*> where alpha is a real scalar, x is an n element vector and A is an
*> n by n hermitian matrix.
*> \endverbatim
*/
/*
 * @overload <%=name%>( alpha, x, [,a ,uplo, order] )
 * @param [Numeric]      alpha
 * @param [Numo::DFloat] x  vector (>=1-dimentional NArray)
 * @param [Numo::DFloat] a  n-by-n symmetric matrix [in/out] (>=2-dimentional NArray)
 * @param [option] uplo  (default='upper')
 * @param [option] order (default='rowmajor')
 * @return [Numo::DFloat] return a
 * @raise
 *
 * DSYR   performs the symmetric rank 1 operation
 *
 *    A := alpha*x*x**T + A,
 *
 * where alpha is a real scalar, x is an n element vector and A is an
 * n by n symmetric matrix.

 * ZHER   performs the hermitian rank 1 operation
 *
 *    A := alpha*x*x**H + A,
 *
 * where alpha is a real scalar, x is an n element vector and A is an
 * n by n hermitian matrix.
 */
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE mod)
{
    int       i;
    VALUE     ans, arg_a = Qnil;
    VALUE     x, a, alpha, uplo, order;
    narray_t *na1, *na3;
    blasint   nx, na;
    size_t    shape[2];
    ndfunc_arg_in_t ain[3] = {{cT,1},{OVERWRITE,2},{sym_init,0}};
    ndfunc_arg_out_t aout[1] = {{cT,2,shape}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, 2, 0, ain, aout};
    args_t g = {CblasRowMajor, CblasUpper, 1.0};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    i = rb_scan_args(argc, argv, "23",
                     &alpha, &x, &a, &uplo, &order);
    switch (i) {
    case 5: g.order = option_order(order);
    case 4: g.uplo  = option_uplo(uplo);
    case 3: arg_a = a;
    }
    if (RTEST(alpha)) {g.alpha = NUM2DBL(alpha);}

    GetNArray(x,na1);
    CHECK_DIM_GE(na1,1);
    nx = na1->shape[na1->ndim-1]; // n
    g.n = nx;

    if (arg_a == Qnil) { // c is not given.
        ndf.nout = 1;
        ain[1] = ain[2];
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

    ans = na_ndloop3(&ndf, &g, 2, x, a);

    if (arg_a == Qnil) { // a is not given.
        return ans;
    } else {
        return a;
    }
}

#undef func_p
#undef args_t
