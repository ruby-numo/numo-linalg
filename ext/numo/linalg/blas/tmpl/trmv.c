#define args_t <%=name%>_args_t

typedef struct {
  enum CBLAS_ORDER order;
  enum CBLAS_UPLO uplo;
  enum CBLAS_TRANSPOSE trans;
  enum CBLAS_DIAG diag;
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

    a = (dtype*)NDL_PTR(lp,0);
    INIT_PTR(lp,1,p1,s1);
    g = (args_t*)(lp->opt_ptr);

    (*func_p)(g->order, g->uplo, g->trans, g->diag, g->n,
              a, g->lda, (dtype*)p1, s1/sizeof(dtype));
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,LDA,N
*       CHARACTER DIAG,TRANS,UPLO
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
*> DTRMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular matrix.
*> \endverbatim
*
*/
/*
 * @overload <%=name%>( a, x [,uplo, trans, diag, order] )
 * @param [Numo::DFloat] a  n-by-n matric (>=2-dimentional NArray)
 * @param [Numo::DFloat] x  n-size vector (>=1-dimentional NArray, allow inplace)
 * @return [Numo::DFloat] return x
 * @raise
 *
 * performs one of the matrix-vector operations
 *
 *    x := A*x,   or   x := A**T*x,
 *
 * where x is an n element vector and  A is an n by n unit, or non-unit,
 * upper or lower triangular matrix.
 */
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE mod)
{
    int       i;
    VALUE     a, x, uplo, trans, diag, order;
    narray_t *na1, *na2;
    blasint   na, nx;
    ndfunc_arg_in_t ain[2] = {{cT,2},{OVERWRITE,1}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, 2, 0, ain, 0};
    args_t g = {CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    i = rb_scan_args(argc, argv, "24",
                     &a, &x, &uplo, &trans, &diag, &order);
    switch (i) {
    case 6: g.order = option_order(order);
    case 5: g.diag = option_diag(diag);
    case 4: g.trans = option_trans(trans);
    case 3: g.uplo = option_uplo(uplo);
    }

    GetNArray(a,na1);
    GetNArray(x,na2);
    CHECK_DIM_GE(na1,2);
    CHECK_DIM_GE(na2,1);

    CHECK_SQUARE("a",na1);

    na = na1->shape[na1->ndim-1]; // n (lda)
    g.lda = na;

    nx = na2->shape[na2->ndim-1];
    CHECK_INT_EQ("na",na,"nx",nx);
    g.n = nx;

    if (!TEST_INPLACE(x)) {
        x = na_copy(x);
    }

    na_ndloop3(&ndf, &g, 2, a, x);
    return x;
}

#undef func_p
#undef args_t
