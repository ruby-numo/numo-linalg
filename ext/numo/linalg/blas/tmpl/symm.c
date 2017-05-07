#define args_t <%=name%>_args_t

typedef struct {
  enum CBLAS_ORDER order;
  enum CBLAS_SIDE side;
  enum CBLAS_UPLO uplo;
  dtype alpha, beta;
  blasint m, n, lda, ldb, ldc;
} args_t;

#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t *const lp)
{
    dtype *a, *b, *c;
    args_t *g;

    a = (dtype*)NDL_PTR(lp,0);
    b = (dtype*)NDL_PTR(lp,1);
    c = (dtype*)NDL_PTR(lp,2);
    g = (args_t*)(lp->opt_ptr);

    (*func_p)(g->order, g->side, g->uplo, g->m, g->n,
              DP(g->alpha), a, g->lda, b, g->ldb, DP(g->beta), c, g->ldc);
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER LDA,LDB,LDC,M,N
*       CHARACTER SIDE,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*A*B + beta*C,
*>
*> or
*>
*>    C := alpha*B*A + beta*C,
*>
*> where alpha and beta are scalars,  A is a symmetric matrix and  B and
*> C are  m by n matrices.
*> \endverbatim
*
*/
/*
  @overload dgemm( alpha, a, b [,beta, c, side, uplo, order] )
  @param [Numo::DFloat] a  >=2-dimentional NArray.
  @param [Numo::DFloat] b  >=2-dimentional NArray.
  @return [Numo::DFloat]
  @raise

  dgemm - performs matrix-matrix multiplication : C = A B
  where A, B and C are matrices,
  with A an m by k matrix, B a k by n matrix and
  C an m by n matrix.
*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE mod)
{
    int       i;
    VALUE     ans, arg_c = Qnil;
    VALUE     a, b, c, alpha, beta, side, uplo, order;
    narray_t *na1, *na2, *na3;
    int       mc, tmp;
    size_t    shape[2];
    ndfunc_arg_in_t ain[4] = {{cT,2},{cT,2},{OVERWRITE,2},{sym_init,0}};
    ndfunc_arg_out_t aout[1] = {{cT,2,shape}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, 3, 0, ain, aout};
    args_t g = {CblasRowMajor, CblasLeft, CblasUpper, m_one, m_zero};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    i = rb_scan_args(argc, argv, "35",
                     &alpha, &a, &b, &beta, &c, &side, &uplo, &order);
    switch (i) {
    case 8: g.order = option_order(order);
    case 7: g.uplo = option_uplo(uplo);
    case 6: g.side = option_side(side);
    case 5: arg_c = c;
    case 4: if (RTEST(beta)) {g.beta = m_num_to_data(beta);}
    }
    if (RTEST(alpha)) {g.alpha = m_num_to_data(alpha);}

    GetNArray(a,na1);
    GetNArray(b,na2);
    CHECK_DIM_GE(na1,2);
    CHECK_DIM_GE(na2,2);

    CHECK_SQUARE("a",na1);

    // row major                        L    R
    //ma    = na1->shape[na1->ndim-2]; // m or n
    g.lda = na1->shape[na1->ndim-1]; // m or n (lda)

    g.m   = na2->shape[na2->ndim-2]; // m
    g.n   = na2->shape[na2->ndim-1]; // n (ldb)
    g.ldb = g.n;

    if (g.side == CblasLeft) {
        CHECK_SIZE_EQ(g.lda,g.m);
    } else {
        CHECK_SIZE_EQ(g.lda,g.n);
    }

    if (arg_c == Qnil) { // c is not given.
        ndf.nout = 1;
        ain[2] = ain[3];
        c = INT2FIX(0);
        shape[0] = g.m;
        shape[1] = g.ldc = g.n;
    } else {
        c = rb_funcall(cT,rb_intern("cast"),1,arg_c);
        if (!TEST_INPLACE(c)) {
            c = na_copy(c);
        }
        GetNArray(c,na3);
        CHECK_DIM_GE(na3,2);
        mc    = na3->shape[na3->ndim-2]; // m
        g.ldc = na3->shape[na3->ndim-1]; // n (ldc)
        if (mc < g.m) {
            rb_raise(nary_eShapeError,"mc=%d must be >= m=%d",mc,g.n);
        }
        CHECK_LEADING_GE("ldc",g.ldc,"m",g.m);
    }

    if ( g.order != CblasRowMajor ) {
        swap(g.m,g.n);
    }

    ans = na_ndloop3(&ndf, &g, 3, a, b, c);

    if (arg_c == Qnil) { // c is not given.
        return ans;
    } else {
        return c;
    }
}

#undef func_p
#undef args_t
