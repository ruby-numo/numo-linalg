#define args_t <%=name%>_args_t

typedef struct {
  enum CBLAS_ORDER order;
  enum CBLAS_TRANSPOSE trans_a, trans_b;
  dtype alpha, beta;
  blasint m, n, k, lda, ldb, ldc;
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

    (*func_p)(g->order, g->trans_a, g->trans_b, g->m, g->n, g->k,
              DP(g->alpha), a, g->lda, b, g->ldb, DP(g->beta), c, g->ldc);
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*               C := alpha*op( A )*op( B ) + beta*C,

*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*/
/*
  @overload dgemm( alpha, a, b [,beta, c, transa, transb, order] )
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
    VALUE     a, b, c, alpha, beta, trans_a, trans_b, order;
    narray_t *na1, *na2, *na3;
    blasint   ma, ka, kb, nb, nc, tmp;
    size_t    shape[2];
    ndfunc_arg_in_t ain[4] = {{cT,2},{cT,2},{OVERWRITE,2},{sym_init,0}};
    ndfunc_arg_out_t aout[1] = {{cT,2,shape}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, 3, 0, ain, aout};
    args_t g = {CblasRowMajor, CblasNoTrans, CblasNoTrans, m_one, m_zero};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    i = rb_scan_args(argc, argv, "35",
                     &alpha, &a, &b, &beta, &c, &trans_a, &trans_b, &order);
    switch (i) {
    case 8: g.order   = option_order(order);
    case 7: g.trans_b = option_trans(trans_b);
    case 6: g.trans_a = option_trans(trans_a);
    case 5: arg_c = c;
    case 4: if (RTEST(beta)) {g.beta = m_num_to_data(beta);}
    }
    if (RTEST(alpha)) {g.alpha = m_num_to_data(alpha);}

    GetNArray(a,na1);
    GetNArray(b,na2);
    CHECK_DIM_GE(na1,2);
    CHECK_DIM_GE(na2,2);

    ma = na1->shape[na1->ndim-2]; // m
    ka = na1->shape[na1->ndim-1]; // k (lda)
    g.lda = ka;
    if ( (g.order==CblasRowMajor && g.trans_a!=CblasNoTrans) ||
         (g.order!=CblasRowMajor && g.trans_a==CblasNoTrans) ) {
        swap(ma,ka);
    }
    g.m = ma;

    kb = na2->shape[na2->ndim-2]; // k
    nb = na2->shape[na2->ndim-1]; // n (ldb)
    g.ldb = nb;
    if ( (g.order==CblasRowMajor && g.trans_b!=CblasNoTrans) ||
         (g.order!=CblasRowMajor && g.trans_b==CblasNoTrans) ) {
        swap(kb,nb);
    }
    g.n = nb;
    CHECK_INT_EQ("ka",ka,"kb",kb);
    //g.k = (ka < kb) ? ka : kb;

    if ( g.order==CblasRowMajor ) {
        swap(ma,nb);
    }

    if (arg_c == Qnil) { // c is not given.
        ndf.nout = 1;
        ain[2] = ain[3];
        c = INT2FIX(0);
        shape[0] = nb;
        shape[1] = g.ldc = ma;
    } else {
        c = rb_funcall(cT,rb_intern("cast"),1,arg_c);
        if (!TEST_INPLACE(c)) {
            c = na_copy(c);
        }
        GetNArray(c,na3);
        CHECK_DIM_GE(na3,2);
        nc = na3->shape[na3->ndim-2];    // n
        g.ldc = na3->shape[na3->ndim-1]; // m (ldc)
        if (nc < nb) {
            rb_raise(nary_eShapeError,"nc=%d must be >= nb=%d",nc,nb);
        }
        CHECK_LEADING_GE("ldc",g.ldc,"m",ma);
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
