#define args_t <%=name%>_args_t

typedef struct {
    enum CBLAS_ORDER order;
    enum CBLAS_UPLO uplo;
    enum CBLAS_TRANSPOSE trans;
    dtype alpha, beta;
    blasint n, k, lda, ldb, ldc;
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

    (*func_p)(g->order, g->uplo, g->trans, g->n, g->k,
              DP(g->alpha), a, g->lda, b, g->ldb, DP(g->beta), c, g->ldc);
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,N
*       CHARACTER TRANS,UPLO
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
*> DSYR2K  performs one of the symmetric rank 2k operations
*>
*>    C := alpha*A*B**T + alpha*B*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*B + alpha*B**T*A + beta*C,
*>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
*> matrices in the second case.
*> \endverbatim
*/
/*
 * @overload <%=name%>( alpha, a, b [,beta, c, uplo, trans, order] )
 * @param [Numeric]      alpha
 * @param [Numo::DFloat] a  n-by-k matrix (>=2-dimentional NArray)
 * @param [Numo::DFloat] b  n-by-k matrix (>=2-dimentional NArray)
 * @param [Numeric]      beta (default=0)
 * @param [Numo::DFloat] c  n-by-n matrix [in/out] (>=2-dimentional NArray)
 * @param [option] uplo  (default='upper')
 * @param [option] trans (default='notrans')
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
    VALUE     ans, arg_c = Qnil;
    VALUE     a, b, c, alpha, beta, uplo, trans, order;
    narray_t *na1, *na2, *na3;
    blasint   na, ka, kb, nb, nc, tmp;
    size_t    shape[2];
    ndfunc_arg_in_t ain[4] = {{cT,2},{cT,2},{OVERWRITE,2},{sym_init,0}};
    ndfunc_arg_out_t aout[1] = {{cT,2,shape}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, 3, 0, ain, aout};
    args_t g = {CblasRowMajor, CblasUpper, CblasNoTrans, m_one, m_zero};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    i = rb_scan_args(argc, argv, "35",
                     &alpha, &a, &b, &beta, &c, &uplo, &trans, &order);
    switch (i) {
    case 8: g.order = option_order(order);
    case 7: g.trans = option_trans(trans);
    case 6: g.uplo  = option_uplo(uplo);
    case 5: arg_c = c;
    case 4: if (RTEST(beta)) {g.beta = m_num_to_data(beta);}
    }
    if (RTEST(alpha)) {g.alpha = m_num_to_data(alpha);}

    GetNArray(a,na1);
    GetNArray(b,na2);
    CHECK_DIM_GE(na1,2);
    CHECK_DIM_GE(na2,2);

    na = na1->shape[na1->ndim-2]; // n
    ka = na1->shape[na1->ndim-1]; // k (lda)
    g.lda = ka;
    if ( (g.order==CblasRowMajor && g.trans!=CblasNoTrans) ||
         (g.order!=CblasRowMajor && g.trans==CblasNoTrans) ) {
        swap(na,ka);
    }

    nb = na2->shape[na2->ndim-2]; // n
    kb = na2->shape[na2->ndim-1]; // k (ldb)
    g.ldb = kb;
    if ( (g.order==CblasRowMajor && g.trans!=CblasNoTrans) ||
         (g.order!=CblasRowMajor && g.trans==CblasNoTrans) ) {
        swap(kb,nb);
    }
    CHECK_INT_EQ("na",na,"nb",nb);
    CHECK_INT_EQ("ka",ka,"kb",kb);
    g.n = nb;
    g.k = kb;

    if ( g.order==CblasRowMajor ) {
        swap(na,nb);
    }

    if (arg_c == Qnil) { // c is not given.
        ndf.nout = 1;
        ain[2] = ain[3];
        c = INT2FIX(0);
        shape[0] = nb;
        shape[1] = g.ldc = na;
    } else {
        c = rb_funcall(cT,rb_intern("cast"),1,arg_c);
        if (!TEST_INPLACE(c)) {
            c = na_copy(c);
        }
        GetNArray(c,na3);
        CHECK_DIM_GE(na3,2);
        nc = na3->shape[na3->ndim-2];    // n
        g.ldc = na3->shape[na3->ndim-1]; // n (ldc)
        if (nc < nb) {
            rb_raise(nary_eShapeError,"nc=%d must be >= nb=%d",nc,nb);
        }
        CHECK_LEADING_GE("ldc",g.ldc,"n",na);
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
