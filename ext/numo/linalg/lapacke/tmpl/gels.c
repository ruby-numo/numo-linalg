/*
* DGELS solves overdetermined or underdetermined systems for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
* ZGELS solves overdetermined or underdetermined systems for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*/
<%
 is_lss = (/gels(s|d|y)/ =~ name)
 is_lsd = (/gelsd/ =~ name)
 is_lsy = (/gelsy/ =~ name)
 if is_complex
   if is_lsd
     rwork = ", &rwork_q"
   else
     rwork = ", 0"
   end
   g_rwork = ", g->rwork"
 else
   rwork = ""
   g_rwork = ""
 end
 if is_lsy
   out_t = "int"    # jpvt
   cOut = "cInt"
 else
   out_t = "rtype"  # s
   cOut = "cRT"
 end
 %>
#define IS_LSS <%=is_lss ? "1":"0"%>
#define IS_LSD <%=is_lsd ? "1":"0"%>
#define IS_LSY <%=is_lsy ? "1":"0"%>
#define args_t <%=func_name%>_args_t
#define func_p <%=func_name%>_p

typedef struct {
    int   order;
    char  trans;
    rtype rcond;
    int   lwork;
    dtype *work;
    rtype *rwork;
    int   *iwork;
} args_t;

static <%=func_name%>_work_t func_p = 0;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a, *b;
    int   *info;
    int    m, n, nb, nrhs, lda, ldb, tmp;
    args_t *g;
#if IS_LSS
    <%=out_t%> *s; // or jpvt
    int   *rank;
#endif

    a = (dtype*)NDL_PTR(lp,0);
    b = (dtype*)NDL_PTR(lp,1);
#if IS_LSS
    s = (<%=out_t%>*)NDL_PTR(lp,2);
    rank = (int*)NDL_PTR(lp,3);
    info = (int*)NDL_PTR(lp,4);
#else
    info = (int*)NDL_PTR(lp,2);
#endif
    g = (args_t*)(lp->opt_ptr);

    m = lp->args[0].shape[0];
    n = lp->args[0].shape[1];
    SWAP_IFCOL(g->order,m,n);
    lda = lp->args[0].iter[0].step / sizeof(dtype);

    if (lp->args[1].ndim == 1) {
        nrhs = 1;
        nb = lp->args[1].shape[0];
        ldb = (g->order==LAPACK_COL_MAJOR) ? nb : 1;
    } else {
        nrhs = lp->args[1].shape[0];
        nb = lp->args[1].shape[1];
        ldb = nrhs;
        SWAP_IFCOL(g->order,nb,nrhs);
    }

    //printf("order=%d trans=%c m=%d n=%d nb=%d nrhs=%d lda=%d ldb=%d\n",g->order,g->trans,m,n,nb,nrhs,lda,ldb);

#if IS_LSD
    *info = (*func_p)(g->order, m, n, nrhs, a, lda, b, ldb, s, g->rcond, rank,
                      g->work, g->lwork <%=g_rwork%>, g->iwork);
#elif IS_LSS
    *info = (*func_p)(g->order, m, n, nrhs, a, lda, b, ldb, s, g->rcond, rank,
                      g->work, g->lwork <%=g_rwork%> );
#else
    *info = (*func_p)(g->order, g->trans, m, n, nrhs, a, lda, b, ldb,
                      g->work, g->lwork);
#endif
    CHECK_ERROR(*info);
}

/*
  @overload <%=name%>(a, b)
  @param [Numo::NArray] a  >=2-dimentional NArray.
  @param [Numo::NArray] b  >=1-dimentional NArray.
  @return [Numo::NArray]

  <%=name%> - solves overdetermined or underdetermined real linear systems
  involving an M-by-N matrix A, or its transpose, using a QR or LQ factorization of A.
  It is assumed that A has full rank.

  1. If m >= n:  find the least squares solution of an overdetermined
     system, i.e., solve the least squares problem minimize || B - A*X ||.

  2. If m < n:  find the minimum norm solution of an underdetermined system A * X = B.
*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE const mod)
{
    VALUE order, tmpwork;
#if IS_LSS && IS_COMPLEX
    VALUE tmprwork;
    int   lrwork;
#endif
#if IS_LSD
    VALUE tmpiwork;
#endif
    VALUE a, b, ans;
    int   m, n, nb, nrhs, lda, ldb, lwork, info, tmp;
    int   max_mn;
    dtype work_q;
    narray_t *na1, *na2;
    VALUE opt1;
#if IS_LSY
    narray_t *na3;
    VALUE jpvt;
#endif
#if IS_LSS
    size_t shape_s[1];
    ndfunc_arg_in_t ain[3] = {{OVERWRITE,2},{OVERWRITE,2},{cInt,1}};
    ndfunc_arg_out_t aout[3] = {{cT,1,shape_s},{cInt,0},{cInt,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 2, 3, ain, aout};
#else
    ndfunc_arg_in_t ain[2] = {{OVERWRITE,2}, {OVERWRITE,2}};
    ndfunc_arg_out_t aout[1] = {{cInt,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 2, 1, ain, aout};
#endif
    args_t g = {0,0,-1};
    int i;

#if IS_LSY
    i = rb_scan_args(argc, argv, "32", &a, &b, &jpvt, &opt1, &order);
#define N 5
#else
    i = rb_scan_args(argc, argv, "22", &a, &b, &opt1, &order);
#define N 4
#endif
    switch (i) {
    case N: g.order = option_order(order);
#if IS_LSS
    case N-1: if (RTEST(opt1)) {g.rcond = NUM2DBL(opt1);}
#else
    case N-1: g.trans = option_trans(opt1);
#endif
#undef N
    }
    if (g.trans==0) { g.trans='N'; }
    if (g.order==0) { g.order=LAPACK_ROW_MAJOR; }

    CHECK_FUNC(func_p,"<%=func_name%>_work");

    COPY_OR_CAST_TO(a,cT);
    COPY_OR_CAST_TO(b,cT);

    //A is DOUBLE PRECISION array, dimension (LDA,N)
    //On entry, the M-by-N matrix A.
    GetNArray(a, na1);
    CHECK_DIM_GE(na1, 2);

    //B is DOUBLE PRECISION array, dimension (LDB,NRHS)
    //B is M-by-NRHS if TRANS = 'N'
    //     N-by-NRHS if TRANS = 'T'
    GetNArray(b, na2);
    CHECK_DIM_GE(na2, 1);

    //The number of rows of the matrix A.
    m = na1->shape[na1->ndim-2];
    //The number of columns of the matrix A.
    n = na1->shape[na1->ndim-1];
    lda = n;
    max_mn = (m > n) ? m : n;
    SWAP_IFCOL(g.order,m,n);

#if IS_LSS
# if IS_LSY
    ndf.nin++;
    ndf.nout--;
    ndf.aout++;
    COPY_OR_CAST_TO(jpvt,cInt);
    GetNArray(jpvt, na3);
    CHECK_DIM_GE(na3, 1);
    { int jpvt_sz = na3->shape[na3->ndim-1];
      CHECK_INT_EQ("jpvt_size",jpvt_sz,"n",n);
    }
# else
    shape_s[0] = (m < n) ? m : n;
# endif
#endif

    //The number of columns of the matrix B.
    nb = na2->shape[na2->ndim-1];
    if (na2->ndim == 1) {
        ain[1].dim = 1; // reduce dimension
        nrhs = 1;
        ldb = (g.order==LAPACK_COL_MAJOR) ? nb : 1;
    } else {
        //The number of rows of the matrix B.
        nrhs = na2->shape[na2->ndim-2];
        ldb = nrhs;
        SWAP_IFCOL(g.order,nb,nrhs);
    }

    if (nb < max_mn) {
        rb_raise(nary_eShapeError,
                 "ldb should be >= max(m,n): ldb=%d m=%d n=%d",nb,m,n);
    }
    //printf("order=%d trans=%c m=%d n=%d nb=%d nrhs=%d lda=%d ldb=%d\n",g.order,g.trans,m,n,nb,nrhs,lda,ldb);
    lwork = -1;
#if IS_LSD
    {
     int iwork_q;
# if IS_COMPLEX
     rtype rwork_q;
# endif
     info = (*func_p)(g.order, m, n, nrhs, 0, lda, 0, ldb, 0, g.rcond, 0,
                      &work_q, lwork <%=rwork%>, &iwork_q);
     CHECK_ERROR(info);
# if IS_COMPLEX
     lrwork = rwork_q;
     g.rwork = (rtype*)rb_alloc_tmp_buffer(&tmprwork, lrwork*sizeof(rtype));
# endif
     g.iwork = (int*)rb_alloc_tmp_buffer(&tmpiwork, iwork_q*sizeof(int));
    }
// end of IS_LSD

#elif IS_LSS
    info = (*func_p)(g.order, m, n, nrhs, 0, lda, 0, ldb, 0, g.rcond, 0,
                     &work_q, lwork <%=rwork%>);
    CHECK_ERROR(info);

# if IS_COMPLEX
#  if IS_LSY
    lrwork = 2 * n;
#  else
    lrwork = 5 * ((m < n) ? m : n);
#  endif
    g.rwork = (rtype*)rb_alloc_tmp_buffer(&tmprwork, lrwork*sizeof(rtype));
# endif
// end of IS_LSS

#else // when gels
    info = (*func_p)(g.order, g.trans, m, n, nrhs, 0, lda, 0, ldb,
                     &work_q, lwork);
    CHECK_ERROR(info);
#endif
    g.lwork = m_real(work_q);
    g.work = (dtype*)rb_alloc_tmp_buffer(&tmpwork, g.lwork*sizeof(dtype));

    // ndloop
#if IS_LSY
    ans = na_ndloop3(&ndf, &g, 3, a, b, jpvt);
#else
    ans = na_ndloop3(&ndf, &g, 2, a, b);
#endif

    // free_tmp_buffer
    rb_free_tmp_buffer(&tmpwork);
#if IS_LSD
    rb_free_tmp_buffer(&tmpiwork);
#endif
#if IS_LSS
# if IS_COMPLEX
    rb_free_tmp_buffer(&tmprwork);
# endif

    // return
# if IS_LSY
    return rb_ary_new3(5,a,b,jpvt,RARRAY_AREF(ans,0),RARRAY_AREF(ans,1));
# else
    return rb_ary_new3(5,a,b,RARRAY_AREF(ans,0),RARRAY_AREF(ans,1),RARRAY_AREF(ans,2));
# endif
#else
    return rb_ary_new3(3,a,b,ans);
#endif
}

#undef func_p
#undef args_t
#undef IS_LSS
#undef IS_LSD
#undef IS_LSY
