/*
*  Definition:
*  ===========
*
*       SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
*                          WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBU, JOBVT
*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
*      $                   VT( LDVT, * ), WORK( * )
*       ..
*

lapack_int LAPACKE_dgesvd_work( int matrix_layout, char jobu, char jobvt,
                                lapack_int m, lapack_int n, double* a,
                                lapack_int lda, double* s, double* u,
                                lapack_int ldu, double* vt, lapack_int ldvt,
                                double* work, lapack_int lwork )

lapack_int LAPACKE_dgesdd_work( int matrix_layout, char jobz, lapack_int m,
                                lapack_int n, double* a, lapack_int lda,
                                double* s, double* u, lapack_int ldu,
                                double* vt, lapack_int ldvt, double* work,
                                lapack_int lwork, lapack_int* iwork )

lapack_int LAPACKE_zgesvd_work( int matrix_layout, char jobu, char jobvt,
                                lapack_int m, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                double* s, lapack_complex_double* u,
                                lapack_int ldu, lapack_complex_double* vt,
                                lapack_int ldvt, lapack_complex_double* work,
                                lapack_int lwork, double* rwork )
<% is_sdd = (/gesdd/ =~ name) %>
*/
#define args_t <%=func_name%>_args_t
typedef struct {
    int order;
    char jobu, jobvt, jobz;
    int  lwork;
    dtype *work;
    rtype *rwork;
    int   *iwork;
} args_t;

#define func_p <%=func_name%>_p
static <%=func_name%>_work_t func_p = 0;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a, *u=0, *vt=0;
    rtype *s;
    int   *info;
    int    m, n, lda, ldu, ldvt, tmp;
    args_t *g;

    a = (dtype*)NDL_PTR(lp,0);
    s = (rtype*)NDL_PTR(lp,1);
    u = (dtype*)NDL_PTR(lp,2);
    vt = (dtype*)NDL_PTR(lp,3);
    info = (int*)NDL_PTR(lp,4);
    g = (args_t*)(lp->opt_ptr);

    m = lp->args[0].shape[0];
    n = lp->args[0].shape[1];
    SWAP_IFCOL(g->order,m,n);
    lda = lp->args[0].iter[0].step / sizeof(dtype);
    ldu = lp->args[2].iter[0].step / sizeof(dtype);
    if (ldu == 0) { ldu = m; } // jobu == 'O' or 'N'
    ldvt = lp->args[3].iter[0].step / sizeof(dtype);
    if (ldvt == 0) { ldvt = n; } // jobvt == 'O' or 'N'

    //printf("order=%d jobu=%c jobvt=%c jobz=%c m=%d n=%d lda=%d ldu=%d ldvt=%d\n",g->order,g->jobu, g->jobvt,g->jobz, m,n,lda,ldu,ldvt);

    <% job = (is_sdd) ? "g->jobz" : "g->jobu, g->jobvt" %>
    <% rwork = (is_complex) ? ", g->rwork" : "" %>
    <% iwork = (is_sdd) ? ", g->iwork" : "" %>
    *info = (*func_p)( g->order, <%=job%>,
                       m, n, a, lda, s, u, ldu, vt, ldvt,
                       g->work, g->lwork <%=rwork%> <%=iwork%> );
    CHECK_ERROR(*info);
}

/*
*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE const mod)
{
    VALUE order, tmpwork=Qnil;
  <% if is_complex %>
    VALUE tmprwork=Qnil;
    int lrwork;
  <% end %>
  <% if is_sdd %>
    VALUE jobz, tmpiwork=Qnil;
  <% else %>
    VALUE jobu, jobvt;
  <% end %>
    VALUE a, ans;
    int   m, n, min_mn, lda, ldu, ldvt, info, tmp;
    dtype work_q;
    narray_t *na1;
    size_t shape_s[1], shape_u[2], shape_vt[2];
    ndfunc_arg_in_t ain[1] = {{OVERWRITE,2}};
    ndfunc_arg_out_t aout[4] = {{cRT,1,shape_s},{cT,2,shape_u},
                                {cT,2,shape_vt},{cInt,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 1, 4, ain, aout};
    int i;
    args_t g = {0,0,0,0};

    check_func((void*)(&func_p),"<%=func_name%>_work");

  <% if is_sdd %>
    i = rb_scan_args(argc, argv, "12", &a, &jobz, &order);
    switch (i) {
    case 3: g.order = option_order(order);
    case 2: g.jobz  = option_job(jobz);
    }
    // set default
    if (g.jobz==0)  {g.jobz='A';}
    g.jobu = g.jobvt = g.jobz;
  <% else %>
    i = rb_scan_args(argc, argv, "13", &a, &jobu, &jobvt, &order);
    switch (i) {
    case 4: g.order = option_order(order);
    case 3: g.jobvt = option_job(jobvt);
    case 2: g.jobu  = option_job(jobu);
    }
    // set default
    if (g.jobu==0)  {g.jobu='A';}
    if (g.jobvt==0) {g.jobvt='A';}
    if (g.jobu=='O' && g.jobvt=='O') {
        rb_raise(rb_eArgError,"JOBVT and JOBU cannot both be 'O'");
    }
  <% end %>
    if (g.order==0) {g.order=LAPACK_ROW_MAJOR;}

    if (g.jobu=='O' || g.jobvt=='O') {
        if (CLASS_OF(a) != cT) {
            rb_raise(rb_eTypeError,"type of matrix a is invalid for overwrite");
        }
    } else {
        if (CLASS_OF(a) == cT) {
            a = na_copy(a);
        } else {
            a = rb_funcall(cT,rb_intern("cast"),1,a);
        }
    }

    GetNArray(a, na1);
    CHECK_DIM_GE(na1, 2);
    m = na1->shape[na1->ndim-2];
    n = na1->shape[na1->ndim-1];
    lda = n;
    SWAP_IFCOL(g.order,m,n);

  <% if is_sdd %>
    if (g.jobz=='O') {
        if (m >= n) { g.jobvt='A';} else { g.jobu='A';}
    }
  <% end %>

    // output S
    shape_s[0] = min_mn = (m < n) ? m : n;

    // output U
    switch(g.jobu){
    case 'A':
        shape_u[0] = m;
        shape_u[1] = m;
        ldu = m;
        break;
    case 'S':
        shape_u[0] = m;
        shape_u[1] = min_mn;
        SWAP_IFCOL(g.order,shape_u[0],shape_u[1]);
        ldu = shape_u[1];
        break;
    case 'O':
    case 'N':
        aout[1].dim = 0; // dummy
        ldu = m;
        break;
    default:
        rb_raise(rb_eArgError,"invalid option: jobu='%c'",g.jobu);
    }
    // output VT
    switch(g.jobvt){
    case 'A':
        shape_vt[0] = n;
        shape_vt[1] = n;
        ldvt = n;
        break;
    case 'S':
        shape_vt[0] = min_mn;
        shape_vt[1] = n;
        SWAP_IFCOL(g.order, shape_vt[0], shape_vt[1]);
        ldvt = shape_vt[1];
        break;
    case 'O':
    case 'N':
        aout[2].dim = 0; // dummy
        ldvt = n;
        break;
    default:
        rb_raise(rb_eArgError,"invalid option: jobvt='%c'",g.jobvt);
    }

    //printf("order=%d jobu=%c jobvt=%c jobz=%c m=%d n=%d lda=%d ldu=%d ldvt=%d\n",g.order,g.jobu, g.jobvt,g.jobz, m,n,lda,ldu,ldvt);

    <% job = (is_sdd) ? "g.jobz" : "g.jobu, g.jobvt" %>
    <% rwork = (is_complex) ? ", 0" : "" %>
    <% iwork = (is_sdd) ? ", 0" : "" %>
    info = (*func_p)( g.order, <%=job%>,
                      m, n, 0, lda, 0, 0, ldu, 0, ldvt, &work_q, -1
                      <%=rwork%> <%=iwork%> );
    //printf("info=%d\n",info);
    CHECK_ERROR(info);
    g.lwork = m_real(work_q);
    //printf("lwork=%d\n",g.lwork);
    g.work = (dtype*)rb_alloc_tmp_buffer(&tmpwork, g.lwork*sizeof(dtype));
   <% if is_complex %>
    <% if is_sdd %>
    if (g.jobz=='N') {
        lrwork = 7*min_mn;
    } else {
        int mx = (m > n) ? m : n;
        int mn = min_mn;
        int lm = 5*mn*mn + 5*mn;
        int ln = 2*mx*mn + 2*mn*mn + mn;
        lrwork = (lm > ln) ? lm : ln;
    }
    <% else %>
    lrwork = 5*min_mn;
    <% end %>
    g.rwork = (rtype*)rb_alloc_tmp_buffer(&tmprwork, lrwork*sizeof(rtype));
   <% end %>
   <% if is_sdd %>
    g.iwork = (int*)rb_alloc_tmp_buffer(&tmpiwork, 8*min_mn*sizeof(int));
   <% end %>

    ans = na_ndloop3(&ndf, &g, 1, a);

   <% if is_sdd %>
    rb_free_tmp_buffer(&tmpiwork);
   <% end %>
   <% if is_complex %>
    rb_free_tmp_buffer(&tmprwork);
   <% end %>
    rb_free_tmp_buffer(&tmpwork);
    if (aout[2].dim == 0) { rb_ary_delete_at(ans,2); }
    if (aout[1].dim == 0) { rb_ary_delete_at(ans,1); }
    return ans;
}

#undef args_t
#undef func_p
