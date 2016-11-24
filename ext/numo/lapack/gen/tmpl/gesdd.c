/*
* DGESDD computes the singular value decomposition (SVD) of a real
* M-by-N matrix A, optionally computing the left and right singular
* vectors.  If singular vectors are desired, it uses a
* divide-and-conquer algorithm.
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK,
*                          LWORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ
*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
*      $                   VT( LDVT, * ), WORK( * )
*       ..
*
* ZGESDD computes the singular value decomposition (SVD) of a complex
* M-by-N matrix A, optionally computing the left and/or right singular
* vectors, by using divide-and-conquer method.
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK,
*                          LWORK, RWORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ
*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   RWORK( * ), S( * )
*       COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
*      $                   WORK( * )
*       ..
*/

#define gesdd FFUNC(<%=blas_char%>gesdd)

void gesdd(
  char const * /*JOBZ*/,
  fortran_integer * /*M*/, fortran_integer * /*N*/, dtype * /*A*/, fortran_integer * /*LDA*/,
  rtype * /*S*/,
  dtype * /*U*/, fortran_integer * /*LDU*/,
  dtype * /*VT*/, fortran_integer * /*LDVT*/,
  dtype * /*WORK*/, fortran_integer * /*LWORK*/,
  <% if is_complex %>
  rtype * /*RWORK*/,
  <% end %>
  fortran_integer * /*IWORK*/,
  fortran_integer * /*INFO*/);

typedef struct {
    int with_uv;
    fortran_integer lwork;
} gesdd_opt_t;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const job_a="A", * const job_n="N";
    char const *job;
    volatile VALUE vopt;  // !!! DONT REMOVE: volatile qualification !!!
    gesdd_opt_t *opt;
    dtype *a;
    dtype *u;
    rtype *s;
    dtype *vt;
    fortran_integer m, n, lda, ldu, ldvt, lwork, info=0;
    dtype *work;
    <% if is_complex %>
    rtype *rwork;
    <% end %>
    fortran_integer *iwork;

    // a[n,lda], u[m,m], s[min(m,n)], vt[n,n]
    a     = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    u     = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    s     = (rtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    vt    = (dtype *)(lp->args[3].ptr + lp->args[3].iter[0].pos);
    m     =           lp->args[1].shape[0];
    n     =           lp->args[3].shape[0];
    opt   =           lp->opt_ptr;
    lwork = opt->lwork;
    job   = ( opt->with_uv ? job_a : job_n ) ;
    {
        char *ptr;
        size_t min_mn;
        <% if is_complex %>
        size_t max_mn, lrwork;
        size_t ofs[4];
        <% else %>
        size_t ofs[3];
        <% end %>

        min_mn = lp->args[2].shape[0];
        <% if is_complex %>
        max_mn = max_(m, n);
        lrwork = min_mn * max_(5*min_mn+7, 2*max_mn+2*min_mn+1);
        <% end %>

        <% if is_complex %>
        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, lwork);     // work[lwork]
        SET_POS(ofs, 2, rtype, lrwork);  // rwork[lrwork]
        SET_POS(ofs, 3, fortran_integer, 8*min_mn);  // iwork[8*min_mn]
        ptr = rb_alloc_tmp_buffer(&vopt, ofs[3]);
        work  =           (dtype *) ptr;
        rwork =           (rtype *)(ptr + ofs[1]);
        iwork = (fortran_integer *)(ptr + ofs[2]);

        <% else %>

        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, lwork);     // work[lwork]
        SET_POS(ofs, 2, fortran_integer, 8*min_mn);  // iwork[8*min_mn]
        ptr = rb_alloc_tmp_buffer(&vopt, ofs[2]);
        work  =           (dtype *) ptr;
        iwork = (fortran_integer *)(ptr + ofs[1]);
        <% end %>
    }
    lda   = m;
    ldu   = m;
    ldvt  = n;

    <% if is_complex %>
    gesdd(job, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, rwork, iwork, &info);
    <% else %>
    gesdd(job, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, iwork, &info);
    <% end %>

    rb_free_tmp_buffer(&vopt);
    RB_GC_GUARD(vopt);
}

#define sub_func_name(f, args) f##_sub args

/*
  @overload gesdd(a)
  @param [Numo::<%=class_name%>] a >=2-dimentional NArray.
  @return ***TBD***
  @raise

  TBD
*/
static VALUE
sub_func_name(<%=c_func%>, (VALUE UNUSED(mod), VALUE a, int with_uv))
{
    char const * const chr = "N";
    dtype wk[1];
    fortran_integer m, n, min_mn, lwork, info=0;
    volatile VALUE ans;  // !!! DONT REMOVE: volatile qualification !!!

    narray_t *na;
    size_t u_shape[2];
    size_t s_shape[1];
    size_t vt_shape[2];
    ndfunc_arg_in_t ain[] = {{cT,2}};
    ndfunc_arg_out_t aout[] = {
      {cT,COUNT_OF_(u_shape),u_shape},
      {cRT,COUNT_OF_(s_shape),s_shape},
      {cT,COUNT_OF_(vt_shape),vt_shape}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};

    GetNArray(a, na);
    CHECK_DIM_GE(na, 2);
    m = na->shape[na->ndim-1];
    n = na->shape[na->ndim-2];
    min_mn = min_(m, n);
    u_shape[0] = m;
    u_shape[1] = m;
    s_shape[0] = min_mn;
    vt_shape[0] = n;
    vt_shape[1] = n;

    lwork = -1;
    <% if is_complex %>
    gesdd(chr, &m, &n, 0, &m, 0, 0, &m, 0, &n, wk, &lwork, 0, 0, &info);
    lwork = REAL(wk[0]);
    <% else %>
    gesdd(chr, &m, &n, 0, &m, 0, 0, &m, 0, &n, wk, &lwork, 0, &info);
    lwork = wk[0];
    <% end %>

    {
        volatile VALUE vopt;
        gesdd_opt_t *opt;

        opt = rb_alloc_tmp_buffer(&vopt, sizeof(gesdd_opt_t));
        opt->with_uv = with_uv;
        opt->lwork = lwork;
        ans = na_ndloop3(&ndf, opt, 1, a);
        rb_free_tmp_buffer(&vopt);
        RB_GC_GUARD(vopt);
    }

    return ans;
}

static VALUE
<%=c_func%>(VALUE mod, VALUE a)
{
    return sub_func_name(<%=c_func%>, (mod, a, 1));
}

#undef sub_func_name
