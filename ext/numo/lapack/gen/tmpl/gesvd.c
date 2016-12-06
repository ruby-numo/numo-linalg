/*
* DGESVD computes the singular value decomposition (SVD) for GE
* matrices
*
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
* ZGESVD computes the singular value decomposition (SVD) for GE
* matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
*                          WORK, LWORK, RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBU, JOBVT
*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * ), S( * )
*       COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
*      $                   WORK( * )
*       ..
*/

#define gesvd FFUNC(<%=blas_char%>gesvd)

void gesvd(
  char const * /*JOBU*/, char const * /*JOBVT*/,
  fortran_integer * /*M*/, fortran_integer * /*N*/, dtype * /*A*/, fortran_integer * /*LDA*/,
  rtype * /*S*/,
  dtype * /*U*/, fortran_integer * /*LDU*/,
  dtype * /*VT*/, fortran_integer * /*LDVT*/,
  dtype * /*WORK*/, fortran_integer * /*LWORK*/,
  <% if is_complex %>
  rtype * /*RWORK*/,
  <% end %>
  fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const jobu="A", * const jobvt="A";
    volatile VALUE vopt;  // !!! DONT REMOVE: volatile qualification !!!
    dtype *a;
    dtype *u;
    rtype *s;
    dtype *vt;
    fortran_integer m, n, lda, ldu, ldvt, lwork, info=0;
    dtype *work;
    <% if is_complex %>
    rtype *rwork;
    <% end %>

    // a[n,lda], u[m,m], s[min(m,n)], vt[n,n]
    a     = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    u     = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    s     = (rtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    vt    = (dtype *)(lp->args[3].ptr + lp->args[3].iter[0].pos);
    lwork = *(fortran_integer *)lp->opt_ptr;
    m      = lp->args[1].shape[0];
    n      = lp->args[3].shape[0];
    {
        char *ptr;
        size_t wksize = (size_t)lwork;
        <% if is_complex %>
        size_t ofs[3];
        size_t min_mn;
        min_mn = lp->args[2].shape[0];
        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, wksize);    // work[lwork]
        SET_POS(ofs, 2, rtype, 5*min_mn);  // rwork[5*min_mn]
        ptr = rb_alloc_tmp_buffer(&vopt, ofs[2]);
        rwork = (rtype *)(ptr + ofs[1]);

        <% else %>

        size_t ofs[2];
        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, wksize);  // work[lwork]
        ptr = rb_alloc_tmp_buffer(&vopt, ofs[1]);
        <% end %>
        work = (dtype *)ptr;
    }
    lda   = m;
    ldu   = m;
    ldvt  = n;

    <% if is_complex %>
    gesvd(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, rwork, &info);
    <% else %>
    gesvd(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, &info);
    <% end %>

    rb_free_tmp_buffer(&vopt);
    RB_GC_GUARD(vopt);
}

#define sub_func_name(f, args) f##_sub args

/*
  @overload gesvd(a)
  @param [Numo::<%=class_name%>] a >=2-dimentional NArray.
  @return ***TBD***
  @raise

  <%=blas_char%>gesvd - computes the singular value decomposition (SVD) of a
  M-by-N Matrix A, optionally computing the left and/or right singular
  vectors. The SVD is written

       A = U * SIGMA * (conjugate-)transpose(V)

  where SIGMA is an M-by-N matrix which is zero except for its
  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
  are the singular values of A; they are real and non-negative, and
  are returned in descending order.  The first min(m,n) columns of
  U and V are the left and right singular vectors of A.
*/
static VALUE
sub_func_name(<%=c_func%>, (VALUE UNUSED(mod), VALUE a, int UNUSED(full), int UNUSED(with_uv)))
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
    gesvd(chr, chr, &m, &n, 0, &m, 0, 0, &m, 0, &n, wk, &lwork, 0, &info);
    lwork = (fortran_integer)REAL(wk[0]);
    <% else %>
    gesvd(chr, chr, &m, &n, 0, &m, 0, 0, &m, 0, &n, wk, &lwork, &info);
    lwork = (fortran_integer)wk[0];
    <% end %>

    ans = na_ndloop3(&ndf, &lwork, 1, a);

    return ans;
}

static VALUE
<%=c_func%>(VALUE mod, VALUE a)
{
    return sub_func_name(<%=c_func%>, (mod, a, 1, 1));
}

#undef sub_func_name
