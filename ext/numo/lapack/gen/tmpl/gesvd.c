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

typedef struct {
    int vals_only;
    size_t lwork;
} gesvd_opt_t;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    volatile VALUE tmp_ptr;  // !!! DONT REMOVE: volatile qualification !!!
    dtype *a, *a2, *u;
    rtype *s;
    dtype *vt;
    size_t n1, n2, min_n;
    dtype *work;
    <% if is_complex %>
    rtype *rwork;
    <% end %>
    fortran_integer lwork;

    gesvd_opt_t *opt = lp->opt_ptr;

    if (opt->vals_only) {
        // a[n,lda], u[m,m], s[min(m,n)], vt[n,n]
        a     = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
        u     = 0;
        s     = (rtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
        vt    = 0;
        min_n =           lp->args[1].shape[0];
    } else {
        // a[n,lda], u[m,m], s[min(m,n)], vt[n,n]
        a     = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
        u     = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
        s     = (rtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
        vt    = (dtype *)(lp->args[3].ptr + lp->args[3].iter[0].pos);
        min_n =           lp->args[2].shape[0];
    }
    n1 = lp->args[0].shape[1];
    n2 = lp->args[0].shape[0];
    {
        char *ptr;
        <% if is_complex %>
        size_t ofs[4];
        <% else %>
        size_t ofs[3];
        <% end %>
        size_t const wksize = opt->lwork;
        lwork = (fortran_integer)wksize;

        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, n1*n2);  // a2[n2,n1]
        SET_POS(ofs, 2, dtype, wksize);  // work[lwork]
        <% if is_complex %>
        SET_POS(ofs, 3, rtype, 5*min_n);  // rwork[5*min_mn]
        ptr = rb_alloc_tmp_buffer(&tmp_ptr, ofs[3]);
        rwork = (rtype *)(ptr + ofs[2]);
        <% else %>
        (void)min_n;  // DUMMY
        ptr = rb_alloc_tmp_buffer(&tmp_ptr, ofs[2]);
        <% end %>
        work  = (dtype *)(ptr + ofs[1]);
        a2    = (dtype *) ptr;

        memcpy(a2, a, sizeof(dtype)*n1*n2);
    }
    {
        char const * const job_a="A", * const job_n="N";
        char const *jobu, *jobvt;
        fortran_integer m = (fortran_integer)n1, n = (fortran_integer)n2;
        fortran_integer lda = m, ldu, ldvt, info=0;
        if (opt->vals_only) {
            jobu  = job_n;
            jobvt = job_n;
            ldu   = 1;
            ldvt  = 1;
        } else {
            jobu  = job_a;
            jobvt = job_a;
            ldu   = m;
            ldvt  = n;
        }
        <% if is_complex %>
        gesvd(jobu, jobvt, &m, &n, a2, &lda, s, u, &ldu, vt, &ldvt,
              work, &lwork, rwork, &info);
        <% else %>
        gesvd(jobu, jobvt, &m, &n, a2, &lda, s, u, &ldu, vt, &ldvt,
              work, &lwork, &info);
        <% end %>
    }
    rb_free_tmp_buffer(&tmp_ptr);
    RB_GC_GUARD(tmp_ptr);
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
sub_func_name(<%=c_func%>, (VALUE a, int const vals_only))
{
    volatile VALUE ans;  // !!! DONT REMOVE: volatile qualification !!!
    gesvd_opt_t opt;
    size_t n1, n2;
    {
        narray_t *na;

        GetNArray(a, na);
        CHECK_DIM_GE(na, 2);
        n1 = na->shape[na->ndim-1];
        n2 = na->shape[na->ndim-2];
    }
    {
        char const * const chr_a = "A", * const chr_n = "N";
        char const *chr;
        fortran_integer m = (fortran_integer)n1, n = (fortran_integer)n2;
        fortran_integer ldu, ldvt;
        dtype wk[1];
        fortran_integer lwork, info=0;

        if (vals_only) {
            chr = chr_n;
            ldu = 1;
            ldvt = 1;
        } else {
            chr = chr_a;
            ldu = m;
            ldvt = n;
        }
        lwork = -1;
        <% if is_complex %>
        gesvd(chr, chr, &m, &n, 0, &m, 0, 0, &ldu, 0, &ldvt, wk, &lwork, 0, &info);
        opt.lwork = (size_t)REAL(wk[0]);
        <% else %>
        gesvd(chr, chr, &m, &n, 0, &m, 0, 0, &ldu, 0, &ldvt, wk, &lwork, &info);
        opt.lwork = (size_t)wk[0];
        <% end %>
    }
    opt.vals_only = vals_only;

    if (vals_only) {
        size_t s_shape[1] = { min_(n1, n2) };
        ndfunc_arg_in_t ain[1] = {{cT,2}};
        ndfunc_arg_out_t aout[1] = {
            {cRT,COUNT_OF_(s_shape),s_shape}};
        ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};

        ans = na_ndloop3(&ndf, &opt, 1, a);
    } else {
        size_t u_shape[2] = { n1, n1 };
        size_t s_shape[1] = { min_(n1, n2) };
        size_t vt_shape[2] = { n2, n2 };
        ndfunc_arg_in_t ain[1] = {{cT,2}};
        ndfunc_arg_out_t aout[3] = {
            {cT,COUNT_OF_(u_shape),u_shape},
            {cRT,COUNT_OF_(s_shape),s_shape},
            {cT,COUNT_OF_(vt_shape),vt_shape}};
        ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};

        ans = na_ndloop3(&ndf, &opt, 1, a);
    }
    return ans;
}

static VALUE
<%=c_func%>(int const argc, VALUE const argv[], VALUE UNUSED(mod))
{
    int flg_vals_only=0;

    rb_check_arity(argc, 1, 2);

    if (argc == 2) {
        ID tbl;
        VALUE v;
        tbl = rb_intern("vals_only");
        rb_get_kwargs(argv[1], &tbl, 0, 1, &v);
        if (v != Qundef) {
            flg_vals_only = RTEST(v);
        }
    }
    return sub_func_name(<%=c_func%>, (argv[0], flg_vals_only));
}

#undef sub_func_name
