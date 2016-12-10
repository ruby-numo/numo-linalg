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
    svd_job job;
    size_t lwork;
} gesdd_opt_t;

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
    fortran_integer lwork, *iwork;

    gesdd_opt_t *opt = lp->opt_ptr;

    switch (opt->job) {
    case SVD_FULL:
        /*FALLTHROUGH*/
    case SVD_THIN:
        // a[n,lda], u[m,m], s[min(m,n)], vt[n,n]
        a     = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
        u     = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
        s     = (rtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
        vt    = (dtype *)(lp->args[3].ptr + lp->args[3].iter[0].pos);
        min_n =           lp->args[2].shape[0];
        break;
    case SVD_VALS_ONLY:
        // a[n,lda], s[min(m,n)] (m=lda)
        a     = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
        u     = 0;
        s     = (rtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
        vt    = 0;
        min_n =           lp->args[1].shape[0];
        break;
    default:
        assert(0);
    }
    n1 = lp->args[0].shape[1];
    n2 = lp->args[0].shape[0];
    {
        char *ptr;
        <% if is_complex %>
        size_t max_n, lrwork;
        size_t ofs[5];
        <% else %>
        size_t ofs[4];
        <% end %>
        size_t const wksize = opt->lwork;
        lwork = (fortran_integer)wksize;

        <% if is_complex %>
        max_n = max_(n1, n2);
        lrwork = min_n * max_(5*min_n+7, 2*max_n+2*min_n+1);
        <% end %>

        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, n1*n2);  // a2[n2,n1]
        SET_POS(ofs, 2, dtype, wksize);  // work[lwork]
        <% if is_complex %>
        SET_POS(ofs, 3, rtype, lrwork);  // rwork[lrwork]
        SET_POS(ofs, 4, fortran_integer, 8*min_n);  // iwork[8*min_n]
        ptr = rb_alloc_tmp_buffer(&tmp_ptr, ofs[4]);
        iwork = (fortran_integer *)(ptr + ofs[3]);
        rwork =           (rtype *)(ptr + ofs[2]);
        <% else %>
        SET_POS(ofs, 3, fortran_integer, 8*min_n);  // iwork[8*min_n]
        ptr = rb_alloc_tmp_buffer(&tmp_ptr, ofs[3]);
        iwork = (fortran_integer *)(ptr + ofs[2]);
        <% end %>
        work  =           (dtype *)(ptr + ofs[1]);
        a2    =           (dtype *) ptr;

        memcpy(a2, a, sizeof(dtype)*n1*n2);
    }
    {
        char const * const job_a="A", * const job_s="S", * const job_n="N";
        char const *job;
        fortran_integer m = (fortran_integer)n1, n = (fortran_integer)n2;
        fortran_integer k = (fortran_integer)min_n;
        fortran_integer lda = m, ldu, ldvt, info=0;
        switch (opt->job) {
        case SVD_FULL:
            job  = job_a;
            ldu  = m;
            ldvt = n;
            break;
        case SVD_THIN:
            job  = job_s;
            ldu  = m;
            ldvt = k;
            break;
        case SVD_VALS_ONLY:
            job  = job_n;
            ldu  = 1;
            ldvt = 1;
            break;
        default:
            assert(0);
        }
        <% if is_complex %>
        gesdd(job, &m, &n, a2, &lda, s, u, &ldu, vt, &ldvt,
              work, &lwork, rwork, iwork, &info);
        <% else %>
        gesdd(job, &m, &n, a2, &lda, s, u, &ldu, vt, &ldvt,
              work, &lwork, iwork, &info);
        <% end %>
    }
    rb_free_tmp_buffer(&tmp_ptr);
    RB_GC_GUARD(tmp_ptr);
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
sub_func_name(<%=c_func%>, (VALUE a, svd_job job))
{
    volatile VALUE ans;  // !!! DONT REMOVE: volatile qualification !!!
    gesdd_opt_t opt;
    size_t n1, n2, min_n;
    {
        narray_t *na;

        GetNArray(a, na);
        CHECK_DIM_GE(na, 2);
        n1 = na->shape[na->ndim-1];
        n2 = na->shape[na->ndim-2];
        min_n = min_(n1, n2);
    }
    {
        char const * const chr_a = "A", * const chr_s = "S", * const chr_n = "N";
        char const *chr;
        fortran_integer m = (fortran_integer)n1, n = (fortran_integer)n2;
        fortran_integer k = (fortran_integer)min_n;
        fortran_integer ldu, ldvt;
        dtype wk[1];
        fortran_integer lwork, info=0;

        switch (job) {
        case SVD_FULL:
            chr = chr_a;
            ldu = m;
            ldvt = n;
            break;
        case SVD_THIN:
            chr = chr_s;
            ldu = m;
            ldvt = k;
            break;
        case SVD_VALS_ONLY:
            chr = chr_n;
            ldu = 1;
            ldvt = 1;
            break;
        default:
            assert(0);
        }
        lwork = -1;
        <% if is_complex %>
        gesdd(chr, &m, &n, 0, &m, 0, 0, &ldu, 0, &ldvt, wk, &lwork, 0, 0, &info);
        opt.lwork = (size_t)REAL(wk[0]);
        <% else %>
        gesdd(chr, &m, &n, 0, &m, 0, 0, &ldu, 0, &ldvt, wk, &lwork, 0, &info);
        opt.lwork = (size_t)wk[0];
        <% end %>
    }
    opt.job = job;

    switch (job) {
    case SVD_FULL:
        {
            size_t u_shape[2] = { n1, n1 };
            size_t s_shape[1] = { min_n };
            size_t vt_shape[2] = { n2, n2 };
            ndfunc_arg_in_t ain[1] = {{cT,2}};
            ndfunc_arg_out_t aout[3] = {
                {cT,COUNT_OF_(u_shape),u_shape},
                {cRT,COUNT_OF_(s_shape),s_shape},
                {cT,COUNT_OF_(vt_shape),vt_shape}};
            ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};

            ans = na_ndloop3(&ndf, &opt, 1, a);
        }
        break;
    case SVD_THIN:
        {
            size_t u_shape[2] = { min_n, n1 };
            size_t s_shape[1] = { min_n };
            size_t vt_shape[2] = { n2, min_n };
            ndfunc_arg_in_t ain[1] = {{cT,2}};
            ndfunc_arg_out_t aout[3] = {
                {cT,COUNT_OF_(u_shape),u_shape},
                {cRT,COUNT_OF_(s_shape),s_shape},
                {cT,COUNT_OF_(vt_shape),vt_shape}};
            ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};

            ans = na_ndloop3(&ndf, &opt, 1, a);
        }
        break;
    case SVD_VALS_ONLY:
        {
            size_t s_shape[1] = { min_(n1, n2) };
            ndfunc_arg_in_t ain[1] = {{cT,2}};
            ndfunc_arg_out_t aout[1] = {
                {cRT,COUNT_OF_(s_shape),s_shape}};
            ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, COUNT_OF_(ain), COUNT_OF_(aout), ain, aout};

            ans = na_ndloop3(&ndf, &opt, 1, a);
        }
        break;
    default:
        assert(0);
    }
    return ans;
}

static VALUE
<%=c_func%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    svd_job job=SVD_FULL;
    {
        VALUE h;
        if ( ! NIL_P(h = rb_check_hash_type(argv[argc-1]))) {
            ID tbl;
            VALUE v;
            tbl = rb_intern("job");
            rb_get_kwargs(h, &tbl, 0, 1, &v);
            if (v != Qundef) {
                ID job_id = rb_to_id(v);
                if (job_id == rb_intern("full")) {
                    job = SVD_FULL;
                } else if (job_id == rb_intern("thin")) {
                    job = SVD_THIN;
                } else if (job_id == rb_intern("vals_only")) {
                    job = SVD_VALS_ONLY;
                } else {
                    rb_raise(rb_eArgError, "not valid argument for job:");
                }
            }
            --argc;
        }
    }
    rb_check_arity(argc, 1, 1);
    return sub_func_name(<%=c_func%>, (argv[0], job));
}

#undef sub_func_name
