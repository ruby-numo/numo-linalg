/*
* DGEEV computes the eigenvalues and, optionally, the left and/or
* right eigenvectors for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
*                         LDVR, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WI( * ), WORK( * ), WR( * )
*       ..
*
* ZGEEV computes the eigenvalues and, optionally, the left and/or
* right eigenvectors for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
*                         WORK, LWORK, RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   W( * ), WORK( * )
*       ..
*/

#define geev FFUNC(<%=blas_char%>geev)

typedef struct {
    int vals_only;
    size_t lwork;
} geev_opt_t;

void geev(
  char const * /*JOBVL*/, char const * /*JOBVR*/,
  fortran_integer * /*N*/, dtype * /*A*/, fortran_integer * /*LDA*/,
<% if is_complex %>
  dtype * /*W*/,
<% else %>
  dtype * /*WR*/, dtype * /*WI*/,
<% end %>
  dtype * /*VL*/,   fortran_integer * /*LDVL*/,
  dtype * /*VR*/,   fortran_integer * /*LDVR*/,
  dtype * /*WORK*/, fortran_integer * /*LWORK*/,
<% if is_complex %>
  rtype * /*RWORK*/,
<% end %>
  fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const jobvl="N";
    char const * const jobvr_n="N"; char const * const jobvr_v="V";
    char const *jobvr;
    size_t n0;
    dtype *a;
    ctype *w, *vr;
    fortran_integer n, lda, ldvl=1, ldvr, info=0;
    dtype *work;
    <% if is_complex %>
    rtype *rwork;
    <% else %>
    dtype *wr, *wi, *vrr;
    <% end %>
    geev_opt_t const * const opt = lp->opt_ptr;

    // a[n, lda], w[n], vr[n, ldvr]
    a  = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    w  = (ctype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    n0 = lp->args[1].shape[0];
    n  = n0;
    lda  = n;
    if (opt->vals_only) {
        jobvr = jobvr_n;
        vr = 0;
        ldvr = 1;
    } else {
        jobvr = jobvr_v;
        vr = (ctype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
        ldvr = n;
    }
    {
        volatile VALUE tmp_wk;
        char *ptr;
        size_t wksize = opt->lwork;
        fortran_integer lwork = (fortran_integer)wksize;

        <% if is_complex %>

        size_t ofs[3];
        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, wksize);
        SET_POS(ofs, 2, rtype, 2*n);
        ptr = rb_alloc_tmp_buffer(&tmp_wk, ofs[2]);
        work  = (dtype *) ptr;
        rwork = (rtype *)(ptr + ofs[1]);

        geev(jobvl, jobvr, &n, a, &lda, w, 0, &ldvl, vr, &ldvr,
             work, &lwork, rwork, &info);

        <% else %>

        if (opt->vals_only) {
            size_t ofs[4];
            ofs[0] = 0;
            SET_POS(ofs, 1, dtype, wksize);
            SET_POS(ofs, 2, dtype, n);
            SET_POS(ofs, 3, dtype, n);
            ptr = rb_alloc_tmp_buffer(&tmp_wk, ofs[3]);
            work = (dtype *) ptr;
            wr   = (dtype *)(ptr + ofs[1]);
            wi   = (dtype *)(ptr + ofs[2]);
            vrr  = 0;
        } else {
            size_t ofs[5];
            ofs[0] = 0;
            SET_POS(ofs, 1, dtype, wksize);
            SET_POS(ofs, 2, dtype, n);
            SET_POS(ofs, 3, dtype, n);
            SET_POS(ofs, 4, dtype, n*n);
            ptr = rb_alloc_tmp_buffer(&tmp_wk, ofs[4]);
            work = (dtype *) ptr;
            wr   = (dtype *)(ptr + ofs[1]);
            wi   = (dtype *)(ptr + ofs[2]);
            vrr  = (dtype *)(ptr + ofs[3]);
        }
        geev(jobvl, jobvr, &n, a, &lda, wr, wi, 0, &ldvl, vrr, &ldvr,
             work, &lwork, &info);
        {
            size_t j;
            for (j=0; j<n0; ++j) {
                REAL(w[j]) = wr[j];
                IMAG(w[j]) = wi[j];
            }
            if ( ! opt->vals_only) {
                for (j=0; j<n0; ++j) {
                    dtype * const p1 = &(vrr[j*ldvr]);
                    ctype * const p2 = &(vr[j*ldvr]);
                    if (IMAG(w[j]) == 0.0) { /* real eigenvalue */
                        size_t i;
                        for (i=0; i<n0; ++i) {
                            REAL(p2[i]) =  p1[i];
                            IMAG(p2[i]) =  0;
                        }
                    } else { /* complex eigenvalue */
                        size_t i;
                        /* v(j)   = VR(:,j) + i*VR(:,j+1) and
                           v(j+1) = VR(:,j) - i*VR(:,j+1) */
                        for (i=0; i<n0; ++i) {
                            REAL(p2[i])      =  p1[i];
                            IMAG(p2[i])      =  p1[i+ldvr];
                            REAL(p2[i+ldvr]) =  p1[i];
                            IMAG(p2[i+ldvr]) = -p1[i+ldvr];
                        }
                        ++j;
                    }
                }
            }
        }

        <% end %>

        rb_free_tmp_buffer(&tmp_wk);
        RB_GC_GUARD(tmp_wk);
    }
}

#define sub_func_name(f, args) f##_sub args

/*
  @overload geev(a)
  @param ***TODO*** [Numo::<%=class_name%>] a >=2-dimentional NArray.
  @return ***TODO*** [[Numo::<%=complex_class_name%>,Numo::<%=complex_class_name%>]] pair of eigenvalue and right eigenvector
  @raise

  ***TODO***
  <%=blas_char%>geev - computes the eigenvalues and the right eigenvectors
  for an N-by-N real nonsymmetric matrix A.
  The right eigenvector v(j) of A satisfies
                        A * v(j) = lambda(j) * v(j)
  where lambda(j) is its eigenvalue.
  The computed eigenvectors are normalized to have
  Euclidean norm equal to 1 and largest component real.
*/
static VALUE
sub_func_name(<%=c_func%>, (VALUE const a, int const vals_only))
{
    volatile VALUE ans;
    geev_opt_t opt;

    narray_t *na;
    size_t n;

    ndfunc_arg_in_t ain[1] = {{cT, 2}};

    GetNArray(a, na);
    CHECK_DIM_GE(na, 2);
    n = na->shape[na->ndim-1];
    if (n != na->shape[na->ndim-2]) {
        rb_raise(nary_eShapeError, "not square-matrix");
    }
    {
        char const * const chr_n="N", * const chr_v="V";
        char const * chr;
        fortran_integer m=(fortran_integer)n, one=1, ldvr, lwork, info=0;
        dtype wk[1];
        if (vals_only) {
            chr = chr_n;
            ldvr = 1;
        } else {
            chr = chr_v;
            ldvr = m;
        }
        lwork = -1;
        <% if is_complex %>
        geev(chr_n, chr, &m, 0, &m, 0, 0, &one,  0, &ldvr, wk, &lwork, 0, &info);
        lwork = (fortran_integer)REAL(wk[0]);
        <% else %>
        geev(chr_n, chr, &m, 0, &m, 0, 0, 0, &one, 0, &ldvr, wk, &lwork, &info);
        lwork = (fortran_integer)wk[0];
        <% end %>
        opt.lwork = (size_t)lwork;
    }
    opt.vals_only = vals_only;
    if (vals_only) {
        size_t eval_shape[1] = { n };
        ndfunc_arg_out_t aout[1] = {
            {cCT, COUNT_OF_(eval_shape), eval_shape}};
        ndfunc_t ndf = {<%=c_iter%>, NO_LOOP,
            COUNT_OF_(ain),
            COUNT_OF_(aout),
            ain, aout};
        ans = na_ndloop3(&ndf, &opt, 1, a);
    } else {
        size_t eval_shape[1] = { n };
        size_t evec_shape[2] = { n, n };
        ndfunc_arg_out_t aout[2] = {
            {cCT, COUNT_OF_(eval_shape), eval_shape},
            {cCT, COUNT_OF_(evec_shape), evec_shape}};
        ndfunc_t ndf = {<%=c_iter%>, NO_LOOP,
            COUNT_OF_(ain),
            COUNT_OF_(aout),
            ain, aout};
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
