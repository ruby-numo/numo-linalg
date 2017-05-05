/*
* DSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*       ..
*
* ZHEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * ), W( * )
*       COMPLEX*16         A( LDA, * ), WORK( * )
*       ..
*/

<% if is_complex %>
#define heev FFUNC(<%=blas_char%>heev)
<% else %>
#define syev FFUNC(<%=blas_char%>syev)
<% end %>

typedef struct {
    int vals_only;
    int upper;
    int overwrite;
    size_t lwork;
} heev_opt_t;

<% if is_complex %>
void heev(
<% else %>
void syev(
<% end %>
  char const * /*JOBZ*/, char const * /*UPLO*/,
  fortran_integer * /*N*/, dtype * /*A*/, fortran_integer * /*LDA*/,
  rtype * /*W*/, dtype * /*WORK*/, fortran_integer * /*LWORK*/,
  <% if is_complex %>
  rtype * /*RWORK*/,
  <% end %>
  fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const jobz_n="N", * const jobz_v="V";
    char const * const uplo_u="U", * const uplo_l="L";
    char const *jobz, *uplo;
    heev_opt_t *opt;
    size_t n1, wksize;
    dtype *ai, *a;
    rtype *w;
    fortran_integer n, lda, lwork, info=0;
    dtype *work;
    <% if is_complex %>
    rtype *rwork;
    <% end %>
    volatile VALUE tmp_work;
    <% if is_complex %>
    volatile VALUE tmp_rwork;
    <% end %>

    opt = lp->opt_ptr;

    jobz = ( opt->vals_only ? jobz_n : jobz_v ) ;
    uplo = ( opt->upper     ? uplo_u : uplo_l ) ;

    n1 = lp->args[0].shape[1];
    if (opt->vals_only || opt->overwrite) {
        a  = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
        w  = (rtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    } else {// !opt->vals_only && !opt->overwrite
        ai = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
        w  = (rtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
        a  = (dtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
        memcpy(a, ai, n1*n1*sizeof(dtype));
    }
    wksize = opt->lwork;
    lwork = (fortran_integer)wksize;

    n = n1;
    lda = n;

    work = rb_alloc_tmp_buffer(&tmp_work, wksize*sizeof(dtype));
    <% if is_complex %>
    rwork = rb_alloc_tmp_buffer(&tmp_rwork, (3*n1-2)*sizeof(rtype));
    <% end %>

    <% if is_complex %>
    heev(jobz, uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
    <% else %>
    syev(jobz, uplo, &n, a, &lda, w, work, &lwork, &info);
    <% end %>

    rb_free_tmp_buffer(&tmp_work);
    <% if is_complex %>
    rb_free_tmp_buffer(&tmp_rwork);
    <% end %>
    RB_GC_GUARD(tmp_work);
    <% if is_complex %>
    RB_GC_GUARD(tmp_rwork);
    <% end %>
}

#define sub_func_name(f, args) f##_sub args

/*
  @overload heev(a)
  @param [Numo::<%=class_name%>] a  >=2-dimentional Hamiltonian (symmetric) NArray.
  @return [Numo::<%=class_name%>] TBD
  @raise

  TBD
*/
static VALUE
sub_func_name(<%=c_func%>, (VALUE const a, int const vals_only, int const upper, int const overwrite))
{
    volatile VALUE ans;
    narray_t *na;
    size_t n1, n2;
    heev_opt_t opt;
    fortran_integer lwork;

    GetNArray(a, na);
    CHECK_DIM_GE(na, 2);
    n1 = na->shape[na->ndim-1];
    n2 = na->shape[na->ndim-2];

    if (n1 != n2) {
        rb_raise(nary_eShapeError, "not square-matrix");
    }

    opt.vals_only = vals_only;
    opt.upper = upper;
    opt.overwrite = overwrite;
    {
        char const * const jobz_n="N", * const jobz_v="V";
        char const * const uplo_u="U", * const uplo_l="L";
        char const *jobz, *uplo;
        fortran_integer info=0;
        fortran_integer n, lda;
        dtype wk[1];

        jobz = ( vals_only ? jobz_n : jobz_v ) ;
        uplo = ( upper     ? uplo_u : uplo_l ) ;

        n = n1;
        lda = n;
        lwork = -1;
        <% if is_complex %>
        heev(jobz, uplo, &n, 0, &lda, 0, wk, &lwork, 0, &info);
        lwork = (fortran_integer)REAL(wk[0]);
        <% else %>
        syev(jobz, uplo, &n, 0, &lda, 0, wk, &lwork, &info);
        lwork = (fortran_integer)wk[0];
        <% end %>
    }
    opt.lwork = (size_t)lwork;

    if (vals_only || overwrite) {
        size_t shape[1];
        ndfunc_arg_in_t ain[1] = {{cT, 2}};
        ndfunc_arg_out_t aout[1] = {
          {cRT, COUNT_OF_(shape), shape}};
        ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
          COUNT_OF_(ain), COUNT_OF_(aout),
          ain, aout};
        shape[0] = n1;
        ans = na_ndloop3(&ndf, &opt, 1, a);
    } else {// !vals_only && !overwrite
        size_t shape1[1], shape2[2];
        ndfunc_arg_in_t ain[1] = {{cT, 2}};
        ndfunc_arg_out_t aout[2] = {
          {cRT, COUNT_OF_(shape1), shape1},
          {cT, COUNT_OF_(shape2), shape2}};
        ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
          COUNT_OF_(ain), COUNT_OF_(aout),
          ain, aout};
        shape1[0] = n1;
        shape2[1] = n1;
        shape2[0] = n1;
        ans = na_ndloop3(&ndf, &opt, 1, a);
    }
    return ans;
}

static VALUE
<%=c_func%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    int const flg_overwrite = 0;
    int flg_vals_only=0, flg_upper=0;
    {
        VALUE const h = rb_check_hash_type(argv[argc-1]);
        if ( ! NIL_P(h)) {
            --argc;
        }
        rb_check_arity(argc, 1, 1);
        if ( ! NIL_P(h)) {
            ID table[2];
            VALUE values[COUNT_OF_(table)];
            table[0] = rb_intern("vals_only");
            table[1] = rb_intern("upper");
            rb_get_kwargs(h, table, 0, COUNT_OF_(table), values);
            if (values[0] != Qundef) {
                flg_vals_only = RTEST(values[0]);
            }
            if (values[1] != Qundef) {
                flg_upper = RTEST(values[1]);
            }
        }
    }
    return sub_func_name(<%=c_func%>, (argv[0], flg_vals_only, flg_upper, flg_overwrite));
}

#undef sub_func_name
