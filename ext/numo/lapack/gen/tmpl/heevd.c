/*
* DSYEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
*                          LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LIWORK, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*       ..
*
* ZHEEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
*                          LRWORK, IWORK, LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   RWORK( * ), W( * )
*       COMPLEX*16         A( LDA, * ), WORK( * )
*       ..
*/

<% if is_complex %>
#define heevd FFUNC(<%=blas_char%>heevd)
<% else %>
#define syevd FFUNC(<%=blas_char%>syevd)
<% end %>

typedef struct {
    int vals_only;
    int upper;
    int overwrite;
    size_t lwork;
    <% if is_complex %>
    size_t lrwork;
    <% end %>
    size_t liwork;
} heevd_opt_t;

<% if is_complex %>
void heevd(
<% else %>
void syevd(
<% end %>
  char const * /*JOBZ*/, char const * /*UPLO*/,
  fortran_integer * /*N*/, dtype * /*A*/, fortran_integer * /*LDA*/,
  rtype * /*W*/, dtype * /*WORK*/, fortran_integer * /*LWORK*/,
  <% if is_complex %>
  rtype * /*RWORK*/, fortran_integer * /*LRWORK*/,
  <% end %>
  fortran_integer * /*IWORK*/, fortran_integer * /*LIWORK*/,
  fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const jobz_n="N", * const jobz_v="V";
    char const * const uplo_u="U", * const uplo_l="L";
    char const *jobz, *uplo;
    heevd_opt_t *opt;
    size_t n1;
    dtype *ai, *a;
    rtype *w;
    fortran_integer n, lda, info=0;
    fortran_integer lwork;
    <% if is_complex %>
    fortran_integer lrwork;
    <% end %>
    fortran_integer liwork;
    dtype *work;
    <% if is_complex %>
    rtype *rwork;
    <% end %>
    fortran_integer *iwork;
    volatile VALUE tmp_ptr;

    opt = lp->opt_ptr;

    jobz = ( opt->vals_only ? jobz_n : jobz_v ) ;
    uplo = ( opt->upper    ? uplo_u : uplo_l ) ;

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
    n = n1;
    lda = n;
    {
        char *ptr;
        <% if is_complex %>
        size_t ofs[4];
        <% else %>
        size_t ofs[3];
        <% end %>
        size_t wksize  = opt->lwork;
        <% if is_complex %>
        size_t rwksize = opt->lrwork;
        <% end %>
        size_t iwksize = opt->liwork;

        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, wksize);  // work[lwork]
        <% if is_complex %>
        SET_POS(ofs, 2, rtype, rwksize);  // rwork[lrwork]
        SET_POS(ofs, 3, fortran_integer, iwksize);  // iwork[liwork]
        ptr = rb_alloc_tmp_buffer(&tmp_ptr, ofs[3]);
        <% else %>
        SET_POS(ofs, 2, fortran_integer, iwksize);  // rwork[liwork]
        ptr = rb_alloc_tmp_buffer(&tmp_ptr, ofs[2]);
        <% end %>

        work  =           (dtype *) ptr;
        <% if is_complex %>
        rwork =           (rtype *)(ptr + ofs[1]);
        iwork = (fortran_integer *)(ptr + ofs[2]);
        <% else %>
        iwork = (fortran_integer *)(ptr + ofs[1]);
        <% end %>

        lwork = (fortran_integer)wksize;
        <% if is_complex %>
        lrwork = (fortran_integer)rwksize;
        <% end %>
        liwork = (fortran_integer)iwksize;
    }
    <% if is_complex %>
    heevd(
    <% else %>
    syevd(
    <% end %>
      jobz, uplo, &n, a, &lda, w, work, &lwork,
      <% if is_complex %>
      rwork, &lrwork,
      <% end %>
      iwork, &liwork, &info);

    rb_free_tmp_buffer(&tmp_ptr);
    RB_GC_GUARD(tmp_ptr);
}

#define sub_func_name(f, args) f##_sub args

/*
  @overload heevd(a)
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
    heevd_opt_t opt;
    fortran_integer lwork;
    <% if is_complex %>
    fortran_integer lrwork;
    <% end %>
    fortran_integer liwork;

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
        <% if is_complex %>
        rtype rwk[1];
        <% end %>
        fortran_integer iwk[1];

        jobz = ( vals_only ? jobz_n : jobz_v ) ;
        uplo = ( upper    ? uplo_u : uplo_l ) ;

        n = n1;
        lda = n;
        lwork = -1;
        liwork = -1;
        <% if is_complex %>
        lrwork = -1;
        heevd(jobz, uplo, &n, 0, &lda, 0, wk, &lwork, rwk, &lrwork, iwk, &liwork, &info);
        lwork = (fortran_integer)REAL(wk[0]);
        lrwork = (fortran_integer)rwk[0];
        <% else %>
        syevd(jobz, uplo, &n, 0, &lda, 0, wk, &lwork, iwk, &liwork, &info);
        lwork = (fortran_integer)wk[0];
        <% end %>
        liwork = iwk[0];
    }
    opt.lwork = (size_t)lwork;
    <% if is_complex %>
    opt.lrwork = (size_t)lrwork;
    <% end %>
    opt.liwork = (size_t)liwork;

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
