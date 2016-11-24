/*
* DGEQRF computes a QR factorization of a real matrix
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
* ZGEQRF computes a QR factorization of a complex matrix
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*       ..
*/

#define geqrf FFUNC(<%=blas_char%>geqrf)

void geqrf(
  fortran_integer * /*M*/, fortran_integer * /*N*/,
  dtype * /*A*/, fortran_integer * /*LDA*/, dtype * /*TAU*/,
  dtype * /*WORK*/, fortran_integer * /*LWORK*/, fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    size_t n1, n2, wksize;
    dtype *ai, *ao, *tau;
    fortran_integer m, n, lda, lwork, info=0;
    dtype *work;
    volatile VALUE tmp_work;

    ai  = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    ao  = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    tau = (dtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    n1 = lp->args[0].shape[1];  // m
    n2 = lp->args[0].shape[0];  // n
    lwork = *(fortran_integer *)lp->opt_ptr;
    wksize = lwork;
    memcpy(ao, ai, n1*n2*sizeof(dtype));

    m = n1; n = n2;
    lda = m;
    work = rb_alloc_tmp_buffer(&tmp_work, wksize*sizeof(dtype));

    geqrf(&m, &n, ao, &lda, tau, work, &lwork, &info);

    rb_free_tmp_buffer(&tmp_work);
    RB_GC_GUARD(tmp_work);
}

/*
  @overload geqrf(a, b)
  @param [Numo::<%=class_name%>] a  >=2-dimentional NArray.
  @return [Numo::<%=class_name%>] TBD
  @raise

  TBD
*/
static VALUE
<%=c_func%>(VALUE const UNUSED(mod), VALUE const a)
{
    volatile VALUE ans;

    narray_t *na;
    size_t n1, n2;
    size_t shape1[2], shape2[1];
    ndfunc_arg_in_t ain[1] = {{cT, 2}};
    ndfunc_arg_out_t aout[2] = {
      {cT, COUNT_OF_(shape1), shape1},
      {cT, COUNT_OF_(shape2), shape2}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
      COUNT_OF_(ain), COUNT_OF_(aout),
      ain, aout};
    fortran_integer lwork;

    GetNArray(a, na);
    CHECK_DIM_GE(na, 2);
    n1 = na->shape[na->ndim-1];  // m
    n2 = na->shape[na->ndim-2];  // n

    shape1[1] = n1;  // m
    shape1[0] = n2;  // n
    shape2[0] = min_(n1, n2);  // min(m, n)
    {
        fortran_integer info=0;
        fortran_integer m, n, lda;
        dtype wk[1];
        m = n1;
        n = n2;
        lda = m;
        lwork = -1;
        geqrf(&m, &n, 0, &lda, 0, wk, &lwork, &info);
        <% if is_complex %>
        lwork = (fortran_integer)REAL(wk[0]);
        <% else %>
        lwork = (fortran_integer)wk[0];
        <% end %>
    }
    ans = na_ndloop3(&ndf, &lwork, 1, a);

    return ans;
}
