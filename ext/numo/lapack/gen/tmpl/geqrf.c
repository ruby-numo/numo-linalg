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
<% if is_complex %>
#define ungqr FFUNC(<%=blas_char%>ungqr)
<% else %>
#define orgqr FFUNC(<%=blas_char%>orgqr)
<% end %>

void geqrf(
    fortran_integer * /*M*/, fortran_integer * /*N*/,
    dtype * /*A*/, fortran_integer * /*LDA*/, dtype * /*TAU*/,
    dtype * /*WORK*/, fortran_integer * /*LWORK*/, fortran_integer * /*INFO*/);

<% if is_complex %>
void ungqr(
<% else %>
void orgqr(
<% end %>
    fortran_integer * /*M*/, fortran_integer * /*N*/, fortran_integer * /*K*/,
    dtype * /*A*/, fortran_integer * /*LDA*/, dtype * /*TAU*/,
    dtype * /*WORK*/, fortran_integer * /*LWORK*/, fortran_integer * /*INFO*/);

typedef struct {
    size_t lwork;
    size_t lwork2;
} geqrf_opt_t;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    geqrf_opt_t *opt;
    size_t n1, n2, n3;
    dtype *a, *q, *r, *tau, *p1, *p2;
    fortran_integer m, n, k, lda, lwork, lwork2, info=0;
    dtype *work, *work2;
    volatile VALUE tmp_ptr;
    <% if is_complex %>
    dtype zero;
    REAL(zero) = 0.0;
    IMAG(zero) = 0.0;
    <% else %>
    dtype const zero = 0.0;
    <% end %>

    a = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    q = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    r = (dtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    n1 = lp->args[0].shape[1];  // m
    n2 = lp->args[0].shape[0];  // n
    n3 = lp->args[1].shape[0];  // MIN(m, n) == k
    if (n3 == n2) {
        p1 = q; p2 = r;
    } else {
        assert(n3 == n1);
        p1 = r; p2 = q;
    }
    memcpy(p1, a, n1*n2*sizeof(dtype));

    opt = lp->opt_ptr;
    {
        char *ptr;
        size_t ofs[4];
        size_t wksize  = opt->lwork;
        size_t wksize2 = opt->lwork2;

        lwork  = (fortran_integer)wksize;
        lwork2 = (fortran_integer)wksize2;

        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, n3);       // tau[MIN(m, n)]
        SET_POS(ofs, 2, dtype, wksize);   // work[lwork]
        SET_POS(ofs, 3, dtype, wksize2);  // work2[lwork2]

        ptr = rb_alloc_tmp_buffer(&tmp_ptr, ofs[3]);

        tau   = (dtype *) ptr;
        work  = (dtype *)(ptr + ofs[1]);
        work2 = (dtype *)(ptr + ofs[2]);
    }
    m = n1; n = n2;
    lda = m;

    geqrf(&m, &n, p1, &lda, tau, work, &lwork, &info);

    {
        size_t i, j;
        for (i = 0; i < n3; ++i) {
            for (j = 0; j <= i; ++j) {
                p2[n3*i + j] = p1[lda*i + j];
            }
            for (; j < n3; ++j) {
                p2[n3*i + j] = zero;
            }
        }
    }
    if (n3 == n2) {
        //NOP
    } else {
        assert(n3 == n1);
        n = m;
    }
    k = n3;

    <% if is_complex %>
    ungqr(
    <% else %>
    orgqr(
    <% end %>
        &m, &n, &k, p1, &lda, tau, work2, &lwork2, &info);

    if (n3 == n2) {
        //NOP
    } else {
        size_t i, j;
        assert(n3 == n1);
        for (i = 0; i < n3; ++i) {
            for (j = 0; j < n3; ++j) {
                dtype * const pp1 = &p1[lda*i + j];
                dtype * const pp2 = &p2[n3*i + j];
                dtype const tmp = *pp1;
                *pp1 = *pp2;
                *pp2 = tmp;
            }
        }
    }

    rb_free_tmp_buffer(&tmp_ptr);
    RB_GC_GUARD(tmp_ptr);
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
    geqrf_opt_t opt;

    narray_t *na;
    size_t n1, n2, n3;
    size_t shape1[2], shape2[2];
    ndfunc_arg_in_t ain[1] = {{cT, 2}};
    ndfunc_arg_out_t aout[2] = {
        {cT, COUNT_OF_(shape1), shape1},   // Q
        {cT, COUNT_OF_(shape2), shape2}};  // R
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
        COUNT_OF_(ain), COUNT_OF_(aout),
        ain, aout};

    GetNArray(a, na);
    CHECK_DIM_GE(na, 2);
    n1 = na->shape[na->ndim-1];  // m
    n2 = na->shape[na->ndim-2];  // n

    if (n1 < n2) {
        shape1[1] = n1;
        shape1[0] = n1;
        shape2[1] = n1;
        shape2[0] = n2;
        n3 = n1;
    } else if (n1 == n2) {
        n3 = n1;
        shape1[1] = n3;
        shape1[0] = n3;
        shape2[1] = n3;
        shape2[0] = n3;
    } else if (n1 > n2) {
        shape1[1] = n1;
        shape1[0] = n2;
        shape2[1] = n2;
        shape2[0] = n2;
        n3 = n2;
    } else {
        n3 = 0;
        assert(0);
    }

    {
        fortran_integer info=0;
        fortran_integer m, n, k, lda, lwork;
        dtype wk[1];
        m = n1;
        n = n2;
        k = n3;
        lda = m;
        lwork = -1;
        geqrf(&m, &n, 0, &lda, 0, wk, &lwork, &info);
        <% if is_complex %>
        opt.lwork = (size_t)REAL(wk[0]);
        <% else %>
        opt.lwork = (size_t)wk[0];
        <% end %>

        if (m > n) {
            m = k;
        } else if (m < n) {
            n = k;
        }
        <% if is_complex %>
        ungqr(
        <% else %>
        orgqr(
        <% end %>
            &m, &n, &k, 0, &lda, 0, wk, &lwork, &info);
        <% if is_complex %>
        opt.lwork2 = (size_t)REAL(wk[0]);
        <% else %>
        opt.lwork2 = (size_t)wk[0];
        <% end %>
    }
    ans = na_ndloop3(&ndf, &opt, 1, a);

    return ans;
}
