/*
* DGELS solves overdetermined or underdetermined systems for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
* ZGELS solves overdetermined or underdetermined systems for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*/

#define gels FFUNC(<%=blas_char%>gels)

void gels(
  char const * /*TRANS*/, fortran_integer * /*M*/,
  fortran_integer * /*N*/, fortran_integer * /*NRHS*/,
  dtype * /*A*/, fortran_integer * /*LDA*/,
  dtype * /*B*/, fortran_integer * /*LDB*/,
  dtype * /*WORK*/, fortran_integer * /*LWORK*/,
  fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const trans = "N";
    size_t n11, n12, n2, n3, wksize;
    dtype *a, *bi, *bo;
    fortran_integer m, n, nrhs, ldb, lwork, info=0;
    dtype *work;
    volatile VALUE tmp_work;

    a  = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    bi = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    bo = (dtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    n11 = lp->args[0].shape[1];  // m
    n12 = lp->args[0].shape[0];  // n
    n2  = lp->args[1].shape[0];  // m (of b)
    n3  = lp->args[2].shape[0];  // ldb
    lwork = *(fortran_integer *)lp->opt_ptr;
    wksize = lwork;
    memcpy(bo, bi, n2*sizeof(dtype));

    m = n11; n = n12;
    nrhs = 1;
    ldb = n3;
    work = rb_alloc_tmp_buffer(&tmp_work, wksize*sizeof(dtype));

    gels(trans, &m, &n, &nrhs, a, &m, bo, &ldb, work, &lwork, &info);

    rb_free_tmp_buffer(&tmp_work);
    RB_GC_GUARD(tmp_work);
}

/*
  @overload gels(a, b)
  @param [Numo::<%=class_name%>] a  >=2-dimentional NArray.
  @param [Numo::<%=class_name%>] b  >=1-dimentional NArray.
  @return [Numo::<%=class_name%>] TBD
  @raise

  TBD
*/
static VALUE
<%=c_func%>(VALUE const UNUSED(mod), VALUE const a1, VALUE const a2)
{
    char const * const chr = "N";
    volatile VALUE ans;

    narray_t *na1, *na2;
    size_t n11, n12, n2, n3;
    size_t shape[1];
    ndfunc_arg_in_t ain[2] = {{cT, 2}, {cT, 1}};
    ndfunc_arg_out_t aout[1] = {{cT, 1, shape}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP, 2, 1, ain, aout};
    fortran_integer lwork;

    GetNArray(a1, na1);
    GetNArray(a2, na2);
    CHECK_DIM_GE(na1, 2);
    CHECK_DIM_GE(na2, 1);
    n11 = na1->shape[na1->ndim-1]; // m
    n12 = na1->shape[na1->ndim-2]; // n
    n2  = na2->shape[na2->ndim-1];
    n3  = max_(n11, n12);

    if (n2 != n11) {
        rb_raise(nary_eShapeError, "matrix dimension mismatch: "
                                   "TBD");
    }

    shape[0] = n3;
    {
        fortran_integer nrhs, info=0;
        fortran_integer m, n, ldb;
        dtype wk[1];
        m = n11;
        n = n12;
        nrhs = 1;
        ldb = n3;
        lwork = -1;
        gels(chr, &m, &n, &nrhs, 0, &m, 0, &ldb, wk, &lwork, &info);
        <% if is_complex %>
        lwork = (fortran_integer)REAL(wk[0]);
        <% else %>
        lwork = (fortran_integer)wk[0];
        <% end %>
    }
    ans = na_ndloop3(&ndf, &lwork, 2, a1, a2);

    return ans;
}
