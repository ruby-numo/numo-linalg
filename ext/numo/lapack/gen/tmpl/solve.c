/*
* DGESV computes the solution to system of linear equations A * X = B
* for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*
* ZGESV computes the solution to system of linear equations A * X = B
* for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * )
*       ..
*/

#define gesv FFUNC(<%=blas_char%>gesv)

void gesv(
  fortran_integer * /*N*/, fortran_integer * /*NRHS*/,
  dtype * /*A*/, fortran_integer * /*LDA*/,
  fortran_integer * /*IPIV*/,
  dtype * /*B*/, fortran_integer * /*LDB*/,
  fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a, *bi, *bo;
    size_t n11;
    fortran_integer *ipiv;
    fortran_integer n, nrhs, info=0;
    volatile VALUE tmp_ipiv;

    a  = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    bi = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    bo = (dtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    n11 = lp->args[0].shape[0]; // n
    if (lp->args[1].ndim > 1) {
        nrhs = (fortran_integer)lp->args[1].shape[1];
    } else {
        // ndim == 1
        nrhs = 1;
    }
    memcpy(bo, bi, n11*nrhs*sizeof(dtype));
    n = n11;
    ipiv = rb_alloc_tmp_buffer(&tmp_ipiv, n11*sizeof(fortran_integer));

    gesv(&n, &nrhs, a, &n, ipiv, bo, &n, &info);

    rb_free_tmp_buffer(&tmp_ipiv);
    RB_GC_GUARD(tmp_ipiv);
}

/*
  @overload solve(a, b)
  @param [Numo::<%=class_name%>] a  >=2-dimentional NArray.
  @param [Numo::<%=class_name%>] b  >=1-dimentional NArray.
  @return [Numo::<%=class_name%>]
  @raise

  <%=blas_char%>gesv - computes the solution to a complex system of linear equations
     A  *  X = B,
  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
  The LU decomposition with partial pivoting and row interchanges is used to factor A as
     A = P * L * U,
  where P is a permutation matrix, L is unit lower triangular, and U is upper triangular.
  The factored form of A is then used to solve the system of equations A * X = B.
*/
static VALUE
<%=c_func%>(VALUE const UNUSED(mod), VALUE const a1, VALUE const a2)
{
    volatile VALUE ans;

    narray_t *na1, *na2;
    size_t n11, n12, n21, n22;
    size_t shape[2];
    ndfunc_arg_in_t ain[2] = {{cT, 2}, {cT, 2}};
    ndfunc_arg_out_t aout[1] = {{cT, 2, shape}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP, 2, 1, ain, aout};

    GetNArray(a1, na1);
    GetNArray(a2, na2);
    CHECK_DIM_GE(na1, 2);
    CHECK_DIM_GE(na2, 1);
    n11 = na1->shape[na1->ndim-2]; // n
    n12 = na1->shape[na1->ndim-1]; // n
    if (NA_NDIM(na2) == 1) {
        ain[1].dim = 1;
        aout[0].dim = 1;
        n21 = na2->shape[na2->ndim-1]; // n
        n22 = 1;                       // nrhs
    } else {
        n21 = na2->shape[na2->ndim-2]; // n
        n22 = na2->shape[na2->ndim-1]; // nrhs
    }
    if ((n11 != n12) || (n11 != n21)) {
        rb_raise(nary_eShapeError, "matrix dimension mismatch: "
                 "n11=%"SZF"u n12=%"SZF"u n21=%"SZF"u", n11, n12, n21);
    }
    shape[0] = n11; // n
    shape[1] = n22; // nrhs

    ans = na_ndloop3(&ndf, 0, 2, a1, a2);

    return ans;
}
