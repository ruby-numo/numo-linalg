/*
* DGETRF computes an LU factorization
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
* ZGETRF computes an LU factorization
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*/

#define getrf FFUNC(<%=blas_char%>getrf)

void getrf(
  fortran_integer * /*M*/, fortran_integer * /*N*/,
  dtype * /*A*/, fortran_integer * /*LDA*/,
  fortran_integer * /*IPIV*/, fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    size_t n1, n2, n3;
    dtype *ai, *ao;
    fortran_integer *ipiv;

    ai   =           (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    ao   =           (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    ipiv = (fortran_integer *)(lp->args[2].ptr + lp->args[2].iter[0].pos);

    n1 = lp->args[0].shape[1];
    n2 = lp->args[0].shape[0];
    n3 = lp->args[2].shape[0];

    memcpy(ao, ai, n1*n2*sizeof(dtype));

    {
        fortran_integer m = (fortran_integer)n1, n = (fortran_integer)n2;
        fortran_integer lda = m, info=0;
        getrf(&m, &n, ao, &lda, ipiv, &info);
        {
            size_t i;
            for (i = 0; i < n3; ++i) {
                --ipiv[i];
            }
        }
    }
}

/*
  @overload getrf(a)
  @param [Numo::<%=class_name%>] a  >=2-dimentional NArray.
  @return [Numo::<%=class_name%>] TBD
  @raise

  TBD
*/
static VALUE
<%=c_func%>(VALUE const UNUSED(mod), VALUE const a)
{
    volatile VALUE ans;
    size_t n1, n2;
    {
        narray_t *na;
        GetNArray(a, na);
        CHECK_DIM_GE(na, 2);
        n1 = na->shape[na->ndim-1];
        n2 = na->shape[na->ndim-2];
    }
    {
        size_t shape1[2] = { n2, n1 };
        size_t shape2[1] = { min_(n1, n2) };
        ndfunc_arg_in_t ain[1] = {{cT, 2}};
        ndfunc_arg_out_t aout[2] = {
            {cT, COUNT_OF_(shape1), shape1},
            {cIT, COUNT_OF_(shape2), shape2}};
        ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
            COUNT_OF_(ain), COUNT_OF_(aout),
            ain, aout};
        ans = na_ndloop3(&ndf, 0, 1, a);
    }
    return ans;
}
