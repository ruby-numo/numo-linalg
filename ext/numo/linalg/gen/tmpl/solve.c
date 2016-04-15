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

int gesv( fortran_integer *n, fortran_integer *nrhs,
          dtype *a, fortran_integer *lda,
          fortran_integer *ipiv,
          dtype *b, fortran_integer *ldb,
          fortran_integer *info );

static void
<%=c_iter%>(na_loop_t *lp, int rhdim)
{
    void  **ptr;
    dtype  *a, *b;
    fortran_integer *ipiv;
    fortran_integer n, nrhs, info=0;

    ptr  = (void**)(lp->opt_ptr);
    a    = (dtype*)(ptr[0]);
    b    = (dtype*)(ptr[1]);
    ipiv = (fortran_integer*)(ptr[2]);

    transpose_set_2d(a,lp,0,3);
    if (rhdim==1) {
        transpose_set_1d(b,lp,1,3);
        nrhs = 1;
    } else {
        transpose_set_2d(b,lp,1,3);
        nrhs = lp->args[1].shape[1]; //n[1+4*1] = nrhs
    }
    n = lp->args[1].shape[0];        //n[1+4*0] = n

    // call GESV
    gesv(&n, &nrhs, a, &n, ipiv, b, &n, &info);
    //printf("info=%d\n", info);

    // copy result X
    if (rhdim==1) {
        transpose_get_1d(b,lp,2,3);
    } else {
        transpose_get_2d(b,lp,2,3);
    }
}

static void
<%=c_iter%>_1d(na_loop_t *lp)
{
    <%=c_iter%>(lp,1);
}

static void
<%=c_iter%>_2d(na_loop_t *lp)
{
    <%=c_iter%>(lp,2);
}

/*
  @overload solve(a,b)
  @param [Numo::<%=class_name%>] a  >=2-dimentional NArray.
  @param [Numo::<%=class_name%>] b  >=2-dimentional NArray.
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
<%=c_func%>(VALUE mod, VALUE a1, VALUE a2)
{
    volatile VALUE ans, tmp_a, tmp_b, tmp_ipiv;
    narray_t *na1, *na2;
    size_t    n11, n12, n21, n22;
    size_t    shape1[2];
    void     *ptr[3];
    ndfunc_arg_in_t ain[2] = {{cT,2},{cT,2}};
    ndfunc_arg_out_t aout[1] = {{cT,2,shape1}};
    ndfunc_t  ndf = {<%=c_iter%>_2d, NO_LOOP, 2, 1, ain, aout};

    if (sizeof(fortran_integer)==8) {
        aout[1].type = numo_cInt64;
    }
    GetNArray(a1,na1);
    CHECK_DIM_GE(na1,2);
    GetNArray(a2,na2);
    CHECK_DIM_GE(na2,1);
    n11 = na1->shape[na1->ndim-2]; // n
    n12 = na1->shape[na1->ndim-1]; // n
    if (NA_NDIM(na2)==1) {
        ain[1].dim = 1;
        aout[0].dim = 1;
        ndf.func = <%=c_iter%>_1d;
        n21 = na2->shape[na2->ndim-1]; // n
        n22 = 1;                       // nrhs
    } else {
        n21 = na2->shape[na2->ndim-2]; // n
        n22 = na2->shape[na2->ndim-1]; // nrhs
    }
    if (n11!=n12 || n11!=n21) {
        rb_raise(nary_eShapeError,"matrix dimension mismatch: "
                 "n11=%"SZF"u n12=%"SZF"u n21=%"SZF"u",n11,n12,n21);
    }
    shape1[0] = n21; // n
    shape1[1] = n22; // nrhs
    printf("n11=%ld,n12=%ld,n21=%ld,n22=%ld\n",n11,n12,n21,n22);
    //printf("n11=%ld\n",n11);
    //printf("n12=%ld\n",n12);
    //printf("n21=%ld\n",n21);
    //printf("n22=%ld\n",n22);
    // Work memory
    ptr[0] = rb_alloc_tmp_buffer(&tmp_a, n11*n12*sizeof(dtype));
    ptr[1] = rb_alloc_tmp_buffer(&tmp_b, n21*n22*sizeof(dtype));
    ptr[2] = rb_alloc_tmp_buffer(&tmp_ipiv, n11*sizeof(fortran_integer));
    ans = na_ndloop3(&ndf, ptr, 2, a1, a2);
    rb_free_tmp_buffer(&tmp_ipiv);
    rb_free_tmp_buffer(&tmp_b);
    rb_free_tmp_buffer(&tmp_a);
    return ans;
}
