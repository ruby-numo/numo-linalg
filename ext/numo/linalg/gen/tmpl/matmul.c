/*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*/

#define gemm FFUNC(<%=blas_char%>gemm)

fortran_integer
gemm(const char*, const char*,
       fortran_integer*, fortran_integer*, fortran_integer*,
       dtype*, dtype*, fortran_integer*,
       dtype*, fortran_integer*,
       dtype*, dtype*, fortran_integer*);

static void
<%=c_iter%>(na_loop_t *const lp)
{
    static const char *trans = "t";
    <% if is_complex %>
    dtype alpha = c_one();
    dtype beta = c_zero();
    <% else %>
    dtype alpha = 1;
    dtype beta = 0;
    <% end %>
    dtype *a, *b, *c;
    size_t n11, n22, n21;
    fortran_integer m, n, k;

    // a[k,m], b[n,k], c[n,m]
    a = (dtype*)(lp->args[0].ptr + lp->iter[0].pos);
    b = (dtype*)(lp->args[1].ptr + lp->iter[1].pos);
    //c = (dtype*)(lp->args[2].ptr + lp->iter[2].pos);
    n11 = lp->args[0].shape[0]; // m
    n22 = lp->args[1].shape[1]; // n
    n21 = lp->args[1].shape[0]; // k

    //lda = lp->iter[0].step / sizeof(dtype);
    //ldb = lp->iter[1].step / sizeof(dtype);
    //ldc = lp->iter[2].step / sizeof(dtype);
    //printf("lda=%ld, ldb=%ld, ldc=%ld\n", lda, ldb, ldc);
    //cblas_<%=blas_char%>gemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
    //          m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
    c = (dtype*)(lp->opt_ptr);
    m = n11;
    n = n22;
    k = n21;
    //printf("m=%d, n=%d, k=%d\n", m, n, k);
    gemm(trans,trans,&m,&n,&k,&alpha,a,&k,b,&n,&beta,c,&m);
    transpose_get_2d(c,lp,/*i=*/2,/*n=*/3);
}

// a[m,lda>=k]
// b[k,ldb>=n]
// c[n,ldc>=m].trans

/*
  @overload matmul(a, b)
  @param [Numo::<%=class_name%>] a  >=2-dimentional NArray.
  @param [Numo::<%=class_name%>] b  >=2-dimentional NArray.
  @return [Numo::<%=class_name%>]
  @raise

  <%=blas_char%>gemm - performs matrix-matrix multiplication : C = A B
  where A, B and C are matrices,
  with A an m by k matrix, B a k by n matrix and
  C an m by n matrix.
*/
static VALUE
<%=c_func%>(VALUE mod, VALUE a1, VALUE a2)
{
    volatile VALUE ans, tmp_c;
    narray_t *na1, *na2;
    size_t    n11, n12, n21, n22;
    size_t    shape[2];
    void     *ptr;
    ndfunc_arg_in_t ain[2] = {{cT,2},{cT,2}};
    ndfunc_arg_out_t aout[1] = {{cT,2,shape}};
    ndfunc_t  ndf = {<%=c_iter%>, NO_LOOP, 2, 1, ain, aout};

    GetNArray(a1,na1);
    GetNArray(a2,na2);
    CHECK_DIM_GE(na1,2);
    CHECK_DIM_GE(na2,2);
    n11 = na1->shape[na1->ndim-2]; // m
    n12 = na1->shape[na1->ndim-1]; // k
    n21 = na2->shape[na2->ndim-2]; // k
    n22 = na2->shape[na2->ndim-1]; // n
    //printf("n11=%ld\n",n11);
    //printf("n12=%ld\n",n12);
    //printf("n21=%ld\n",n21);
    //printf("n22=%ld\n",n22);
    if (n12!=n21) {
        rb_raise(nary_eShapeError,"matrix dimension mismatch: "
                 "n11=%"SZF"u n12=%"SZF"u n21=%"SZF"u",n11,n12,n21);
    }
    shape[0] = n11; // m
    shape[1] = n22; // n
    ptr = rb_alloc_tmp_buffer(&tmp_c, n11*n22*sizeof(dtype));
    ans = na_ndloop3(&ndf, ptr, 2, a1, a2);
    rb_free_tmp_buffer(&tmp_c);
    return ans;
}
