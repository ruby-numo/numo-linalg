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

void gemm(
  char const * /*TRANSA*/, char const * /*TRANSB*/,
  fortran_integer * /*M*/, fortran_integer * /*N*/, fortran_integer * /*K*/,
  dtype * /*ALPHA*/, dtype * /*A*/, fortran_integer * /*LDA*/,
  dtype * /*B*/, fortran_integer * /*LDB*/,
  dtype * /*BETA*/, dtype * /*C*/, fortran_integer * /*LDC*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const trans="N";
    <% if is_complex %>
    dtype alpha = c_one();
    dtype beta = c_zero();
    <% else %>
    dtype alpha = 1;
    dtype beta = 0;
    <% end %>
    dtype *a, *b, *c;
    size_t n11, n12, n21;
    fortran_integer m, n, k;

    // a[k,m], b[n,k], c[n,m]
    a = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    b = (dtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    c = (dtype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    n11 = lp->args[0].shape[0]; // k
    n12 = lp->args[0].shape[1]; // m
    n21 = lp->args[1].shape[0]; // n
    //n22 = lp->args[1].shape[1]; // k
    m = n12;
    n = n21;
    k = n11;

    gemm(trans, trans, &m, &n, &k, &alpha, a, &m, b, &k, &beta, c, &m);
}

/*
  @overload gemm(a, b)
  @param [Numo::<%=class_name%>] a  >=2-dimentional NArray.
  @param [Numo::<%=class_name%>] b  >=2-dimentional NArray.
  @return [Numo::<%=class_name%>]
  @raise

  ***TBD***<%=blas_char%>gemm - performs matrix-matrix multiplication : C = A B
  where A, B and C are matrices,
  with A an m by k matrix, B a k by n matrix and
  C an m by n matrix.
*/
static VALUE
<%=c_func%>(VALUE const UNUSED(mod), VALUE a1, VALUE a2)
{
    volatile VALUE ans;

    narray_t *na1, *na2, *na;
    size_t n11, n12, n21, n22;
    size_t m, n, k;
    size_t shape[2];
    int trim1=0, trim2=0;
    ndfunc_arg_in_t ain[2] = {{cT, 2}, {cT, 2}};
    ndfunc_arg_out_t aout[1] = {{cT, 2, shape}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP, 2, 1, ain, aout};

    GetNArray(a1, na1);
    GetNArray(a2, na2);
    CHECK_DIM_GE(na1, 1);
    CHECK_DIM_GE(na2, 1);
    if (na1->ndim == 1) {
        trim1 = 1;
        a1 = na_expand_dims(a1, INT2FIX(1));
        GetNArray(a1, na1);
    }
    if (na2->ndim == 1) {
        trim2 = 1;
        a2 = na_expand_dims(a2, INT2FIX(0));
        GetNArray(a2, na2);
    }
    n11 = na1->shape[na1->ndim-2]; // k
    n12 = na1->shape[na1->ndim-1]; // m
    n21 = na2->shape[na2->ndim-2]; // n
    n22 = na2->shape[na2->ndim-1]; // k
    if (n11 != n22) {
        rb_raise(nary_eShapeError, "matrix dimension mismatch: TBD");
        // "n11=%"SZF"u n12=%"SZF"u n21=%"SZF"u", n11, n12, n21
    }
    m = n12;
    n = n21;
    k = n11;
    shape[0] = m;
    shape[1] = n;

    ans = na_ndloop3(&ndf, 0, 2, a1, a2);

    GetNArray(ans, na);
    assert(na->type == NARRAY_DATA_T);
    if (trim2) {
        if (trim1) {
            return rb_funcall(ans, rb_intern("[]"), 1, INT2FIX(0));
        }
        assert(na->shape[na->ndim-1] == 1);
        --na->ndim;
    } else if (trim1) {
        assert(na->shape[na->ndim-2] == 1);
        na->shape[na->ndim-2] = na->shape[na->ndim-1];
        na->shape[na->ndim-1] = 1;
        --na->ndim;
    }
    return ans;
}
