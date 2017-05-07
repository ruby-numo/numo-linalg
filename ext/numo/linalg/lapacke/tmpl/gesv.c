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

<% uplo = (/^?ge/=~name) ? nil : "g->uplo," %>
<% ipiv = (/^?po/=~name) ? nil : "ipiv," %>
*/
#define args_t <%=func_name%>_args_t
typedef struct {
    int order; char uplo;
} args_t;

static <%=func_name%>_t <%=func_name%>_p = 0;

static void
<%=c_iter%>(na_loop_t * const lp)
{
    dtype *a, *b;
    int   *info;
    int    n, nrhs;
    int    lda, ldb;
    args_t *g;
 <% if ipiv %>
    int *ipiv;
    ipiv = (int*)NDL_PTR(lp,2);
    info = (int*)NDL_PTR(lp,3);
 <% else %>
    info = (int*)NDL_PTR(lp,2);
 <% end %>
    a = (dtype*)NDL_PTR(lp,0);
    b = (dtype*)NDL_PTR(lp,1);
    g = (args_t*)(lp->opt_ptr);

    n = lp->args[0].shape[0];
    lda = lp->args[0].iter[0].step / sizeof(dtype);
    if (lp->args[1].ndim == 1) {
        nrhs = 1;
        ldb = (g->order==LAPACK_COL_MAJOR) ? n : 1;
    } else {
        nrhs = lp->args[1].shape[1];
        ldb = lp->args[1].iter[0].step / sizeof(dtype);
    }
    //printf("order=%d n=%d nrhs=%d lda=%d ldb=%d b.ndim=%d\n",
    //       g->order,n,nrhs,lda,ldb,lp->args[1].ndim);
    *info = (*<%=func_name%>_p)(g->order, <%=uplo%>
                                n, nrhs, a, lda, <%=ipiv%> b, ldb);
    CHECK_ERROR(*info);
}

/*
<% if uplo %>
  @overload <%=name%>(a, b [,uplo, order])
<% else %>
  @overload <%=name%>(a, b [,order])
<% end %>
  @param [Numo::NArray] a  >=2-dimentional NArray.
  @param [Numo::NArray] b  >=1-dimentional NArray.
<% if ipiv %>
  @return [Array] array of [lu,x,piv,info]
<% else %>
  @return [Array] array of [lu,x,info]
<% end %>

  <%=name%> - computes the solution to a complex system of linear equations
     A  *  X = B,
  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
  The LU decomposition with partial pivoting and row interchanges is used to factor A as
     A = P * L * U,
  where P is a permutation matrix, L is unit lower triangular, and U is upper triangular.
  The factored form of A is then used to solve the system of equations A * X = B.
*/
static VALUE
<%=c_func(-1)%>(int argc, VALUE *argv, VALUE const mod)
{
    VALUE a, b, order, ans;
    narray_t *na1, *na2;
    size_t n11, n12, n21, n22;
    ndfunc_arg_in_t ain[2] = {{OVERWRITE,2},{OVERWRITE,2}};
    size_t shape[2];
    ndfunc_arg_out_t aout[2] = {{cInt,1,shape},{cInt,0}};
    ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP|NDF_EXTRACT, 2, 2, ain, aout};
    int i;
    args_t g = {LAPACK_ROW_MAJOR, 'U'};

 <% if uplo %>
    VALUE uplo;
    i = rb_scan_args(argc, argv, "22", &a, &b, &uplo, &order);
    switch (i) {
    case 4: g.order = option_order(order);
    case 3: g.uplo = option_uplo(uplo);
    }
 <% else %>
    i = rb_scan_args(argc, argv, "21", &a, &b, &order);
    switch (i) {
    case 3: g.order = option_order(order);
    }
 <% end %>

    check_func((void*)(&<%=func_name%>_p),"<%=func_name%>");

    a = rb_funcall(cT,rb_intern("cast"),1,a);
    if (!TEST_INPLACE(a)) {a = na_copy(a);}
    b = rb_funcall(cT,rb_intern("cast"),1,b);
    if (!TEST_INPLACE(b)) {b = na_copy(b);}

    GetNArray(a, na1);
    GetNArray(b, na2);
    CHECK_DIM_GE(na1, 2);
    CHECK_DIM_GE(na2, 1);
    n11 = na1->shape[na1->ndim-2]; // n
    n12 = na1->shape[na1->ndim-1]; // n
    if (NA_NDIM(na2) == 1) {
        ain[1].dim = 1;
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
 <% if !ipiv %>
    aout[0] = aout[1];
    ndf.nout--;
 <% end %>

    ans = na_ndloop3(&ndf, &g, 2, a, b);
 <% if ipiv %>
    return rb_ary_concat(rb_assoc_new(a,b),ans);
 <% else %>
    return rb_ary_push(rb_assoc_new(a,b),ans);
 <% end %>
}

#undef args_t
