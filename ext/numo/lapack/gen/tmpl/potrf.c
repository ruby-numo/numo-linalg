/*
* DPOTRF computes the Cholesky factorization
*
*  Definition:
*  ===========
*
*       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
* ZPOTRF computes the Cholesky factorization
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*/

#define potrf FFUNC(<%=blas_char%>potrf)

typedef struct {
    int upper;
    int overwrite;
} potrf_opt_t;

void potrf(
  char const * /*UPLO*/,
  fortran_integer * /*N*/, dtype * /*A*/, fortran_integer * /*LDA*/,
  fortran_integer * /*INFO*/);

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const uplo_u = "U", * const uplo_l = "L";
    char const *uplo;
    potrf_opt_t *opt;
    size_t n1;
    dtype *ai, *ao, *a;
    fortran_integer n, lda, info=0;

    opt = lp->opt_ptr;

    uplo = ( opt->upper ? uplo_u : uplo_l ) ;

    ai = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    ao = (rtype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    n1    = lp->args[0].shape[0];
    //n2 = lp->args[0].shape[1];  // (lda)
    if (opt->overwrite) {
        a = ai;
    } else {
        memcpy(ao, ai, n1*n1*sizeof(dtype));
        a = ao;
    }
    n = (fortran_integer)n1;
    lda = n;

    potrf(uplo, &n, a, &lda, &info);

    if (opt->overwrite) {
        memcpy(ao, a, n1*n1*sizeof(dtype));
    }
    {
        size_t i, j;
        dtype *p;
        <% if is_complex %>
        dtype zero;
        REAL(zero) = 0.0;
        IMAG(zero) = 0.0;
        <% else %>
        dtype const zero = 0.0;
        <% end %>
        if (opt->upper) {
            for (i = 0; i < n1; ++i) {
                p = &(ao[i*n1]);
                for (j = i+1; j < n1; ++j) {
                    *(p+j) = zero;
                }
            }
        } else {
            for (i = 0; i < n1; ++i) {
                p = &(ao[i*n1]);
                for (j = 0; j < i; ++j) {
                    *(p+j) = zero;
                }
            }
        }
    }
}

#define sub_func_name(f, args) f##_sub args

/*
  @overload potrf(a)
  @param [Numo::<%=class_name%>] a  >=2-dimentional symmetric NArray.
  @return [Numo::<%=class_name%>] TBD
  @raise

  TBD
*/
static VALUE
sub_func_name(<%=c_func%>, (VALUE const a, int const upper, int const overwrite))
{
    volatile VALUE ans;
    narray_t *na;
    size_t n1, n2;
    potrf_opt_t opt;

    GetNArray(a, na);
    CHECK_DIM_GE(na, 2);
    n1 = na->shape[na->ndim-1];
    n2 = na->shape[na->ndim-2];

    if (n1 != n2) {
        rb_raise(nary_eShapeError, "not square-matrix");
    }
    opt.upper = upper;
    opt.overwrite = overwrite;
    {
        size_t shape[2];
        ndfunc_arg_in_t ain[1] = {{cT, 2}};
        ndfunc_arg_out_t aout[1] = {
            {cRT, COUNT_OF_(shape), shape}};
        ndfunc_t ndf = {&<%=c_iter%>, NO_LOOP,
            COUNT_OF_(ain), COUNT_OF_(aout),
            ain, aout};
        shape[1] = n1;
        shape[0] = n1;
        ans = na_ndloop3(&ndf, &opt, 1, a);
    }
    return ans;
}

static VALUE
<%=c_func%>(int argc, VALUE const argv[], VALUE UNUSED(mod))
{
    int const flg_overwrite = 0;
    int flg_upper=0;
    if (argc == 2) {
        VALUE h;
        if ( ! NIL_P(h = rb_check_hash_type(argv[1]))) {
            ID tbl;
            VALUE v;
            tbl = rb_intern("upper");
            rb_get_kwargs(h, &tbl, 0, 1, &v);
            if (v != Qundef) {
                flg_upper = RTEST(v);
            }
            --argc;
        }
        
    }
    rb_check_arity(argc, 1, 1);
    return sub_func_name(<%=c_func%>, (argv[0], flg_upper, flg_overwrite));
}

#undef sub_func_name
