/*
* DGEEV computes the eigenvalues and, optionally, the left and/or
* right eigenvectors for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
*                         LDVR, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WI( * ), WORK( * ), WR( * )
*       ..
*
* ZGEEV computes the eigenvalues and, optionally, the left and/or
* right eigenvectors for GE matrices
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
*                         WORK, LWORK, RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   W( * ), WORK( * )
*       ..
*/

#define geev FFUNC(<%=blas_char%>geev)

void geev(
  const char *jobvl, const char *jobvr,
  fortran_integer *n, dtype *a, fortran_integer *lda,
<% if is_complex %>
  dtype *w,
<% else %>
  dtype *wr, dtype *wi,
<% end %>
  dtype *vl, fortran_integer *ldvl,
  dtype *vr, fortran_integer *ldvr,
  dtype *work, fortran_integer *lwork,
<% if is_complex %>
  rtype *rwork,
<% end %>
  fortran_integer *info
);

typedef struct {
    dtype *work;
    rtype *rwork;
    dtype *wr, *wi, *vrr;
    fortran_integer lwork;
} geev_work_t;


static void
<%=c_iter%>(na_loop_t *const lp)
{
    const char *jobvl="N", *jobvr="V";
    dtype *a;
    ctype *w, *vr;
    fortran_integer n, lda, ldvr, lwork, info=0;
    geev_work_t *opt;
    dtype *work;
    size_t m;
    <% if is_complex %>
    rtype *rwork;
    <% else %>
    size_t i, j;
    dtype *wr, *wi, *vrr;
    dtype *p1;
    ctype *p2;
    <% end %>

    opt = (geev_work_t*)(lp->opt_ptr);

    // a[n,lda], w[n], vr[n,ldvl]
    a   = (dtype*)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    w   = (ctype*)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    vr  = (ctype*)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    n   = m = lp->args[1].shape[0];
    lda  = lp->args[0].iter[0].step / sizeof(dtype);
    ldvr = lp->args[2].iter[0].step / sizeof(ctype);
    lwork = opt->lwork;

    work = opt->work;
    <% if is_complex %>
    rwork = opt->rwork;
    geev(jobvl, jobvr, &n, a, &lda, w, 0, &ldvr, vr, &ldvr,
         work, &lwork, rwork, &info);
    <% else %>
    wr = opt->wr;
    wi = opt->wi;
    vrr = opt->vrr;
    ldvr = n;
    geev(jobvl, jobvr, &n, a, &lda, wr, wi, 0, &ldvr, vrr, &ldvr,
         work, &lwork, &info);
    for (j=0; j<m; j++) {
        REAL(w[j]) = wr[j];
        IMAG(w[j]) = wi[j];
    }
    for (j=0; j<m; j++) {
        if (IMAG(w[j]) == 0.0) { /* real eigenvalue */
            p1 = &(vrr[j*ldvr]);
            p2 = &(vr[j*ldvr]);
            for (i=0; i<m; i++) {
                REAL(p2[i]) =  p1[i];
                IMAG(p2[i]) =  0;
            }
        } else { /* complex eigenvalue */
            /* v(j)   = VR(:,j) + i*VR(:,j+1) and
               v(j+1) = VR(:,j) - i*VR(:,j+1) */
            p1 = &(vrr[j*ldvr]);
            p2 = &(vr[j*ldvr]);
            for (i=0; i<m; i++) {
                REAL(p2[i]) = p1[i];
                IMAG(p2[i]) = p1[i+ldvr];
                REAL(p2[i+ldvr]) = p1[i];
                IMAG(p2[i+ldvr]) = -p1[i+ldvr];
            }
            j++;
        }
    }
    <% end %>
}


#define SET_POS(pos,i,type,n) {(pos)[i] = (pos)[(i)-1] + ((sizeof(type)*(n)-1)/16+1)*16;}

static inline geev_work_t *
set_opt_complex(size_t n, fortran_integer lwork, VALUE *vopt)
{
    size_t *sz;
    geev_work_t *opt;
    char *ptr;

    sz = ALLOCA_N(size_t,4);
    sz[0] = 0;
    SET_POS(sz,1,geev_work_t,1);
    SET_POS(sz,2,dtype,lwork); // work
    SET_POS(sz,3,rtype,n*2);   // rwork
    ptr = rb_alloc_tmp_buffer(vopt, sz[3]);
    opt = (geev_work_t*)ptr;
    opt->work  = (dtype*)(ptr+sz[1]);
    opt->rwork = (rtype*)(ptr+sz[2]);
    opt->lwork = lwork;
    return opt;
}

static inline geev_work_t *
set_opt_real(size_t n, fortran_integer lwork, VALUE *vopt)
{
    size_t *sz;
    geev_work_t *opt;
    char *ptr;

    sz = ALLOCA_N(size_t,6);
    sz[0] = 0;
    SET_POS(sz,1,geev_work_t,1);
    SET_POS(sz,2,dtype,lwork); // work
    SET_POS(sz,3,dtype,n);     // wr
    SET_POS(sz,4,dtype,n);     // wi
    SET_POS(sz,5,dtype,n*n);   // vrr
    ptr = rb_alloc_tmp_buffer(vopt, sz[5]);
    opt = (geev_work_t*)ptr;
    opt->work = (dtype*)(ptr+sz[1]);
    opt->wr   = (dtype*)(ptr+sz[2]);
    opt->wi   = (dtype*)(ptr+sz[3]);
    opt->vrr  = (dtype*)(ptr+sz[4]);
    opt->lwork = lwork;
    return opt;
}

/*
  @overload eigen(a)
  @param [Numo::<%=class_name%>] a >=2-dimentional NArray.
  @return [[Numo::<%=complex_class_name%>,Numo::<%=complex_class_name%>]] pair of eigenvalue and right eigenvector
  @raise

  <%=blas_char%>geev - computes the eigenvalues and the right eigenvectors
  for an N-by-N real nonsymmetric matrix A.
  The right eigenvector v(j) of A satisfies
                        A * v(j) = lambda(j) * v(j)
  where lambda(j) is its eigenvalue.
  The computed eigenvectors are normalized to have
  Euclidean norm equal to 1 and largest component real.
*/
static VALUE
<%=c_func%>(VALUE mod, VALUE a)
{
    const char *chr = "N";
    dtype wk[1];
    fortran_integer m, lwork, info=0;
    geev_work_t *opt;
    VALUE vopt, ans;

    narray_t *na;
    size_t n;
    size_t eval_shape[1];
    size_t evec_shape[2];
    ndfunc_arg_in_t ain[1] = {{cT,2}};
    ndfunc_arg_out_t aout[2] = {{cCT,1,eval_shape},{cCT,2,evec_shape}};
    ndfunc_t  ndf = {<%=c_iter%>, NO_LOOP, 1, 2, ain, aout};

    GetNArray(a,na);
    CHECK_DIM_GE(na,2);
    n = na->shape[na->ndim-1];
    if (n != na->shape[na->ndim-2]) {
        rb_raise(nary_eShapeError,"not square-matrix");
    }
    eval_shape[0] = n;
    evec_shape[0] = n;
    evec_shape[1] = n;

    m = n;
    lwork = -1;
    info = 0;
    <% if is_complex %>
    geev(chr, chr, &m, 0, &m, 0, 0, &m,  0, &m, wk, &lwork, 0, &info);
    lwork = REAL(wk[0]);
    opt = set_opt_complex(n, lwork, &vopt);
    <% else %>
    geev(chr, chr, &m, 0, &m, 0, 0, 0, &m, 0, &m, wk, &lwork, &info);
    lwork = wk[0];
    opt = set_opt_real(n, lwork, &vopt);
    <% end %>

    ans = na_ndloop3(&ndf, opt, 1, na_copy(a));

    rb_free_tmp_buffer(&vopt);
    return ans;
}
