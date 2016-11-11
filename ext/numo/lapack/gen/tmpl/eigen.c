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
  char const * /*JOBVL*/, char const * /*JOBVR*/,
  fortran_integer * /*N*/, dtype * /*A*/, fortran_integer * /*LDA*/,
<% if is_complex %>
  dtype * /*W*/,
<% else %>
  dtype * /*WR*/, dtype * /*WI*/,
<% end %>
  dtype * /*VL*/,   fortran_integer * /*LDVL*/,
  dtype * /*VR*/,   fortran_integer * /*LDVR*/,
  dtype * /*WORK*/, fortran_integer * /*LWORK*/,
<% if is_complex %>
  rtype * /*RWORK*/,
<% end %>
  fortran_integer * /*INFO*/);

#define SET_POS(pos, i, type, n) do {(pos)[i] = (pos)[(i)-1] + ((sizeof(type)*(n)-1)/16+1)*16;} while (0)

static void
<%=c_iter%>(na_loop_t * const lp)
{
    char const * const jobvl="N", * const jobvr="V";
    size_t n0;
    dtype *a;
    ctype *w, *vr;
    fortran_integer n, lda, ldvr, lwork, info=0;
    dtype *work;
    <% if is_complex %>
    rtype *rwork;
    <% else %>
    dtype *wr, *wi, *vrr;
    <% end %>

    // a[n, lda], w[n], vr[n, ldvl]
    a  = (dtype *)(lp->args[0].ptr + lp->args[0].iter[0].pos);
    w  = (ctype *)(lp->args[1].ptr + lp->args[1].iter[0].pos);
    vr = (ctype *)(lp->args[2].ptr + lp->args[2].iter[0].pos);
    n0 = lp->args[1].shape[0];
    n  = n0;
    lda  = n;
    ldvr = n;
    lwork = *(fortran_integer *)lp->opt_ptr;
    {
        volatile VALUE tmp_wk;
        char *ptr;

        <% if is_complex %>

        size_t ofs[3];
        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, lwork);
        SET_POS(ofs, 2, rtype, 2*n);
        ptr = rb_alloc_tmp_buffer(&tmp_wk, ofs[2]);
        work  = (dtype *) ptr;
        rwork = (rtype *)(ptr + ofs[1]);

        geev(jobvl, jobvr, &n, a, &lda, w, 0, &ldvr, vr, &ldvr,
             work, &lwork, rwork, &info);

        <% else %>

        size_t ofs[5];
        ofs[0] = 0;
        SET_POS(ofs, 1, dtype, lwork);
        SET_POS(ofs, 2, dtype, n);
        SET_POS(ofs, 3, dtype, n);
        SET_POS(ofs, 4, dtype, n*n);
        ptr = rb_alloc_tmp_buffer(&tmp_wk, ofs[4]);
        work = (dtype *) ptr;
        wr   = (dtype *)(ptr + ofs[1]);
        wi   = (dtype *)(ptr + ofs[2]);
        vrr  = (dtype *)(ptr + ofs[3]);

        geev(jobvl, jobvr, &n, a, &lda, wr, wi, 0, &ldvr, vrr, &ldvr,
             work, &lwork, &info);
        {
            size_t j;
            for (j=0; j<n0; ++j) {
                REAL(w[j]) = wr[j];
                IMAG(w[j]) = wi[j];
            }
            for (j=0; j<n0; ++j) {
                dtype * const p1 = &(vrr[j*ldvr]);
                ctype * const p2 = &(vr[j*ldvr]);
                if (IMAG(w[j]) == 0.0) { /* real eigenvalue */
                    size_t i;
                    for (i=0; i<n0; ++i) {
                        REAL(p2[i]) =  p1[i];
                        IMAG(p2[i]) =  0;
                    }
                } else { /* complex eigenvalue */
                    size_t i;
                    /* v(j)   = VR(:,j) + i*VR(:,j+1) and
                       v(j+1) = VR(:,j) - i*VR(:,j+1) */
                    for (i=0; i<n0; ++i) {
                        REAL(p2[i])      =  p1[i];
                        IMAG(p2[i])      =  p1[i+ldvr];
                        REAL(p2[i+ldvr]) =  p1[i];
                        IMAG(p2[i+ldvr]) = -p1[i+ldvr];
                    }
                    ++j;
                }
            }
        }

        <% end %>

        rb_free_tmp_buffer(&tmp_wk);
        RB_GC_GUARD(tmp_wk);
    }
}

#undef SET_POS

/*
  @overload eigen(a)
  @param ***TODO*** [Numo::<%=class_name%>] a >=2-dimentional NArray.
  @return ***TODO*** [[Numo::<%=complex_class_name%>,Numo::<%=complex_class_name%>]] pair of eigenvalue and right eigenvector
  @raise

  ***TODO***
  <%=blas_char%>geev - computes the eigenvalues and the right eigenvectors
  for an N-by-N real nonsymmetric matrix A.
  The right eigenvector v(j) of A satisfies
                        A * v(j) = lambda(j) * v(j)
  where lambda(j) is its eigenvalue.
  The computed eigenvectors are normalized to have
  Euclidean norm equal to 1 and largest component real.
*/
static VALUE
<%=c_func%>(VALUE const mod, VALUE const a)
{
    volatile VALUE ans;
    fortran_integer lwork;

    narray_t *na;
    size_t n;
    size_t eval_shape[1];
    size_t evec_shape[2];
    ndfunc_arg_in_t ain[1] = {{cT, 2}};
    ndfunc_arg_out_t aout[2] = {{cCT, 1, eval_shape}, {cCT, 2, evec_shape}};
    ndfunc_t ndf = {<%=c_iter%>, NO_LOOP, 1, 2, ain, aout};

    GetNArray(a, na);
    CHECK_DIM_GE(na, 2);
    n = na->shape[na->ndim-1];
    if (n != na->shape[na->ndim-2]) {
        rb_raise(nary_eShapeError, "not square-matrix");
    }
    eval_shape[0] = n;
    evec_shape[0] = n;
    evec_shape[1] = n;
    {
        char const * const chr="N";
        fortran_integer m=n;
        fortran_integer info=0;
        dtype wk[1];
        lwork = -1;
        <% if is_complex %>
        geev(chr, chr, &m, 0, &m, 0, 0, &m,  0, &m, wk, &lwork, 0, &info);
        lwork = REAL(wk[0]);
        <% else %>
        geev(chr, chr, &m, 0, &m, 0, 0, 0, &m, 0, &m, wk, &lwork, &info);
        lwork = wk[0];
        <% end %>
    }
    ans = na_ndloop3(&ndf, &lwork, 1, a);
    return ans;
}
