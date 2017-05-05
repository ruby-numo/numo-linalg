#define func_p <%=cblas_func%>_p

static <%=cblas_func%>_t func_p = 0;

static void
<%=c_iter%>(na_loop_t *const lp)
{
    char *p1, *p2;
    size_t n;
    ssize_t s1, s2;
    dtype  *g;

    INIT_COUNTER(lp,n);
    INIT_PTR(lp,0,p1,s1);
    INIT_PTR(lp,1,p2,s2);
    g = (dtype*)(lp->opt_ptr);

    (*func_p)(n, (dtype*)p1, s1/sizeof(dtype),
                 (dtype*)p2, s2/sizeof(dtype), g);
}

/*
*  Definition:
*  ===========
*
*       SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       REAL SPARAM(5),SX(*),SY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
*>
*>    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
*>    (SX**T)
*>
*>    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
*>    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
*>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*>
*>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
*>
*>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
*>    H=(          )    (          )    (          )    (          )
*>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
*>    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in,out] SX
*> \verbatim
*>          SX is REAL array, dimension N
*>         double precision vector with N elements
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of SX
*> \endverbatim
*>
*> \param[in,out] SY
*> \verbatim
*>          SY is REAL array, dimension N
*>         double precision vector with N elements
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of SY
*> \endverbatim
*>
*> \param[in,out] SPARAM
*> \verbatim
*>          SPARAM is REAL array, dimension 5
*>     SPARAM(1)=SFLAG
*>     SPARAM(2)=SH11
*>     SPARAM(3)=SH21
*>     SPARAM(4)=SH12
*>     SPARAM(5)=SH22
*> \endverbatim
*/
/*
 *  @overload <%=name%>( x, y, param )
 *  @param [Numo::NArray] x  1-dimentional NArray. [in/out]
 *  @param [Numo::NArray] y  1-dimentional NArray. [in/out]
 *  @param [Array,NArray] param
 *  @return [nil]
 *  @raise
 *
 *    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
 *
 *    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
 *    (SX**T)
 *
 *    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
 *    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
 *    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
 *
 *    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
 *
 *      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
 *    H=(          )    (          )    (          )    (          )
 *      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
 *    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
 */
static VALUE
<%=c_func(3)%>(VALUE mod, VALUE x, VALUE y, VALUE param)
{
    dtype *g;
    narray_t *na1, *na2, *nap;
    ndfunc_arg_in_t ain[2] = {{OVERWRITE,0},{OVERWRITE,0}};
    ndfunc_t ndf = {<%=c_iter%>, STRIDE_LOOP, 2,0, ain,0};

    check_func((void*)(&func_p),"<%=cblas_func%>");

    CHECK_NARRAY_TYPE(x,cT);
    CHECK_NARRAY_TYPE(y,cT);
    GetNArray(x,na1);
    GetNArray(y,na2);
    CHECK_DIM_GE(na1,1);
    CHECK_DIM_GE(na2,1);
    CHECK_NON_EMPTY(na1);
    CHECK_NON_EMPTY(na2);
    CHECK_SAME_SHAPE(na1,na2);

    param = rb_funcall(cT,rb_intern("cast"),1,param);
    GetNArray(param,nap);
    CHECK_DIM_EQ(nap,1);
    CHECK_SIZE_GE(nap,5);
    g = (dtype*)nary_get_pointer_for_read(param);

    na_ndloop3(&ndf, g, 2, x, y);

    RB_GC_GUARD(param);
    return Qnil;
}

#undef func_p
