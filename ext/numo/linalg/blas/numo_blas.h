#include "cblas.h"
#include "cblas_t.h"

typedef int blasint;

extern VALUE na_expand_dims(VALUE self, VALUE vdim);

#define CHECK_NARRAY_TYPE(x,t)                                 \
    if (CLASS_OF(x)!=(t)) {                                    \
        rb_raise(rb_eTypeError,"invalid NArray type (class)"); \
    }

// Error Class ??
#define CHECK_DIM_GE(na,nd)                                     \
    if ((na)->ndim<(nd)) {                                      \
        rb_raise(nary_eShapeError,                              \
                 "n-dimension=%d, but >=%d is expected",        \
                 (na)->ndim, (nd));                             \
    }

#define CHECK_DIM_EQ(na1,nd)                                    \
    if ((na1)->ndim != (nd)) {                                  \
        rb_raise(nary_eShapeError,                              \
                 "dimention mismatch: %d != %d",                \
                 (na1)->ndim, (nd));                            \
    }

#define CHECK_SQUARE(name,na)                                           \
    if ((na)->shape[(na)->ndim-1] != (na)->shape[(na)->ndim-2]) {       \
        rb_raise(nary_eShapeError,"%s is not square matrix",name);      \
    }

#define CHECK_SIZE_GE(na,sz)                                    \
    if ((na)->size < (size_t)(sz)) {                            \
        rb_raise(nary_eShapeError,                              \
                 "NArray size must be >= %"SZF"u",(size_t)(sz));\
    }

#define CHECK_NON_EMPTY(na)                                     \
    if ((na)->size==0) {                                        \
        rb_raise(nary_eShapeError,"empty NArray");              \
    }

#define CHECK_SIZE_EQ(n,m)                                      \
    if ((n)!=(m)) {                                             \
        rb_raise(nary_eShapeError,                              \
                 "size mismatch: %"SZF"d != %"SZF"d",           \
                 (size_t)(n),(size_t)(m));                      \
    }

#define CHECK_SAME_SHAPE(na1,na2)                                \
    {   int i;                                                   \
        CHECK_DIM_EQ(na1,na2->ndim);                             \
        for (i=0; i<na1->ndim; i++) {                            \
            CHECK_SIZE_EQ(na1->shape[i],na2->shape[i]);          \
        }                                                        \
    }

#define CHECK_INT_EQ(sm,m,sn,n)                          \
    if ((m) < (n)) {                                     \
        rb_raise(nary_eShapeError,                       \
                 "%s must be == %s: %s=%d %s=%d",        \
                 sm,sn,sm,m,sn,n);                       \
    }

// Error Class ??
#define CHECK_LEADING_GE(sld,ld,sn,n)                    \
    if ((ld) < (n)) {                                    \
        rb_raise(nary_eShapeError,                       \
                 "%s must be >= max(%s,1): %s=%d %s=%d", \
                 sld,sn,sld,ld,sn,n);                    \
    }

#define swap(a,b) {tmp=a;a=b;b=tmp;}

#define option_order numo_cblas_option_order
extern enum CBLAS_ORDER numo_cblas_option_order(VALUE trans);

#define option_trans numo_cblas_option_trans
extern enum CBLAS_TRANSPOSE numo_cblas_option_trans(VALUE trans);

#define option_uplo numo_cblas_option_uplo
extern enum CBLAS_UPLO numo_cblas_option_uplo(VALUE uplo);

#define option_diag numo_cblas_option_diag
extern enum CBLAS_DIAG numo_cblas_option_diag(VALUE diag);

#define option_side numo_cblas_option_side
extern enum CBLAS_SIDE numo_cblas_option_side(VALUE side);

#define check_func numo_cblas_check_func
extern void numo_cblas_check_func(void **func, const char *name);
