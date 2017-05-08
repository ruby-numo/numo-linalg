#define lapack_complex_float scomplex
#define lapack_complex_double dcomplex
#include "lapacke.h"
#include "lapacke_t.h"
#include "../numo_linalg.h"

typedef int blasint;

#define option_order numo_lapacke_option_order
extern int numo_lapacke_option_order(VALUE order);

#define option_job numo_lapacke_option_job
extern char numo_lapacke_option_job(VALUE job);

#define option_trans numo_lapacke_option_trans
extern char numo_lapacke_option_trans(VALUE trans);

#define option_uplo numo_lapacke_option_uplo
extern char numo_lapacke_option_uplo(VALUE uplo);

#define option_diag numo_lapacke_option_diag
extern char numo_lapacke_option_diag(VALUE diag);

#define option_side numo_lapacke_option_side
extern char numo_lapacke_option_side(VALUE side);

#define check_func numo_lapacke_check_func
extern void numo_lapacke_check_func(void **func, const char *name);

#define CHECK_ERROR(info)                                           \
    { if ((info)<0) {                                               \
            rb_raise(eLapackError,"LAPACK error, info=%d",(info)); \
        }}

#define SWAP_IFCOL(order,a,b)                                           \
    { if ((order)==LAPACK_COL_MAJOR) {tmp=(a);(a)=(b);(b)=tmp;} }
