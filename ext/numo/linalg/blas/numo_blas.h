#include "cblas.h"
#include "cblas_t.h"
#include "../numo_linalg.h"

typedef int blasint;

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
