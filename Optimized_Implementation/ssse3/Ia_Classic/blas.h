#ifndef _BLAS_H_
#define _BLAS_H_


#include <stdint.h>

#include <stdio.h>




#include "blas_config.h"




#ifdef _BLAS_SSE_

#include "blas_sse.h"

#define gf16v_mul_scalar   gf16v_mul_scalar_sse
#define gf16v_madd         gf16v_madd_sse
#define gf16v_madd_multab  gf16v_madd_multab_sse

#define gf256v_add         gf256v_add_sse
#define gf256v_mul_scalar  gf256v_mul_scalar_sse
#define gf256v_madd        gf256v_madd_sse
#define gf256v_madd_multab gf256v_madd_multab_sse

#define gf16mat_prod              gf16mat_prod_sse
#define gf16mat_gauss_elim        gf16mat_gauss_elim_sse
#define gf16mat_solve_linear_eq   gf16mat_solve_linear_eq_sse

#define gf256mat_prod             gf256mat_prod_sse
#define gf256mat_gauss_elim       gf256mat_gauss_elim_sse
#define gf256mat_solve_linear_eq  gf256mat_solve_linear_eq_sse


#include "blas_u64.h"
#define gf256v_predicated_add	_gf256v_predicated_add_u64

/// faster
#define gf16v_dot	gf16v_dot_sse
//#define gf16v_dot	gf16v_dot_avx2

#else
error here.
#endif



#include "blas_comm.h"




#endif

