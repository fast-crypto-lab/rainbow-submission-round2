#ifndef _BLAS_H_
#define _BLAS_H_

#include "blas_comm.h"
#include "blas_u32.h"

#define gf16v_mul_scalar _gf16v_mul_scalar_u32
#define gf16v_madd _gf16v_madd_u32

#define gf256v_add _gf256v_add_u32
#define gf256v_mul_scalar _gf256v_mul_scalar_u32
#define gf256v_madd _gf256v_madd_u32

#define gf256v_predicated_add _gf256v_predicated_add_u32
#define gf16v_dot _gf16v_dot_u32

#define gf16mat_prod _gf16mat_prod
#define gf16mat_gauss_elim _gf16mat_gauss_elim
#define gf16mat_solve_linear_eq _gf16mat_solve_linear_eq

#define gf256mat_prod _gf256mat_prod
#define gf256mat_gauss_elim _gf256mat_gauss_elim
#define gf256mat_solve_linear_eq _gf256mat_solve_linear_eq

#endif
