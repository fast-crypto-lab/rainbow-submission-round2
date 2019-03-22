#ifndef _BLAS_H_
#define _BLAS_H_


#include <stdint.h>

#include <stdio.h>




#include "blas_config.h"



#ifdef _BLAS_UINT64_

#include "blas_u64.h"


#define gf16v_mul_scalar  _gf16v_mul_scalar_u64
#define gf16v_madd        _gf16v_madd_u64

#define gf256v_add        _gf256v_add_u64
#define gf256v_mul_scalar  _gf256v_mul_scalar_u64
#define gf256v_madd        _gf256v_madd_u64

#define gf256v_predicated_add      _gf256v_predicated_add_u64

#else
error here.
#endif


#include "blas_comm.h"




#endif

