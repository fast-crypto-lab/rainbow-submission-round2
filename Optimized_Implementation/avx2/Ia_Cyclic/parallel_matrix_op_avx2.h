
#ifndef _P_MATRIX_OP_AVX2_H_
#define _P_MATRIX_OP_AVX2_H_


#include "blas_sse.h"
#include "blas_avx2.h"

#include "parallel_matrix_op.h"


#ifdef  __cplusplus
extern  "C" {
#endif




////////////////////  "quadratric" matrix evaluation  ///////////////////////////////


///  y =  x^Tr * trimat * x

void batch_quad_trimat_eval_multab_gf16_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * multab_x, unsigned dim , unsigned size_batch );

void batch_quad_trimat_eval_multab_gf256_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * multab_x, unsigned dim , unsigned size_batch );






#ifdef  __cplusplus
}
#endif


#endif
