
#ifndef _P_MATRIX_OP_SSE_H_
#define _P_MATRIX_OP_SSE_H_


#include "blas_sse.h"

#include "parallel_matrix_op.h"


#ifdef  __cplusplus
extern  "C" {
#endif




////////////////////  "quadratric" matrix evaluation  ///////////////////////////////


///  y =  x^Tr * trimat * x

void batch_quad_trimat_eval_multab_gf16_sse( unsigned char * y, const unsigned char * trimat, const unsigned char * multab_x, unsigned dim , unsigned size_batch );

void batch_quad_trimat_eval_multab_gf256_sse( unsigned char * y, const unsigned char * trimat, const unsigned char * multab_x, unsigned dim , unsigned size_batch );






#ifdef  __cplusplus
}
#endif


#endif
