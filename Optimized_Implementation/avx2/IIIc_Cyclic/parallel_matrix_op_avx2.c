
#include "blas_comm.h"
#include "blas.h"

#include "parallel_matrix_op_avx2.h"



////////////////////  "quadratric" matrix evaluation  ///////////////////////////////


///  y =  x^Tr * trimat * x

void batch_quad_trimat_eval_multab_gf16_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * x_multab, unsigned dim , unsigned size_batch )
{
///    assert( dim <= 128 );
///    assert( size_batch <= 128 );
    unsigned char tmp[256] __attribute__((aligned(32)));

    gf256v_set_zero( y , size_batch );
    for(unsigned i=0;i<dim;i++) {
        gf256v_set_zero( tmp , size_batch );
        for(unsigned j=i;j<dim;j++) {
           gf16v_madd_multab( tmp , trimat , x_multab+16*j , size_batch );
           trimat += size_batch;
        }
        gf16v_madd_multab( y , tmp , x_multab+16*i , size_batch );
    }
}

void batch_quad_trimat_eval_multab_gf256_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * x_multab, unsigned dim , unsigned size_batch )
{
///    assert( dim <= 128 );
///    assert( size_batch <= 128 );
    unsigned char tmp[256] __attribute__((aligned(32)));

    gf256v_set_zero( y , size_batch );
    for(unsigned i=0;i<dim;i++) {
        gf256v_set_zero( tmp , size_batch );
        for(unsigned j=i;j<dim;j++) {
           gf256v_madd_multab( tmp , trimat , x_multab+32*j , size_batch );
           trimat += size_batch;
        }
        gf256v_madd_multab( y , tmp , x_multab+32*i , size_batch );
    }
}





