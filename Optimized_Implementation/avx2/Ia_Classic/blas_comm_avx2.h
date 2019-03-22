#ifndef _BLAS_COMM_AVX2_H_
#define _BLAS_COMM_AVX2_H_


#include "stdint.h"



#ifdef  __cplusplus
extern  "C" {
#endif


///////////////////////////////  GF( 16 ) ////////////////////////////////////////////////////


void gf16mat_prod_multab_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * multab );

void gf16mat_prod_avx2( uint8_t * c , const uint8_t * mat_a , unsigned a_h_byte , unsigned a_w , const uint8_t * b );

unsigned gf16mat_gauss_elim_avx2( uint8_t * mat , unsigned h , unsigned w );

unsigned gf16mat_solve_linear_eq_avx2( uint8_t * sol , const uint8_t * inp_mat , const uint8_t * c_terms , unsigned n );


///////////////////////////////  GF( 256 ) ////////////////////////////////////////////////////


//void gf256mat_prod_add_multab_avx2( __m128i * r , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * multab );

void gf256mat_prod_multab_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * multab );

//void gf256mat_prod_add_avx2( __m128i * r , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * b );

void gf256mat_prod_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * b );

unsigned gf256mat_gauss_elim_avx2( uint8_t * mat , unsigned h , unsigned w );

unsigned gf256mat_solve_linear_eq_avx2( uint8_t * sol , const uint8_t * inp_mat , const uint8_t * c_terms, unsigned n );



#ifdef  __cplusplus
}
#endif



#endif
