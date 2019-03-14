
#include "gf16.h"

#include <immintrin.h>

#include "blas_config.h"

#include "gf16_avx2.h"

#include "blas_sse.h"

#include "blas_avx2.h"

#include "assert.h"





void gf16mat_prod_multab_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * multab ) {
	assert( n_A_width <= 256 );
	assert( n_A_vec_byte <= 128 );
	if( 16 >= n_A_vec_byte ) { return gf16mat_prod_multab_sse(c,matA,n_A_vec_byte,n_A_width,multab); }

	__m256i mask_f = _mm256_load_si256( (__m256i*)__mask_low);

	__m256i r0[4];
	__m256i r1[4];
	unsigned n_ymm = ((n_A_vec_byte + 31)>>5);
	for(unsigned i=0;i<n_ymm;i++) r0[i] = _mm256_setzero_si256();
	for(unsigned i=0;i<n_ymm;i++) r1[i] = _mm256_setzero_si256();

	for(unsigned i=0;i<n_A_width;i++) {
		__m128i ml = _mm_load_si128( (__m128i*)( multab + i*16) );
		__m256i mt = _mm256_inserti128_si256(_mm256_castsi128_si256(ml),ml,1);
		for(unsigned j=0;j<n_ymm;j++) {
			__m256i inp = _mm256_loadu_si256( (__m256i*)(matA+j*32) );
			r0[j] ^= _mm256_shuffle_epi8( mt , inp&mask_f );
			r1[j] ^= _mm256_shuffle_epi8( mt , _mm256_srli_epi16(inp,4)&mask_f );
		}
		matA += n_A_vec_byte;
	}
	uint8_t temp[128] __attribute__((aligned(32)));
	for(unsigned i=0;i<n_ymm;i++) _mm256_store_si256( (__m256i*)(temp + i*32) , r0[i]^_mm256_slli_epi16(r1[i],4) );
	for(unsigned i=0;i<n_A_vec_byte;i++) c[i] = temp[i];
}

#if 0
void gf16mat_prod_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * b ) {
	assert( n_A_width <= 128 );
	assert( n_A_vec_byte <= 64 );

	uint8_t multab[128*16] __attribute__((aligned(32)));
	gf16v_generate_multab_sse( multab , b , n_A_width );

	gf16mat_prod_multab_avx2( c , matA , n_A_vec_byte , n_A_width , multab );
}
#else
void gf16mat_prod_avx2( uint8_t * c , const uint8_t * mat_a , unsigned a_h_byte , unsigned a_w , const uint8_t * b ) {
	assert( a_w <= 256 );
	assert( a_h_byte <= 512 );
	if( 16 >= a_h_byte ) { return gf16mat_prod_sse(c,mat_a,a_h_byte,a_w,b); }

	__m256i mask_f = _mm256_load_si256( (__m256i*)__mask_low);

	__m256i r0[16];
	__m256i r1[16];
	unsigned n_ymm = ((a_h_byte+31)>>5);
	for(unsigned i=0;i<n_ymm;i++) r0[i] = _mm256_setzero_si256();
	for(unsigned i=0;i<n_ymm;i++) r1[i] = _mm256_setzero_si256();

	uint8_t _x[256] __attribute__((aligned(32)));
	gf16v_split_sse( _x , b , a_w );
	for(unsigned i=0;i < ((a_w+31)>>5); i++) {
		__m256i lxi = tbl32_gf16_log( _mm256_load_si256( (__m256i*)(_x+i*32) ) );
		_mm256_store_si256( (__m256i*)(_x+i*32) , lxi );
	}

	for(unsigned i=0;i< a_w -1;i++) {
		_x[0] = _x[i];
		__m256i ml = _mm256_broadcastb_epi8( _mm_load_si128((__m128i*)_x) );
		for(unsigned j=0;j<n_ymm;j++) {
			__m256i inp = _mm256_loadu_si256( (__m256i*)(mat_a+j*32) );
			r0[j] ^= tbl32_gf16_mul_log( inp&mask_f , ml , mask_f );
			r1[j] ^= tbl32_gf16_mul_log( _mm256_srli_epi16(inp,4)&mask_f , ml , mask_f );
		}
		mat_a += a_h_byte;
	}
	unsigned n_32 = (a_h_byte>>5);
	unsigned n_32_rem = a_h_byte&0x1f;
	{
		/// last column
		unsigned i=a_w-1;
		_x[0] = _x[i];
		__m256i ml = _mm256_broadcastb_epi8( _mm_load_si128((__m128i*)_x) );
		for(unsigned j=0;j<n_32;j++) {
			__m256i inp = _mm256_loadu_si256( (__m256i*)(mat_a+j*32) );
			r0[j] ^= tbl32_gf16_mul_log( inp&mask_f , ml , mask_f );
			r1[j] ^= tbl32_gf16_mul_log( _mm256_srli_epi16(inp,4)&mask_f , ml , mask_f );
		}
		if( n_32_rem ) {
			unsigned j = n_32;
			__m256i inp = _load_ymm( mat_a+j*32 , n_32_rem );
			r0[j] ^= tbl32_gf16_mul_log( inp&mask_f , ml , mask_f );
			r1[j] ^= tbl32_gf16_mul_log( _mm256_srli_epi16(inp,4)&mask_f , ml , mask_f );
		}
	}

	for(unsigned i=0;i<n_32;i++) _mm256_storeu_si256( (__m256i*)(c + i*32) , r0[i]^_mm256_slli_epi16(r1[i],4) );
	if( n_32_rem ) _store_ymm( c + n_32*32 , n_32_rem , r0[n_32]^_mm256_slli_epi16(r1[n_32],4) );
}
#endif






static inline
uint8_t _if_zero_then_0xf(uint8_t p ) {
	return (p-1)>>4;
}

static inline
unsigned _linear_solver_32x32_avx2( uint8_t * r , const uint8_t * mat_32x32 , const uint8_t * cc )
{

	uint8_t mat[32*32] __attribute__((aligned(32)));
	for(unsigned i=0;i<32;i++) gf16v_split_sse( mat + i*32 , mat_32x32 + i*16 , 32 );

	__m256i mask_f = _mm256_load_si256((__m256i const *) __mask_low);

	uint8_t temp[32] __attribute__((aligned(32)));
	uint8_t pivots[32] __attribute__((aligned(32)));

	uint8_t rr8 = 1;
	for(unsigned i=0;i<32;i++) {
		for(unsigned j=0;j<32;j++) pivots[j] = mat[j*32+i];
			if( 0 == i ) {
				gf16v_split_sse( temp , cc , 32 );
				for(unsigned j=0;j<32;j++) mat[j*32] = temp[j];
			}
		__m256i rowi = _mm256_load_si256( (__m256i*)(mat+i*32) );
		for(unsigned j=i+1;j<32;j++) {
			temp[0] = _if_zero_then_0xf( pivots[i] );
			__m256i mask_zero = _mm256_broadcastb_epi8(_mm_load_si128((__m128i*)temp));

			__m256i rowj = _mm256_load_si256( (__m256i*)(mat+j*32) );
			rowi ^= mask_zero&rowj;
			//rowi ^= predicate_zero&(*(__m256i*)(mat+j*32));
			pivots[i] ^= temp[0]&pivots[j];
		}
		uint8_t is_pi_nz = _if_zero_then_0xf(pivots[i]);
		is_pi_nz = ~is_pi_nz;
		rr8 &= is_pi_nz;

		temp[0] = pivots[i];
		__m128i inv_rowi = tbl_gf16_inv( _mm_load_si128((__m128i*)temp) );
		pivots[i] = _mm_extract_epi8( inv_rowi , 0 );

		__m256i log_pivots = tbl32_gf16_log( _mm256_load_si256( (__m256i*)pivots ) );
		_mm256_store_si256( (__m256i*)pivots , log_pivots );

		temp[0] = pivots[i];
		__m256i logpi = _mm256_broadcastb_epi8( _mm_load_si128((__m128i*)temp) );
		rowi = tbl32_gf16_mul_log( rowi , logpi , mask_f );
		__m256i log_rowi = tbl32_gf16_log( rowi );
		for(unsigned j=0;j<32;j++) {
			if(i==j) {
				_mm256_store_si256( (__m256i*)(mat+j*32) , rowi );
				continue;
			}
			__m256i rowj = _mm256_load_si256( (__m256i*)(mat+j*32) );
			temp[0] = pivots[j];
			__m256i logpj = _mm256_broadcastb_epi8( _mm_load_si128((__m128i*)temp) );
			rowj ^= tbl32_gf16_mul_log_log( log_rowi , logpj , mask_f );
			_mm256_store_si256( (__m256i*)(mat+j*32) , rowj );
		}
	}

	for(unsigned i=0;i<32;i+=2) {
		//gf16v_set_ele( r , i , mat[i*32] );
		r[i>>1] = mat[i*32]| (mat[(i+1)*32]<<4);
	}
	return rr8;
}



unsigned gf16mat_solve_linear_eq_avx2( uint8_t * sol , const uint8_t * inp_mat , const uint8_t * c_terms , unsigned n )
{
	if( 32 == n ) {
		return _linear_solver_32x32_avx2( sol , inp_mat , c_terms );
	} else  {
		return gf16mat_solve_linear_eq_sse( sol , inp_mat , c_terms , n );
	}
}





///////////////////////////////  GF( 256 ) ////////////////////////////////////////////////////





void gf256mat_prod_multab_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * multab ) {
	assert( n_A_width <= 256 );
	assert( n_A_vec_byte <= 48*48 );
	if( 16 >= n_A_vec_byte ) { return gf256mat_prod_multab_sse(c,matA,n_A_vec_byte,n_A_width,multab); }

	__m256i mask_f = _mm256_load_si256((__m256i const *) __mask_low);

	__m256i r[72];   /// 48*48/32
	unsigned n_ymm = ((n_A_vec_byte + 31)>>5);
	for(unsigned i=0;i<n_ymm;i++) r[i] = _mm256_setzero_si256();

	for(unsigned i=0;i<n_A_width;i++) {
		__m256i mt = _mm256_load_si256( (__m256i*)( multab + i*32) );
		__m256i ml = _mm256_permute2x128_si256(mt,mt,0x00 );
		__m256i mh = _mm256_permute2x128_si256(mt,mt,0x11 );
		for(unsigned j=0;j<n_ymm;j++) {
			__m256i inp = _mm256_loadu_si256( (__m256i*)(matA+j*32) );
			r[j] ^= _mm256_shuffle_epi8( ml , inp&mask_f );
			r[j] ^= _mm256_shuffle_epi8( mh , _mm256_srli_epi16(inp,4)&mask_f );
		}
		matA += n_A_vec_byte;
	}
	uint8_t r8[256] __attribute__((aligned(32)));
	for(unsigned i=0;i<n_ymm;i++) _mm256_store_si256( (__m256i*)(r8 + i*32) , r[i] );
	for(unsigned i=0;i<n_A_vec_byte;i++) c[i] = r8[i];
}


#if 0
// slower
void gf256mat_prod_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * b ) {
	assert( n_A_width <= 128 );
	assert( n_A_vec_byte <= 80 );

	uint8_t multab[256*16] __attribute__((aligned(32)));
	gf256v_generate_multab_sse( multab , b , n_A_width );

	gf256mat_prod_multab_avx2( c , matA , n_A_vec_byte , n_A_width , multab );
}
#else
void gf256mat_prod_avx2( uint8_t * c , const uint8_t * matA , unsigned n_A_vec_byte , unsigned n_A_width , const uint8_t * b ) {
	assert( n_A_width <= 256 );
	assert( n_A_vec_byte <= 48*48 );
	//if( 256 < n_A_vec_byte ) return gf256mat_prod_no_buffer_avx2(c,matA,n_A_vec_byte,n_A_width,b);
	if( 16 >= n_A_vec_byte ) { return gf256mat_prod_sse(c,matA,n_A_vec_byte,n_A_width,b); }

	__m256i mask_f = _mm256_load_si256( (__m256i*)__mask_low);

	__m256i r[48*48/32];
	unsigned n_ymm = ((n_A_vec_byte + 31)>>5);
	for(unsigned i=0;i<n_ymm;i++) r[i] = _mm256_setzero_si256();

	uint8_t x0[256] __attribute__((aligned(32)));
	uint8_t x1[256] __attribute__((aligned(32)));
	for(unsigned i=0;i<n_A_width;i++) x0[i] = b[i];
	for(unsigned i=0;i< ((n_A_width+31)>>5);i++) {
		__m256i inp = _mm256_load_si256((__m256i*)(x0+i*32));
		__m256i i0 = inp&mask_f;
		__m256i i1 = _mm256_srli_epi16(inp,4)&mask_f;
		_mm256_store_si256((__m256i*)(x0+i*32),tbl32_gf16_log(i0));
		_mm256_store_si256((__m256i*)(x1+i*32),tbl32_gf16_log(i1));
	}

	for(unsigned i=0;i<n_A_width;i++) {
		x0[0] = x0[i]; __m256i m0 = _mm256_broadcastb_epi8( _mm_load_si128((__m128i*)x0) );
		x1[0] = x1[i]; __m256i m1 = _mm256_broadcastb_epi8( _mm_load_si128((__m128i*)x1) );
		//__m128i ml = _mm_set1_epi8(x[i]);
		for(unsigned j=0;j<n_ymm;j++) {
			__m256i inp = _mm256_loadu_si256( (__m256i*)(matA+j*32) );
			__m256i l_i0 = tbl32_gf16_log(inp&mask_f);
			__m256i l_i1 = tbl32_gf16_log(_mm256_srli_epi16(inp,4)&mask_f);

			__m256i ab0 = tbl32_gf16_mul_log_log( l_i0 , m0 , mask_f );
			__m256i ab1 = tbl32_gf16_mul_log_log( l_i1 , m0 , mask_f )^tbl32_gf16_mul_log_log( l_i0 , m1 , mask_f );
			__m256i ab2 = tbl32_gf16_mul_log_log( l_i1 , m1 , mask_f );
			__m256i ab2x8 = tbl32_gf16_mul_0x8( ab2 );

			r[j] ^= ab0 ^ ab2x8 ^ _mm256_slli_epi16( ab1^ab2 , 4 );
		}
		matA += n_A_vec_byte;
	}

        if( 0 == (n_A_vec_byte & 31 ) ) {
		for(unsigned i=0;i<n_ymm;i++) _mm256_storeu_si256( (__m256i*)(c + i*32) , r[i] );
	} else {
		__m256i *rr = r;
		n_ymm--;
		while( n_ymm ) {
			_mm256_storeu_si256( (__m256i*)c , *rr );
			rr += 1;
			c += 32;
			n_ymm--;
		}
		_mm256_store_si256( (__m256i*)x0, *rr );
		for(unsigned i=0;i<(n_A_vec_byte&31);i++) { c[i] = x0[i]; }
	}
}
#endif




