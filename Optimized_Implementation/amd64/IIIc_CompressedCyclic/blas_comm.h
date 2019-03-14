#ifndef _BLAS_COMM_H_
#define _BLAS_COMM_H_

#include <stdint.h>

static inline uint8_t gf16v_get_ele(const uint8_t *a, unsigned i) {
    uint8_t r = a[i >> 1];
    uint8_t r0 = r&0xf;
    uint8_t r1 = r>>4;
    uint8_t m = (uint8_t)(-(i&1));
    return (r1&m)|((~m)&r0);
}

static inline uint8_t gf16v_set_ele(uint8_t *a, unsigned i, uint8_t v) {
    uint8_t m = 0xf ^ (-(i&1));   ///  1--> 0xf0 , 0--> 0x0f
    uint8_t ai_remaining = a[i>>1] & (~m);   /// erase
    a[i>>1] = ai_remaining | (m&(v<<4))|(m&v&0xf);  /// set
    return v;
}


static inline uint8_t gf256v_get_ele(const uint8_t *a, unsigned i) { return a[i]; }

static inline uint8_t gf256v_set_ele(uint8_t *a, unsigned i, uint8_t v) { a[i]=v; return v; }


#ifdef  __cplusplus
extern  "C" {
#endif


/////////////////////////////////////

void gf256v_set_zero(uint8_t *b, unsigned _num_byte);

unsigned gf256v_is_zero(const uint8_t *a, unsigned _num_byte);

///////////////// multiplications  ////////////////////////////////

/// polynomial multiplication
/// School boook
void gf256v_polymul(uint8_t *c, const uint8_t *a, const uint8_t *b, unsigned _num);

/// matrix-vector

void _gf16mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

void _gf256mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

/// matrix-matrix

void gf16mat_mul(uint8_t *c, const uint8_t *a, const uint8_t *b, unsigned len_vec);

void gf256mat_mul(uint8_t *c, const uint8_t *a, const uint8_t *b, unsigned len_vec);

/////////////////   algorithms:  gaussian elim  //////////////////

unsigned _gf16mat_gauss_elim(uint8_t *mat, unsigned h, unsigned w);

unsigned _gf16mat_solve_linear_eq(uint8_t *sol, const uint8_t *inp_mat, const uint8_t *c_terms, unsigned n);

unsigned _gf256mat_gauss_elim(uint8_t *mat, unsigned h, unsigned w);

unsigned _gf256mat_solve_linear_eq(uint8_t *sol, const uint8_t *inp_mat, const uint8_t *c_terms, unsigned n);

////////////////  rand for matrices   //////////////////////////

/// buffer has to be as large as the input matrix
unsigned gf16mat_inv(uint8_t *inv_a, const uint8_t *a, unsigned H, uint8_t *buffer);

unsigned gf256mat_inv(uint8_t *inv_a, const uint8_t *a, unsigned H, uint8_t *buffer);


#ifdef  __cplusplus
}
#endif

#endif
