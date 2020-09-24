/// @file blas_matrix_vartime.h
/// @brief linear algebra functions for matrix op.
///
#ifndef _BLAS_MATRIX_VARTIME_H_
#define _BLAS_MATRIX_VARTIME_H_

#include <stdint.h>



#ifdef  __cplusplus
extern  "C" {
#endif


///////////////// Section: multiplications  ////////////////////////////////


/// @brief matrix-vector multiplication:  c = matA * b , in GF(16)
///
/// @param[out]  c         - the output vector c
/// @param[in]   matA      - a column-major matrix A.
/// @param[in]   n_A_vec_byte  - the size of column vectors in bytes.
/// @param[in]   n_A_width   - the width of matrix A.
/// @param[in]   b          - the vector b.
///
void gf16mat_prod_vartime(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);


/// @brief matrix-vector multiplication:  c = matA * b , in GF(256)
///
/// @param[out]  c         - the output vector c
/// @param[in]   matA      - a column-major matrix A.
/// @param[in]   n_A_vec_byte  - the size of column vectors in bytes.
/// @param[in]   n_A_width   - the width of matrix A.
/// @param[in]   b          - the vector b.
///
void gf256mat_prod_vartime(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);


#ifdef  __cplusplus
}
#endif

#endif  // _BLAS_MATRIX_VARTIME_H_

