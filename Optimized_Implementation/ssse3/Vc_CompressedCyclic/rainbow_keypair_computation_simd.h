#ifndef _RAINBOW_KEYPAIR_COMP_SIMD_H_
#define _RAINBOW_KEYPAIR_COMP_SIMD_H_


#include "rainbow_config.h"

#include "rainbow_keypair.h"

#include "rainbow_keypair_computation.h"

#ifdef  __cplusplus
extern  "C" {
#endif



void calculate_Q_from_F_simd( ext_cpk_t * Qs, const sk_t * Fs , const sk_t * Ts );

void calculate_F_from_Q_simd( sk_t * Fs , const sk_t * Qs , sk_t * Ts );

void calculate_Q_from_F_cyclic_simd( cpk_t * Qs, const sk_t * Fs , const sk_t * Ts );



#ifdef  __cplusplus
}
#endif

#endif
