#ifndef _RAINBOW_KEYPAIR_COMP_H_
#define _RAINBOW_KEYPAIR_COMP_H_


#include "rainbow_config.h"


#include "rainbow_keypair.h"


#ifdef  __cplusplus
extern  "C" {
#endif




typedef
struct rainbow_extend_publickey {
    unsigned char l1_Q1[_O1_BYTE * N_TRIANGLE_TERMS(_V1)];
    unsigned char l1_Q2[_O1_BYTE * _V1*_O1];
    unsigned char l1_Q3[_O1_BYTE * _V1*_O2];
    unsigned char l1_Q5[_O1_BYTE * N_TRIANGLE_TERMS(_O1)];
    unsigned char l1_Q6[_O1_BYTE * _O1*_O2];
    unsigned char l1_Q9[_O1_BYTE * N_TRIANGLE_TERMS(_O2)];

    unsigned char l2_Q1[_O2_BYTE * N_TRIANGLE_TERMS(_V1)];
    unsigned char l2_Q2[_O2_BYTE * _V1*_O1];
    unsigned char l2_Q3[_O2_BYTE * _V1*_O2];
    unsigned char l2_Q5[_O2_BYTE * N_TRIANGLE_TERMS(_O1)];
    unsigned char l2_Q6[_O2_BYTE * _O1*_O2];
    unsigned char l2_Q9[_O2_BYTE * N_TRIANGLE_TERMS(_O2)];
} ext_cpk_t;



void extcpk_to_pk( pk_t * pk , const ext_cpk_t * cpk );


/////////////////////////////////////////////////


void calculate_Q_from_F( ext_cpk_t * Qs, const sk_t * Fs , const sk_t * Ts );

void calculate_F_from_Q( sk_t * Fs , const sk_t * Qs , sk_t * Ts );

void calculate_Q_from_F_cyclic( cpk_t * Qs, const sk_t * Fs , const sk_t * Ts );






#ifdef  __cplusplus
}
#endif

#endif
