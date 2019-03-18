#ifndef _RAINBOW_KEYPAIR_H_
#define _RAINBOW_KEYPAIR_H_


#include "rainbow_config.h"


#define N_TRIANGLE_TERMS(n_var) (n_var*(n_var+1)/2)


#ifdef  __cplusplus
extern  "C" {
#endif




typedef
struct rainbow_publickey {
    unsigned char pk[(_PUB_M_BYTE) * N_TRIANGLE_TERMS(_PUB_N)];
} pk_t;


typedef
struct rainbow_secretkey {
    unsigned char sk_seed[LEN_SKSEED];

    unsigned char s1[_O1_BYTE*_O2];
    unsigned char t1[_V1_BYTE*_O1];
    unsigned char t4[_V1_BYTE*_O2];
    unsigned char t3[_O1_BYTE*_O2];

    unsigned char l1_F1[_O1_BYTE * N_TRIANGLE_TERMS(_V1)];
    unsigned char l1_F2[_O1_BYTE * _V1*_O1];

    unsigned char l2_F1[_O2_BYTE * N_TRIANGLE_TERMS(_V1)];
    unsigned char l2_F2[_O2_BYTE * _V1*_O1];

    unsigned char l2_F3[_O2_BYTE * _V1*_O2];
    unsigned char l2_F5[_O2_BYTE * N_TRIANGLE_TERMS(_O1)];
    unsigned char l2_F6[_O2_BYTE * _O1*_O2];
} sk_t;




typedef
struct rainbow_publickey_cyclic {
    unsigned char pk_seed[LEN_PKSEED];

    unsigned char l1_Q3[_O1_BYTE * _V1*_O2];
    unsigned char l1_Q5[_O1_BYTE * N_TRIANGLE_TERMS(_O1)];
    unsigned char l1_Q6[_O1_BYTE * _O1*_O2];
    unsigned char l1_Q9[_O1_BYTE * N_TRIANGLE_TERMS(_O2)];

    unsigned char l2_Q9[_O2_BYTE * N_TRIANGLE_TERMS(_O2)];
} cpk_t;



typedef
struct rainbow_secretkey_cyclic {
    unsigned char pk_seed[LEN_PKSEED];
    unsigned char sk_seed[LEN_SKSEED];
} csk_t;



/////////////////////////////////////

void generate_keypair( pk_t * pk, sk_t* sk, const unsigned char *sk_seed );

void generate_keypair_cyclic( cpk_t * pk, sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed );

void generate_compact_keypair_cyclic( cpk_t * pk, csk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed );

////////////////////////////////////

void generate_secretkey( sk_t* sk, const unsigned char *sk_seed );

void generate_secretkey_cyclic( sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed );

////////////////////////////////////


/// converting formats of public keys : from cyclic version to classic key

void cpk_to_pk( pk_t * pk , const cpk_t * cpk );










#ifdef  __cplusplus
}
#endif

#endif
