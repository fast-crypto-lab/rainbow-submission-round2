

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "rainbow_config.h"

#include "rainbow_keypair.h"

#include "rainbow.h"

#include "prng_utils.h"

#include "blas.h"

#include "rainbow_blas.h"


#define MAX_ATTEMPT_FRMAT  128

/////////////////////////////

#include "blas_sse.h"

#include "rainbow_blas_simd.h"


#include "hash_utils.h"


#define _MAX_O  ((_O1>_O2)?_O1:_O2)
#define _MAX_O_BYTE  ((_O1_BYTE>_O2_BYTE)?_O1_BYTE:_O2_BYTE)



/// algorithm 7
int rainbow_sign( uint8_t * signature , const sk_t * sk , const uint8_t * _digest )
{
    //const sk_t * sk = (const sk_t *)_sk;
//// line 1 - 5
    uint8_t * mat_l1 = aligned_alloc( 32 , _O1*_O1_BYTE + 32 );
    uint8_t * mat_l2 = aligned_alloc( 32 , _O2*_O2_BYTE + 32 );
    uint8_t * mat_buffer = aligned_alloc( 32 , 2*_MAX_O*_MAX_O_BYTE + 32 );
    //uint8_t mat_l1[_O1*_O1_BYTE];
    //uint8_t mat_l2[_O2*_O2_BYTE];
    //uint8_t mat_buffer[_MAX_O*_MAX_O_BYTE*2];

    /// setup PRNG
    prng_t prng_sign;
    uint8_t prng_preseed[LEN_SKSEED+_HASH_LEN];
    memcpy( prng_preseed , sk->sk_seed , LEN_SKSEED );
    memcpy( prng_preseed + LEN_SKSEED , _digest , _HASH_LEN );
    uint8_t prng_seed[_HASH_LEN];
    hash_msg( prng_seed , _HASH_LEN , prng_preseed , _HASH_LEN+LEN_SKSEED );
    prng_set( &prng_sign , prng_seed , _HASH_LEN );  /// seed = H( sk_seed || digest )
    for(unsigned i=0;i<LEN_SKSEED+_HASH_LEN;i++) prng_preseed[i] ^= prng_preseed[i];   /// clean
    for(unsigned i=0;i<_HASH_LEN;i++) prng_seed[i] ^= prng_seed[i];   /// clean

    /// roll vinegars.
    uint8_t vinegar[_V1_BYTE]  __attribute__((aligned(32)));
    unsigned n_attempt = 0;
    unsigned l1_succ = 0;
    while( !l1_succ ) {
        if( MAX_ATTEMPT_FRMAT <= n_attempt ) break;
        prng_gen( &prng_sign , vinegar , _V1_BYTE );
        gfmat_prod( mat_l1 , sk->l1_F2 , _O1*_O1_BYTE , _V1 , vinegar );
        l1_succ = gfmat_inv( mat_l1 , mat_l1 , _O1 , mat_buffer );
        n_attempt ++;
    }
    uint8_t multab[_V1*2*16] __attribute__((aligned(32)));
    gfv_generate_multab( multab , vinegar , _V1 );

    /// pre-compute variables needed for layer 2
    uint8_t r_l1_F1[_O1_BYTE]  __attribute__((aligned(32))) = {0};
    //batch_quad_trimat_eval( r_l1_F1, sk->l1_F1, vinegar, _V1, _O1_BYTE );
    batch_quad_trimat_eval_multab( r_l1_F1, sk->l1_F1, multab , _V1, _O1_BYTE );
    uint8_t r_l2_F1[_O2_BYTE]  __attribute__((aligned(32))) = {0};
    //batch_quad_trimat_eval( r_l2_F1, sk->l2_F1, vinegar, _V1, _O2_BYTE );
    batch_quad_trimat_eval_multab( r_l2_F1, sk->l2_F1, multab, _V1, _O2_BYTE );

    uint8_t * mat_l2_F3 = aligned_alloc( 32 , _O2*_O2_BYTE + 32 );
    //gfmat_prod( mat_l2_F3 , sk->l2_F3 , _O2*_O2_BYTE , _V1 , vinegar );
    gfmat_prod_multab( mat_l2_F3 , sk->l2_F3 , _O2*_O2_BYTE , _V1 , multab );
    uint8_t * mat_l2_F2 = aligned_alloc( 32 , _O1*_O2_BYTE + 32 );
    //gfmat_prod( mat_l2_F2 , sk->l2_F2 , _O1*_O2_BYTE , _V1 , vinegar );
    gfmat_prod_multab( mat_l2_F2 , sk->l2_F2 , _O1*_O2_BYTE , _V1 , multab );

    //// line 7 - 14
    uint8_t _z[_PUB_M_BYTE];
    uint8_t * x_v1 = vinegar;
    uint8_t x_o1[_O1_BYTE];
    uint8_t x_o2[_O1_BYTE];
    uint8_t y[_PUB_M_BYTE];
    uint8_t digest_salt[_HASH_LEN + _SALT_BYTE];
    memcpy( digest_salt , _digest , _HASH_LEN );
    uint8_t * salt = digest_salt + _HASH_LEN;

    uint8_t temp_o[_MAX_O_BYTE + 32]  = {0};
    unsigned succ = 0;
    while( !succ ) {
        if( MAX_ATTEMPT_FRMAT <= n_attempt ) break;
        /// roll salt
        prng_gen( &prng_sign , salt , _SALT_BYTE );
        hash_msg( _z , _PUB_M_BYTE , digest_salt , _HASH_LEN+_SALT_BYTE );

        ///  y = S^-1 * z
        memcpy(y, _z, _PUB_M_BYTE);  /// identity part of S
        gfmat_prod(temp_o, sk->s1, _O1_BYTE, _O2, _z+_O1_BYTE);
        gf256v_add(y, temp_o, _O1_BYTE);

        /// layer 1. calculate o1
        memcpy( temp_o , r_l1_F1 , _O1_BYTE );
        gf256v_add( temp_o , y , _O1_BYTE );
        gfmat_prod( x_o1 , mat_l1, _O1_BYTE , _O1 , temp_o );

        gfv_generate_multab( multab , x_o1 , _O1 );

        /// layer 2. calculate o2
        gf256v_set_zero( temp_o , _O2_BYTE );
        //gfmat_prod( temp_o , mat_l2_F2, _O2_BYTE , _O1 , x_o1 );  /// F2
        gfmat_prod_multab( temp_o , mat_l2_F2, _O2_BYTE , _O1 , multab );  /// F2
        //batch_quad_trimat_eval( mat_l2 , sk->l2_F5, x_o1 , _O1, _O2_BYTE );  /// F5
        batch_quad_trimat_eval_multab( mat_l2 , sk->l2_F5, multab , _O1, _O2_BYTE );  /// F5
        gf256v_add( temp_o , mat_l2 , _O2_BYTE );
        gf256v_add( temp_o , r_l2_F1 , _O2_BYTE );      /// F1
        gf256v_add( temp_o , y + _O1_BYTE , _O2_BYTE );

        /// generate inv_mat
        //gfmat_prod( mat_l2 , sk->l2_F6 , _O2*_O2_BYTE , _O1 , x_o1 );   /// F6
        gfmat_prod_multab( mat_l2 , sk->l2_F6 , _O2*_O2_BYTE , _O1 , multab );   /// F6
        gf256v_add( mat_l2 , mat_l2_F3 , _O2*_O2_BYTE);    /// F3
        succ = gfmat_inv( mat_l2 , mat_l2 , _O2 , mat_buffer );
        /// solve l2 eqs
        gfmat_prod( x_o2 , mat_l2 , _O2_BYTE , _O2 , temp_o );

        n_attempt ++;
    };
    uint8_t w[_PUB_N_BYTE];
    //memcpy(w, x , _PUB_N_BYTE );   /// identity part of matrix
    memcpy( w , x_v1 , _V1_BYTE );
    memcpy( w + _V1_BYTE , x_o1 , _O1_BYTE );
    memcpy( w + _V2_BYTE , x_o2 , _O2_BYTE );
    gfmat_prod(y, sk->t1, _V1_BYTE , _O1 , x_o1 );
    gf256v_add(w, y, _V1_BYTE );

    gfmat_prod(y, sk->t4, _V1_BYTE , _O2 , x_o2 );
    gf256v_add(w, y, _V1_BYTE );

    gfmat_prod(y, sk->t3, _O1_BYTE , _O2 , x_o2 );
    gf256v_add(w+_V1_BYTE, y, _O1_BYTE );

    memset( signature , 0 , _SIGNATURE_BYTE );

    /// clean
    memset( mat_l1 , 0 , _O1*_O1_BYTE + 32 );  free( mat_l1 );
    memset( mat_l2 , 0 , _O2*_O2_BYTE + 32 );  free( mat_l2 );
    memset( mat_buffer , 0 , 2*_MAX_O*_MAX_O_BYTE + 32 );  free( mat_buffer );
    memset( &prng_sign , 0 , sizeof(prng_t) );
    memset( vinegar , 0 , _V1_BYTE );
    memset( r_l1_F1 , 0 , _O1_BYTE );
    memset( r_l2_F1 , 0 , _O2_BYTE );
    memset( mat_l2_F3 , 0 , _O2*_O2_BYTE + 32 );  free( mat_l2_F3 );
    memset( mat_l2_F2 , 0 , _O1*_O2_BYTE + 32 );  free( mat_l2_F2 );
    memset( _z , 0 , _PUB_M_BYTE );
    memset( y , 0 , _PUB_M_BYTE );
    memset( x_o1 , 0 , _O1_BYTE );
    memset( x_o2 , 0 , _O2_BYTE );
    memset( temp_o , 0 , _MAX_O_BYTE +32 );

    // return time;
    if( MAX_ATTEMPT_FRMAT <= n_attempt ) return -1;
    gf256v_add( signature , w , _PUB_N_BYTE );
    gf256v_add( signature + _PUB_N_BYTE , salt , _SALT_BYTE );
    return 0;
}







/// algorithm 8
int rainbow_verify( const uint8_t * digest , const uint8_t * signature , const pk_t * pk )
{
    unsigned char digest_ck[_PUB_M_BYTE];
    //public_map( digest_ck , pk , signature );
    batch_quad_trimat_eval( digest_ck , pk->pk , signature , _PUB_N , _PUB_M_BYTE );

    unsigned char correct[_PUB_M_BYTE];
    unsigned char digest_salt[_HASH_LEN + _SALT_BYTE];
    memcpy( digest_salt , digest , _HASH_LEN );
    memcpy( digest_salt+_HASH_LEN , signature+_PUB_N_BYTE , _SALT_BYTE );
    hash_msg( correct , _PUB_M_BYTE , digest_salt , _HASH_LEN+_SALT_BYTE );

    unsigned char cc = 0;
    for(unsigned i=0;i<_PUB_M_BYTE;i++) {
        cc |= (digest_ck[i]^correct[i]);
    }
    return (0==cc)? 0: -1;
}



///////////////  cyclic version  ///////////////////////////


int rainbow_sign_cyclic( uint8_t * signature , const csk_t * csk , const uint8_t * digest )
{
    sk_t * sk = aligned_alloc( 32 , sizeof(sk_t) + 32 );
    if( NULL == sk ) return -1;
    generate_secretkey_cyclic( sk, csk->pk_seed , csk->sk_seed );

    int r = rainbow_sign( signature , sk , digest );
    free( sk );
    return r;
}

int rainbow_verify_cyclic( const uint8_t * digest , const uint8_t * signature , const cpk_t * _pk )
{
    pk_t * pk = aligned_alloc( 32 , sizeof(pk_t) + 32 );
    if( NULL == pk ) return -1;
    cpk_to_pk( pk , _pk );

    int r = rainbow_verify( digest , signature , pk );
    free( pk );
    return r;
}


