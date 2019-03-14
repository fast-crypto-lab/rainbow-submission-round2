#include "rainbow_keypair.h"

#include "blas_comm.h"
#include "blas.h"
#include "rainbow_blas.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>



#include "prng_utils.h"


static
void generate_S_T( unsigned char * s_and_t , prng_t * prng0 )
{
    prng_gen( prng0 , s_and_t , _O1_BYTE*_O2 ); // S1
    s_and_t += _O1_BYTE*_O2;
    prng_gen( prng0 , s_and_t , _V1_BYTE*_O1 ); // T1
    s_and_t += _V1_BYTE*_O1;
    prng_gen( prng0 , s_and_t , _V1_BYTE*_O2 ); // T2
    s_and_t += _V1_BYTE*_O2;
    prng_gen( prng0 , s_and_t , _O1_BYTE*_O2 ); // T3
}


static
unsigned generate_l1_F12( unsigned char * sk, prng_t * prng0 )
{
    unsigned n_byte_generated = 0;
    prng_gen( prng0 , sk , _O1_BYTE * N_TRIANGLE_TERMS(_V1) ); // l1_F1
    sk += _O1_BYTE * N_TRIANGLE_TERMS(_V1);
    n_byte_generated += _O1_BYTE * N_TRIANGLE_TERMS(_V1);

    prng_gen( prng0 , sk , _O1_BYTE * _V1*_O1 ); // l1_F2
    sk += _O1_BYTE * _V1*_O1;
    n_byte_generated += _O1_BYTE * _V1*_O1;

    return n_byte_generated;
}


static
unsigned generate_l2_F12356( unsigned char * sk, prng_t * prng0 )
{
    unsigned n_byte_generated = 0;

    prng_gen( prng0 , sk , _O2_BYTE * N_TRIANGLE_TERMS(_V1) ); // l2_F1
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_V1);
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_V1);

    prng_gen( prng0 , sk , _O2_BYTE * _V1*_O1 ); // l2_F2
    sk += _O2_BYTE * _V1*_O1;
    n_byte_generated += _O2_BYTE * _V1*_O1;

    prng_gen( prng0 , sk , _O2_BYTE * _V1*_O2 ); // l2_F3
    sk += _O2_BYTE * _V1*_O1;
    n_byte_generated += _O2_BYTE * _V1*_O1;

    prng_gen( prng0 , sk , _O2_BYTE * N_TRIANGLE_TERMS(_O1) ); // l2_F5
    sk += _O2_BYTE * N_TRIANGLE_TERMS(_O1);
    n_byte_generated += _O2_BYTE * N_TRIANGLE_TERMS(_O1);

    prng_gen( prng0 , sk , _O2_BYTE * _O1*_O2 ); // l2_F6
    n_byte_generated += _O2_BYTE * _O1*_O2;

    return n_byte_generated;
}


static
void generate_B1_B2( unsigned char * sk , prng_t * prng0 )
{
    sk += generate_l1_F12( sk , prng0 );
    generate_l2_F12356( sk , prng0 );
}


//////////////////////////////////////////////////////////


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


/////////////////////////////////////////////////

static
void extcpk_to_pk( pk_t * pk , const ext_cpk_t * cpk )
{
    const unsigned char * idx_l1 = cpk->l1_Q1;
    const unsigned char * idx_l2 = cpk->l2_Q1;
    for(unsigned i=0;i<_V1;i++) {
        for(unsigned j=i;j<_V1;j++) {
            unsigned pub_idx = idx_of_trimat(i,j,_PUB_N);
            memcpy( & pk->pk[ _PUB_M_BYTE*pub_idx ] , idx_l1 , _O1_BYTE );
            memcpy( (&pk->pk[ _PUB_M_BYTE*pub_idx ]) + _O1_BYTE , idx_l2 , _O2_BYTE );
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q2;
    idx_l2 = cpk->l2_Q2;
    for(unsigned i=0;i<_V1;i++) {
        for(unsigned j=_V1;j<_V1+_O1;j++) {
            unsigned pub_idx = idx_of_trimat(i,j,_PUB_N);
            memcpy( & pk->pk[ _PUB_M_BYTE*pub_idx ] , idx_l1 , _O1_BYTE );
            memcpy( (&pk->pk[ _PUB_M_BYTE*pub_idx ]) + _O1_BYTE , idx_l2 , _O2_BYTE );
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q3;
    idx_l2 = cpk->l2_Q3;
    for(unsigned i=0;i<_V1;i++) {
        for(unsigned j=_V1+_O1;j<_PUB_N;j++) {
            unsigned pub_idx = idx_of_trimat(i,j,_PUB_N);
            memcpy( & pk->pk[ _PUB_M_BYTE*pub_idx ]             , idx_l1 , _O1_BYTE );
            memcpy( (&pk->pk[ _PUB_M_BYTE*pub_idx ]) + _O1_BYTE , idx_l2 , _O2_BYTE );
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q5;
    idx_l2 = cpk->l2_Q5;
    for(unsigned i=_V1;i<_V1+_O1;i++) {
        for(unsigned j=i;j<_V1+_O1;j++) {
            unsigned pub_idx = idx_of_trimat(i,j,_PUB_N);
            memcpy( & pk->pk[ _PUB_M_BYTE*pub_idx ]             , idx_l1 , _O1_BYTE );
            memcpy( (&pk->pk[ _PUB_M_BYTE*pub_idx ]) + _O1_BYTE , idx_l2 , _O2_BYTE );
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q6;
    idx_l2 = cpk->l2_Q6;
    for(unsigned i=_V1;i<_V1+_O1;i++) {
        for(unsigned j=_V1+_O1;j<_PUB_N;j++) {
            unsigned pub_idx = idx_of_trimat(i,j,_PUB_N);
            memcpy( & pk->pk[ _PUB_M_BYTE*pub_idx ]             , idx_l1 , _O1_BYTE );
            memcpy( (&pk->pk[ _PUB_M_BYTE*pub_idx ]) + _O1_BYTE , idx_l2 , _O2_BYTE );
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
    idx_l1 = cpk->l1_Q9;
    idx_l2 = cpk->l2_Q9;
    for(unsigned i=_V1+_O1;i<_PUB_N;i++) {
        for(unsigned j=i;j<_PUB_N;j++) {
            unsigned pub_idx = idx_of_trimat(i,j,_PUB_N);
            memcpy( & pk->pk[ _PUB_M_BYTE*pub_idx ]             , idx_l1 , _O1_BYTE );
            memcpy( (&pk->pk[ _PUB_M_BYTE*pub_idx ]) + _O1_BYTE , idx_l2 , _O2_BYTE );
            idx_l1 += _O1_BYTE;
            idx_l2 += _O2_BYTE;
        }
    }
}


void cpk_to_pk( pk_t * rpk, const cpk_t * cpk )
{
    ext_cpk_t * pk = (ext_cpk_t *) aligned_alloc( 32, sizeof(ext_cpk_t)+32 );

    prng_t prng0;
    prng_set( &prng0 , cpk->pk_seed , LEN_SKSEED );

    generate_l1_F12( pk->l1_Q1 , &prng0 );
    memcpy( pk->l1_Q3 , cpk->l1_Q3 , _O1_BYTE*( _V1*_O2 + N_TRIANGLE_TERMS(_O1) + _O1*_O2 + N_TRIANGLE_TERMS(_O2) ) );

    generate_l2_F12356( pk->l2_Q1 , &prng0 );
    memcpy( pk->l2_Q9 , cpk->l2_Q9 , _O2_BYTE* N_TRIANGLE_TERMS(_O2) );

    extcpk_to_pk( rpk , pk );

    free( pk );
}




/////////////////////////////////////////////////////////




static
void calculate_t4( unsigned char * t2_to_t4 , const unsigned char *t1 , const unsigned char *t3 )
{
///  t4 = T_sk.t1 * T_sk.t3 - T_sk.t2
    unsigned char * t4 = t2_to_t4;
    for(unsigned i=0;i<_O2;i++) {  /// t3 width
        for(unsigned j=0;j<_O1;j++) { /// t3 height
            gfv_madd( t4 , &t1[j*_V1_BYTE] , gfv_get_ele( t3 , j ) , _V1_BYTE );
        }
        t4 += _V1_BYTE;
        t3 += _O1_BYTE;
    }
}



static
void obsfucate_l1_polys( unsigned char * l1_polys , const unsigned char * l2_polys , unsigned n_terms , const unsigned char * s1 )
{
    while( n_terms-- ) {
        for(unsigned i=0;i<_O2;i++) {
            gfv_madd( l1_polys , & s1[_O1_BYTE*i] , gfv_get_ele( l2_polys , i ) , _O1_BYTE );
        }
        l1_polys += _O1_BYTE;
        l2_polys += _O2_BYTE;
    }
}




/////////////////////////////////////////////////////



static
void calculate_Q_from_F( ext_cpk_t * Qs, const sk_t * Fs , const sk_t * Ts )
{
/// Layer 1
/*
    Q_pk.l1_F1s[i] = F_sk.l1_F1s[i]

    Q_pk.l1_F2s[i] = (F1* T1 + F2) + F1tr * t1
    Q_pk.l1_F5s[i] = UT( T1tr* (F1 * T1 + F2) )
*/
    const unsigned char * t2 = Ts->t4;

    memcpy( Qs->l1_Q1 , Fs->l1_F1 , _O1_BYTE * N_TRIANGLE_TERMS(_V1) );

    memcpy( Qs->l1_Q2 , Fs->l1_F2 , _O1_BYTE * _V1 * _O1 );
    batch_trimat_madd( Qs->l1_Q2 , Fs->l1_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );    /// F1*T1 + F2

    memset( Qs->l1_Q3 , 0 , _O1_BYTE * _V1 * _O2 );
    memset( Qs->l1_Q5 , 0 , _O1_BYTE * N_TRIANGLE_TERMS(_O1) );
    memset( Qs->l1_Q6 , 0 , _O1_BYTE * _O1 * _O2 );
    memset( Qs->l1_Q9 , 0 , _O1_BYTE * N_TRIANGLE_TERMS(_O2) );

    /// l1_Q5 : _O1_BYTE * _O1 * _O1
    /// l1_Q9 : _O1_BYTE * _O2 * _O2
    /// l2_Q5 : _O2_BYTE * _V1 * _O1
    /// l2_Q9 : _O2_BYTE * _V1 * _O2
    unsigned size_tempQ = _O1_BYTE * _O1 * _O1;
    if( _O1_BYTE*_O2*_O2 > size_tempQ ) size_tempQ = _O1_BYTE*_O2*_O2;
    if( _O2_BYTE*_O1*_O1 > size_tempQ ) size_tempQ = _O2_BYTE*_O1*_O1;
    if( _O2_BYTE*_O2*_O2 > size_tempQ ) size_tempQ = _O2_BYTE*_O2*_O2;
    unsigned char * tempQ = (unsigned char *) aligned_alloc( 32 , size_tempQ + 32 );

    memset( tempQ , 0 , _O1_BYTE * _O1 * _O1 );   /// l1_Q5
    batch_matTr_madd( tempQ , Ts->t1 , _V1, _V1_BYTE, _O1, Qs->l1_Q2, _O1, _O1_BYTE );  //// t1_tr*(F1*T1 + F2)
    UpperTrianglize( Qs->l1_Q5 , tempQ , _O1, _O1_BYTE );    /// UT( ... )   /// Q5

    batch_trimatTr_madd( Qs->l1_Q2 , Fs->l1_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );    /// Q2
/*
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2
    Q_pk.l1_F3s[i] =         F1_F1T_T2 + F2_T3
    Q_pk.l1_F6s[i] = T1tr* ( F1_F1T_T2 + F2_T3 ) + F2tr * t2
    Q_pk.l1_F9s[i] = UT( T2tr* ( F1_T2 + F2_T3 ) )
*/
    batch_trimat_madd( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE , _O2, _O1_BYTE );    /// F1*T2
    batch_mat_madd( Qs->l1_Q3 , Fs->l1_F2 , _V1, Ts->t3 , _O1, _O1_BYTE , _O2, _O1_BYTE );  /// F1_T2 + F2_T3

    memset( tempQ , 0 , _O1_BYTE * _O2 * _O2 );   /// l1_Q9
    batch_matTr_madd( tempQ , t2 , _V1, _V1_BYTE, _O2, Qs->l1_Q3, _O2, _O1_BYTE );  /// T2tr * ( F1_T2 + F2_T3 )
    UpperTrianglize( Qs->l1_Q9 , tempQ , _O2 , _O1_BYTE );   /// Q9

    batch_trimatTr_madd( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE, _O2, _O1_BYTE );  /// F1_F1T_T2 + F2_T3  /// Q3

    batch_bmatTr_madd( Qs->l1_Q6 , Fs->l1_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE ); /// F2tr*T2
    batch_matTr_madd( Qs->l1_Q6 , Ts->t1, _V1, _V1_BYTE, _O1, Qs->l1_Q3, _O2, _O1_BYTE ); /// Q6

/*
  # 2nd layer
    Q1 = F1
    Q2 = F1_F1T*T1 + F2
    Q5 = UT( T1tr( F1*T1 + F2 )  + F5 )
*/
    memcpy( Qs->l2_Q1 , Fs->l2_F1 , _O2_BYTE * N_TRIANGLE_TERMS(_V1) );

    memcpy( Qs->l2_Q2 , Fs->l2_F2 , _O2_BYTE * _V1 * _O1 );
    batch_trimat_madd( Qs->l2_Q2 , Fs->l2_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// F1*T1 + F2

    memcpy( Qs->l2_Q5 , Fs->l2_F5 , _O2_BYTE * N_TRIANGLE_TERMS(_O1) );
    memset( tempQ , 0 , _O2_BYTE * _O1 * _O1 );   /// l2_Q5
    batch_matTr_madd( tempQ , Ts->t1 , _V1, _V1_BYTE, _O1, Qs->l2_Q2, _O1, _O2_BYTE );  //// t1_tr*(F1*T1 + F2)
    UpperTrianglize( Qs->l2_Q5 , tempQ , _O1, _O2_BYTE );    /// UT( ... )   /// Q5

    batch_trimatTr_madd( Qs->l2_Q2 , Fs->l2_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// Q2

/*
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2

    Q3 =        F1_F1T*T2 + F2*T3 + F3
    Q9 = UT( T2tr*( F1*T2 + F2*T3 + F3 )  +      T3tr*( F5*T3 + F6 ) )
    Q6 = T1tr*( F1_F1T*T2 + F2*T3 + F3 )  + F2Tr*T2 + F5_F5T*T3 + F6
*/
    memcpy( Qs->l2_Q3 , Fs->l2_F3 , _O2_BYTE * _V1 * _O2 );
    batch_trimat_madd( Qs->l2_Q3 , Fs->l2_F1 , t2 , _V1, _V1_BYTE , _O2, _O2_BYTE );    /// F1*T2 + F3
    batch_mat_madd( Qs->l2_Q3 , Fs->l2_F2 , _V1, Ts->t3 , _O1, _O1_BYTE , _O2, _O2_BYTE );  /// F1_T2 + F2_T3 + F3

    memset( tempQ , 0 , _O2_BYTE * _O2 * _O2 );   /// l2_Q9
    batch_matTr_madd( tempQ , t2 , _V1, _V1_BYTE, _O2, Qs->l2_Q3, _O2, _O2_BYTE );  /// T2tr * ( ..... )

    memcpy( Qs->l2_Q6 , Fs->l2_F6 , _O2_BYTE * _O1 *_O2 );

    batch_trimat_madd( Qs->l2_Q6 , Fs->l2_F5 , Ts->t3 , _O1, _O1_BYTE, _O2, _O2_BYTE ); /// F5*T3 + F6
    batch_matTr_madd( tempQ , Ts->t3 , _O1, _O1_BYTE, _O2, Qs->l2_Q6, _O2, _O2_BYTE );  /// T2tr*( ..... ) + T3tr*( ..... )
    memset( Qs->l2_Q9 , 0 , _O2_BYTE * N_TRIANGLE_TERMS(_O2) );
    UpperTrianglize( Qs->l2_Q9 , tempQ , _O2 , _O2_BYTE );   /// Q9

    batch_trimatTr_madd( Qs->l2_Q3 , Fs->l2_F1 , t2 , _V1, _V1_BYTE, _O2, _O2_BYTE );  /// F1_F1T_T2 + F2_T3 + F3 /// Q3

    batch_bmatTr_madd( Qs->l2_Q6 , Fs->l2_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O2_BYTE ); ///  F5*T3 + F6 +  F2tr*T2
    batch_trimatTr_madd( Qs->l2_Q6 , Fs->l2_F5 , Ts->t3 , _O1, _O1_BYTE, _O2, _O2_BYTE );  ///   F2tr*T2 + F5_F5T*T3 + F6
    batch_matTr_madd( Qs->l2_Q6 , Ts->t1, _V1, _V1_BYTE, _O1, Qs->l2_Q3, _O2, _O2_BYTE ); /// Q6

    memset( tempQ , 0 , size_tempQ + 32 );
    free( tempQ );
}





/////////////////////////////////////////////////////


static
void _generate_secretkey( sk_t* sk, const unsigned char *sk_seed )
{
    memcpy( sk->sk_seed , sk_seed , LEN_SKSEED );

    prng_t prng0;
    prng_set( &prng0 , sk_seed , LEN_SKSEED );

    generate_S_T( sk->s1 , &prng0 );

    generate_B1_B2( sk->l1_F1 , &prng0 );

    memset( &prng0 , 0 , sizeof(prng_t) );
}


void generate_secretkey( sk_t* sk, const unsigned char *sk_seed )
{
    _generate_secretkey( sk , sk_seed );
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );
}



void generate_keypair( pk_t * rpk, sk_t* sk, const unsigned char *sk_seed )
{
    _generate_secretkey( sk , sk_seed );
    ext_cpk_t * pk = (ext_cpk_t*) aligned_alloc( 32 , sizeof(ext_cpk_t) + 32 );
    calculate_Q_from_F( pk, sk , sk );
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );

    obsfucate_l1_polys( pk->l1_Q1 , pk->l2_Q1 , N_TRIANGLE_TERMS(_V1) , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q2 , pk->l2_Q2 , _V1*_O1 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q3 , pk->l2_Q3 , _V1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q5 , pk->l2_Q5 , N_TRIANGLE_TERMS(_O1) , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q6 , pk->l2_Q6 , _O1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q9 , pk->l2_Q9 , N_TRIANGLE_TERMS(_O2) , sk->s1 );

    extcpk_to_pk( rpk , pk );
    free( pk );
}


//////////////////////////////////////////////////


static
void calculate_F_from_Q( sk_t * Fs , const sk_t * Qs , sk_t * Ts )
{
    /// F_sk.l1_F1s[i] = Q_pk.l1_F1s[i]
    memcpy( Fs->l1_F1 , Qs->l1_F1 , _O1_BYTE * N_TRIANGLE_TERMS(_V1) );

    /// F_sk.l1_F2s[i] = ( Q_pk.l1_F1s[i] + Q_pk.l1_F1s[i].transpose() ) * T_sk.t1 + Q_pk.l1_F2s[i]
    memcpy( Fs->l1_F2 , Qs->l1_F2 , _O1_BYTE * _V1*_O1 );
    batch_2trimat_madd( Fs->l1_F2 , Qs->l1_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );
    //batch_trimat_madd( Fs->l1_F2 , Qs->l1_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );
    //batch_trimatTr_madd( Fs->l1_F2 , Qs->l1_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );

/*
    computations for layer 2:

    F_sk.l2_F1s[i] = Q_pk.l2_F1s[i]

    Q1_T1 = Q_pk.l2_F1s[i]*T_sk.t1
    F_sk.l2_F2s[i] =              Q1_T1 + Q_pk.l2_F2s[i]     + Q_pk.l2_F1s[i].transpose() * T_sk.t1
    F_sk.l2_F5s[i] = UT( t1_tr* ( Q1_T1 + Q_pk.l2_F2s[i] ) ) + Q_pk.l2_F5s[i]

    Q1_Q1T_T4 =  (Q_pk.l2_F1s[i] + Q_pk.l2_F1s[i].transpose()) * t4
    #Q1_Q1T_T4 =  Q1_Q1T * t4
    Q2_T3 = Q_pk.l2_F2s[i]*T_sk.t3
    F_sk.l2_F3s[i] =           Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i]
    F_sk.l2_F6s[i] = t1_tr * ( Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i] )
                    +  Q_pk.l2_F2s[i].transpose() * t4
                    + (Q_pk.l2_F5s[i] + Q_pk.l2_F5s[i].transpose())*T_sk.t3   + Q_pk.l2_F6s[i]

*/
    /// F_sk.l2_F1s[i] = Q_pk.l2_F1s[i]
    memcpy( Fs->l2_F1 , Qs->l2_F1 , _O2_BYTE * N_TRIANGLE_TERMS(_V1) );

    /// F_sk.l2_F2s[i] =              Q1_T1 + Q_pk.l2_F2s[i]     + Q_pk.l2_F1s[i].transpose() * T_sk.t1
    /// F_sk.l2_F5s[i] = UT( t1_tr* ( Q1_T1 + Q_pk.l2_F2s[i] ) ) + Q_pk.l2_F5s[i]
    memcpy( Fs->l2_F2 , Qs->l2_F2 , _O2_BYTE * _V1*_O1 );
    batch_trimat_madd( Fs->l2_F2 , Qs->l2_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// Q1_T1+ Q2

    unsigned char * tempQ = (unsigned char *) aligned_alloc( 32 , _O1*_O1*_O2_BYTE + 32 );
    memset( tempQ , 0 , _O1*_O1*_O2_BYTE );
    //memset( tempQ.l2_F2 , 0 , _V1*_O1 * _O2_BYTE   );  /// the size should be _O1*_O1*_O2_BYTE for the next line. presuming _V1 >= _O1
    batch_matTr_madd( tempQ , Ts->t1 , _V1, _V1_BYTE, _O1, Fs->l2_F2, _O1, _O2_BYTE );  //// t1_tr*(Q1_T1+Q2)
    //gf256v_add( Fs->l2_F5, Qs->l2_F5, _O2_BYTE * N_TRIANGLE_TERMS(_O1) );   /// F5
    memcpy( Fs->l2_F5, Qs->l2_F5, _O2_BYTE * N_TRIANGLE_TERMS(_O1) );   /// F5
    UpperTrianglize( Fs->l2_F5 , tempQ , _O1, _O2_BYTE );    /// UT( ... )
    memset( tempQ , 0 , _O1*_O1*_O2_BYTE + 32);
    free( tempQ );

    batch_trimatTr_madd( Fs->l2_F2 , Qs->l2_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// F2 = Q1_T1 + Q2 + Q1^tr*t1

    /// Q1_Q1T_T4 =  (Q_pk.l2_F1s[i] + Q_pk.l2_F1s[i].transpose()) * t4
    /// Q2_T3 = Q_pk.l2_F2s[i]*T_sk.t3
    /// F_sk.l2_F3s[i] =           Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i]
    memcpy( Fs->l2_F3 , Qs->l2_F3 , _V1*_O2*_O2_BYTE );
    batch_2trimat_madd( Fs->l2_F3 , Qs->l2_F1 , Ts->t4 , _V1, _V1_BYTE , _O2, _O2_BYTE );   /// Q1_Q1T_T4
    batch_mat_madd( Fs->l2_F3 , Qs->l2_F2 , _V1, Ts->t3 , _O1, _O1_BYTE , _O2, _O2_BYTE );  /// Q2_T3

    /// F_sk.l2_F6s[i] = t1_tr * ( Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i] )
    ///                +  Q_pk.l2_F2s[i].transpose() * t4
    ///                + (Q_pk.l2_F5s[i] + Q_pk.l2_F5s[i].transpose())*T_sk.t3   + Q_pk.l2_F6s[i]

    memcpy( Fs->l2_F6 , Qs->l2_F6 , _O1*_O2*_O2_BYTE );
    batch_matTr_madd( Fs->l2_F6 , Ts->t1 , _V1, _V1_BYTE, _O1, Fs->l2_F3, _O2, _O2_BYTE );  /// t1_tr * ( Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i] )
    batch_2trimat_madd( Fs->l2_F6 , Qs->l2_F5 , Ts->t3, _O1, _O1_BYTE, _O2, _O2_BYTE );  /// (Q_pk.l2_F5s[i] + Q_pk.l2_F5s[i].transpose())*T_sk.t3
    batch_bmatTr_madd( Fs->l2_F6 , Qs->l2_F2, _O1, Ts->t4, _V1, _V1_BYTE, _O2, _O2_BYTE );

}


///////////////////////////////////////////////////////


void generate_secretkey_cyclic( sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed )
{
    memcpy( sk->sk_seed , sk_seed , LEN_SKSEED );
    //memcpy( sk->pk_seed , pk_seed , LEN_PKSEED );

    prng_t prng0;
    prng_set( &prng0 , sk_seed , LEN_SKSEED );
    generate_S_T( sk->s1 , &prng0 );
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );

    sk_t inst_Qs;
    sk_t * Qs = &inst_Qs;
    prng_t prng1;
    prng_set( &prng1 , pk_seed , LEN_PKSEED );
    generate_B1_B2( Qs->l1_F1 , &prng1 );

    obsfucate_l1_polys( Qs->l1_F1 , Qs->l2_F1 , N_TRIANGLE_TERMS(_V1) , sk->s1 );
    obsfucate_l1_polys( Qs->l1_F2 , Qs->l2_F2 , _V1*_O1 , sk->s1 );

    calculate_F_from_Q( sk , Qs , sk );
}


static
void calculate_Q_from_F_cyclic( cpk_t * Qs, const sk_t * Fs , const sk_t * Ts )
{
/// Layer 1 , Q5, Q3, Q6, Q9
/*
    Q_pk.l1_F5s[i] = UT( T1tr* (F1 * T1 + F2) )
*/
    const unsigned char * t2 = Ts->t4;

    sk_t tempQ;
    memcpy( tempQ.l1_F2 , Fs->l1_F2 , _O1_BYTE * _V1 * _O1 );
    batch_trimat_madd( tempQ.l1_F2 , Fs->l1_F1 , Ts->t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );    /// F1*T1 + F2
    memset( tempQ.l2_F1 , 0 , _O1_BYTE * _V1 * _O2 );
    batch_matTr_madd( tempQ.l2_F1 , Ts->t1 , _V1, _V1_BYTE, _O1, tempQ.l1_F2, _O1, _O1_BYTE );  //// T1tr*(F1*T1 + F2)
    memset( Qs->l1_Q5 , 0 , _O1_BYTE * N_TRIANGLE_TERMS(_O1) );
    UpperTrianglize( Qs->l1_Q5 , tempQ.l2_F1 , _O1, _O1_BYTE );    /// UT( ... )   /// Q5
/*
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2
    Q_pk.l1_F3s[i] =         F1_F1T_T2 + F2_T3
    Q_pk.l1_F6s[i] = T1tr* ( F1_F1T_T2 + F2_T3 ) + F2tr * t2
    Q_pk.l1_F9s[i] = UT( T2tr* ( F1_T2 + F2_T3 ) )
*/
    memset( Qs->l1_Q3 , 0 , _O1_BYTE * _V1 * _O2 );
    memset( Qs->l1_Q6 , 0 , _O1_BYTE * _O1 * _O2 );
    memset( Qs->l1_Q9 , 0 , _O1_BYTE * N_TRIANGLE_TERMS(_O2) );

    batch_trimat_madd( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE , _O2, _O1_BYTE );    /// F1*T2
    batch_mat_madd( Qs->l1_Q3 , Fs->l1_F2 , _V1, Ts->t3 , _O1, _O1_BYTE , _O2, _O1_BYTE );  /// F1_T2 + F2_T3

    memset( tempQ.l1_F2 , 0 , _O1_BYTE * _V1 * _O2 );   //// should be F3. assuming: _O1 >= _O2
    batch_matTr_madd( tempQ.l1_F2 , t2 , _V1, _V1_BYTE, _O2, Qs->l1_Q3, _O2, _O1_BYTE );  /// T2tr * ( F1_T2 + F2_T3 )
    UpperTrianglize( Qs->l1_Q9 , tempQ.l1_F2 , _O2 , _O1_BYTE );   /// Q9

    batch_trimatTr_madd( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE, _O2, _O1_BYTE );  /// F1_F1T_T2 + F2_T3  /// Q3

    batch_bmatTr_madd( Qs->l1_Q6 , Fs->l1_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE ); /// F2tr*T2
    batch_matTr_madd( Qs->l1_Q6 , Ts->t1, _V1, _V1_BYTE, _O1, Qs->l1_Q3, _O2, _O1_BYTE ); /// Q6
/*
  # 2nd layer , Q9

    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    Q9 = UT( T2tr*( F1*T2 + F2*T3 + F3 )  +  T3tr*( F5*T3 + F6 ) )
*/
    sk_t tempQ2;
    memcpy( tempQ2.l2_F3 , Fs->l2_F3 , _O2_BYTE * _V1 * _O2 );  /// F3 actually.
    batch_trimat_madd( tempQ2.l2_F3 , Fs->l2_F1 , t2 , _V1, _V1_BYTE , _O2, _O2_BYTE );    /// F1*T2 + F3
    batch_mat_madd( tempQ2.l2_F3 , Fs->l2_F2 , _V1, Ts->t3 , _O1, _O1_BYTE , _O2, _O2_BYTE );  /// F1_T2 + F2_T3 + F3

    memset( tempQ.l2_F3 , 0 , _O2_BYTE * _V1 * _O2 );
    batch_matTr_madd( tempQ.l2_F3 , t2 , _V1, _V1_BYTE, _O2, tempQ2.l2_F3, _O2, _O2_BYTE );  /// T2tr * ( ..... )

    memcpy( tempQ.l2_F6 , Fs->l2_F6 , _O2_BYTE * _O1 *_O2 );
    batch_trimat_madd( tempQ.l2_F6 , Fs->l2_F5 , Ts->t3 , _O1, _O1_BYTE, _O2, _O2_BYTE ); /// F5*T3 + F6

    batch_matTr_madd( tempQ.l2_F3 , Ts->t3 , _O1, _O1_BYTE, _O2, tempQ.l2_F6, _O2, _O2_BYTE );  /// T2tr*( ..... ) + T3tr*( ..... )
    memset( Qs->l2_Q9 , 0 , _O2_BYTE * N_TRIANGLE_TERMS(_O2) );
    UpperTrianglize( Qs->l2_Q9 , tempQ.l2_F3 , _O2 , _O2_BYTE );   /// Q9
}



void generate_keypair_cyclic( cpk_t * pk, sk_t* sk, const unsigned char *pk_seed , const unsigned char *sk_seed )
{
    memcpy( pk->pk_seed , pk_seed , LEN_PKSEED );

    prng_t prng;
    prng_t * prng0 = &prng;
    prng_set( prng0 , sk_seed , LEN_SKSEED );
    generate_S_T( sk->s1 , prng0 );

    //unsigned char t2[_V1_BYTE*_O2];
    unsigned char * t2 = (unsigned char *) aligned_alloc( 32, sizeof(sk->t4) );
    memcpy( t2 , sk->t4 , _V1_BYTE*_O2 );   /// save t2
    calculate_t4( sk->t4 , sk->t1 , sk->t3 );  /// t2 --> t4

    sk_t inst_Qs;
    sk_t * Qs = &inst_Qs;
    prng_t * prng1 = &prng;
    prng_set( prng1 , pk_seed , LEN_PKSEED );
    generate_B1_B2( Qs->l1_F1 , prng1 );
    obsfucate_l1_polys( Qs->l1_F1 , Qs->l2_F1 , N_TRIANGLE_TERMS(_V1) , sk->s1 );
    obsfucate_l1_polys( Qs->l1_F2 , Qs->l2_F2 , _V1*_O1 , sk->s1 );

    calculate_F_from_Q( sk , Qs , sk );  /// secret key

    memcpy( sk->t4 , t2 , _V1_BYTE*_O2 );   /// restore t2
    calculate_Q_from_F_cyclic( pk, sk , sk );  /// pubkey

    obsfucate_l1_polys( pk->l1_Q3 , Qs->l2_F3 , _V1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q5 , Qs->l2_F5 , N_TRIANGLE_TERMS(_O1) , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q6 , Qs->l2_F6 , _O1*_O2 , sk->s1 );
    obsfucate_l1_polys( pk->l1_Q9 , pk->l2_Q9 , N_TRIANGLE_TERMS(_O2) , sk->s1 );

    /// clean
    memset( &prng , 0 , sizeof(prng_t) );
    memset( t2 , 0 , sizeof(sk->t4) );
    free( t2 );
}



void generate_compact_keypair_cyclic( cpk_t * pk, csk_t* rsk, const unsigned char *pk_seed , const unsigned char *sk_seed )
{
    memcpy( rsk->pk_seed , pk_seed , LEN_PKSEED );
    memcpy( rsk->sk_seed , sk_seed , LEN_SKSEED );

    sk_t * sk = (sk_t *) aligned_alloc( 32 , sizeof(sk_t) + 32 );
    generate_keypair_cyclic( pk , sk , pk_seed , sk_seed );
    memset( sk , 0 , sizeof(sk_t) + 32 );
    free( sk );
}




