#include "rainbow_keypair.h"
#include "rainbow_keypair_computation.h"
#include "rainbow_keypair_computation_simd.h"

#include "blas_comm.h"
#include "blas.h"
#include "rainbow_blas.h"
#include "rainbow_blas_simd.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "utils_malloc.h"

////////////////////////////////////////////////////////////////





void calculate_Q_from_F_simd( ext_cpk_t * Qs, const sk_t * Fs , const sk_t * Ts )
{
/// Layer 1
/*
    Q_pk.l1_F1s[i] = F_sk.l1_F1s[i]

    Q_pk.l1_F2s[i] = (F1* T1 + F2) + F1tr * t1
    Q_pk.l1_F5s[i] = UT( T1tr* (F1 * T1 + F2) )
*/
    unsigned char * t1 = (unsigned char *) adapted_alloc( 32 , _V1*_O1*32 );
    unsigned char * t2 = (unsigned char *) adapted_alloc( 32 , _V1*_O2*32 );
    unsigned char * t3 = (unsigned char *) adapted_alloc( 32 , _O1*_O2*32 );
    gfv_generate_multab( t1 , Ts->t1 , _V1*_O1 );
    gfv_generate_multab( t2 , Ts->t4 , _V1*_O2 );
    gfv_generate_multab( t3 , Ts->t3 , _O1*_O2 );

    memcpy( Qs->l1_Q1 , Fs->l1_F1 , _O1_BYTE * N_TRIANGLE_TERMS(_V1) );

    memcpy( Qs->l1_Q2 , Fs->l1_F2 , _O1_BYTE * _V1 * _O1 );
    batch_trimat_madd_multab( Qs->l1_Q2 , Fs->l1_F1 , t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );    /// F1*T1 + F2

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
    unsigned char * tempQ = (unsigned char *) adapted_alloc( 32 , size_tempQ + 32 );

    memset( tempQ , 0 , _O1_BYTE * _O1 * _O1 );   /// l1_Q5
    batch_matTr_madd_multab( tempQ , t1 , _V1, _V1_BYTE, _O1, Qs->l1_Q2, _O1, _O1_BYTE );  //// t1_tr*(F1*T1 + F2)
    UpperTrianglize( Qs->l1_Q5 , tempQ , _O1, _O1_BYTE );    /// UT( ... )   /// Q5

    batch_trimatTr_madd_multab( Qs->l1_Q2 , Fs->l1_F1 , t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );    /// Q2
/*
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2
    Q_pk.l1_F3s[i] =         F1_F1T_T2 + F2_T3
    Q_pk.l1_F6s[i] = T1tr* ( F1_F1T_T2 + F2_T3 ) + F2tr * t2
    Q_pk.l1_F9s[i] = UT( T2tr* ( F1_T2 + F2_T3 ) )
*/
    batch_trimat_madd_multab( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE , _O2, _O1_BYTE );    /// F1*T2
    batch_mat_madd_multab( Qs->l1_Q3 , Fs->l1_F2 , _V1, t3 , _O1, _O1_BYTE , _O2, _O1_BYTE );  /// F1_T2 + F2_T3

    memset( tempQ , 0 , _O1_BYTE * _O2 * _O2 );   /// l1_Q9
    batch_matTr_madd_multab( tempQ , t2 , _V1, _V1_BYTE, _O2, Qs->l1_Q3, _O2, _O1_BYTE );  /// T2tr * ( F1_T2 + F2_T3 )
    UpperTrianglize( Qs->l1_Q9 , tempQ , _O2 , _O1_BYTE );   /// Q9

    batch_trimatTr_madd_multab( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE, _O2, _O1_BYTE );  /// F1_F1T_T2 + F2_T3  /// Q3

    batch_bmatTr_madd_multab( Qs->l1_Q6 , Fs->l1_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE ); /// F2tr*T2
    batch_matTr_madd_multab( Qs->l1_Q6 , t1, _V1, _V1_BYTE, _O1, Qs->l1_Q3, _O2, _O1_BYTE ); /// Q6

/*
  # 2nd layer
    Q1 = F1
    Q2 = F1_F1T*T1 + F2
    Q5 = UT( T1tr( F1*T1 + F2 )  + F5 )
*/
    memcpy( Qs->l2_Q1 , Fs->l2_F1 , _O2_BYTE * N_TRIANGLE_TERMS(_V1) );

    memcpy( Qs->l2_Q2 , Fs->l2_F2 , _O2_BYTE * _V1 * _O1 );
    batch_trimat_madd_multab( Qs->l2_Q2 , Fs->l2_F1 , t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// F1*T1 + F2

    memcpy( Qs->l2_Q5 , Fs->l2_F5 , _O2_BYTE * N_TRIANGLE_TERMS(_O1) );
    memset( tempQ , 0 , _O2_BYTE * _O1 * _O1 );   /// l2_Q5
    batch_matTr_madd_multab( tempQ , t1 , _V1, _V1_BYTE, _O1, Qs->l2_Q2, _O1, _O2_BYTE );  //// t1_tr*(F1*T1 + F2)
    UpperTrianglize( Qs->l2_Q5 , tempQ , _O1, _O2_BYTE );    /// UT( ... )   /// Q5

    batch_trimatTr_madd_multab( Qs->l2_Q2 , Fs->l2_F1 , t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// Q2

/*
    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    F1_F1T_T2 + F2_T3 = F1_T2 + F2_T3 + F1tr * t2

    Q3 =        F1_F1T*T2 + F2*T3 + F3
    Q9 = UT( T2tr*( F1*T2 + F2*T3 + F3 )  +      T3tr*( F5*T3 + F6 ) )
    Q6 = T1tr*( F1_F1T*T2 + F2*T3 + F3 )  + F2Tr*T2 + F5_F5T*T3 + F6
*/
    memcpy( Qs->l2_Q3 , Fs->l2_F3 , _O2_BYTE * _V1 * _O2 );
    batch_trimat_madd_multab( Qs->l2_Q3 , Fs->l2_F1 , t2 , _V1, _V1_BYTE , _O2, _O2_BYTE );    /// F1*T2 + F3
    batch_mat_madd_multab( Qs->l2_Q3 , Fs->l2_F2 , _V1, t3 , _O1, _O1_BYTE , _O2, _O2_BYTE );  /// F1_T2 + F2_T3 + F3

    memset( tempQ , 0 , _O2_BYTE * _O2 * _O2 );   /// l2_Q9
    batch_matTr_madd_multab( tempQ , t2 , _V1, _V1_BYTE, _O2, Qs->l2_Q3, _O2, _O2_BYTE );  /// T2tr * ( ..... )

    memcpy( Qs->l2_Q6 , Fs->l2_F6 , _O2_BYTE * _O1 *_O2 );

    batch_trimat_madd_multab( Qs->l2_Q6 , Fs->l2_F5 , t3 , _O1, _O1_BYTE, _O2, _O2_BYTE ); /// F5*T3 + F6
    batch_matTr_madd_multab( tempQ , t3 , _O1, _O1_BYTE, _O2, Qs->l2_Q6, _O2, _O2_BYTE );  /// T2tr*( ..... ) + T3tr*( ..... )
    memset( Qs->l2_Q9 , 0 , _O2_BYTE * N_TRIANGLE_TERMS(_O2) );
    UpperTrianglize( Qs->l2_Q9 , tempQ , _O2 , _O2_BYTE );   /// Q9

    batch_trimatTr_madd_multab( Qs->l2_Q3 , Fs->l2_F1 , t2 , _V1, _V1_BYTE, _O2, _O2_BYTE );  /// F1_F1T_T2 + F2_T3 + F3 /// Q3

    batch_bmatTr_madd_multab( Qs->l2_Q6 , Fs->l2_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O2_BYTE ); ///  F5*T3 + F6 +  F2tr*T2
    batch_trimatTr_madd_multab( Qs->l2_Q6 , Fs->l2_F5 , t3 , _O1, _O1_BYTE, _O2, _O2_BYTE );  ///   F2tr*T2 + F5_F5T*T3 + F6
    batch_matTr_madd_multab( Qs->l2_Q6 , t1, _V1, _V1_BYTE, _O1, Qs->l2_Q3, _O2, _O2_BYTE ); /// Q6

    memset( tempQ , 0 , size_tempQ + 32 );
    free( tempQ );
    memset( t1 , 0 , _V1*_O1*32 );
    memset( t2 , 0 , _V1*_O2*32 );
    memset( t3 , 0 , _O1*_O2*32 );
    free( t1 );
    free( t2 );
    free( t3 );
}





/////////////////////////////////////////////////////


void calculate_F_from_Q_simd( sk_t * Fs , const sk_t * Qs , sk_t * Ts )
{
    unsigned char * t1 = (unsigned char *) adapted_alloc( 32 , _V1*_O1*32 );
    unsigned char * t4 = (unsigned char *) adapted_alloc( 32 , _V1*_O2*32 );
    unsigned char * t3 = (unsigned char *) adapted_alloc( 32 , _O1*_O2*32 );
    gfv_generate_multab( t1 , Ts->t1 , _V1*_O1 );
    gfv_generate_multab( t4 , Ts->t4 , _V1*_O2 );
    gfv_generate_multab( t3 , Ts->t3 , _O1*_O2 );


    /// F_sk.l1_F1s[i] = Q_pk.l1_F1s[i]
    memcpy( Fs->l1_F1 , Qs->l1_F1 , _O1_BYTE * N_TRIANGLE_TERMS(_V1) );

    /// F_sk.l1_F2s[i] = ( Q_pk.l1_F1s[i] + Q_pk.l1_F1s[i].transpose() ) * T_sk.t1 + Q_pk.l1_F2s[i]
    memcpy( Fs->l1_F2 , Qs->l1_F2 , _O1_BYTE * _V1*_O1 );
    batch_2trimat_madd_multab( Fs->l1_F2 , Qs->l1_F1 , t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );
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
    batch_trimat_madd_multab( Fs->l2_F2 , Qs->l2_F1 , t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// Q1_T1+ Q2

    unsigned char * tempQ = (unsigned char *) adapted_alloc( 32 , _O1*_O1*_O2_BYTE + 32 );
    memset( tempQ , 0 , _O1*_O1*_O2_BYTE );
    //memset( tempQ.l2_F2 , 0 , _V1*_O1 * _O2_BYTE   );  /// the size should be _O1*_O1*_O2_BYTE for the next line. presuming _V1 >= _O1
    batch_matTr_madd_multab( tempQ , t1 , _V1, _V1_BYTE, _O1, Fs->l2_F2, _O1, _O2_BYTE );  //// t1_tr*(Q1_T1+Q2)
    //gf256v_add( Fs->l2_F5, Qs->l2_F5, _O2_BYTE * N_TRIANGLE_TERMS(_O1) );   /// F5
    memcpy( Fs->l2_F5, Qs->l2_F5, _O2_BYTE * N_TRIANGLE_TERMS(_O1) );   /// F5
    UpperTrianglize( Fs->l2_F5 , tempQ , _O1, _O2_BYTE );    /// UT( ... )
    memset( tempQ , 0 , _O1*_O1*_O2_BYTE + 32);
    free( tempQ );

    batch_trimatTr_madd_multab( Fs->l2_F2 , Qs->l2_F1 , t1 , _V1, _V1_BYTE , _O1, _O2_BYTE );    /// F2 = Q1_T1 + Q2 + Q1^tr*t1

    /// Q1_Q1T_T4 =  (Q_pk.l2_F1s[i] + Q_pk.l2_F1s[i].transpose()) * t4
    /// Q2_T3 = Q_pk.l2_F2s[i]*T_sk.t3
    /// F_sk.l2_F3s[i] =           Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i]
    memcpy( Fs->l2_F3 , Qs->l2_F3 , _V1*_O2*_O2_BYTE );
    batch_2trimat_madd_multab( Fs->l2_F3 , Qs->l2_F1 , t4 , _V1, _V1_BYTE , _O2, _O2_BYTE );   /// Q1_Q1T_T4
    batch_mat_madd_multab( Fs->l2_F3 , Qs->l2_F2 , _V1, t3 , _O1, _O1_BYTE , _O2, _O2_BYTE );  /// Q2_T3

    /// F_sk.l2_F6s[i] = t1_tr * ( Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i] )
    ///                +  Q_pk.l2_F2s[i].transpose() * t4
    ///                + (Q_pk.l2_F5s[i] + Q_pk.l2_F5s[i].transpose())*T_sk.t3   + Q_pk.l2_F6s[i]

    memcpy( Fs->l2_F6 , Qs->l2_F6 , _O1*_O2*_O2_BYTE );
    batch_matTr_madd_multab( Fs->l2_F6 , t1 , _V1, _V1_BYTE, _O1, Fs->l2_F3, _O2, _O2_BYTE );  /// t1_tr * ( Q1_Q1T_T4 + Q2_T3 + Q_pk.l2_F3s[i] )
    batch_2trimat_madd_multab( Fs->l2_F6 , Qs->l2_F5 , t3, _O1, _O1_BYTE, _O2, _O2_BYTE );  /// (Q_pk.l2_F5s[i] + Q_pk.l2_F5s[i].transpose())*T_sk.t3
    batch_bmatTr_madd_multab( Fs->l2_F6 , Qs->l2_F2, _O1, t4, _V1, _V1_BYTE, _O2, _O2_BYTE );

    memset( t1 , 0 , _V1*_O1*32 );
    memset( t4 , 0 , _V1*_O2*32 );
    memset( t3 , 0 , _O1*_O2*32 );
    free( t1 );
    free( t4 );
    free( t3 );

}


//////////////////////////////////////////////////////////////////////////////////////////////////


void calculate_Q_from_F_cyclic_simd( cpk_t * Qs, const sk_t * Fs , const sk_t * Ts )
{
/// Layer 1 , Q5, Q3, Q6, Q9
/*
    Q_pk.l1_F5s[i] = UT( T1tr* (F1 * T1 + F2) )
*/
    unsigned char * t1 = (unsigned char *) adapted_alloc( 32 , _V1*_O1*32 );
    unsigned char * t2 = (unsigned char *) adapted_alloc( 32 , _V1*_O2*32 );
    unsigned char * t3 = (unsigned char *) adapted_alloc( 32 , _O1*_O2*32 );
    gfv_generate_multab( t1 , Ts->t1 , _V1*_O1 );
    gfv_generate_multab( t2 , Ts->t4 , _V1*_O2 );
    gfv_generate_multab( t3 , Ts->t3 , _O1*_O2 );

    sk_t tempQ;
    memcpy( tempQ.l1_F2 , Fs->l1_F2 , _O1_BYTE * _V1 * _O1 );
    batch_trimat_madd_multab( tempQ.l1_F2 , Fs->l1_F1 , t1 , _V1, _V1_BYTE , _O1, _O1_BYTE );    /// F1*T1 + F2
    memset( tempQ.l2_F1 , 0 , _O1_BYTE * _V1 * _O2 );
    batch_matTr_madd_multab( tempQ.l2_F1 , t1 , _V1, _V1_BYTE, _O1, tempQ.l1_F2, _O1, _O1_BYTE );  //// T1tr*(F1*T1 + F2)
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

    batch_trimat_madd_multab( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE , _O2, _O1_BYTE );    /// F1*T2
    batch_mat_madd_multab( Qs->l1_Q3 , Fs->l1_F2 , _V1, t3 , _O1, _O1_BYTE , _O2, _O1_BYTE );  /// F1_T2 + F2_T3

    memset( tempQ.l1_F2 , 0 , _O1_BYTE * _V1 * _O2 );   //// should be F3. assuming: _O1 >= _O2
    batch_matTr_madd_multab( tempQ.l1_F2 , t2 , _V1, _V1_BYTE, _O2, Qs->l1_Q3, _O2, _O1_BYTE );  /// T2tr * ( F1_T2 + F2_T3 )
    UpperTrianglize( Qs->l1_Q9 , tempQ.l1_F2 , _O2 , _O1_BYTE );   /// Q9

    batch_trimatTr_madd_multab( Qs->l1_Q3 , Fs->l1_F1 , t2 , _V1, _V1_BYTE, _O2, _O1_BYTE );  /// F1_F1T_T2 + F2_T3  /// Q3

    batch_bmatTr_madd_multab( Qs->l1_Q6 , Fs->l1_F2, _O1, t2, _V1, _V1_BYTE, _O2, _O1_BYTE ); /// F2tr*T2
    batch_matTr_madd_multab( Qs->l1_Q6 , t1, _V1, _V1_BYTE, _O1, Qs->l1_Q3, _O2, _O1_BYTE ); /// Q6
/*
  # 2nd layer , Q9

    F1_T2     = F1 * t2
    F2_T3     = F2 * t3
    Q9 = UT( T2tr*( F1*T2 + F2*T3 + F3 )  +  T3tr*( F5*T3 + F6 ) )
*/
    sk_t tempQ2;
    memcpy( tempQ2.l2_F3 , Fs->l2_F3 , _O2_BYTE * _V1 * _O2 );  /// F3 actually.
    batch_trimat_madd_multab( tempQ2.l2_F3 , Fs->l2_F1 , t2 , _V1, _V1_BYTE , _O2, _O2_BYTE );    /// F1*T2 + F3
    batch_mat_madd_multab( tempQ2.l2_F3 , Fs->l2_F2 , _V1, t3 , _O1, _O1_BYTE , _O2, _O2_BYTE );  /// F1_T2 + F2_T3 + F3

    memset( tempQ.l2_F3 , 0 , _O2_BYTE * _V1 * _O2 );
    batch_matTr_madd_multab( tempQ.l2_F3 , t2 , _V1, _V1_BYTE, _O2, tempQ2.l2_F3, _O2, _O2_BYTE );  /// T2tr * ( ..... )

    memcpy( tempQ.l2_F6 , Fs->l2_F6 , _O2_BYTE * _O1 *_O2 );
    batch_trimat_madd_multab( tempQ.l2_F6 , Fs->l2_F5 , t3 , _O1, _O1_BYTE, _O2, _O2_BYTE ); /// F5*T3 + F6

    batch_matTr_madd_multab( tempQ.l2_F3 , t3 , _O1, _O1_BYTE, _O2, tempQ.l2_F6, _O2, _O2_BYTE );  /// T2tr*( ..... ) + T3tr*( ..... )
    memset( Qs->l2_Q9 , 0 , _O2_BYTE * N_TRIANGLE_TERMS(_O2) );
    UpperTrianglize( Qs->l2_Q9 , tempQ.l2_F3 , _O2 , _O2_BYTE );   /// Q9

    memset( t1 , 0 , _V1*_O1*32 );
    memset( t2 , 0 , _V1*_O2*32 );
    memset( t3 , 0 , _O1*_O2*32 );
    free( t1 );
    free( t2 );
    free( t3 );
}



