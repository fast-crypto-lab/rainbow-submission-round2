
#ifndef _P_MATRIX_OP_AVX2_H_
#define _P_MATRIX_OP_AVX2_H_



#ifdef  __cplusplus
extern  "C" {
#endif





////////////////////  matrix multiplications  ///////////////////////////////



///  bC += btriA * B
void batch_trimat_madd_multab_gf16_avx2( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );

void batch_trimat_madd_multab_gf256_avx2( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );



///  bC += btriA^Tr * B
void batch_trimatTr_madd_multab_gf16_avx2( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );

void batch_trimatTr_madd_multab_gf256_avx2( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );



///  bC +=  (btriA + btriA^Tr) *B
void batch_2trimat_madd_multab_gf16_avx2( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );

void batch_2trimat_madd_multab_gf256_avx2( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );



/// bC += A^Tr * bB
void batch_matTr_madd_multab_gf16_avx2( unsigned char * bC ,
        const unsigned char* A_to_tr , unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
        const unsigned char* bB, unsigned Bwidth, unsigned size_batch );

void batch_matTr_madd_multab_gf256_avx2( unsigned char * bC ,
        const unsigned char* A_to_tr , unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
        const unsigned char* bB, unsigned Bwidth, unsigned size_batch );


///  bC += bA^Tr * B
void batch_bmatTr_madd_multab_gf16_avx2( unsigned char *bC , const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch );

void batch_bmatTr_madd_multab_gf256_avx2( unsigned char *bC , const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch );



///  bC += bA * B
void batch_mat_madd_multab_gf16_avx2( unsigned char * bC , const unsigned char* bA , unsigned Aheight,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );

void batch_mat_madd_multab_gf256_avx2( unsigned char * bC , const unsigned char* bA , unsigned Aheight,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch );





////////////////////  "quadratric" matrix evaluation  ///////////////////////////////


///  y =  x^Tr * trimat * x

void batch_quad_trimat_eval_gf16_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * x, unsigned dim , unsigned size_batch );

void batch_quad_trimat_eval_gf256_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * x, unsigned dim , unsigned size_batch );



///  y =  x^Tr * trimat * x

void batch_quad_trimat_eval_multab_gf16_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * multab_x, unsigned dim , unsigned size_batch );

void batch_quad_trimat_eval_multab_gf256_avx2( unsigned char * y, const unsigned char * trimat, const unsigned char * multab_x, unsigned dim , unsigned size_batch );







#ifdef  __cplusplus
}
#endif


#endif
