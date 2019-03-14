#ifndef _SHA256_UTILS_H_
#define _SHA256_UTILS_H_


/// for the definition of _HASH_LEN.
#include "hash_len_config.h"


#ifdef  __cplusplus
extern  "C" {
#endif


int hash_msg( unsigned char * digest , unsigned len_digest , const unsigned char * m , unsigned long long mlen );



#ifdef  __cplusplus
}
#endif



#endif

