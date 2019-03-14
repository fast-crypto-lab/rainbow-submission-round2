
#ifndef _RAINBOW_H_
#define _RAINBOW_H_

#include "rainbow_config.h"
#include "rainbow_keypair.h"

#include <stdint.h>

#ifdef  __cplusplus
extern  "C" {
#endif


int rainbow_sign( uint8_t * signature , const sk_t * sk , const uint8_t * digest );

int rainbow_verify( const uint8_t * digest , const uint8_t * signature , const pk_t * pk );



int rainbow_sign_cyclic( uint8_t * signature , const csk_t * sk , const uint8_t * digest );

int rainbow_verify_cyclic( const uint8_t * digest , const uint8_t * signature , const cpk_t * pk );



#ifdef  __cplusplus
}
#endif


#endif
