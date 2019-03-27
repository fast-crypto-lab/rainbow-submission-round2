///  @file rainbow-genkey.c
///  @brief A command-line tool for generating key pairs.
///

#include <stdio.h>
#include <stdint.h>

#include "rainbow_config.h"

#include "utils.h"

#include "rng.h"

#include "api.h"

int main( int argc , char ** argv )
{
	printf( "%s\n", CRYPTO_ALGNAME );

        printf("sk size: %lu\n", CRYPTO_SECRETKEYBYTES );
        printf("pk size: %lu\n",  CRYPTO_PUBLICKEYBYTES );
        printf("hash size: %d\n", _HASH_LEN );
        printf("signature size: %d\n\n", CRYPTO_BYTES );

	if( !((3 == argc) || (4 == argc)) ) {
		printf("Usage:\n\n\trainbow-genkey pk_file_name sk_file_name [random_seed_file]\n\n");
		return -1;
	}

	// set random seed
	unsigned char rnd_seed[48] = {0};
	int rr = byte_from_binfile( rnd_seed , 48 , (4==argc)? argv[3] : "/dev/random" );
	if( 0 != rr ) printf("read seed file fail.\n");
	randombytes_init( rnd_seed , NULL , 256 );


	uint8_t *_sk = (uint8_t*)malloc( CRYPTO_SECRETKEYBYTES );
	uint8_t *_pk = (uint8_t*)malloc( CRYPTO_PUBLICKEYBYTES );
	FILE * fp;

	int r = crypto_sign_keypair( _pk, _sk );
	if( 0 != r ) {
		printf("%s genkey fails.\n", CRYPTO_ALGNAME );
		return -1;
	}

	fp = fopen( argv[1] , "w+");
	if( NULL == fp ) {
		printf("fail to open public key file.\n");
		return -1;
	}
	byte_fdump( fp , CRYPTO_ALGNAME " public key" , _pk , CRYPTO_PUBLICKEYBYTES );
	fclose( fp );

	fp = fopen( argv[2] , "w+");
	if( NULL == fp ) {
		printf("fail to open secret key file.\n");
		return -1;
	}
	//ptr = (unsigned char *)&sk;
	//sprintf(msg,"%s secret key", name);
	byte_fdump( fp ,  CRYPTO_ALGNAME " secret key" , _sk , CRYPTO_SECRETKEYBYTES );
	fclose( fp );

	printf("generate %s pk/sk success.\n" , CRYPTO_ALGNAME );

	free( _sk );
	free( _pk );

	return 0;
}

