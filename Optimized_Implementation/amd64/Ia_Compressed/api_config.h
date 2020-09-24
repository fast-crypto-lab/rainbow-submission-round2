/// @file api_config.h
/// @brief Defining which key format used for api.h/sign.c
///
/// Defining which key format used for api.h/sign.c
///

#ifndef _API_CONFIG_H_
#define _API_CONFIG_H_




#if (!defined(_RAINBOW_CLASSIC))&&(!defined(_RAINBOW_CIRCUMZENITHAL))&&(!defined(_RAINBOW_COMPRESSED))
#define _RAINBOW_COMPRESSED
//#define _RAINBOW_CIRCUMZENITHAL
//#define _RAINBOW_COMPRESSED
#endif


#if defined _RAINBOW_CLASSIC
#define _SUFFIX " - classic"
#elif defined _RAINBOW_CIRCUMZENITHAL
#define _SUFFIX " - circumzenithal"
#elif defined _RAINBOW_COMPRESSED
#define _SUFFIX " - compressed"
#else
error here
#endif



#endif  // _API_CONFIG_H_
