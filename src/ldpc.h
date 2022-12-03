#ifndef _HEADER_LDPC_H_
#define _HEADER_LDPC_H_

#include "jabTypes.h"

#define LPDC_METADATA_SEED 	38545
#define LPDC_MESSAGE_SEED 	785465

//#define LDPC_DEFAULT_WC		4	//default error correction level 3
//#define LDPC_DEFAULT_WR		9	//default error correction level 3

//static const jab_vector2d default_ecl = {4, 7};	//default (wc, wr) for LDPC, corresponding to ecc level 5.
//static const jab_vector2d default_ecl = {5, 6};	//This (wc, wr) could be used, if higher robustness is preferred to capacity.

std::vector<jab_char> encodeLDPC(jab_char* data, jab_int32* coderate_params, const jab_int32 size) noexcept; // This is the main function

jab_int32 decodeLDPChd(jab_byte* data, jab_int32 length, jab_int32 wc, jab_int32 wr) noexcept;
jab_int32 decodeLDPC(jab_float* enc, jab_int32 length, jab_int32 wc, jab_int32 wr, jab_byte* dec) noexcept;

#endif // !_HEADER_LDPC_H_
