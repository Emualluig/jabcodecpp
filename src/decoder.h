#ifndef _HEADER_DECODER_H_
#define _HEADER_DECODER_H_

#include "jabTypes.h"
#include "decoderTypes.h"

void getNextMetadataModuleInMaster(jab_int32 matrix_height, jab_int32 matrix_width, jab_int32 next_module_count, jab_int32& x, jab_int32& y) noexcept;
jab_int32 decodeMaster(const Bitmap& matrix, jab_decoded_symbol& symbol) noexcept;
jab_int32 decodeSlave(Bitmap& matrix, jab_decoded_symbol& symbol) noexcept;
std::optional<std::vector<jab_byte>> decodeData(const std::vector<jab_char>& bits) noexcept;

#endif // !_HEADER_DECODER_H_
