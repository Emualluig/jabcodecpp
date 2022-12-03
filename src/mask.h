#ifndef _HEADER_MASK_H_
#define _HEADER_MASK_H_

#include "jabTypes.h"
#include "encoderTypes.h"

constexpr void maskSymbols(Encode& enc, jab_int32 mask_type, std::vector<jab_int32>& masked, const CodeParams& cp) noexcept;
constexpr void maskSymbols(Encode& enc, jab_int32 mask_type) noexcept;

jab_int32 maskCode(Encode& enc, const CodeParams&) noexcept;
void demaskSymbol(std::vector<jab_char>& data, const std::vector<jab_byte>& data_map, jab_vector2d symbol_size, jab_int32 mask_type, jab_int32 color_number) noexcept;


#endif // !_HEADER_MASK_H_
