#ifndef _HEADER_ENCODER_H_
#define _HEADER_ENCODER_H_

#include "jabTypes.h"
#include "encoderTypes.h"

Encode createEncode(ColorPaletteType color_number, jab_int32 symbol_number) noexcept;
std::optional<Bitmap> generateJABCode(Encode& enc, const std::vector<jab_char>& data, jab_int32& error) noexcept;

#endif // !_HEADER_ENCODER_H_
