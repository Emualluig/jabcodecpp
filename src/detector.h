#ifndef _HEADER_DETECTOR_H_
#define _HEADER_DETECTOR_H_

#include "jabTypes.h"
#include "detectorTypes.h"

std::optional<std::vector<jab_byte>> decodeJABCode(Bitmap& bitmap, const jab_int32 mode, jab_int32& status) noexcept;

#endif // !_HEADER_DETECTOR_H_
