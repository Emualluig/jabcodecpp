#ifndef _HEADER_SAMPLE_H_
#define _HEADER_SAMPLE_H_

#include "jabTypes.h"
#include "detectorTypes.h"

std::optional<Bitmap> sampleSymbol(const Bitmap& bitmap, const jab_perspective_transform& pt, const jab_vector2d side_size) noexcept;
constexpr std::optional<Bitmap> sampleCrossArea(const Bitmap& bitmap, const jab_perspective_transform& pt) noexcept;

#endif // !_HEADER_SAMPLE_H_
