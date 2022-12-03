#ifndef _HEADER_IMAGE_H_
#define _HEADER_IMAGE_H_

#include "jabTypes.h"

jab_boolean saveImage(const Bitmap& bitmap, const char* filename) noexcept;
std::optional<Bitmap> readImage(const char* filename) noexcept;

std::pair<unsigned char*, unsigned int> saveImageMemory(const Bitmap& bitmap) noexcept;
std::optional<Bitmap> readImageMemory(const unsigned char* pointer, int length) noexcept;

#endif // !_HEADER_IMAGE_H_

