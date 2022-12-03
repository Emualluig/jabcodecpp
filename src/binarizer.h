#ifndef _HEADER_BINARZIER_H_
#define _HEADER_BINARZIER_H_

#include "jabTypes.h"
#include <utility>

void balanceRGB(Bitmap& bitmap) noexcept;

constexpr std::pair<jab_double, jab_double> getAveVar(const jab_byte* rgb) noexcept;

void getMinMax(const jab_byte* rgb, jab_byte* min, jab_byte* mid, jab_byte* max, jab_int32* index_min, jab_int32* index_mid, jab_int32* index_max) noexcept;
constexpr std::tuple<jab_byte, jab_byte, jab_byte, jab_byte, jab_byte, jab_byte> getMinMidMax(const jab_byte* rgb) noexcept {
	jab_byte index_min = 0;
	jab_byte index_mid = 1;
	jab_byte index_max = 2;
	if (rgb[index_min] > rgb[index_max]) {
		std::swap(index_min, index_max);
	}
	if (rgb[index_min] > rgb[index_mid]) {
		std::swap(index_min, index_mid);
	}
	if (rgb[index_mid] > rgb[index_max]) {
		std::swap(index_mid, index_max);
	}
	jab_byte min = rgb[index_min];
	jab_byte mid = rgb[index_mid];
	jab_byte max = rgb[index_max];
	return { min, mid, max, index_min, index_mid, index_max };
}

jab_boolean binarizerRGB(const Bitmap& bitmap, std::array<Bitmap, 3>& rgb, const std::array<jab_float, 3>& blk_ths) noexcept;
std::array<Bitmap, 3> binarizerRGB(const Bitmap& bitmap) noexcept;

#endif // !_HEADER_BINARZIER_H_
