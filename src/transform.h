#ifndef _HEADER_TRANSFORM_H_
#define _HEADER_TRANSFORM_H_

#include "jabTypes.h"
#include "detectorTypes.h"

#include <vector>

constexpr jab_perspective_transform perspectiveTransform(
	const jab_float x0, const jab_float y0,
	const jab_float x1, const jab_float y1,
	const jab_float x2, const jab_float y2,
	const jab_float x3, const jab_float y3,
	const jab_float x0p, const jab_float y0p,
	const jab_float x1p, const jab_float y1p,
	const jab_float x2p, const jab_float y2p,
	const jab_float x3p, const jab_float y3p
) noexcept;
jab_perspective_transform getPerspectiveTransform(
	const jab_point p0,
	const jab_point p1,
	const jab_point p2,
	const jab_point p3,
	const jab_vector2d side_size
) noexcept;

/**
 * @brief Warp points from source image to destination image in place
 * @param pt the transformation matrix
 * @param points the source points
 * @param length the number of source points
*/
void warpPoints(const jab_perspective_transform& pt, std::vector<jab_point>& points) noexcept;
template<std::size_t SIZE> constexpr void warpPoints(const jab_perspective_transform& pt, std::array<jab_point, SIZE>& points) noexcept {
	for (jab_point& element : points) {
		jab_float x = element.x;
		jab_float y = element.y;
		jab_float denominator = pt.a13 * x + pt.a23 * y + pt.a33;
		element = {
			.x = (pt.a11 * x + pt.a21 * y + pt.a31) / denominator,
			.y = (pt.a12 * x + pt.a22 * y + pt.a32) / denominator,
		};
	}
}

#endif // !_HEADER_TRANSFORM_H_
