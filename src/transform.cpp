#include "transform.h"
#include "reportError.h"

/**
 * @brief Calculate transformation matrix of square to quadrilateral
 * @param x0 the x coordinate of the 1st destination point
 * @param y0 the y coordinate of the 1st destination point
 * @param x1 the x coordinate of the 2nd destination point
 * @param y1 the y coordinate of the 2nd destination point
 * @param x2 the x coordinate of the 3rd destination point
 * @param y2 the y coordinate of the 3rd destination point
 * @param x3 the x coordinate of the 4th destination point
 * @param y3 the y coordinate of the 4th destination point
 * @return the transformation matrix
*/
static constexpr jab_perspective_transform square2Quad(
	const jab_float x0, const jab_float y0,
	const jab_float x1, const jab_float y1,
	const jab_float x2, const jab_float y2,
	const jab_float x3, const jab_float y3
) noexcept {
	jab_float dx3 = x0 - x1 + x2 - x3;
	jab_float dy3 = y0 - y1 + y2 - y3;
	if (dx3 == 0 && dy3 == 0) {
		return {
			.a11 = x1 - x0,
			.a12 = y1 - y0,
			.a13 = 0,
			.a21 = x2 - x1,
			.a22 = y2 - y1,
			.a23 = 0,
			.a31 = x0,
			.a32 = y0,
			.a33 = 1,
		};
	}
	else
	{
		jab_float dx1 = x1 - x2;
		jab_float dx2 = x3 - x2;
		jab_float dy1 = y1 - y2;
		jab_float dy2 = y3 - y2;
		jab_float denominator = dx1 * dy2 - dx2 * dy1;
		jab_float a13 = (dx3 * dy2 - dx2 * dy3) / denominator;
		jab_float a23 = (dx1 * dy3 - dx3 * dy1) / denominator;
		return {
			.a11 = x1 - x0 + a13 * x1,
			.a12 = y1 - y0 + a13 * y1,
			.a13 = a13,
			.a21 = x3 - x0 + a23 * x3,
			.a22 = y3 - y0 + a23 * y3,
			.a23 = a23,
			.a31 = x0,
			.a32 = y0,
			.a33 = 1,
		};
	}
}

/**
 * @brief Calculate transformation matrix of quadrilateral to square
 * @param x0 the x coordinate of the 1st source point
 * @param y0 the y coordinate of the 1st source point
 * @param x1 the x coordinate of the 2nd source point
 * @param y1 the y coordinate of the 2nd source point
 * @param x2 the x coordinate of the 3rd source point
 * @param y2 the y coordinate of the 3rd source point
 * @param x3 the x coordinate of the 4th source point
 * @param y3 the y coordinate of the 4th source point
 * @return the transformation matrix
*/
static constexpr jab_perspective_transform quad2Square(
	const jab_float x0, const jab_float y0,
	const jab_float x1, const jab_float y1,
	const jab_float x2, const jab_float y2,
	const jab_float x3, const jab_float y3
) noexcept {
	jab_perspective_transform s2q = square2Quad(x0, y0, x1, y1, x2, y2, x3, y3);
	//calculate the adjugate matrix of s2q
	return {
		.a11 = s2q.a22 * s2q.a33 - s2q.a23 * s2q.a32,
		.a12 = s2q.a13 * s2q.a32 - s2q.a12 * s2q.a33,
		.a13 = s2q.a12 * s2q.a23 - s2q.a13 * s2q.a22,
		.a21 = s2q.a23 * s2q.a31 - s2q.a21 * s2q.a33,
		.a22 = s2q.a11 * s2q.a33 - s2q.a13 * s2q.a31,
		.a23 = s2q.a13 * s2q.a21 - s2q.a11 * s2q.a23,
		.a31 = s2q.a21 * s2q.a32 - s2q.a22 * s2q.a31,
		.a32 = s2q.a12 * s2q.a31 - s2q.a11 * s2q.a32,
		.a33 = s2q.a11 * s2q.a22 - s2q.a12 * s2q.a21,
	};
}

/**
 * @brief Calculate matrix multiplication
 * @param m1 the multiplicand
 * @param m2 the multiplier
 * @return m1 x m2
*/
static constexpr jab_perspective_transform multiply(const jab_perspective_transform& m1, const jab_perspective_transform& m2) noexcept {
	return {
		.a11 = m1.a11 * m2.a11 + m1.a12 * m2.a21 + m1.a13 * m2.a31,
		.a12 = m1.a11 * m2.a12 + m1.a12 * m2.a22 + m1.a13 * m2.a32,
		.a13 = m1.a11 * m2.a13 + m1.a12 * m2.a23 + m1.a13 * m2.a33,
		.a21 = m1.a21 * m2.a11 + m1.a22 * m2.a21 + m1.a23 * m2.a31,
		.a22 = m1.a21 * m2.a12 + m1.a22 * m2.a22 + m1.a23 * m2.a32,
		.a23 = m1.a21 * m2.a13 + m1.a22 * m2.a23 + m1.a23 * m2.a33,
		.a31 = m1.a31 * m2.a11 + m1.a32 * m2.a21 + m1.a33 * m2.a31,
		.a32 = m1.a31 * m2.a12 + m1.a32 * m2.a22 + m1.a33 * m2.a32,
		.a33 = m1.a31 * m2.a13 + m1.a32 * m2.a23 + m1.a33 * m2.a33,
	};
}

/**
 * @brief Calculate transformation matrix of quadrilateral to quadrilateral
 * @param x0 the x coordinate of the 1st source point
 * @param y0 the y coordinate of the 1st source point
 * @param x1 the x coordinate of the 2nd source point
 * @param y1 the y coordinate of the 2nd source point
 * @param x2 the x coordinate of the 3rd source point
 * @param y2 the y coordinate of the 3rd source point
 * @param x3 the x coordinate of the 4th source point
 * @param y3 the y coordinate of the 4th source point
 * @param x0p the x coordinate of the 1st destination point
 * @param y0p the y coordinate of the 1st destination point
 * @param x1p the x coordinate of the 2nd destination point
 * @param y1p the y coordinate of the 2nd destination point
 * @param x2p the x coordinate of the 3rd destination point
 * @param y2p the y coordinate of the 3rd destination point
 * @param x3p the x coordinate of the 4th destination point
 * @param y3p the y coordinate of the 4th destination point
 * @return the transformation matrix
*/
constexpr jab_perspective_transform perspectiveTransform(
	const jab_float x0, const jab_float y0,
	const jab_float x1, const jab_float y1,
	const jab_float x2, const jab_float y2,
	const jab_float x3, const jab_float y3,
	const jab_float x0p, const jab_float y0p,
	const jab_float x1p, const jab_float y1p,
	const jab_float x2p, const jab_float y2p,
	const jab_float x3p, const jab_float y3p
) noexcept {
	return multiply(
			quad2Square(x0, y0, x1, y1, x2, y2, x3, y3), 
			square2Quad(x0p, y0p, x1p, y1p, x2p, y2p, x3p, y3p));
}

constexpr static jab_perspective_transform tmpTransform(
	const jab_float x1, 
	const jab_float x2, const jab_float y2,
	const jab_float y3,
	const jab_float x0p, const jab_float y0p,
	const jab_float x1p, const jab_float y1p,
	const jab_float x2p, const jab_float y2p,
	const jab_float x3p, const jab_float y3p
) noexcept {
	return multiply(
		quad2Square(3.5f, 3.5f, (jab_float)x1 - 3.5f, 3.5f, (jab_float)x2 - 3.5f, (jab_float)y2 - 3.5f, 3.5f, y3 - 3.5f),
		square2Quad(x0p, y0p, x1p, y1p, x2p, y2p, x3p, y3p));
}

/**
 * @brief Get perspetive transformation matrix
 * @param p0 the coordinate of the 1st finder/alignment pattern
 * @param p1 the coordinate of the 2nd finder/alignment pattern
 * @param p2 the coordinate of the 3rd finder/alignment pattern
 * @param p3 the coordinate of the 4th finder/alignment pattern
 * @param side_size the side size of the symbol
 * @return the transformation matrix
*/
jab_perspective_transform getPerspectiveTransform(
	const jab_point p0,
	const jab_point p1,
	const jab_point p2,
	const jab_point p3,
	const jab_vector2d side_size
) noexcept {
	return perspectiveTransform(3.5f, 3.5f,
		(jab_float)side_size.x - 3.5f, 3.5f,
		(jab_float)side_size.x - 3.5f, (jab_float)side_size.y - 3.5f,
		3.5f, (jab_float)side_size.y - 3.5f,
		p0.x, p0.y,
		p1.x, p1.y,
		p2.x, p2.y,
		p3.x, p3.y
	);
}

void warpPoints(const jab_perspective_transform& pt, std::vector<jab_point>& points) noexcept {
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