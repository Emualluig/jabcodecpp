#ifndef _HEADER_DETECTOR_TYPES_H_
#define _HEADER_DETECTOR_TYPES_H_

#include "jabTypes.h"

#define TEST_MODE			0
#if TEST_MODE
jab_bitmap* test_mode_bitmap;
jab_int32	test_mode_color;
#endif

#define MAX_MODULES 		145	//the number of modules in side-version 32
#define MAX_SYMBOL_ROWS		3
#define MAX_SYMBOL_COLUMNS	3
#define MAX_FINDER_PATTERNS 500
#define PI 					3.14159265
#define CROSS_AREA_WIDTH	14	//the width of the area across the host and slave symbols

#define DIST(x1, y1, x2, y2) (jab_float)(sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)))

/**
 * @brief Detection modes
*/
typedef enum
{
	QUICK_DETECT = 0,
	NORMAL_DETECT,
	INTENSIVE_DETECT
}jab_detect_mode;

/**
 * @brief Finder pattern, alignment pattern
*/
typedef struct {
	jab_int32		type;
	jab_float		module_size;
	jab_point		center;			//coordinates of the center
	jab_int32		found_count;
	jab_int32 		direction;
}jab_finder_pattern, jab_alignment_pattern;

/**
 * @brief Perspective transform
*/
typedef struct {
	jab_float a11;
	jab_float a12;
	jab_float a13;
	jab_float a21;
	jab_float a22;
	jab_float a23;
	jab_float a31;
	jab_float a32;
	jab_float a33;
}jab_perspective_transform;

#endif // !_HEADER_DETECTOR_TYPES_H_
