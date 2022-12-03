#ifndef _HEADER_JAB_TYPES_H_
#define _HEADER_JAB_TYPES_H_

typedef unsigned char 		jab_byte;
typedef char 				jab_char;
typedef unsigned char 		jab_boolean;
typedef int 				jab_int32;
typedef unsigned int 		jab_uint32;
typedef short 				jab_int16;
typedef unsigned short 		jab_uint16;
typedef long long 			jab_int64;
typedef unsigned long long	jab_uint64;
typedef float				jab_float;
typedef double              jab_double;

#define MAX_SYMBOL_NUMBER       61
#define MAX_COLOR_NUMBER        256
#define MAX_SIZE_ENCODING_MODE  256
#define JAB_ENCODING_MODES      6
#define ENC_MAX                 1000000
#define NUMBER_OF_MASK_PATTERNS	8

#define DEFAULT_SYMBOL_NUMBER 			1
#define DEFAULT_MODULE_SIZE				12
#define DEFAULT_COLOR_NUMBER 			8
#define DEFAULT_MODULE_COLOR_MODE 		2
#define DEFAULT_ECC_LEVEL				3
#define DEFAULT_MASKING_REFERENCE 		7

#define DISTANCE_TO_BORDER      4
#define MAX_ALIGNMENT_NUMBER    9
constexpr const jab_uint32 COLOR_PALETTE_NUMBER = 4;

#define BITMAP_BITS_PER_PIXEL	32
#define BITMAP_BITS_PER_CHANNEL	8
#define BITMAP_CHANNEL_COUNT	4

#define	JAB_SUCCESS		1
#define	JAB_FAILURE		0

#define NORMAL_DECODE		0
#define COMPATIBLE_DECODE	1

#define VERSION2SIZE(x)		(x * 4 + 17)
#define SIZE2VERSION(x)		((x - 17) / 4)

/**
 * @brief 2-dimensional integer vector
*/
typedef struct {
	jab_int32	x;
	jab_int32	y;
}jab_vector2d;

/**
 * @brief 2-dimensional float vector
*/
typedef struct {
	jab_float	x;
	jab_float	y;
}jab_point;

#include <optional>
#include <vector>
#include <memory>
#include <array>
class Bitmap {
public:
	jab_int32 width;
	jab_int32 height;
	jab_int32 bits_per_pixel;
	jab_int32 bits_per_channel;
	jab_int32 channel_count;
	std::vector<jab_byte> pixel;
	constexpr explicit Bitmap(
		jab_int32 width, jab_int32 height, 
		jab_int32 bits_per_pixel, jab_int32 bits_per_channel, jab_int32 channel_count,
		jab_uint64 size, jab_byte initializer
	) noexcept : 
		width{ width }, height{ height },  bits_per_pixel{ bits_per_pixel }, 
		bits_per_channel{ bits_per_channel }, channel_count{ channel_count }, 
		pixel{ std::vector<jab_byte>(size, initializer) } {}
	constexpr explicit Bitmap(
		jab_int32 width, jab_int32 height,
		jab_int32 bits_per_pixel, jab_int32 bits_per_channel, jab_int32 channel_count,
		jab_uint64 size
	) noexcept :
		width{ width }, height{ height }, bits_per_pixel{ bits_per_pixel },
		bits_per_channel{ bits_per_channel }, channel_count{ channel_count },
		pixel{ std::vector<jab_byte>(size) } {}
};

enum ColorPaletteType {
	FULL = 8,
	CMYK = 4,
};

/**
 * @brief Symbol parameters
*/
struct Symbol {
	jab_int32		index;
	jab_vector2d	side_size;
	jab_int32		host;
	jab_int32		slaves[4];
	jab_int32 		wcwr[2];
	std::vector<jab_char> data;
	std::vector<jab_byte> data_map;
	std::vector<jab_char> metadata;
	std::unique_ptr<jab_byte[]> matrix;
};

/**
 * @brief Encode parameters
*/
struct Encode {
	jab_int32		color_number;
	jab_int32		symbol_number;
	jab_int32		module_size;
	jab_int32		master_symbol_width;
	jab_int32		master_symbol_height;
	std::unique_ptr<jab_byte[]> palette;				///< Palette holding used module colors in format RGB
	std::unique_ptr<jab_vector2d[]> symbol_versions;
	std::unique_ptr<jab_byte[]> symbol_ecc_levels;
	std::unique_ptr<jab_int32[]> symbol_positions;
	std::unique_ptr<Symbol[]> symbols;				///< Pointer to internal representation of JAB Code symbols
};

/**
 * @brief Decoded metadata
*/
typedef struct {
	jab_boolean default_mode;
	jab_byte Nc;
	jab_byte mask_type;
	jab_byte docked_position;
	jab_vector2d side_version;
	jab_vector2d ecl;
} jab_metadata;

/**
 * @brief Decoded symbol
*/
typedef struct {
	jab_int32 index;
	jab_int32 host_index;
	jab_int32 host_position;
	jab_vector2d side_size;
	jab_float module_size;
	std::array<jab_point, 4> pattern_positions;
	jab_metadata metadata;
	std::array<jab_metadata,4 > slave_metadata;
	std::unique_ptr<jab_byte[]> palette;
	std::vector<jab_char> data;
} jab_decoded_symbol;


#endif // !_HEADER_JAB_TYPES_H_
