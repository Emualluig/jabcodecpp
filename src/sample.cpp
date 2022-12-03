#include "sample.h"

#include "reportError.h"
#include "transform.h"
#include "detector.h"
#include "decoder.h"

#include <array>

#define SAMPLE_AREA_WIDTH	(CROSS_AREA_WIDTH / 2 - 2) //width of the columns where the metadata and palette in slave symbol are located
#define SAMPLE_AREA_HEIGHT	20	//height of the metadata rows including the first row, though it does not contain metadata

/**
 * @brief Sample a symbol
 * @param bitmap the image bitmap
 * @param pt the transformation matrix
 * @param side_size the symbol size in module
 * @return the sampled symbol matrix
*/
std::optional<Bitmap> sampleSymbol(const Bitmap& bitmap, const jab_perspective_transform& pt, const jab_vector2d side_size) noexcept {
	jab_int32 mtx_bytes_per_pixel = bitmap.bits_per_pixel / 8;
	jab_int32 mtx_bytes_per_row = side_size.x * mtx_bytes_per_pixel;

	Bitmap matrix = Bitmap(
		side_size.x, side_size.y, 
		bitmap.bits_per_channel * bitmap.channel_count, 
		bitmap.bits_per_channel, bitmap.channel_count,
		side_size.x * side_size.y * mtx_bytes_per_pixel
	);

	jab_int32 bmp_bytes_per_pixel = bitmap.bits_per_pixel / 8;
	jab_int32 bmp_bytes_per_row = bitmap.width * bmp_bytes_per_pixel;

	std::vector<jab_point> points = std::vector<jab_point>(side_size.x);
    for(jab_int32 i=0; i < side_size.y; i++)
    {
		for(jab_int32 j=0; j < side_size.x; j++)
		{
			points[j] = {
				.x = (jab_float)j + 0.5f,
				.y = (jab_float)i + 0.5f,
			};
		}
		warpPoints(pt, points);
		for(jab_int32 j = 0; j < side_size.x; j++) {
			jab_int32 mapped_x = (jab_int32)points[j].x;
			jab_int32 mapped_y = (jab_int32)points[j].y;
			if(mapped_x < 0 || mapped_x > bitmap.width - 1) {
				if (mapped_x == -1) { 
					mapped_x = 0; 
				}
				else if (mapped_x == bitmap.width) { 
					mapped_x = bitmap.width - 1; 
				}
				else { 
					return std::nullopt; 
				}
			}
			if(mapped_y < 0 || mapped_y > bitmap.height-1) {
				if (mapped_y == -1) { 
					mapped_y = 0; 
				}
				else if (mapped_y == bitmap.height) { 
					mapped_y = bitmap.height - 1; 
				}
				else {
					return std::nullopt;
				}
			}
			for(jab_int32 c=0; c< matrix.channel_count; c++)
			{
				//get the average of pixel values in 3x3 neighborhood as the sampled value
				jab_float sum = 0;
				for(jab_int32 dx=-1; dx<=1; dx++)
				{
					for(jab_int32 dy=-1; dy<=1; dy++)
					{
						jab_int32 px = mapped_x + dx;
						jab_int32 py = mapped_y + dy;
						if(px < 0 || px > bitmap.width - 1)  px = mapped_x;
						if(py < 0 || py > bitmap.height - 1) py = mapped_y;
						sum += bitmap.pixel[py*bmp_bytes_per_row + px*bmp_bytes_per_pixel + c];
					}
				}
				jab_byte ave = (jab_byte)(sum / 9.0f + 0.5f);
				matrix.pixel[i*mtx_bytes_per_row + j*mtx_bytes_per_pixel + c] = ave;
			}
		}
    }
	return matrix;
}

/**
 * @brief Sample a cross area between the host and slave symbols
 * @param bitmap the image bitmap
 * @param pt the transformation matrix
 * @return the sampled area matrix
*/
constexpr std::optional<Bitmap> sampleCrossArea(const Bitmap& bitmap, const jab_perspective_transform& pt) noexcept {
	jab_int32 mtx_bytes_per_pixel = bitmap.bits_per_pixel / 8;
	jab_int32 mtx_bytes_per_row = SAMPLE_AREA_WIDTH * mtx_bytes_per_pixel;

	Bitmap matrix = Bitmap(
		SAMPLE_AREA_WIDTH, SAMPLE_AREA_HEIGHT, bitmap.bits_per_channel * bitmap.channel_count, 
		bitmap.bits_per_channel, bitmap.channel_count,
		SAMPLE_AREA_WIDTH * SAMPLE_AREA_HEIGHT * mtx_bytes_per_pixel
	);

	jab_int32 bmp_bytes_per_pixel = bitmap.bits_per_pixel / 8;
	jab_int32 bmp_bytes_per_row = bitmap.width * bmp_bytes_per_pixel;

	//only sample the area where the metadata and palette are located
	std::array<jab_point, SAMPLE_AREA_WIDTH> points;
    for(jab_int32 i=0; i<SAMPLE_AREA_HEIGHT; i++)
    {
		for(jab_int32 j=0; j<SAMPLE_AREA_WIDTH; j++)
		{
			points[j] = {
				.x = (jab_float)j + CROSS_AREA_WIDTH / 2 + 0.5f,
				.y = (jab_float)i + 0.5f,
			};
		}
		warpPoints<points.size()>(pt, points);
		for(jab_int32 j=0; j<SAMPLE_AREA_WIDTH; j++)
		{
			jab_int32 mapped_x = (jab_int32)points[j].x;
			jab_int32 mapped_y = (jab_int32)points[j].y;
			if(mapped_x < 0 || mapped_x > bitmap.width-1)
			{
				if (mapped_x == -1) { 
					mapped_x = 0; 
				}
				else if (mapped_x == bitmap.width) {
					mapped_x = bitmap.width - 1;
				}
				else {
					return std::nullopt;
				}
			}
			if(mapped_y < 0 || mapped_y > bitmap.height-1)
			{
				if (mapped_y == -1) {
					mapped_y = 0;
				}
				else if (mapped_y == bitmap.height) {
					mapped_y = bitmap.height - 1;
				}
				else {
					return std::nullopt;
				}
			}
			for(jab_int32 c=0; c < matrix.channel_count; c++)
			{
				//get the average of pixel values in 3x3 neighborhood as the sampled value
				jab_float sum = 0;
				for(jab_int32 dx=-1; dx<=1; dx++)
				{
					for(jab_int32 dy=-1; dy<=1; dy++)
					{
						jab_int32 px = mapped_x + dx;
						jab_int32 py = mapped_y + dy;
						if (px < 0 || px > bitmap.width - 1) { px = mapped_x; }
						if (py < 0 || py > bitmap.height - 1) { py = mapped_y; }
						sum += bitmap.pixel[py*bmp_bytes_per_row + px*bmp_bytes_per_pixel + c];
					}
				}
				jab_byte ave = (jab_byte)(sum / 9.0f + 0.5f);
				matrix.pixel[i*mtx_bytes_per_row + j*mtx_bytes_per_pixel + c] = ave;
			}
		}
    }
	return matrix;
}

