#include "binarizer.h"

#include "reportError.h"

#include <cmath>
#include <array>

#define BLOCK_SIZE_POWER	5
#define BLOCK_SIZE 			(1 << BLOCK_SIZE_POWER)
#define BLOCK_SIZE_MASK 	(BLOCK_SIZE - 1)
#define MINIMUM_DIMENSION 	(BLOCK_SIZE * 5)
#define CAP(val, min, max)	(val < min ? min : (val > max ? max : val))

/**
 * @brief Filter out noises in binary bitmap
 * @param binary the binarized bitmap
*/
constexpr static void filterBinary(Bitmap& binary) noexcept {
	const jab_int32 width = binary.width;
	const jab_int32 height= binary.height;

	constexpr jab_int32 filter_size = 5;
	constexpr jab_int32 half_size = (filter_size - 1)/2; // This has a value of 2

	for(jab_int32 i = half_size; i < height-half_size; i++) {
		for(jab_int32 j = half_size; j < width-half_size; j++) {
			jab_int32 hsum = 0;
			jab_int32 vsum = 0;
			// TODO: investigate data dependence, could move "binary.pixel[i * width + j] =" to outer loop and reuse some memory
			for (jab_int32 k = -half_size; k <= half_size; k++) {
				vsum += binary.pixel[(i + k) * width + j] != 0; // Replaced > 0 ? 1 : 0 with != 0 since the lowest value of jab_byte is 0
				hsum += binary.pixel[i * width + (j + k)] != 0;
			}
			binary.pixel[i * width + j] = std::max(hsum, vsum) > half_size ? 255 : 0;
		}
	}
	/*
	The resulting binary.pixel should look this this:
	xx............xx
	xx............xx
	................
	................
	................
	................
	xx............xx
	xx............xx
	Where . is either 255 or 0 and x is unmodified data
	*/
}

/**
 * @brief Get the histogram of a color channel
 * @param bitmap the image
 * @param channel the channel
 * @param hist the histogram
*/
/**
 * @brief Get the min and max index in the histogram whose value is larger than the threshold
 * @param hist the histogram
 * @param max the max index
 * @param min the min index
 * @param ths the threshold
*/
template <jab_byte channel>
constexpr static std::pair<jab_int32, jab_int32> histMaxMin(const Bitmap& bitmap) noexcept {
	// Get histogram of color channel
	std::array<uint32_t, 256> hist = {};
	jab_int32 bytes_per_pixel = bitmap.bits_per_pixel / 8;
	for (jab_int32 i = 0; i < bitmap.width * bitmap.height; i++)
	{
		hist[bitmap.pixel[i * bytes_per_pixel + channel]]++;
	}

	// Calculate min and max based on threshold
	jab_int32 min = 0;
	jab_int32 max = 255;

	//threshold for the number of pixels having the max or min values
	constexpr const int threshold = 20;
	//get min
	for (jab_int32 i = 0; i < 256; i++)
	{
		if (hist[i] > threshold)
		{
			min = i;
			break;
		}
	}

	//get max
	for (jab_int32 i = 255; i >= 0; i--)
	{
		if (hist[i] > threshold)
		{
			max = i;
			break;
		}
	}
	return { max, min };
}

/**
 * @brief Stretch the histograms of R, G and B channels
 * @param bitmap the image
*/
void balanceRGB(Bitmap& bitmap) noexcept {
	jab_int32 bytes_per_pixel = bitmap.bits_per_pixel / 8;
    jab_int32 bytes_per_row = bitmap.width * bytes_per_pixel;

	//calculate max and min for each channel
	const auto [max_r, min_r] = histMaxMin<0>(bitmap);
	const auto [max_g, min_g] = histMaxMin<1>(bitmap);
	const auto [max_b, min_b] = histMaxMin<2>(bitmap);

	//normalize each channel
	for(jab_int32 i=0; i<bitmap.height; i++)
	{
		for(jab_int32 j=0; j<bitmap.width; j++)
		{
			jab_int32 offset = i * bytes_per_row + j * bytes_per_pixel;
			//R channel
			if (bitmap.pixel[offset + 0] < min_r) { 
				bitmap.pixel[offset + 0] = 0; 
			}
			else if (bitmap.pixel[offset + 0] > max_r) { 
				bitmap.pixel[offset + 0] = 255; 
			}
			else {
				bitmap.pixel[offset + 0] = (jab_byte)((jab_double)(bitmap.pixel[offset + 0] - min_r) / (jab_double)(max_r - min_r) * 255.0);
			}
			//G channel
			if (bitmap.pixel[offset + 1] < min_g) { 
				bitmap.pixel[offset + 1] = 0; 
			}
			else if (bitmap.pixel[offset + 1] > max_g) { 
				bitmap.pixel[offset + 1] = 255; 
			}
			else { 
				bitmap.pixel[offset + 1] = (jab_byte)((jab_double)(bitmap.pixel[offset + 1] - min_g) / (jab_double)(max_g - min_g) * 255.0);		
			}
			//B channel
			if (bitmap.pixel[offset + 2] < min_b) { 
				bitmap.pixel[offset + 2] = 0; 
			}
			else if (bitmap.pixel[offset + 2] > max_b) { 
				bitmap.pixel[offset + 2] = 255; 
			}
			else { 
				bitmap.pixel[offset + 2] = (jab_byte)((jab_double)(bitmap.pixel[offset + 2] - min_b) / (jab_double)(max_b - min_b) * 255.0);
			}
		}
	}
}

/**
 * @brief Get the average and variance of RGB values
 * @param rgb the pixel with RGB values
 * @param ave the average value
 * @param var the variance value
*/
constexpr std::pair<jab_double, jab_double> getAveVar(const jab_byte* rgb) noexcept {
	const jab_double ave = (rgb[0] + rgb[1] + rgb[2]) / 3;

	//calculate variance
	const jab_double sum = (rgb[0] - ave) * (rgb[0] - ave) + 
						   (rgb[1] - ave) * (rgb[1] - ave) + 
		                   (rgb[2] - ave) * (rgb[2] - ave);
	const jab_double var = sum / 3;
	return { ave, var };
}

/**
 * @brief Swap two variables
*/
static void swap(jab_int32* a, jab_int32* b)
{
	jab_int32 tmp = *a;
	*a = *b;
	*b = tmp;
}

/**
 * @brief Get the min, middle and max value of three values and the corresponding indexes
 * @param rgb the pixel with RGB values
 * @param min the min value
 * @param mid the middle value
 * @param max the max value
*/
void getMinMax(const jab_byte* rgb, jab_byte* min, jab_byte* mid, jab_byte* max, jab_int32* index_min, jab_int32* index_mid, jab_int32* index_max) noexcept {
	*index_min = 0;
	*index_mid = 1;
	*index_max = 2;
	if (rgb[*index_min] > rgb[*index_max]) {
		swap(index_min, index_max);
	}
	if (rgb[*index_min] > rgb[*index_mid]) {
		swap(index_min, index_mid);
	}
	if (rgb[*index_mid] > rgb[*index_max]) {
		swap(index_mid, index_max);
	}
	*min = rgb[*index_min];
	*mid = rgb[*index_mid];
	*max = rgb[*index_max];
}

/**
 * @brief Binarize a color channel of a bitmap using local binarization algorithm
 * @param bitmap the input bitmap
 * @param rgb the binarized RGB channels
 * @param blk_ths the black color thresholds for RGB channels
 * @return JAB_SUCCESS | JAB_FAILURE
*/
jab_boolean binarizerRGB(const Bitmap& bitmap, std::array<Bitmap, 3>& rgb, const std::array<jab_float, 3>& blk_ths) noexcept {
#if false
	// TODO: Investigate disabling this block.
	//	     Disabling this so far has not failed any tests
	for (auto& color : rgb) {
		color.pixel = std::vector<jab_byte>(color.width * color.height, 0);
	}
#endif

	jab_int32 bytes_per_pixel = bitmap.bits_per_pixel / 8;
    jab_int32 bytes_per_row = bitmap.width * bytes_per_pixel;

    //calculate the average pixel value, block-wise
    jab_int32 max_block_size = std::max(bitmap.width, bitmap.height) / 2;
    jab_int32 block_num_x = (bitmap.width % max_block_size) != 0 ? (bitmap.width / max_block_size) + 1 : (bitmap.width / max_block_size);
    jab_int32 block_num_y = (bitmap.height% max_block_size) != 0 ? (bitmap.height/ max_block_size) + 1 : (bitmap.height/ max_block_size);
    jab_int32 block_size_x = bitmap.width / block_num_x;
    jab_int32 block_size_y = bitmap.height/ block_num_y;

	//binarize each pixel in each channel
	jab_double ths_std = 0.08;
	jab_float rgb_ths[3] = {0, 0, 0};
    for(jab_int32 i=0; i<bitmap.height; i++)
	{
		for(jab_int32 j=0; j<bitmap.width; j++)
		{
			jab_int32 offset = i * bytes_per_row + j * bytes_per_pixel;
			//check black pixel
			rgb_ths[0] = blk_ths[0];
			rgb_ths[1] = blk_ths[1];
			rgb_ths[2] = blk_ths[2];

			if(bitmap.pixel[offset + 0] < rgb_ths[0] && bitmap.pixel[offset + 1] < rgb_ths[1] && bitmap.pixel[offset + 2] < rgb_ths[2])
            {
                rgb[0].pixel[i*bitmap.width + j] = 0;
                rgb[1].pixel[i*bitmap.width + j] = 0;
                rgb[2].pixel[i*bitmap.width + j] = 0;
                continue;
            }

			const auto [ ave, var ] = getAveVar(&bitmap.pixel[offset]);
			jab_double std = sqrt(var);	//standard deviation
			jab_byte min, mid, max;
			jab_int32 index_min, index_mid, index_max;
			getMinMax(&bitmap.pixel[offset], &min, &mid, &max, &index_min, &index_mid, &index_max);
			std /= (jab_double)max;	//normalize std

			if(std < ths_std && (bitmap.pixel[offset + 0] > rgb_ths[0] && bitmap.pixel[offset + 1] > rgb_ths[1] && bitmap.pixel[offset + 2] > rgb_ths[2]))
			{
				rgb[0].pixel[i*bitmap.width + j] = 255;
				rgb[1].pixel[i*bitmap.width + j] = 255;
				rgb[2].pixel[i*bitmap.width + j] = 255;
			}
			else
			{
				rgb[index_max].pixel[i*bitmap.width + j] = 255;
				rgb[index_min].pixel[i*bitmap.width + j] = 0;
				jab_double r1 = (jab_double)bitmap.pixel[offset + index_mid] / (jab_double)bitmap.pixel[offset + index_min];
				jab_double r2 = (jab_double)bitmap.pixel[offset + index_max] / (jab_double)bitmap.pixel[offset + index_mid];
				if(r1 > r2)
					rgb[index_mid].pixel[i*bitmap.width + j] = 255;
				else
					rgb[index_mid].pixel[i*bitmap.width + j] = 0;
			}
		}
	}
	filterBinary(rgb[0]);
	filterBinary(rgb[1]);
	filterBinary(rgb[2]);
	return JAB_SUCCESS;
}
std::array<Bitmap, 3> binarizerRGB(const Bitmap& bitmap) noexcept {
	jab_int32 bytes_per_pixel = bitmap.bits_per_pixel / 8;
	jab_int32 bytes_per_row = bitmap.width * bytes_per_pixel;

	//calculate the average pixel value, block-wise
	jab_int32 max_block_size = std::max(bitmap.width, bitmap.height) / 2;
	jab_int32 block_num_x = (bitmap.width % max_block_size) != 0 ? (bitmap.width / max_block_size) + 1 : (bitmap.width / max_block_size);
	jab_int32 block_num_y = (bitmap.height % max_block_size) != 0 ? (bitmap.height / max_block_size) + 1 : (bitmap.height / max_block_size);
	jab_int32 block_size_x = bitmap.width / block_num_x;
	jab_int32 block_size_y = bitmap.height / block_num_y;

	std::array<jab_float, 3>* pixel_ave = new std::array<jab_float, 3>[block_num_x * block_num_y]{};
	for (jab_int32 i = 0; i < block_num_y; i++)
	{
		for (jab_int32 j = 0; j < block_num_x; j++)
		{
			jab_int32 block_index = i * block_num_x + j;

			jab_int32 sx = j * block_size_x;
			jab_int32 ex = (j == block_num_x - 1) ? bitmap.width : (sx + block_size_x);
			jab_int32 sy = i * block_size_y;
			jab_int32 ey = (i == block_num_y - 1) ? bitmap.height : (sy + block_size_y);
			jab_int32 counter = 0;
			for (jab_int32 y = sy; y < ey; y++)
			{
				for (jab_int32 x = sx; x < ex; x++)
				{
					jab_int32 offset = y * bytes_per_row + x * bytes_per_pixel;
					pixel_ave[block_index][0] += bitmap.pixel[offset + 0];
					pixel_ave[block_index][1] += bitmap.pixel[offset + 1];
					pixel_ave[block_index][2] += bitmap.pixel[offset + 2];
					counter++;
				}
			}
			pixel_ave[block_index][0] /= (jab_float)counter;
			pixel_ave[block_index][1] /= (jab_float)counter;
			pixel_ave[block_index][2] /= (jab_float)counter;
		}
	}

	std::array<Bitmap, 3> rgb = {
		Bitmap(bitmap.width, bitmap.height, 8, 8, 1, bitmap.width * bitmap.height, 0),
		Bitmap(bitmap.width, bitmap.height, 8, 8, 1, bitmap.width * bitmap.height, 0),
		Bitmap(bitmap.width, bitmap.height, 8, 8, 1, bitmap.width * bitmap.height, 0),
	};

	//binarize each pixel in each channel
	jab_double ths_std = 0.08;
	jab_float rgb_ths[3] = { 0, 0, 0 };
	for (jab_int32 i = 0; i < bitmap.height; i++)
	{
		for (jab_int32 j = 0; j < bitmap.width; j++)
		{
			jab_int32 offset = i * bytes_per_row + j * bytes_per_pixel;
			//check black pixel
			jab_int32 block_index = std::min(i / block_size_y, block_num_y - 1) * block_num_x + std::min(j / block_size_x, block_num_x - 1);
			rgb_ths[0] = pixel_ave[block_index][0];
			rgb_ths[1] = pixel_ave[block_index][1];
			rgb_ths[2] = pixel_ave[block_index][2];

			if (bitmap.pixel[offset + 0] < rgb_ths[0] && bitmap.pixel[offset + 1] < rgb_ths[1] && bitmap.pixel[offset + 2] < rgb_ths[2])
			{
				rgb[0].pixel[i * bitmap.width + j] = 0;
				rgb[1].pixel[i * bitmap.width + j] = 0;
				rgb[2].pixel[i * bitmap.width + j] = 0;
				continue;
			}

			const auto [ ave, var ] = getAveVar(&bitmap.pixel[offset]);
			jab_double std = std::sqrt(var);	//standard deviation
			jab_byte min, mid, max;
			jab_int32 index_min, index_mid, index_max;
			getMinMax(&bitmap.pixel[offset], &min, &mid, &max, &index_min, &index_mid, &index_max);
			std /= (jab_double)bitmap.pixel[offset + index_max];	//normalize std

			if (std < ths_std && (bitmap.pixel[offset + 0] > rgb_ths[0] && bitmap.pixel[offset + 1] > rgb_ths[1] && bitmap.pixel[offset + 2] > rgb_ths[2]))
			{
				rgb[0].pixel[i * bitmap.width + j] = 255;
				rgb[1].pixel[i * bitmap.width + j] = 255;
				rgb[2].pixel[i * bitmap.width + j] = 255;
			}
			else
			{
				rgb[index_max].pixel[i * bitmap.width + j] = 255;
				rgb[index_min].pixel[i * bitmap.width + j] = 0;
				jab_double r1 = (jab_double)bitmap.pixel[offset + index_mid] / (jab_double)bitmap.pixel[offset + index_min];
				jab_double r2 = (jab_double)bitmap.pixel[offset + index_max] / (jab_double)bitmap.pixel[offset + index_mid];
				if (r1 > r2)
					rgb[index_mid].pixel[i * bitmap.width + j] = 255;
				else
					rgb[index_mid].pixel[i * bitmap.width + j] = 0;
			}
		}
	}
	delete[] pixel_ave;
	filterBinary(rgb[0]);
	filterBinary(rgb[1]);
	filterBinary(rgb[2]);
	return rgb;
}
