

#include "image.h"

#include "reportError.h"

extern "C" {
#include "png.h"
}

/**
 * @brief Save code bitmap in RGB as png image
 * @param bitmap the code bitmap
 * @param filename the image filename
 * @return JAB_SUCCESS | JAB_FAILURE
*/
jab_boolean saveImage(const Bitmap& bitmap, const char* filename) noexcept { // TODO:
	png_image image;
    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;

    if(bitmap.channel_count == 4) {
		image.format = PNG_FORMAT_RGBA;
		image.flags  = PNG_FORMAT_FLAG_ALPHA | PNG_FORMAT_FLAG_COLOR;
    }
    else {
		image.format = PNG_FORMAT_GRAY;
    }
    image.width  = bitmap.width;
    image.height = bitmap.height;

    if (png_image_write_to_file(&image,
								filename,
								0/*convert_to_8bit*/,
								bitmap.pixel.data(),
								0/*row_stride*/,
								NULL/*colormap*/) == 0)
	{
		reportError(image.message);
		reportError("Saving png image failed");
		return JAB_FAILURE;
	}
	return JAB_SUCCESS;
}

/**
 * @brief Read image into code bitmap
 * @param filename the image filename
 * @return Pointer to the code bitmap read from image | NULL
*/
std::optional<Bitmap> readImage(const char* filename) noexcept { // TODO:
	png_image image;
    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;

	if (png_image_begin_read_from_file(&image, filename) == 0) {
		reportError(image.message);
		reportError("Opening png image failed");
		return std::nullopt;
	}
	image.format = PNG_FORMAT_RGBA;

	Bitmap bitmap = Bitmap(
		image.width,
		image.height,
		BITMAP_BITS_PER_PIXEL,
		BITMAP_BITS_PER_CHANNEL,
		BITMAP_CHANNEL_COUNT,
		PNG_IMAGE_SIZE(image), 0
	);

	if (png_image_finish_read(&image, NULL/*background*/, bitmap.pixel.data(), 0/*row_stride*/, NULL/*colormap*/) == 0) {
		reportError(image.message);
		reportError("Reading png image failed");
		return std::nullopt;
	}
	return bitmap;
}

std::pair<unsigned char*, unsigned int> saveImageMemory(const Bitmap& bitmap) noexcept {
	png_image image{};
	std::memset(&image, 0, sizeof(image));
	image.version = PNG_IMAGE_VERSION;
	if (bitmap.channel_count == 4) {
		image.format = PNG_FORMAT_RGBA;
		image.flags = PNG_FORMAT_FLAG_ALPHA | PNG_FORMAT_FLAG_COLOR;
	}
	else {
		image.format = PNG_FORMAT_GRAY;
	}
	image.width = bitmap.width;
	image.height = bitmap.height;

	png_alloc_size_t png_size = 0;
	bool result = png_image_write_get_memory_size(image, png_size, false, bitmap.pixel.data(), 0, NULL);
	if (!result) {
		std::cout << "Error writing to memory, 1.\n";
		std::cout << image.message << "\n";
		return { nullptr, (unsigned int)png_size };
	}
	unsigned char* ptr = new unsigned char[png_size];
	if (png_image_write_to_memory(&image, ptr, &png_size, false, bitmap.pixel.data(), 0, NULL) == 0) {
		std::cout << "Error writing to memory, 2.\n";
		std::cout << image.message << "\n";
		return { nullptr, (unsigned int)png_size };
	}

	return { ptr, (unsigned int)png_size };
}
/**
 * @brief Read image into code bitmap
 * @return Pointer to the code bitmap read from image | NULL
*/
std::optional<Bitmap> readImageMemory(const unsigned char* pointer, int length) noexcept {
	png_image image{}; // Default initialize (Will be all zeroes since this is a POD C struct)
    image.version = PNG_IMAGE_VERSION;

	if (png_image_begin_read_from_memory(&image, pointer, length) == 0) {
		reportError("Error reading PNG from memory");
		reportError(image.message);
		return std::nullopt;
	}
	image.format = PNG_FORMAT_RGBA; // This needs to be here

	Bitmap bitmap = Bitmap(
		image.width, image.height, BITMAP_BITS_PER_PIXEL, 
		BITMAP_BITS_PER_CHANNEL, BITMAP_CHANNEL_COUNT,
		PNG_IMAGE_SIZE(image), 0
	);

	if (png_image_finish_read(&image, NULL, bitmap.pixel.data(), 0, NULL) == 0) {
		reportError("Reading png image failed");
		reportError(image.message);
		return std::nullopt;
	}

	return bitmap;
}
