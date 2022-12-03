#include "JabCPP.h"

#include "src/jabTypes.h"
#include "src/image.h"
#include "src/encoder.h"
#include "src/detector.h"
#include "src/reportError.h"

#include <iostream>

std::optional<std::vector<unsigned char>> JabCPP::decodeFromMemory(const unsigned char* pointer, unsigned int size) noexcept {
	std::optional<Bitmap> jbo = readImageMemory(pointer, size);
	if (!jbo.has_value()) {
		return std::nullopt;
	}
	// TODO: do something more with the decode_status
	jab_int32 decode_status = 0;
	std::optional<std::vector<jab_byte>> decoded_dataOp = decodeJABCode(jbo.value(), NORMAL_DECODE, decode_status);
	if (!decoded_dataOp.has_value()) {
		reportError("Decoding JABCode failed\n");
		return std::nullopt;
	}
	if (decode_status == 2) {
		reportError("The code is only partly decoded. Some slave symbols have not been decoded and are ignored.\n");
	}
	return decoded_dataOp.value();
}

constexpr static ColorPaletteType getColorPaletteType(JabCPP::ColorPalette palette) noexcept {
	if (palette == JabCPP::ColorPalette::FULL) {
		return ColorPaletteType::FULL;
	}
	else {
		return ColorPaletteType::CMYK;
	}
}

std::pair<unsigned char*, unsigned int> JabCPP::encodeToMemory(const std::vector<char>& bytes, const JabCPP::Settings& settings) noexcept {
	const jab_int32 symbol_number = settings.numberOfSymbols;
	const ColorPaletteType palette = getColorPaletteType(settings.palette);

	Encode enc = createEncode(palette, symbol_number);
	enc.module_size = settings.moduleSize;
	const int pixelsPerBlock = enc.module_size; // Typically 12 (minimum 1) (unable to decode if under 4)
	const int numberOfModules = 21; // Typically 21
	enc.master_symbol_width = numberOfModules * pixelsPerBlock; // Must be at least 84 for image to be decodable
	enc.master_symbol_height = numberOfModules * pixelsPerBlock;
	for (jab_int32 loop = 0; loop < symbol_number; loop++) {
		if (settings.symbolEccLevels.size() > 0) {
			enc.symbol_ecc_levels[loop] = settings.symbolEccLevels[loop];
		}
		if (settings.symbolVersions.size() > 0) {
			enc.symbol_versions[loop] = {
				.x = settings.symbolVersions[loop].x,
				.y = settings.symbolVersions[loop].y,
			};
		}
		if (settings.symbolPositions.size() > 0) {
			enc.symbol_positions[loop] = settings.symbolPositions[loop];
		}
	}
	// TODO: return error codes. Example: error = 4 means too much data for the current jab code size
	jab_int32 error = 0;
	std::optional<Bitmap> bitmapOp = generateJABCode(enc, bytes, error);
	if (!bitmapOp.has_value()) {
		JAB_REPORT_ERROR("Error code: ", error, "\n");
		return { nullptr, 0 };
	}
	return saveImageMemory(bitmapOp.value());
};