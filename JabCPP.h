#ifndef _HEADER_JABCODE_H_
#define _HEADER_JABCODE_H_

#include <optional>
#include <vector>
namespace JabCPP {
	enum ColorPalette {
		FULL = 8,
		CMYK = 4,
	};
	struct Point2d {
		int32_t x;
		int32_t y;
	};
	static const int32_t DefaultNumberOfSymbols = 1;
	static const int32_t DefaultEccLevel = 1;
	// The Symbol version essentially is the x and y size of a symbol, higher = more data but larger jabcode.
	static const Point2d DefaultSymbolVersion = Point2d{ 1, 1 };
	// The position 0 is for the master symbol.
	static const int32_t DefaultSymbolPosition = 0;
	// Each module (square) in the jabcode will be 8x8 pixels. Note that the minimum decodable size of image is 84x84.
	static const int32_t DefaultModuleSize = 8;

	// These settings strive to replicate the options available at https://jabcode.org/create
	// The settings can adjust many aspects of the created code.
	// Settings play no role in decoding.
	struct Settings {
		// The color palette that will be used to generate the jabcode. FULL is 8 colors and CMYK is 4 colors.
		ColorPalette palette;
		// This is the total number of symbols in the jabcode. numberOfSymbols > 0.
		int numberOfSymbols;
		// The Error Correction Level for each symbol.
		// The size of the vector must match with the number of symbols. Minimum ECC is 1, maximum is 10.
		// The first value is for the master symbol, then the secondary symbols.
		std::vector<int32_t> symbolEccLevels;
		// The symbol version of each symbol.
		// The first value is for the master symbol, then the secondary symbols.
		// The version is essentially the size. For example if the version.x = 5, there will be 17 + version.x * 4 modules as width of the symbol.
		// Minimum 1, maximum 32.
		std::vector<Point2d> symbolVersions;
		// The symbol position of each symbol. 
		// View https://jabcode.org/create for a information about symbol position.
		std::vector<int32_t> symbolPositions;
		// The Module size. A module is an individual square in the jabcode.
		// If moduleSize = 1, then each individual module would be 1x1 pixels in the png.
		// Minimum 1
		int32_t moduleSize;
	};

	// These are the default encoding settings. They are those found at https://jabcode.org/create
	static const Settings DefaultSettings = {
		// Default is all 8 colors.
		.palette = ColorPalette::FULL,
		// Default is 1 symbol
		.numberOfSymbols = JabCPP::DefaultNumberOfSymbols,
		// The default is 1 on every symbol
		.symbolEccLevels = std::vector<int32_t>(JabCPP::DefaultNumberOfSymbols, JabCPP::DefaultEccLevel),
		// Every symbol has a version of {1, 1}. For a total of 21 block in each dimension.
		.symbolVersions = std::vector<Point2d>(JabCPP::DefaultNumberOfSymbols, JabCPP::DefaultSymbolVersion),
		// By default we only have 1 symbol, the master symbol.
		.symbolPositions = std::vector<int32_t>(JabCPP::DefaultNumberOfSymbols, JabCPP::DefaultSymbolPosition),
		// By default we have 8x8 modules.
		.moduleSize = JabCPP::DefaultModuleSize,
	};
	
	// decodeFromMemory takes a pointer to a png in memory and the length of the memory and attempts to decode it.
	// returns: If the decoding is sucessful, a vector of bytes is returned. Else std::nullopt.
	// The png should be RGB and the size is the number of bytes long it is.
	// Note: There will be error when decoding if the height or width is less than 84 pixels.
	std::optional<std::vector<unsigned char>> decodeFromMemory(const unsigned char* pointer, unsigned int size) noexcept;

	// encodeToMemory takes bytes and settings and writes a png to memory.
	// returns: the pointer to the memory and number of bytes
	// The returned pointer points to a rgb png. Caller must free.
	std::pair<unsigned char*, unsigned int> encodeToMemory(const std::vector<char>& bytes, const JabCPP::Settings& settings = JabCPP::DefaultSettings) noexcept;
};

#endif // !_HEADER_JABCODE_H_
