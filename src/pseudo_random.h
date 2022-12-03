#ifndef _HEADER_PSEUDO_RANDOM_H_
#define _HEADER_PSEUDO_RANDOM_H_

#include <cstdint>

class LCG {
	uint64_t _seed;
	static constexpr uint64_t defaultSeed = 42;
	constexpr uint32_t _temper(uint32_t x) const noexcept {
		x ^= x >> 11;
		x ^= x << 7 & 0x9D2C5680;
		x ^= x << 15 & 0xEFC60000;
		x ^= x >> 18;
		return x;
	}
public:
	constexpr explicit LCG(uint64_t seed = LCG::defaultSeed) noexcept : _seed{ seed }  {}
	constexpr uint32_t temper() noexcept {
		_seed = 6364136223846793005ULL * _seed + 1;
		return _temper(_seed >> 32);
	}
	constexpr void next() noexcept {
		_seed = 6364136223846793005ULL * _seed + 1;
	}
	constexpr uint32_t rtemper() noexcept {
		// Calculates the previous _seed (can't go infinitely backwards)
		_seed = 13877824140714322085ULL * (_seed - 1); // The number is the modular multiplicative inverse, found with the extented euclidean algorithm.
		return _temper(_seed >> 32);
	}
	constexpr uint64_t getSeed() const noexcept {
		return _seed;
	}
};

#endif // !_HEADER_PSEUDO_RANDOM_H_

