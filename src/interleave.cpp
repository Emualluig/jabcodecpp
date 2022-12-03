#include "interleave.h"

#include "reportError.h"

#include "encoder.h"
#include "pseudo_random.h"

#define INTERLEAVE_SEED 226759

// Interleaving then deinterleaving is about 2x faster using these functions
static std::array<uint64_t, 2> recent = { 0, 0 }; // A small cache, most images will be the same size, so it is likely to have same data.size()
constexpr void FastInterleaveData(std::vector<jab_char>& data) noexcept {
    LCG lcg = LCG(INTERLEAVE_SEED);
    const std::size_t size = data.size();
    for (jab_int32 i = 0; i < size; i++) {
        jab_int32 pos = (jab_int32)((jab_float)lcg.temper() / (jab_float)UINT32_MAX * (size - i));
        std::swap(data[size - 1 - i], data[pos]);
    }
    if (recent[0] != size) {
        // Update cache
        lcg.temper();
        recent = { size, lcg.getSeed() };
    }
}
constexpr void FastDeinterleaveData(std::vector<jab_char>& data) noexcept {
    const std::size_t size = data.size();
    LCG lcg = LCG(recent[1]);
    if (recent[0] != size) {
        // Even if there is a recent miss, it is still faster to calculate the LCG
        lcg = LCG(INTERLEAVE_SEED);
        for (int i = 0; i <= size; i++) {
            lcg.next();
        }
        recent = { size, lcg.getSeed() }; // Update cache
    }
    for (int i = size - 1; i >= 0; i--) {
        jab_int32 pos = (jab_int32)((jab_float)lcg.rtemper() / (jab_float)UINT32_MAX * (size - i));
        std::swap(data[size - 1 - i], data[pos]);
    }
}

/**
 * @brief In-place interleaving
 * @param data the input data to be interleaved
*/
void interleaveData(std::vector<jab_char>& data) noexcept {
    FastInterleaveData(data);
}

/**
 * @brief In-place deinterleaving
 * @param data the input data to be deinterleaved
*/
void deinterleaveData(std::vector<jab_char>& data) noexcept {
    FastDeinterleaveData(data);
}