#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace fpm {

/**
 * @brief Stable SplitMix64-based RNG used by claim-bearing simulation code.
 *
 * The integer-to-random mapping, rejection sampling and Fisher-Yates shuffle are
 * part of the reproducibility contract. Replacing this class with standard-library
 * distributions can change seed-to-trajectory mapping across platforms.
 */
class StableRng {
    uint64_t state;
public:
    explicit StableRng(uint64_t seed = 0) : state(seed) {}

    /** Return the next reproducible 64-bit value. */
    uint64_t next_u64() {
        uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    /** Return a uniform double in [0,1) using the high 53 bits. */
    double uniform01() { return (next_u64() >> 11) * 0x1.0p-53; }

    /** Uniform integer in [0,bound), using rejection sampling to avoid modulo bias. */
    uint64_t bounded(uint64_t bound) {
        if (bound == 0) throw std::runtime_error("RNG bound is zero");
        const uint64_t threshold = uint64_t(-bound) % bound;
        while (true) {
            const uint64_t value = next_u64();
            if (value >= threshold) return value % bound;
        }
    }

    /** In-place reproducible Fisher-Yates shuffle. */
    template <typename T>
    void shuffle(std::vector<T> &items) {
        for (size_t i = items.size(); i > 1; --i) {
            const size_t j = size_t(bounded(i));
            std::swap(items[i - 1], items[j]);
        }
    }
};

} // namespace fpm
