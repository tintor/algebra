#pragma once
#include <catch2/catch_test_macros.hpp>
#include "algebra/types.h"
#include <limits>
#include <random>
#include <print>

namespace algebra { namespace literals { } }

using std::format;
using std::print;
using namespace algebra;
using namespace algebra::literals;

using ucent = unsigned __int128;
using cent = __int128;
using ulong = unsigned long;
using uint = unsigned int;

// Workaround for a libc++ std::print() issue on Apple platforms. __sFILE is a BSD type, so
// this does not compile with libstdc++ or with libc++ on Linux.
#if defined(_LIBCPP_VERSION) && defined(__APPLE__)
namespace std {
inline namespace __1 {
bool __is_posix_terminal(__sFILE*) { return true; }
}
}
#endif

class Random {
public:
    Random() : _rng(std::random_device{}()) {}
    Random(unsigned long seed) : _rng(seed) {}
    operator std::mt19937_64&() { return _rng; }
    std::mt19937_64& get() { return _rng; }

    template<typename T>
    T Uniform(T min, T max) {
        // note: std::is_integral_v is false for __int128 in strict standard mode, and
        // std::uniform_int_distribution does not accept it either way
        static_assert(algebra::std_int<T> || std::is_floating_point_v<T>);

        if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min, max);
            return dist(_rng);
        } else if constexpr (sizeof(T) <= 8) {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(_rng);
        } else {
            using U = algebra::make_unsigned_t<T>;
            const U max_u = algebra::UINT128_MAX;
            const U span = static_cast<U>(max) - static_cast<U>(min);
            if (span == 0)
                return min;

            auto draw = [this] { return (static_cast<U>(_rng()) << 64) | _rng(); };
            if (span == max_u)
                return static_cast<T>(draw()); // every draw is in range, and min is the lowest value

            // Rejection sampling: taking the remainder of a raw draw would favour the values below
            // (2**128 mod count), so the draws that fall in that last partial block are discarded.
            const U count = span + 1;
            const U keep = (max_u / count) * count; // the largest multiple of count that fits
            U r;
            do {
                r = draw();
            } while (r >= keep);
            return static_cast<T>(static_cast<U>(min) + r % count);
        }
    }

private:
    std::mt19937_64 _rng;
};
