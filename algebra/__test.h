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

// Workaround for a libc++ std::print() issue; not present (and not compilable) with libstdc++.
#ifdef _LIBCPP_VERSION
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
            const U span = static_cast<U>(max) - static_cast<U>(min);
            if (span == 0)
                return min;
            U r;
            const U limit = (span == std::numeric_limits<U>::max()) ? span : (span + 1);
            do { // rejection sampling, so the result is uniform
                r = (static_cast<U>(_rng()) << 64) | _rng();
                if (limit != 0)
                    r %= limit;
            } while (false);
            return static_cast<T>(static_cast<U>(min) + r);
        }
    }

private:
    std::mt19937_64 _rng;
};
