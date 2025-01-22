#pragma once
#include <catch2/catch_test_macros.hpp>
#include <random>
#include <print>

using std::format;
using std::print;
using namespace algebra;
using namespace algebra::literals;

using ucent = unsigned __int128;
using cent = __int128;
using ulong = unsigned long;
using uint = unsigned int;

namespace std {
inline namespace __1 {
bool __is_posix_terminal(__sFILE*) { return true; }
}
}

class Random {
public:
    Random() : _rng(std::random_device{}()) {}
    Random(unsigned long seed) : _rng(seed) {}
    operator std::mt19937_64&() { return _rng; }
    std::mt19937_64& get() { return _rng; }

    template<typename T>
    T Uniform(T min, T max) {
        static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>);

        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(_rng);
        }
        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(_rng);
        }
        if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min, max);
            return dist(_rng);
        }
    }

    rational SampleQ() { return rational{Uniform<int>(-1000, 1000), Uniform<int>(1, 1000)}; }
    Vec3<rational> SampleVec3Q() { return {SampleQ(), SampleQ(), SampleQ()}; }

private:
    std::mt19937_64 _rng;
};
