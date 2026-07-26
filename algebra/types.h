#pragma once
#include <bit>
#include <concepts>
#include <cstdint>
#include <limits>
#include <cmath>
#include <type_traits>

namespace algebra {

using int128_t = __int128;
using uint128_t = unsigned __int128;
static const auto UINT128_MAX = std::numeric_limits<uint128_t>::max();

// gcc only treats __int128 as an integral type outside of strict standard mode
// (-std=gnu++23 but not -std=c++23), and libstdc++ leaves std::make_unsigned and
// std::countr_zero undefined for it either way, so both are named explicitly here.
template<typename T> concept __int128_like =
    std::same_as<std::remove_cv_t<T>, int128_t> || std::same_as<std::remove_cv_t<T>, uint128_t>;

template<typename T> concept std_int = (std::integral<T> || __int128_like<T>) && !std::same_as<std::remove_cv_t<T>, bool>;
template<typename T> concept std_signed_int =
    (std::signed_integral<T> || std::same_as<std::remove_cv_t<T>, int128_t>) && !std::same_as<std::remove_cv_t<T>, bool>;
template<typename T> concept std_unsigned_int =
    (std::unsigned_integral<T> || std::same_as<std::remove_cv_t<T>, uint128_t>) && !std::same_as<std::remove_cv_t<T>, bool>;

template<typename T> struct __make_unsigned { using type = std::make_unsigned_t<T>; };
template<> struct __make_unsigned<int128_t> { using type = uint128_t; };
template<> struct __make_unsigned<uint128_t> { using type = uint128_t; };
template<std_int T> using make_unsigned_t = typename __make_unsigned<std::remove_cv_t<T>>::type;

template<std_int T> constexpr make_unsigned_t<T> make_unsigned(T a) { return a; }
constexpr auto abs_unsigned(std_int auto a) { return (a < 0) ? (~make_unsigned(a) + 1) : make_unsigned(a); }

// std::countr_zero() does not accept __int128
template<std_unsigned_int T>
constexpr int countr_zero(T a) {
    if constexpr (sizeof(T) <= 8) {
        return std::countr_zero(a);
    } else {
        const uint64_t low = static_cast<uint64_t>(a);
        return low ? std::countr_zero(low) : 64 + std::countr_zero(static_cast<uint64_t>(a >> 64));
    }
}

}
