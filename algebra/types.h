#pragma once
#include <concepts>
#include <limits>
#include <cmath>
#include <type_traits>

namespace algebra {

template<typename T> concept std_int = std::integral<T> && !std::same_as<T, bool>;
template<typename T> concept std_signed_int = std::signed_integral<T> && !std::same_as<T, bool>;
template<typename T> concept std_unsigned_int = std::unsigned_integral<T> && !std::same_as<T, bool>;

template<std_int T> constexpr std::make_unsigned_t<T> make_unsigned(T a) { return a; }
constexpr auto abs_unsigned(std_int auto a) { return (a < 0) ? (~make_unsigned(a) + 1) : make_unsigned(a); }

using int128_t = __int128;
using uint128_t = unsigned __int128;
static const auto UINT128_MAX = std::numeric_limits<uint128_t>::max();

}
