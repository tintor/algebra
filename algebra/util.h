#pragma once
#include <source_location>
#include <format>
#include <print>
#include <stdexcept>
#include <bit>
#include <algorithm>
#include "algebra/types.h"

namespace algebra {

constexpr void Check(bool value, std::source_location loc = std::source_location::current()) {
    if (!value)
        throw std::runtime_error(std::format("Check failed at {}:{} in {}", loc.file_name(), loc.line(), loc.function_name()));
}

constexpr void Check(bool value, const char* message, std::source_location loc = std::source_location::current()) {
    if (!value)
        throw std::runtime_error(std::format("Check failed at {}:{} in {} with message: {}", loc.file_name(), loc.line(), loc.function_name(), message));
}

[[noreturn]] constexpr void Fail(const char* message, std::source_location loc = std::source_location::current()) {
    throw std::runtime_error(std::format("Failed at {}:{} in {} with message: {}", loc.file_name(), loc.line(), loc.function_name(), message));
}

template<typename T> struct IsNumberClass : std::false_type {};
template<typename T> concept __ncsi = IsNumberClass<T>::value || std_int<T>;

constexpr bool operator>(const __ncsi auto& a, const __ncsi auto& b) { return b < a; }
constexpr bool operator>=(const __ncsi auto& a, const __ncsi auto& b) { return !(a < b); }
constexpr bool operator<=(const __ncsi auto& a, const __ncsi auto& b) { return !(b < a); }
constexpr bool operator!=(const __ncsi auto& a, const __ncsi auto& b) { return !(a == b); }

constexpr auto operator+(std_int auto a, const __ncsi auto& b) { return b + a; }
constexpr auto operator*(std_int auto a, const __ncsi auto& b) { return b * a; }
constexpr bool operator==(std_int auto a, const __ncsi auto& b) { return b == a; }

constexpr uint128_t __mulq(uint64_t a, uint64_t b) { return uint128_t(a) * b; }

using int128_t = __int128;
using uint128_t = unsigned __int128;

#if defined(__x86_64__) && !defined(__ILP32__)
constexpr void __divq(uint128_t a, uint64_t b, uint64_t& q, uint64_t& r) {
    uint64_t hi = a >> 64;
    uint64_t lo = a;
    __asm__ (
        "divq %[divisor]"
        : "+a" (lo),  // Input: lo in RAX, Output: quotient in RAX
            "=d" (r)    // Input: hi in RDX, Output: remainder in RDX
        : [divisor] "r" (b),
            "d" (hi)
        : "cc"
    );
    q = lo;
}
#else
constexpr void __divq(uint128_t a, uint64_t b, uint64_t& q, uint64_t& r) {
    q = a / b;
    m = a % b;
}
#endif

#if defined(__x86_64__) && !defined(__ILP32__)
constexpr uint64_t __divq(uint128_t a, uint64_t b) {
    uint64_t hi = a >> 64;
    uint64_t lo = a;
    // Assembly clobbers hi (RDX) but we don't need it anymore
    __asm__ (
        "divq %[divisor]"
        : "+a" (lo)      // Input: lo in RAX, Output: quotient in RAX
        : [divisor] "r" (b),
          "d" (hi)       // High bits in RDX (implicitly clobbered)
        : "cc"           // Clobbers condition codes
    );
    return lo;
}
#else
constexpr uint64_t __divq(uint128_t a, uint64_t b) { return a / b; }
#endif

#if defined(__x86_64__) && !defined(__ILP32__)
constexpr uint64_t __divq_mod(uint128_t a, uint64_t b) {
    uint64_t hi = a >> 64;
    uint64_t lo = a;
    uint64_t m;
    __asm__ (
        "divq %[divisor]"
        : "=d" (m), "+a" (lo)
        : [divisor] "r" (b), "d" (hi)
    );
    return m;
}
#else
constexpr uint64_t __divq_mod(uint128_t a, uint64_t b) { return a % b; }
#endif

constexpr uint64_t pow(uint64_t base, unsigned exp) {
    if (base == 2)
        return 1 << exp;
    if (exp == 0)
        return 1;

    uint64_t result = 1;
    if (exp & 1)
        result = base;
    exp >>= 1;

    while (exp) {
        base *= base;
        if (exp & 1)
            result *= base;
        exp >>= 1;
    }
    return result;
}

// assumes both A and B are in [0, M) range
constexpr uint128_t add_mod(uint128_t a, uint128_t b, uint128_t m) {
    return (b >= m - a) ? (a + b - m) : (a + b);
}

// TODO this seems to be buggy!
// assumes both A and B are in [0, M) range
constexpr uint128_t mul_mod(uint128_t a, uint128_t b, uint128_t m) {
    if (a == 0 || b == 0)
        return 0;
    if (a == 1)
        return b;
    if (b == 1)
        return a;
    if (a < UINT128_MAX / b)
        return (a * b) % m;

    uint128_t result = 0;
    while (a && b) {
        if (a < b)
            std::swap(a, b);
        if (b & 1)
            result = add_mod(result, a, m);
        a = add_mod(a, a, m);
        b >>= 1;
    }
    return result;
}

// returns (A ** N) mod M
// assumes A is in [0, M) range
constexpr uint64_t pow_mod(uint64_t a, uint64_t n, uint64_t m) {
    uint64_t b = 1;
    while (n) {
        if (n & 1)
            b = __divq_mod(__mulq(b, a), m);
        n >>= 1;
        a = __divq_mod(__mulq(a, a), m);
    }
    return b;
}

constexpr int num_bits(std::unsigned_integral auto a) { return sizeof(a) * 8 - std::countl_zero(a); }

constexpr bool can_overflow_with_carry(const uint64_t a, const uint64_t b) {
    return (b == UINT64_MAX) || (a > UINT64_MAX - b - 1);
}

#define ALGEBRA_SHIFT_OP(CLASS) \
constexpr CLASS& operator>>=(CLASS& a, int64_t b) { a <<= -b; return a; } \
constexpr CLASS operator>>(CLASS a, int64_t b) { a <<= -b; return a; } \
constexpr CLASS operator<<(CLASS a, int64_t b) { a <<= b; return a; } \
template<typename T> requires (std_int<T> && !std::same_as<T, int64_t>) \
constexpr CLASS& operator<<=(CLASS& a, T b) { a <<= (int64_t)b; return a; } \
constexpr CLASS& operator>>=(CLASS& a, std_int auto b) { a >>= (int64_t)b; return a; } \
constexpr CLASS operator<<(const CLASS& a, std_int auto b) { return a << (int64_t)b; } \
constexpr CLASS operator>>(const CLASS& a, std_int auto b) { return a >> (int64_t)b; }

struct bit_range {
    uint64_t min, max;
    constexpr bit_range(uint64_t min, uint64_t max) : min(min), max(max) { }
    constexpr bit_range(uint64_t a) { min = max = a; }
};

constexpr bit_range operator*(bit_range a, bit_range b) { return {a.min + b.min - 1, a.max + b.max}; }
constexpr bit_range operator+(bit_range a, bit_range b) { return {std::min(a.min, b.min), std::max(a.max, b.max) + 1}; }
constexpr bool operator<(bit_range a, bit_range b) { return a.max < b.min; }

template<typename A, typename B>
using larger_type = typename std::conditional_t<sizeof(A) >= sizeof(B), A, B>;

template<std::unsigned_integral A, std::unsigned_integral B>
constexpr auto min(const A& a, const B& b) -> larger_type<A, B> { return (a < b) ? a : b; }
template<std::unsigned_integral A, std::unsigned_integral B>
constexpr auto max(const A& a, const B& b) -> larger_type<A, B> { return (a > b) ? a : b; }

template<std::signed_integral A, std::unsigned_integral B>
constexpr auto min(const A& a, const B& b) -> larger_type<A, B> { return (a < b) ? a : b; }
template<std::signed_integral A, std::signed_integral B>
constexpr auto max(const A& a, const B& b) -> larger_type<A, B> { return (a > b) ? a : b; }

}
