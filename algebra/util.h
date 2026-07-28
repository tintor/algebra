#pragma once
#include <source_location>
#include <format>
#include <print>
#include <stdexcept>
#include <bit>
#include <cstdint>
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

// int128_t and uint128_t come from types.h

// x86-64 divq is a 128/64 -> 64 division, which the compiler does not generate on its own.
// Define ALGEBRA_NO_ASM to use the portable implementation instead.
#if defined(__x86_64__) && !defined(__ILP32__) && !defined(ALGEBRA_NO_ASM)
#define ALGEBRA_X86_DIVQ 1
#endif

#ifdef ALGEBRA_X86_DIVQ
constexpr void __divq(uint128_t a, uint64_t b, uint64_t& q, uint64_t& r) {
    if consteval { // inline assembly is not allowed in a constant expression
        q = a / b;
        r = a % b;
        return;
    }
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

constexpr uint64_t __divq(uint128_t a, uint64_t b) {
    if consteval { // inline assembly is not allowed in a constant expression
        return a / b;
    }
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

constexpr uint64_t __divq_mod(uint128_t a, uint64_t b) {
    if consteval { // inline assembly is not allowed in a constant expression
        return a % b;
    }
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
// Note: unlike divq, these do not trap when the quotient does not fit into 64 bits.
constexpr void __divq(uint128_t a, uint64_t b, uint64_t& q, uint64_t& r) {
    q = a / b;
    r = a % b;
}

constexpr uint64_t __divq(uint128_t a, uint64_t b) { return a / b; }

constexpr uint64_t __divq_mod(uint128_t a, uint64_t b) { return a % b; }
#endif

// Constrained so a bignum class never converts into this by accident; without it a missing
// pow(integer, ...) overload resolves here and truncates the operand to 64 bits.
template<std_unsigned_int T>
constexpr uint64_t pow(T base, unsigned exp) {
    if (base == 2)
        return (exp < 64) ? (static_cast<uint64_t>(1) << exp) : 0; // 0 matches the wrap-around of the loop below
    if (exp == 0)
        return 1;

    // the squaring below has to happen at the width of the result: in a narrower T it would wrap
    // at that width instead, and a wider T is truncated to the same value modulo 2**64 anyway
    uint64_t b = static_cast<uint64_t>(base);
    uint64_t result = 1;
    if (exp & 1)
        result = b;
    exp >>= 1;

    while (exp) {
        b *= b;
        if (exp & 1)
            result *= b;
        exp >>= 1;
    }
    return result;
}

// assumes both A and B are in [0, M) range
constexpr uint128_t add_mod(uint128_t a, uint128_t b, uint128_t m) {
    return (b >= m - a) ? (a + b - m) : (a + b);
}

// assumes both A and B are in [0, M) range
// (swapping a and b inside the loop is safe: the invariant is result + a*b, which is symmetric)
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

constexpr uint64_t hash_fn_64bit(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdllu;
    return k;
}

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

// min() and max() over two builtin integers of the same signedness, widening to whichever type
// holds both. Mixed signedness is deliberately not offered: the comparison would go through the
// usual arithmetic conversions, where a negative value turns into a huge unsigned one, and there
// is no larger_type that can hold both operands either. Use std::cmp_less and pick the type at
// the call site for that case.
template<std::unsigned_integral A, std::unsigned_integral B>
constexpr auto min(const A& a, const B& b) -> larger_type<A, B> { return (a < b) ? a : b; }
template<std::unsigned_integral A, std::unsigned_integral B>
constexpr auto max(const A& a, const B& b) -> larger_type<A, B> { return (a > b) ? a : b; }

template<std::signed_integral A, std::signed_integral B>
constexpr auto min(const A& a, const B& b) -> larger_type<A, B> { return (a < b) ? a : b; }
template<std::signed_integral A, std::signed_integral B>
constexpr auto max(const A& a, const B& b) -> larger_type<A, B> { return (a > b) ? a : b; }

}
