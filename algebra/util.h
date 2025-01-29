#pragma once
#include <source_location>
#include <format>
#include <print>
#include <stdexcept>
#include <bit>

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

constexpr bool __less(const uint64_t* a, const int A, const uint64_t* b, const int B) {
    if (A > B)
        return false;
    if (A < B)
        return true;
    for (auto i = A; i-- > 0;) {
        if (a[i] > b[i])
            return false;
        if (a[i] < b[i])
            return true;
    }
    return false;
}

constexpr bool can_overflow_with_carry(const uint64_t a, const uint64_t b) {
    return (b == UINT64_MAX) || (a > UINT64_MAX - b - 1);
}

constexpr size_t add_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    if (A > B)
        return (a[A - 1] == UINT64_MAX) ? A + 1 : A;
    if (A < B)
        return (b[B - 1] == UINT64_MAX) ? B + 1 : B;
    return can_overflow_with_carry(a[A - 1], b[A - 1]) ? A + 1 : A;
}

constexpr size_t mul_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    return (A && B) ? A + B - 1 + (127 - std::countl_zero(a[A - 1]) - std::countl_zero(b[B - 1])) / 64 : 0;
}

constexpr size_t div_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    return (A >= B) ? 1 + 64 * (A - B) + std::countl_zero(b[B - 1]) - std::countl_zero(a[A - 1]) : 0;
}

constexpr size_t num_bits(const uint64_t* a, const size_t A) { return A ? 64 * A - std::countl_zero(a[A - 1]) : 0; }

// `q` needs to have capacity of at least A + B!
// supports q == a
// b != q
constexpr void __mul(const uint64_t* a, const int A, const uint64_t* b, const int B, uint64_t* q, int& Q, bool init = true) {
    // TODO if a != q it might be possible to swap a and b depending on their sizes! maybe have A be smaller than B

    Q = A + B;
    if (init)
        for (int i = A; i < Q; i++)
            q[i] = 0;
    for (int i = A; i-- > 0;) {
        const uint64_t w = a[i];
        if (w == 0)
            continue;

        uint64_t* qi = q + i;
        auto acc = __mulq(w, *b);
        *qi = acc;
        acc >>= 64;
        int j = 1;

        while (j < B) {
            acc += __mulq(w, b[j]);
            acc += qi[j];
            qi[j] = acc;
            acc >>= 64;
            j += 1;
        }

        if (acc) {
            acc += qi[j];
            qi[j] = acc;
            acc >>= 64;
            if (acc) {
                j += 1;
                acc += qi[j];
                qi[j] = acc;
            }
        }
    }
    while (Q && !q[Q - 1])
        Q -= 1;
}

// shift must be >= 0
constexpr void __add(uint64_t* a, int& A, const uint64_t* b, const int B, int shift = 0) {
    while (A < B + shift)
        a[A++] = 0;

    uint128_t acc = 0;
    for (int i = 0; i < B; ++i) {
        acc += a[shift + i];
        acc += b[i];
        a[shift + i] = acc;
        acc >>= 64;
    }
    for (int i = shift + B; i < A; ++i) {
        acc += a[i];
        a[i] = acc;
        acc >>= 64;
    }
    if (acc)
        a[A++] = acc;
}

// assuming a >= b
constexpr void __sub(uint64_t* a, int& A, const uint64_t* b, const int B) {
    uint64_t borrow = 0;
    for (int i = 0; i < B; ++i) {
        int128_t diff = (int128_t)a[i] - b[i] - borrow;
        if (diff < 0) {
            a[i] = diff + ((uint128_t)1 << 64);
            borrow = 1;
        } else {
            a[i] = diff;
            borrow = 0;
        }
    }
    for (int i = B; borrow && i < A; ++i) {
        int128_t diff = (int128_t)a[i] - borrow;
        if (diff < 0) {
            a[i] = diff + ((uint128_t)1 << 64);
            borrow = 1;
        } else {
            a[i] = diff;
            borrow = 0;
        }
    }
    while (A && !a[A - 1])
        A -= 1;
}

constexpr int KARATSUBA_LIMIT = 32;

// long mul  16384               takes 235/230ms*
// karatsuba 16384 with 8 limit  takes ****ms (fails, likely due to small W buffer size)
// karatsuba 16384 with 16 limit takes  33ms* (regresses size 16 compared to long mul)
// karatsuba 16384 with 32 limit takes  26ms*
// karatsuba 16384 with 64 limit takes  26ms*

// assuming a.size >= b.size
constexpr void __mul_karatsuba_rec(const uint64_t* a, int A, const uint64_t* b, int B, uint64_t* q, int& Q, uint64_t* w, const uint64_t* we) {
    if (A < B) {
        std::swap(a, b);
        std::swap(A, B);
    }
    if (A < KARATSUBA_LIMIT) {
        __mul(a, A, b, B, q, Q);
        return;
    }

    const int m = (A + 1) / 2; // m >= 2

    // AA and BB are stored in Q, while R and P are stored in W
    const uint64_t* a1 = a + m;
    const int A1 = A - m;

    int R;
    if (B <= m) {
        auto r = w;
        w += A1 + B;
        __mul_karatsuba_rec(a1, A1, b, B, r, R, w, we); // r = a1 * b
        __mul_karatsuba_rec(a, m, b, B, q, Q, w, we); // q = a0 * b
        __add(q, Q, r, R, m); // q = a * b
    } else {
        uint64_t* aa = q;
        std::copy(a, a + m, aa);
        int AA = m;
        __add(aa, AA, a1, A1, 0); // aa = a0 + a1

        const uint64_t* b1 = b + m;
        const int B1 = B - m;

        uint64_t* bb = aa + AA;
        int BB = m;
        std::copy(b, b + m, bb);
        __add(bb, BB, b1, B1, 0); // bb = b0 + b1

        auto r = w;
        w += AA + BB + 1;
        __mul_karatsuba_rec(aa, AA, bb, BB, r, R, w, we); // r = aa * bb
        // TODO ^ How is this working? AA and BB are stored in Q and nested __mul_karatsuba_rec call will overwrite them? Tests are pasing with 256 words.

        uint64_t* p = w;
        w += A1 + B1;
        int P;
        __mul_karatsuba_rec(a1, A1, b1, B1, p, P, w, we);
        __mul_karatsuba_rec(a, m, b, m, q, Q, w, we);

        __sub(r, R, p, P);
        __sub(r, R, q, Q);
        __add(q, Q, p, P, m * 2);
        __add(q, Q, r, R, m);
    }
}

// `Q` must be initialized to buffer capacity
// `q` must be all zero
constexpr void __mini_mul(const uint64_t* a, const int A, const uint64_t* b, const int B, uint64_t* q, int& Q) {
    for (int i = A; i-- > 0;) {
        const uint64_t w = a[i];
        if (w == 0)
            continue;

        uint64_t* qi = q + i;
        auto acc = __mulq(w, *b);
        *qi = acc;
        acc >>= 64;
        int j = 1;

        while (j < B) {
            acc += __mulq(w, b[j]);
            acc += qi[j];
            qi[j] = acc;
            acc >>= 64;
            j += 1;
        }

        if (acc) {
            acc += qi[j];
            qi[j] = acc;
            acc >>= 64;
            if (acc) {
                j += 1;
                acc += qi[j];
                qi[j] = acc;
            }
        }
    }
    while (Q && !q[Q - 1])
        Q -= 1;
}

constexpr bool is_power_of_two(const uint64_t* a, const int A) {
    if (A == 0)
        return false;
    auto e = a + A - 1;
    if (*e & (*e - 1))
        return false;
    while (--e >= a)
        if (*e)
            return false;
    return true;
}

// supports a == q
constexpr uint64_t __div(const uint64_t* a, const int A, uint64_t b, uint64_t* q, int& Q) {
    Check(b != 0, "division by zero");
    Q = A;
    uint128_t acc = 0;
    for (auto i = A; i-- > 0;) {
        acc <<= 64;
        acc |= a[i];
        uint64_t r;
        __divq(acc, b, q[i], r);
        acc = r;
    }
    while (Q > 0 && q[Q - 1] == 0)
        Q -= 1;
    return acc;
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

}
