#pragma once
#include "algebra/integer_backend.h"
#include <format>
#include <string_view>
#include <bit>
#include <stdexcept>
#include <algorithm>
#include <print>
#include <vector>
#include <source_location>

namespace algebra {

template<typename T> struct IsNumberClass : std::false_type {};
template<typename T> concept __ncsi = IsNumberClass<T>::value || std_int<T>;

constexpr bool operator>(const __ncsi auto& a, const __ncsi auto& b) { return b < a; }
constexpr bool operator>=(const __ncsi auto& a, const __ncsi auto& b) { return !(a < b); }
constexpr bool operator<=(const __ncsi auto& a, const __ncsi auto& b) { return !(b < a); }
constexpr bool operator!=(const __ncsi auto& a, const __ncsi auto& b) { return !(a == b); }

constexpr auto operator+(std_int auto a, const __ncsi auto& b) { return b + a; }
constexpr auto operator*(std_int auto a, const __ncsi auto& b) { return b * a; }
constexpr bool operator==(std_int auto a, const __ncsi auto& b) { return b == a; }

struct natural;
template<> struct IsNumberClass<natural> : std::true_type {};

constexpr uint128_t __mulq(uint64_t a, uint64_t b) { return uint128_t(a) * b; }

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

// returns (a ** n) mod p
constexpr uint64_t pow_mod(uint64_t a, uint64_t n, uint64_t p) {
    uint64_t b = 1;
    a %= p;
        while (n) {
        if (n & 1)
            b = (static_cast<__uint128_t>(b) * a) % p;
        n >>= 1;
        a = (static_cast<__uint128_t>(a) * a) % p;
    }
    return b;
}

// TODO test cases for operator float()
struct natural {
    using size_type = integer_backend::size_type;

    integer_backend words;

    constexpr natural() {}
    constexpr natural(std_int auto a) : words(a) {
        if (a < 0)
            throw std::runtime_error("assigning negative number to natural");
    }
    constexpr natural(natural&& o) : words(std::move(o.words)) { }
    constexpr natural(const natural& o) : words(o.words) { }

    constexpr void set_zero() { words.set_zero(); }
    constexpr void operator=(std_int auto a) { words = a; }
    constexpr void operator=(natural&& o) { words = std::move(o.words); }
    constexpr void operator=(const natural& o) { words = o.words; }

    constexpr bool is_one() const { return words[0] == 1 && words.size() == 1; }
    constexpr bool is_even() const { return (words[0] % 2) == 0; }
    constexpr bool is_odd() const { return words[0] % 2; }

    constexpr bool is_uint8() const;
    constexpr bool is_uint16() const;
    constexpr bool is_uint32() const { return words.size() <= 1 && words[0] <= UINT32_MAX; }
    constexpr bool is_uint64() const { return words.size() <= 1; }
    constexpr bool is_uint128() const { return words.size() <= 2; }

    constexpr operator uint8_t() const {
        if (!is_uint8())
            throw std::runtime_error("cast overflow");
        return words[0];
    }
    constexpr operator uint16_t() const {
        if (!is_uint16())
            throw std::runtime_error("cast overflow");
        return words[0];
    }
    constexpr operator uint32_t() const {
        if (!is_uint32())
            throw std::runtime_error("cast overflow");
        return words[0];
    }
    constexpr operator unsigned long() const {
        static_assert(sizeof(unsigned long) == 8);
        if (!is_uint64())
            throw std::runtime_error("cast overflow");
        return words[0];
    }
    constexpr operator unsigned long long() const {
        static_assert(sizeof(unsigned long long) == 8);
        if (!is_uint64())
            throw std::runtime_error("cast overflow");
        return words[0];
    }

    constexpr operator uint128_t() const {
        if (!is_uint128())
            throw std::runtime_error("cast overflow");
        uint128_t a = words[0];
        if (words.size() >= 2)
            a |= static_cast<uint128_t>(words[1]) << 64;
        return a;
    }

    constexpr size_t num_trailing_zeros() const {
        size_t a = 0;
        for (size_type i = 0; i < words.size(); i++) {
            if (words[i])
                return a + std::countr_zero(words[i]);
            a += 64;
        }
        return a;
    }

    constexpr void __add(uint64_t b, size_type i) {
        while (true) {
            if (b == 0)
                return;
            if (i == words.size()) {
                words.push_back(b);
                return;
            }
            uint128_t acc = (uint128_t)words[i] + b;
            words[i] = acc;
            b = acc >> 64;
            i += 1;
        }
    }

    constexpr natural& operator+=(uint64_t b) {
        __add(b, 0);
        return *this;
    }

    constexpr natural& operator+=(uint128_t b) {
        __add(static_cast<uint64_t>(b), 0);
        if (words.size() == 0)
            words.push_back(0);
        __add(static_cast<uint64_t>(b >> 64), 1);
        return *this;
    }

    constexpr natural& operator+=(const natural& b) {
        // TODO predict max possible n taking into account carry!
        if (b.words.size() > words.size()) {
            words.reserve(b.words.size());
            while (words.size() < b.words.size())
                words.push_back(0);
        }
        uint128_t acc = 0;
        for (size_type i = 0; i < words.size(); ++i) {
            acc += words[i];
            if (i < b.words.size())
                acc += b.words[i];
            words[i] = acc;
            acc >>= 64;
        }
        if (acc)
            words.push_back(acc);
        return *this;
    }

    constexpr natural& operator-=(uint64_t b) {
        if (!words.allocated()) {
            if (words[0] < b)
                throw std::runtime_error("natural can't be negative");
            words[0] -= b;
            if (words[0] == 0)
                words.pop_back();
            return *this;
        }

        if (words[0] < b && words.size() < 2)
            throw std::runtime_error("natural can't be negative");

        if (words[0] > b) {
            words[0] -= b;
            return *this;
        }

        words[0] -= b;
        for (size_type i = 1; i < words.size(); ++i)
            if (words[i]--)
                break;
        if (words.back() == 0)
            words.pop_back();
        return *this;
    }

    // TODO fix same issue as above
    constexpr natural& operator-=(uint128_t b) {
        if (b <= UINT64_MAX)
            return operator-=(static_cast<uint64_t>(b));

        // subtract low part of b
        if (words[0] < uint64_t(b) && words.size() < 2)
            throw std::runtime_error("natural can't be negative");

        words[0] -= uint64_t(b);
        for (size_type i = 1; i < words.size(); ++i)
            if (words[i]--)
                break;
        if (words.back() == 0)
            words.pop_back();

        // subtract high part of b
        if (words.size() < 3 && words[1] < uint64_t(b >> 64))
            throw std::runtime_error("natural can't be negative");

        words[1] -= uint64_t(b >> 64);
        for (size_type i = 2; i < words.size(); ++i)
            if (words[i]--)
                break;
        if (words.back() == 0)
            words.pop_back();

        return *this;
    }

    constexpr natural& operator-=(const natural& b) {
        uint64_t borrow = 0;
        for (size_type i = 0; i < b.words.size(); ++i) {
            int128_t diff = (int128_t)words[i] - b.words[i] - borrow;
            if (diff < 0) {
                words[i] = diff + ((uint128_t)1 << 64);
                borrow = 1;
            } else {
                words[i] = diff;
                borrow = 0;
            }
        }
        for (size_t i = b.words.size(); borrow && i < words.size(); ++i) {
            int128_t diff = (int128_t)words[i] - borrow;
            if (diff < 0) {
                words[i] = diff + ((uint128_t)1 << 64);
                borrow = 1;
            } else {
                words[i] = diff;
                borrow = 0;
            }
        }
        if (borrow)
            throw std::format_error("natural subtraction out of domain");
        words.normalize();
        return *this;
    }

    constexpr void mul_add(uint64_t a, uint64_t carry) {
        if (a == 0) {
            words.set_zero();
            if (carry)
                words.push_back(carry);
            return;
        }
        for (size_type i = 0; i < words.size(); ++i) {
            uint128_t acc = (uint128_t)words[i] * a + carry;
            words[i] = acc;
            carry = acc >> 64;
        }
        if (carry)
            words.push_back(carry);
    }
    constexpr natural& operator*=(uint64_t b) { mul_add(b, 0); return *this; }

    constexpr uint64_t operator%(std_int auto b) const {
        static_assert(sizeof(b) <= 8);
        if (b <= 0)
            throw std::runtime_error((b == 0) ? "division by zero" : "division of natural by negative number");
        uint128_t acc = 0;
        for (size_type i = words.size(); i-- > 0;) {
            acc <<= 64;
            acc |= words[i];
            acc %= static_cast<uint64_t>(b);
        }
        return acc;
    }

    constexpr uint128_t operator%(const uint128_t b) const {
        if (b <= UINT64_MAX)
            return operator%(static_cast<uint64_t>(b));
        if (b == 0)
            throw std::runtime_error("division by zero");

        // compute m = (2**128) % b
        uint128_t m = UINT128_MAX % b;
        m += 1;
        if (m == b)
            m = 0;

        uint128_t res = 0;
        for (size_type i = 0; i < words.size(); i += 2) {
            res = mul_mod(res, m, b);
            uint128_t acc = words[i];
            if (i + 1 < words.size())
                acc |= static_cast<uint128_t>(words[i + 1]) << 64;
            res = add_mod(res, acc % b, b);
        }
        return res;
    }

    constexpr uint64_t mod2() const { return words[0] % 2; }

    constexpr uint64_t mod3() const {
        uint64_t acc = 0;
        for (size_type i = words.size(); i-- > 0;)
            acc += words[i] % 3;
        return acc % 3;
    }

    // 1180591620717411302909 = 30101795611 * 39219973319 in 294820 ms
    constexpr uint64_t mod4() const { return words[0] % 4; }

    constexpr uint64_t mod5() const {
        uint64_t acc = 0;
        if (words.size() > UINT64_MAX / 4) {
            for (size_type i = words.size(); i-- > 0;) {
                acc += words[i] % 5;
                if ((i % 65536) == 0)
                    acc %= 5;
            }
        } else {
            // (2**64) mod 5 == 1
            for (size_type i = words.size(); i-- > 0;)
                acc += words[i] % 5;
        }
        return acc % 5;
    }

    constexpr natural& operator%=(std_int auto b) { *this = operator%(b); return *this; }

    constexpr natural(std::string_view s, unsigned base = 10);
    constexpr natural(const char* s, uint32_t base = 10) : natural(std::string_view(s), base) {}

    constexpr void swap(natural& o) {
        words.swap(o.words);
    }

    constexpr size_type str_size_upper_bound(uint32_t base = 10) const;
    constexpr size_type str(char* buffer, int buffer_size, uint32_t base = 10, bool upper = true) const;
    constexpr std::string str(uint32_t base = 10, bool upper = true) const {
        std::string s;
        s.resize(str_size_upper_bound(base));
        s.resize(str(s.data(), s.size(), base, upper));
        return s;
    }
    constexpr std::string hex() const { return str(16); }

    constexpr size_t num_bits() const { return words.size() ? words.size() * 64 - std::countl_zero(words.back()) : 0; }
    constexpr bool bit(size_t i) const {
        size_t w = i / 64;
        size_t b = i % 64;
        return w < words.size() && (words[w] & (uint64_t(1) << b));
    }

    constexpr size_t popcount() const {
        size_t c = 0;
        for (size_type i = 0; i < words.size(); i++)
            c += std::popcount(words[i]);
        return c;
    }

    constexpr size_t size_of() const { return words.size() * 8; }

    constexpr operator bool() const { return words.size(); }

    constexpr natural& operator++() { return operator+=(uint64_t(1)); }
    constexpr natural& operator--() { return operator-=(uint64_t(1)); }

    constexpr natural operator++(int) { natural a = *this; operator++(); return a; }
    constexpr natural operator--(int) { natural a = *this; operator--(); return a; }

    template<std::floating_point T>
    constexpr operator T() const;
};

constexpr int num_bits(std::unsigned_integral auto a) { return sizeof(a) * 8 - std::countl_zero(a); }

constexpr natural operator-(const natural& a) {
    if (a.words.size() == 0)
        return 0;
    throw std::runtime_error("natural can't be negative");
}

constexpr natural operator+(const natural& a, const natural& b) { natural c = a; return c += b; }
constexpr natural operator+(natural a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a += static_cast<uint64_t>(b);
    static_assert(sizeof(b) <= 16);
    return a += static_cast<uint128_t>(b);
}
constexpr natural operator+(natural a, const std_signed_int auto b) {
    if (b < 0)
        return a - abs_unsigned(b);
    return a + make_unsigned(b);
}
constexpr natural operator+(const std_int auto a, natural b) { return std::move(b) + a; }

constexpr natural& operator+=(natural& a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a += static_cast<uint64_t>(b);
    static_assert(sizeof(b) <= 16);
    return a += static_cast<uint128_t>(b);
}

constexpr natural& operator+=(natural& a, const std_signed_int auto b) {
    if (b < 0)
        return a -= abs_unsigned(b);
    return a += make_unsigned(b);
}

constexpr natural operator-(const natural& a, const natural& b) { natural c = a; return c -= b; }

constexpr natural operator-(natural a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a -= uint64_t(b);
    if (a < b)
        throw std::runtime_error("natural can't be negative");
    static_assert(sizeof(b) <= 16);
    return a -= static_cast<uint128_t>(b);
}

constexpr natural operator-(const std_unsigned_int auto a, natural b) {
    if (a <= UINT64_MAX) {
        if (b.words.size() > 1 || uint64_t(a) < b.words[0])
            throw std::runtime_error("natural can't be negative");
        return uint64_t(a) - b.words[0];
    }
    if (a > b)
        throw std::runtime_error("natural can't be negative");
    return a - static_cast<uint128_t>(b);
}

constexpr natural operator-(natural a, const std_signed_int auto b) {
    if (b < 0)
        return a + (~make_unsigned(b) + 1);
    return a - make_unsigned(b);
}

constexpr natural operator-(const std_signed_int auto a, natural b) {
    if (a < 0)
        throw std::runtime_error("natural can't be negative");
    return make_unsigned(a) - b;
}

constexpr natural& operator-=(natural& a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a -= static_cast<uint64_t>(b);
    static_assert(sizeof(b) <= 16);
    return a -= static_cast<uint128_t>(b);
}

constexpr natural& operator-=(natural& a, const std_signed_int auto b) {
    if (b < 0)
        return a += ~make_unsigned(b) + 1;
    return a -= make_unsigned(b);
}

template<std_int Size>
constexpr bool __less(const uint64_t* a, const Size A, const uint64_t* b, const Size B) {
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

constexpr bool operator<(const natural& a, const natural& b) { return __less(a.words.data(), a.words.size(), b.words.data(), b.words.size()); }

constexpr bool operator<(const natural& a, const std_unsigned_int auto b) { return a.words[0] < b && a.words.size() <= 1; }
constexpr bool operator<(const std_unsigned_int auto a, const natural& b) { return a < b.words[0] || b.words.size() > 1; }

constexpr bool operator<(const natural& a, const std_signed_int auto b) { return b >= 0 && a < static_cast<uint64_t>(b); }
constexpr bool operator<(const std_signed_int auto a, const natural b) { return a < 0 || static_cast<uint64_t>(a) < b; }

constexpr bool operator<(const natural& a, const uint128_t b) {
    if (b <= UINT64_MAX)
        return a < static_cast<uint64_t>(b);
    if (a.words.size() < 2)
        return true;
    if (a.words.size() > 2)
        return false;
    if (a.words[1] < static_cast<uint64_t>(b >> 64))
        return true;
    return a.words[1] == static_cast<uint64_t>(b >> 64) && a.words[0] < static_cast<uint64_t>(b);
}

constexpr bool operator<(const uint128_t a, const natural& b) {
    if (a <= UINT64_MAX)
        return static_cast<uint64_t>(b) < b;
    if (b.words.size() < 2)
        return false;
    if (b.words.size() > 2)
        return true;
    if (static_cast<uint64_t>(a >> 64) < b.words[1])
        return true;
    return static_cast<uint64_t>(a >> 64) == b.words[1] && static_cast<uint64_t>(a) < b.words[0];
}

constexpr bool operator==(const natural& a, const natural& b) {
    if (a.words.size() != b.words.size())
        return false;
    for (auto i = a.words.size(); i-- > 0;)
        if (a.words[i] != b.words[i])
            return false;
    return true;
}

constexpr bool operator==(const natural& a, const std_unsigned_int auto b) {
    if constexpr (sizeof(b) <= 8)
        return a.words[0] == b && a.words.size() <= 1;
    if (b <= UINT64_MAX)
        return a.words[0] == b && a.words.size() <= 1;
    return a.words.size() == 2 && a.words[1] == uint64_t(b >> 64) && a.words[0] == uint64_t(b);
}
constexpr bool operator==(const natural& a, const std_signed_int auto b) { return b >= 0 && a == make_unsigned(b); }

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
constexpr void __add(uint64_t* a, int& A, const uint64_t* b, const int B, int shift) {
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

constexpr void Check(bool value, std::source_location loc = std::source_location::current()) {
    if (!value)
        throw std::runtime_error(std::format("Check failed at {}:{} in {}", loc.file_name(), loc.line(), loc.function_name()));
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

constexpr size_t num_bits(const uint64_t* a, const size_t A) { return A ? 64 * A - std::countl_zero(a[A - 1]) : 0; }

constexpr size_t mul_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    return (A && B) ? A + B - 1 + (127 - std::countl_zero(a[A - 1]) - std::countl_zero(b[B - 1])) / 64 : 0;
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

constexpr bool is_power_of_two(const natural& a) {
    if (a.words.empty())
        return false;
    auto w = a.words.data();
    auto e = w + a.words.size() - 1;
    if (*e & (*e - 1))
        return false;
    while (e >= w)
        if (*e--)
            return false;
    return true;
}

constexpr natural& operator<<=(natural& a, int64_t b);

// no support for &a == &q
constexpr void mul_karatsuba(const natural& a, const natural& b, natural& q) {
    auto A = a.words.size();
    auto B = b.words.size();
    // TODO this case will disappear when a.words.empty() is removed!
    if (A == 0 || B == 0) {
        q.set_zero();
        return;
    }
    if (A == 1) {
        if (B == 1) {
            q = __mulq(a.words[0], b.words[0]);
            return;
        }
        if (a.words[0] == 1) {
            q = b;
            return;
        }
    }
    if (b.words[0] == 1 && B == 1) {
        q = a;
        return;
    }

    if (is_power_of_two(a)) {
        const size_t z = (A - 1) * 64 + std::countr_zero(a.words[A - 1]); // = a.num_trailing_zeros() but O(1)
        const size_t bits = b.num_bits() + z;
        const size_t words = (bits + 63) / 64;
        q.words.reset(words, /*init*/false); // preallocates memory!
        // TODO this can be done directly without moving
        q = b;
        q <<= z;
        return;
    }
    if (is_power_of_two(b)) {
        const size_t z = (B - 1) * 64 + std::countr_zero(a.words[B - 1]); // = b.num_trailing_zeros() but O(1)
        const size_t bits = a.num_bits() + z;
        const size_t words = (bits + 63) / 64;
        q.words.reset(words, /*init*/false); // preallocates memory!
        // TODO this can be done directly without moving
        q = a;
        q <<= z;
        return;
    }

    int Q = mul_max_size(a.words.data(), A, b.words.data(), B);
    q.words.reset(Q);
    if (std::min(A, B) <= 2 || std::max(A, B) < KARATSUBA_LIMIT) {
        __mini_mul(a.words.data(), A, b.words.data(), B, q.words.data(), Q);
    } else {
        const int W = 4 * std::max(A, B);
        if (W <= 1024) {
            uint64_t w[1024];
            __mul_karatsuba_rec(a.words.data(), A, b.words.data(), B, q.words.data(), Q, w, w + W);
        } else {
            auto w = new uint64_t[W];
            __mul_karatsuba_rec(a.words.data(), A, b.words.data(), B, q.words.data(), Q, w, w + W);
            delete[] w;
        }
    }
    q.words.resize(Q);
}

constexpr natural mul_karatsuba(const natural& a, const natural& b) {
    natural q;
    mul_karatsuba(a, b, q);
    return q;
}

// supports &a == &q
constexpr void __mul(const natural& a, const natural& b, natural& q) {
    natural::size_type Q;
    if (&a != &q)
        q.set_zero();
    auto A = a.words.size();
    q.words.resize(A + b.words.size());
    __mul(a.words.data(), A, b.words.data(), b.words.size(), q.words.data(), Q, /*init*/false);
    q.words.resize(Q);
}

constexpr void __mul(const natural& a, uint128_t b, natural& out) {
    if (&a != &out)
        out.set_zero();
    auto as = a.words.size();
    out.words.resize(a.words.size() + 2);

    for (auto i = as; i-- > 0;) {
        const auto aw = a.words[i];
        if (aw == 0)
            continue;

        uint64_t* ow = &out.words[i];
        *ow = 0;

        auto acc = __mulq(aw, b);
        acc += *ow;
        *ow = acc;
        acc >>= 64;

        acc += __mulq(aw, b >> 64);
        acc += *++ow;
        *ow = acc;
        acc >>= 64;

        if (acc) {
            acc += *++ow;
            *ow = acc;
            acc >>= 64;

            if (acc) {
                acc += *++ow;
                *ow = acc;
            }
        }
    }
    out.words.normalize();
}

// supports &a == &out
constexpr void __mul(const natural& a, uint64_t b, uint64_t carry, natural& out) {
    if (b == 0) {
        out.words.set_zero();
        if (carry)
            out.words.push_back(carry);
        return;
    }
    if (&a != &out)
        out.words.reset(a.words.size(), /*initialize*/false);
    for (natural::size_type i = 0; i < a.words.size(); ++i) {
        auto acc = __mulq(a.words[i], b) + carry;
        out.words[i] = acc;
        carry = acc >> 64;
    }
    if (carry)
        out.words.push_back(carry);
}

// assumes a.words.size() >= 2
constexpr void __square(natural& a) {
    if (a.words.size() == 2) {
        auto carry = __mulq(a.words[0], a.words[0]);
        uint64_t b0 = carry;
        uint64_t b1 = carry >> 64;

        uint128_t pq = __mulq(a.words[0], a.words[1]);
        carry = b1 + pq; // can't use pq*2 here due to dword overflow
        b1 = carry;
        uint64_t b2 = carry >> 64;

        carry = b1 + pq;
        b1 = carry;
        carry >>= 64;
        carry += b2;
        b2 = carry;

        uint64_t b3 = carry >> 64;
        carry = b2 + __mulq(a.words[1], a.words[1]);
        b2 = carry;
        b3 += carry >> 64;

        if (b3) {
            a.words.reset(4, /*initialize*/false);
            a.words[0] = b0;
            a.words[1] = b1;
            a.words[2] = b2;
            a.words[3] = b3;
            return;
        }
        if (b2) {
            a.words.reset(3, /*initialize*/false);
            a.words[0] = b0;
            a.words[1] = b1;
            a.words[2] = b2;
            return;
        }
        a.words.reset(2, /*initialize*/false);
        a.words[0] = b0;
        a.words[1] = b1;
        return;
    }

    auto n = a.words.size();
    a.words.resize(n << 1);
    for (auto k = (n - 1) << 1; k >= 0; k--) {
        auto w = a.words[k];
        a.words[k] = 0;
        auto i_min = std::max<natural::size_type>(0, k - n + 1);
        auto i_max = std::min<natural::size_type>(n - 1, k);
        for (auto i = i_min; i <= i_max; i++) {
            auto j = k - i;
            const auto ai = (i < k) ? a.words[i] : w;
            const auto aj = (j < k) ? a.words[j] : w;
            uint128_t carry = __mulq(ai, aj);
            for (auto p = k; carry; p++) {
                carry += a.words[p];
                a.words[p] = carry;
                carry >>= 64;
            }
        }
    }
    a.words.normalize();
}

constexpr void square(natural& a) {
    if (a.words.size() == 0)
        return;
    if (a.words.size() == 1) {
        uint128_t p = __mulq(a.words[0], a.words[0]);

        a.words.reset_one_without_init();
        a.words[0] = p;

        uint64_t high = p >> 64;
        if (high)
            a.words.push_back(high);
        return;
    }
    __square(a);
}

constexpr void mul(const natural& a, const natural& b, natural& out) {
    if (a.words.size() == 0 || b.words.size() == 0) {
        out.set_zero();
        return;
    }

    if (a.words.size() == 1) {
        if (b.words.size() == 1) {
            uint128_t p = a.words[0];
            p *= b.words[0];

            out.words.reset_one_without_init();
            out.words[0] = p;

            uint64_t high = p >> 64;
            if (high)
                out.words.push_back(high);
            return;
        }
        __mul(b, a.words[0], /*carry*/0, out);
        return;
    }

    if (b.words.size() == 1) {
        __mul(a, b.words[0], /*carry*/0, out);
        return;
    }

    __mul(a, b, out);
}

constexpr void mul(natural& a, const natural& b) {
    if (a.words.size() == 0 || b.words.size() == 0) {
        a.words.set_zero();
        return;
    }

    if (b.words.size() == 1) {
        a *= b.words[0];
        return;
    }

    if (a.words.size() == 1) {
        __mul(b, a.words[0], /*carry*/0, a);
        return;
    }

    if (&a == &b)
        __square(a);
    else
        __mul(a, b, a);
}

constexpr natural operator*(const natural& a, const natural& b) { natural c; mul(a, b, /*out*/c); return c; }
constexpr natural operator*(natural a, std_unsigned_int auto b) { return a *= b; }
constexpr natural operator*(natural a, std_signed_int auto b) {
    if (b < 0)
        throw std::runtime_error("multiplication of natural with negative number");
    return std::move(a) * make_unsigned(b);
}
constexpr natural operator*(std_int auto a, natural b) { return std::move(b) * a; }

constexpr natural& operator*=(natural& a, const natural& b) { mul(a, b); return a; }
constexpr natural& operator*=(natural& a, std_unsigned_int auto b) {
    static_assert(sizeof(b) == 16 || sizeof(b) <= 8);

    if (b <= UINT64_MAX)
        return a *= static_cast<uint64_t>(b);
    __mul(a, static_cast<uint128_t>(b), a);
    return a;
}
constexpr natural& operator*=(natural& a, std_signed_int auto b) {
    if (b < 0)
        throw std::runtime_error("multiplication of natural with negative number");
    return a *= make_unsigned(b);
}

// acc += a * b
constexpr void add_product(natural& acc, const natural& a, const natural& b) {
    acc.words.resize(1 + std::max(acc.words.size(), a.words.size() * b.words.size()));
    // TODO in this algo it doesn't matter if a*b or b*a
    // TODO figure out how to choose A B order based on sizes of A and B
    for (size_t i = 0; i < a.words.size(); i++) {
        uint128_t carry = 0;
        for (natural::size_type j = 0; j < b.words.size(); j++) {
            carry += __mulq(a.words[i], b.words[j]);
            carry += acc.words[i + j];
            acc.words[i + j] = carry;
            carry >>= 64;
        }
        for (auto j = b.words.size(); carry; j++) {
            carry += acc.words[i + j];
            acc.words[i + j] = carry;
            carry >>= 64;
        }
    }
    acc.words.normalize();
}

// acc += a * b
constexpr void add_product(natural& acc, const natural& a, const uint64_t b) {
    // TODO optimize memory allocation
    acc += a * b;
}

// Assumes acc >= a * b
// acc -= a * b
constexpr void sub_product(natural& acc, const natural& a, const natural& b) {
    // TODO optimize memory allocation
    acc -= a * b;
#if 0
    // TODO in this algo it doesn't matter if a*b or b*a
    // TODO figure out how to choose A B order based on sizes of A and B
    for (decltype(a.words.size()) i = 0; i < a.words.size(); i++) {
        uint128_t carry = 0;
        for (decltype(b.words.size()) j = 0; j < b.words.size(); j++) {
            carry += __mulq(a.words[i], b.words[j]);
            carry += acc.words[i + j];
            acc.words[i + j] = carry;
            carry >>= 64;
        }
        for (auto j = b.words.size(); carry; j++) {
            carry += acc.words[i + j];
            acc.words[i + j] = carry;
            carry >>= 64;
        }
    }
    acc.words.normalize();
#endif
}

// Assumes acc >= a * b
// acc -= a * b
constexpr void sub_product(natural& acc, const natural& a, const uint64_t b) {
    // TODO optimize memory allocation
    natural temp;
    mul(a, b, temp);
    acc -= temp;
}

constexpr uint64_t div(const natural& dividend, uint64_t divisor, natural& quotient) {
    if (divisor == 0)
        throw std::runtime_error("division by zero");
    if (&dividend != &quotient)
        quotient.words.reset(dividend.words.size());
    uint128_t acc = 0;
    for (auto i = dividend.words.size(); i-- > 0;) {
        acc <<= 64;
        acc |= dividend.words[i];
        quotient.words[i] = acc / divisor;
        acc %= divisor;
    }
    quotient.words.normalize();
    return acc;
}

// returns static_cast<ucent>(a >> e) - without memory allocation
constexpr uint128_t extract_128bits(const natural& a, int e) {
    const auto word_shift = e / 64;
    const auto bit_shift = e % 64;

    if (word_shift >= a.words.size())
        return 0;

    uint128_t res = a.words[word_shift] >> bit_shift;
    if (word_shift + 1 >= a.words.size())
        return res;
    res |= static_cast<uint128_t>(a.words[word_shift + 1]) << (64 - bit_shift);
    if (word_shift + 2 >= a.words.size())
        return res;
    res |= static_cast<uint128_t>(a.words[word_shift + 2]) << (128 - bit_shift);
    return res;
}

// returns static_cast<ulong>(a >> e) - without memory allocation
constexpr uint64_t extract_64bits(const natural& a, int e) {
    const auto word_shift = e / 64;
    const auto bit_shift = e % 64;

    if (word_shift >= a.words.size())
        return 0;

    uint64_t res = a.words[word_shift] >> bit_shift;
    if (word_shift + 1 >= a.words.size())
        return res;
    res |= static_cast<uint64_t>(a.words[word_shift + 1]) << (64 - bit_shift);
    return res;
}

// returns largest q such that a * q <= b (assuming a != 0)
constexpr uint64_t __word_div(const natural& a, const natural& b) {
    if (a > b)
        return 0;
    if (a == b)
        return 1;
    if (b.words.size() == 1)
        return b.words[0] / a.words[0];
    if (b.words.size() == 2)
        return std::min<uint128_t>(static_cast<uint128_t>(b) / static_cast<uint128_t>(a), UINT64_MAX);

    const int e = b.num_bits() - 128;
    const uint128_t q = extract_128bits(b, e) / std::max<uint128_t>(1, extract_128bits(a, e));
    if (q > UINT64_MAX)
        return UINT64_MAX;

    uint64_t g = q;
    if (g * a <= b)
        return g;
    return ((g - 1) * a <= b) ? (g - 1) : (g - 2);
}

constexpr void div(const natural& dividend, const natural& divisor, natural& quotient, natural& remainder) {
    if (divisor.is_uint64()) {
        remainder = div(dividend, static_cast<uint64_t>(divisor), quotient);
        return;
    }
    if (divisor.words.size() == 0)
        throw std::runtime_error("division by zero");
    if (&dividend != &quotient)
        quotient.words.reset(dividend.words.size());
    remainder.words.set_zero();
    for (auto i = dividend.words.size(); i-- > 0;) {
        if (remainder.words.size() || dividend.words[i])
            remainder.words.insert_first_word(dividend.words[i]);
        const uint64_t q = __word_div(divisor, remainder);
        quotient.words[i] = q;
        sub_product(remainder, divisor, q); // remainder -= divisor * q
    }
    quotient.words.normalize();
}

constexpr void mod(const natural& dividend, const natural& divisor, natural& remainder) {
    if (divisor.is_uint64()) {
        remainder = dividend % static_cast<uint64_t>(divisor);
        return;
    }
    if (divisor.words.size() == 0)
        throw std::runtime_error("division by zero");
    remainder.words.set_zero();
    for (auto i = dividend.words.size(); i-- > 0;) {
        if (remainder.words.size() || dividend.words[i])
            remainder.words.insert_first_word(dividend.words[i]);
        const uint64_t q = __word_div(divisor, remainder);
        sub_product(remainder, divisor, q); // remainder -= divisor * q
    }
}

constexpr int natural::str(char* buffer, int buffer_size, unsigned base, const bool upper) const {
    char* p = buffer;
    const char* end = buffer + buffer_size;

    const char A = (upper ? 'A' : 'a') - 10;
    if (words.size() == 0) {
        if (p >= end)
            throw new std::runtime_error("buffer too small");
        *p++ = '0';
    } else if (words.size() == 1) {
        auto a = operator uint64_t();
        if (base == 10 && buffer_size >= 20 * words.size()) {
            while (a) {
                *p++ = '0' + int(a % 10);
                a /= 10;
            }
        } else if (base == 16 && buffer_size >= 16 * words.size()) {
            auto w = operator uint64_t();
            while (w) {
                int c = w % 16;
                *p++ = (c < 10) ? ('0' + c) : (A + c);
                w /= 16;
            }
        } else {
            while (a) {
                const int c = a % base;
                if (p >= end)
                    throw new std::runtime_error("buffer too small");
                *p++ = (c < 10) ? ('0' + c) : (A + c);
                a /= base;
            }
        }
    } else {
        if (base == 10 && buffer_size >= 20 * words.size()) {
            natural a = *this;
            while (a)
                *p++ = '0' + div(a, 10ul, /*out*/a);
        } else if (base == 16 && buffer_size >= 16 * words.size()) {
            for (int i = 0; i < words.size(); i++) {
                auto w = words[i];
                if (i == words.size() - 1) {
                    while (w) {
                        int c = w % 16;
                        *p++ = (c < 10) ? ('0' + c) : (A + c);
                        w /= 16;
                    }
                } else {
                    for (int j = 0; j < 16; j++) {
                        int c = w % 16;
                        *p++ = (c < 10) ? ('0' + c) : (A + c);
                        w /= 16;
                    }
                }
            }
        } else {
            natural n = *this;
            while (n) {
                const int c = div(n, static_cast<uint64_t>(base), /*out*/n);
                if (p >= end)
                    throw new std::runtime_error("buffer too small");
                *p++ = (c < 10) ? ('0' + c) : (A + c);
            }
        }
    }
    std::reverse(buffer, p);
    return p - buffer;
}

constexpr natural operator/(const natural& a, const natural& b) { natural quot, rem; div(a, b, /*out*/quot, /*out*/rem); return quot; }
constexpr natural operator/(const natural& a, std_int auto b) {
    if (b < 0)
        throw std::runtime_error("division of natural with negative number");
    natural q;
    div(a, static_cast<uint64_t>(b), q);
    return q;
}

constexpr natural& operator/=(natural& a, const natural &b) { natural rem; div(a, b, /*out*/a, /*out*/rem); return a; }
constexpr natural& operator/=(natural& a, std_int auto b) {
    if (b < 0)
        throw std::runtime_error("division of natural with negative number");
    div(a, static_cast<uint64_t>(b), a);
    return a;
}

constexpr natural operator%(const natural& a, const natural& b) { natural rem; mod(a, b, /*out*/rem); return rem; }
constexpr natural& operator%=(natural& a, const natural& b) { natural rem; mod(a, b, /*out*/rem); a = rem; return a; }

constexpr natural& operator<<=(natural& a, int64_t b) {
    if (b > 0) {
        if (a.words.size() == 0)
            return a;
        auto word_shift = b / 64;
        auto bit_shift = b % 64;

        if (bit_shift) {
            uint64_t carry = 0;
            for (natural::size_type i = 0; i < a.words.size(); ++i) {
                auto current = a.words[i];
                a.words[i] = (current << bit_shift) | carry;
                carry = current >> (64 - bit_shift);
            }
            if (carry)
                a.words.push_back(carry);
        }
        a.words.insert_first_n_words(word_shift);
        return a;
    }
    if (b < 0) {
        b = -b;
        auto word_shift = b / 64;
        auto bit_shift = b % 64;

        if (word_shift >= a.words.size()) {
            a.words.set_zero();
            return a;
        }

        a.words.erase_first_n_words(word_shift);
        if (bit_shift != 0) {
            uint64_t carry = 0;
            for (auto idx = a.words.size(); idx-- > 0;) {
                auto current = a.words[idx];
                a.words[idx] = (current >> bit_shift) | carry;
                carry = (current << (64 - bit_shift));
            }
        }
        a.words.normalize();
    }
    return a;
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

ALGEBRA_SHIFT_OP(natural)

constexpr natural operator~(natural a) {
    for (size_t i = 0; i < a.words.size(); i++)
        a.words[i] = ~a.words[i];
    return a;
}

constexpr natural operator|(const natural& a, const natural& b) {
    natural c;
    auto n = std::max(a.words.size(), b.words.size());
    c.words.reset(n);
    for (int i = 0; i < n; i++) {
        auto aw = (i < a.words.size()) ? a.words[i] : 0;
        auto bw = (i < b.words.size()) ? b.words[i] : 0;
        c.words[i] = aw | bw;
    }
    return c;
}

constexpr natural& operator|=(natural& a, const natural& b) {
    auto bs = b.words.size();
    if (bs > a.words.size())
        a.words.resize(bs);
    for (natural::size_type i = 0; i < bs; i++)
        a.words[i] |= b.words[i];
    a.words.normalize();
    return a;
}

constexpr natural& operator|=(natural& a, uint64_t b) {
    if (a.words.size() == 0) {
        a = b;
    } else {
        a.words[0] |= b;
    }
    return a;
}

constexpr natural operator&(const natural& a, const natural& b) {
    natural c;
    int n = std::max(a.words.size(), b.words.size());
    c.words.reset(n);
    for (int i = 0; i < n; i++) {
        auto aw = (i < a.words.size()) ? a.words[i] : 0;
        auto bw = (i < b.words.size()) ? b.words[i] : 0;
        c.words[i] = aw & bw;
    }
    c.words.normalize();
    return c;
}

constexpr natural& operator&=(natural& a, const natural& b) {
    size_t bs = b.words.size();
    if (bs > a.words.size())
        a.words.resize(bs);
    for (size_t i = 0; i < bs; i++)
        a.words[i] &= b.words[i];
    a.words.normalize();
    return a;
}

constexpr natural& operator&=(natural& a, uint64_t b) {
    a.words[0] &= b;
    a.words.normalize();
    return a;
}

constexpr natural operator^(const natural& a, const natural& b) {
    natural c;
    size_t n = std::max(a.words.size(), b.words.size());
    c.words.reset(n);
    for (size_t i = 0; i < n; i++) {
        auto aw = (i < a.words.size()) ? a.words[i] : 0;
        auto bw = (i < b.words.size()) ? b.words[i] : 0;
        c.words[i] = aw ^ bw;
    }
    c.words.normalize();
    return c;
}

constexpr natural& operator^=(natural& a, const natural& b) {
    size_t bs = b.words.size();
    if (bs > a.words.size())
        a.words.resize(bs);
    for (size_t i = 0; i < bs; i++)
        a.words[i] ^= b.words[i];
    a.words.normalize();
    return a;
}

constexpr natural& operator^=(natural& a, uint64_t b) {
    if (a.words.size() == 0) {
        a = b;
    } else {
        a.words[0] ^= b;
        a.words.normalize();
    }
    return a;
}

namespace literals {
constexpr auto operator""_n(const char* s) { return natural(s); }
}

}

template <>
struct std::formatter<algebra::natural, char> {
    int width = 0;
    char fill = ' ';
    char align = '>';
    unsigned base = 10;
    bool upper = false;

    constexpr auto parse(auto& ctx) {
        auto it = ctx.begin();
        if (it != ctx.end() && it + 1 != ctx.end() && (it[1] == '>' || it[1] == '<' || it[1] == '^')) {
            fill = *it++;
            align = *it++;
        }
        while (it != ctx.end() && '0' <= *it && *it <= '9')
            width = width * 10 + *it++ - '0';
        if (it != ctx.end()) {
            if (*it == 'b' || *it == 'B') {
                base = 2;
                it++;
            } else if (*it == 'o') {
                base = 8;
                it++;
            } else if (*it == 'd') {
                base = 10;
                it++;
            } else if (*it == 'x' || *it == 'X') {
                base = 16;
                upper = (*it++ == 'X');
            }
        }
        if (it == ctx.end() || *it != '}')
            throw std::format_error("Invalid format specifier for natural.");
        return it;
    }

    constexpr auto format(const algebra::natural& a, auto& ctx) const {
        auto it = ctx.out();
        int bound = a.str_size_upper_bound(base);

        char c_array[100];
        std::string str;
        char* buffer;
        if (bound <= 100) {
            buffer = c_array;
        } else {
            str.resize(bound);
            buffer = str.data();
        }

        int n = a.str(buffer, bound, base, upper);
        int pre = 0;
        int post = 0;
        if (width > n) {
            if (align == '>')
                pre = width - n;
            if (align == '<')
                post = width - n;
            if (align == '^') {
                pre = (width - n) / 2;
                post = width - pre - n;
            }
        }
        for (int i = 0; i < pre; i++)
            *it++ = fill;
        for (int i = 0; i < n; i++)
            *it++ = buffer[i];
        for (int i = 0; i < post; i++)
            *it++ = fill;
        return it;
    }
};

constexpr std::ostream& operator<<(std::ostream& os, const algebra::natural& a) { return os << a.str(); }

namespace algebra {

constexpr bool natural::is_uint8() const { return is_uint32() && words[0] <= 255; }
constexpr bool natural::is_uint16() const { return is_uint32() && words[0] <= 65535; }

constexpr int natural::str_size_upper_bound(unsigned base) const {
    if (words.size() == 0)
        return 1;
    int m;
    switch (base) {
    case 2: m = 64; break;
    case 3: m = 41; break;
    case 4: m = 32; break;
    case 5: m = 28; break;
    case 6: m = 25; break;
    case 7: m = 23; break;
    case 8: m = 22; break;
    case 9: m = 21; break;
    case 10: m = 20; break;
    case 11: m = 19; break;
    case 12: m = 18; break;
    case 13: m = 18; break;
    case 14: m = 17; break;
    case 15: m = 17; break;
    case 16: m = 16; break;
    default: m = 16; // avoid dependency on log_upper()
    }
    return words.size() * m;
}

template<std::floating_point T>
constexpr natural::operator T() const {
    const int exponent = static_cast<int>(num_bits()) - std::numeric_limits<T>::digits;
    if (exponent >= std::numeric_limits<T>::max_exponent)
        throw std::runtime_error("natural is too large for float");
    if (exponent <= 0)
        return words[0];
    const auto m = extract_64bits(*this, exponent);
    return std::ldexp(static_cast<T>(m), exponent);
}

static_assert(sizeof(natural) == 16);

constexpr natural::natural(std::string_view s, unsigned base) {
    const char* p = s.data();
    const char* end = s.data() + s.size();
    if (p >= end)
        throw std::runtime_error("expecting digit instead of end of string");

    uint64_t acc = 0;
    unsigned count = 0;
    if (base == 10) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            if ('0' > c || c > '9')
                throw std::runtime_error("expecting 0-9 for base 10");
            acc = acc * 10 + c - '0';
            count += 1;
            if (count == 19) {
                const uint64_t m = 10'000'000'000'000'000'000ull;
                mul_add(m, acc);
                acc = 0;
                count = 0;
            }
        }
        if (count)
            mul_add(pow(10ull, count), acc);
        words.normalize();
        return;
    }

    if (base == 2) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            if (c != '0' && c != '1')
                throw std::runtime_error("expecting 0-1 for base 2");
            acc = acc * 2 + c - '0';
            count += 1;
            if (count == 64) {
                words.insert_first_word(acc);
                acc = 0;
                count = 0;
            }
        }
    } else if (base == 8) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            if ('0' > c || c > '7')
                throw std::runtime_error("expecting 0-7 for base 8");
            acc = acc * 8 + c - '0';
            count += 3;
            if (count == 63) {
                words.insert_first_word(acc);
                acc = 0;
                count = 0;
            }
        }
    } else if (base == 16) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            int d;
            if ('0' <= c && c <= '9')
                d = c - '0';
            else if ('a' <= c && c <= 'f')
                d = c - 'a' + 10;
            else if ('A' <= c && c <= 'F')
                d = c - 'A' + 10;
            else
                throw std::runtime_error("expecting 0-9 or A-F for base 16");
            acc = acc * 16 + d;
            count += 4;
            if (count == 64) {
                words.insert_first_word(acc);
                acc = 0;
                count = 0;
            }
        }
    } else
        throw std::runtime_error("unsupported base");
    if (count) {
        *this <<= count;
        if (words.size() == 0)
            words.push_back(acc);
        else
            words[0] |= acc;
    }
    words.normalize();
}

}

template<>
struct std::hash<algebra::natural> {
    constexpr size_t operator()(const algebra::natural& a) const {
        return std::hash<algebra::integer_backend>()(a.words);
    }
};
