#pragma once
#include "algebra/integer_backend.h"
#include <format>
#include <string_view>
#include <bit>
#include <stdexcept>
#include <algorithm>
#include <print>
#include <vector>

namespace algebra {

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

constexpr unsigned long abs_ulong(std::signed_integral auto a) { return (a >= 0) ? a : (~static_cast<unsigned long>(a) + 1); }

// TODO test cases for operator float()
// TODO add support for long double
struct natural {
    using size_type = integer_backend::size_type;
    using word = integer_backend::word;
    using dword = integer_backend::dword;

    integer_backend words;

    constexpr natural() {}
    constexpr natural(std::integral auto a) : words(a) {
        if (a < 0)
            throw std::runtime_error("assigning negative number to natural");
    }
    constexpr natural(natural&& o) : words(std::move(o.words)) { }
    constexpr natural(const natural& o) : words(o.words) { }

    constexpr void set_zero() { words.set_zero(); }
    constexpr void operator=(std::integral auto a) { words = a; }
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

    constexpr operator unsigned __int128() const {
        if (!is_uint128())
            throw std::runtime_error("cast overflow");
        dword a = words[0];
        if (words.size() >= 2)
            a |= dword(words[1]) << 64;
        return a;
    }

    constexpr size_t num_trailing_zeros() const {
        size_t a = 0;
        for (size_type i = 0; i < words.size(); i++) {
            if (words[i])
                return a + __builtin_ctzl(words[i]);
            a += 64;
        }
        return a;
    }

    constexpr natural& operator+=(word b) {
        for (size_type i = 0; b && i < words.size(); ++i) {
            dword acc = (dword)words[i] + b;
            words[i] = acc;
            b = acc >> 64;
        }
        if (b)
            words += b;
        return *this;
    }

    constexpr natural& operator+=(const natural& b) {
        // TODO predict max possible n taking into account carry!
        if (b.words.size() > words.size()) {
            words.reserve(b.words.size());
            while (words.size() < b.words.size())
                words += 0;
        }
        dword acc = 0;
        for (size_type i = 0; i < words.size(); ++i) {
            acc += words[i];
            if (i < b.words.size())
                acc += b.words[i];
            words[i] = acc;
            acc >>= 64;
        }
        if (acc)
            words += acc;
        return *this;
    }

    constexpr natural& operator-=(word b) {
        if (!words.allocated()) {
            words[0] -= b;
            words.normalize();
            return *this;
        }

        word* w = &words[0];
        size_type s = words.size();

        if (*w >= b) {
            *w -= b;
            return *this;
        }
        *w -= b;
        for (size_type i = 1; i < s; ++i)
            if (w[i]--) {
                if (w[s - 1] == 0)
                    words.pop_back();
                break;
            }
        words.normalize();
        return *this;
    }

    constexpr natural& operator-=(const natural& b) {
        word borrow = 0;
        for (size_type i = 0; i < b.words.size(); ++i) {
            __int128 diff = (__int128)words[i] - b.words[i] - borrow;
            if (diff < 0) {
                words[i] = diff + ((dword)1 << 64);
                borrow = 1;
            } else {
                words[i] = diff;
                borrow = 0;
            }
        }
        for (size_t i = b.words.size(); borrow && i < words.size(); ++i) {
            __int128 diff = (__int128)words[i] - borrow;
            if (diff < 0) {
                words[i] = diff + ((dword)1 << 64);
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
                words += carry;
            return;
        }
        for (size_type i = 0; i < words.size(); ++i) {
            dword acc = (dword)words[i] * a + carry;
            words[i] = acc;
            carry = acc >> 64;
        }
        if (carry)
            words += carry;
    }
    constexpr natural& operator*=(uint64_t b) { mul_add(b, 0); return *this; }

    constexpr uint64_t operator%(std::integral auto b) const {
        if (b <= 0)
            throw std::runtime_error((b == 0) ? "division by zero" : "division of natural by negative number");
        dword acc = 0;
        for (size_type i = words.size(); i-- > 0;) {
            acc <<= 64;
            acc |= words[i];
            acc %= static_cast<word>(b);
        }
        return acc;
    }

    constexpr uint64_t mod2() const { return words[0] % 2; }

    constexpr uint64_t mod3() const {
        word acc = 0;
        for (size_type i = words.size(); i-- > 0;)
            acc += words[i] % 3;
        return acc % 3;
    }

    constexpr uint64_t mod4() const { return words[0] % 4; }

    constexpr uint64_t mod5() const {
        word acc = 0;
        if (words.size() > std::numeric_limits<word>::max() / 4) {
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

    constexpr natural& operator%=(std::integral auto b) { *this = operator%(b); return *this; }

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

    constexpr size_t num_bits() const { return words.size() ? words.size() * 64 - __builtin_clzl(words.back()) : 0; }
    constexpr bool bit(size_t i) const {
        size_t bits_per_word = sizeof(word) * 8;
        size_t w = i / bits_per_word;
        size_t b = i % bits_per_word;
        return w < words.size() && (words[w] & (word(1) << b));
    }

    constexpr size_t popcount() const {
        size_t c = 0;
        for (size_type i = 0; i < words.size(); i++)
            c += std::popcount(words[i]);
        return c;
    }

    constexpr size_t size_of() const { return words.size() * 8; }

    constexpr operator bool() const { return words.size(); }

    constexpr natural& operator++() { operator+=(1); return *this; }
    constexpr natural& operator--() { operator-=(1); return *this; }

    constexpr natural operator++(int) { natural a = *this; operator++(); return a; }
    constexpr natural operator--(int) { natural a = *this; operator--(); return a; }

    constexpr operator float() const;
    constexpr operator double() const;
};

constexpr int num_bits(std::unsigned_integral auto a) { return sizeof(a) * 8 - __builtin_clzl(a); }

constexpr natural operator+(const natural& a, const natural& b) { natural c = a; c += b; return c; }
constexpr natural operator+(natural a, const std::integral auto b) {
    if (b < 0)
        throw std::runtime_error("addition of natural and negative number");
    a += static_cast<natural::word>(b);
    return a;
}
constexpr natural operator+(const std::integral auto a, natural b) { return std::move(b) + a; }

constexpr natural operator-(const natural& a, const natural& b) { natural c = a; c -= b; return c; }
constexpr natural operator-(natural a, const std::integral auto b) {
    if (b < 0)
        throw std::runtime_error("subtraction of natural and negative number");
    a -= static_cast<natural::word>(b);
    return a;
}
constexpr natural operator-(const std::integral auto a, natural b) { return natural(a) - b; }

constexpr bool operator<(const natural& a, const natural& b) {
    if (a.words.size() > b.words.size())
        return false;
    if (a.words.size() < b.words.size())
        return true;
    for (auto i = a.words.size(); i-- > 0;) {
        if (a.words[i] > b.words[i])
            return false;
        if (a.words[i] < b.words[i])
            return true;
    }
    return false;
}

constexpr bool operator<(const natural& a, const std::unsigned_integral auto b) { return a.words[0] < b && a.words.size() <= 1; }
constexpr bool operator<(const std::unsigned_integral auto a, const natural& b) { return a < b.words[0] || b.words.size() > 1; }

constexpr bool operator<(const natural& a, const std::signed_integral auto b) { return b >= 0 && a < static_cast<uint64_t>(b); }
constexpr bool operator<(const std::signed_integral auto a, const natural b) { return a < 0 || static_cast<uint64_t>(a) < b; }

constexpr bool operator<(const natural& a, const unsigned __int128 b) {
    if (b <= UINT64_MAX)
        return a < static_cast<uint64_t>(b);

    if (a.words.size() <= 1)
        return true;
    if (a.words.size() > 2)
        return false;
    natural::word bw0 = b;
    natural::word bw1 = b >> 64;
    return a.words[1] < bw1 || (a.words[1] == bw1 && a.words[0] < bw0);
}

constexpr bool operator>(const auto& a, const auto& b) { return b < a; }
constexpr bool operator>=(const auto& a, const auto& b) { return !(a < b); }
constexpr bool operator<=(const auto& a, const auto& b) { return !(a > b); }

constexpr bool operator==(const natural& a, const natural& b) {
    if (a.words.size() != b.words.size())
        return false;
    for (auto i = a.words.size(); i-- > 0;)
        if (a.words[i] != b.words[i])
            return false;
    return true;
}

constexpr bool operator==(const natural& a, const uint32_t b) { return a.words[0] == b && a.words.size() <= 1; }
constexpr bool operator==(const uint32_t a, const natural b) { return b == a; }

constexpr bool operator==(const natural& a, const int b) { return b >= 0 && a == static_cast<uint32_t>(b); }
constexpr bool operator==(const int a, const natural b) { return a >= 0 && static_cast<uint32_t>(a) == b; }

constexpr bool operator==(const natural& a, const uint64_t b) { return a.words[0] == b && a.words.size() <= 1; }
constexpr bool operator==(const uint64_t a, const natural& b) { return b == a; }

constexpr bool operator==(const natural& a, const unsigned __int128 b) {
    if (b <= UINT64_MAX)
        return a == static_cast<uint64_t>(b);

    if (a.words.size() != 2)
        return false;
    natural::word bw0 = b;
    natural::word bw1 = b >> 64;
    return a.words[1] == bw1 && a.words[0] == bw0;
}
constexpr bool operator==(const unsigned __int128 a, const natural& b) { return b == a; }

constexpr bool operator!=(const auto& a, const auto& b) { return !(a == b); }

// TODO can this be optimized for mul(a, a, out)?

// supports &a == &out
constexpr void __mul(const natural& a, const natural& b, natural& out) {
    if (&a != &out)
        out.set_zero();
    auto as = a.words.size();
    out.words.resize(a.words.size() + b.words.size());
    for (auto i = as; i-- > 0;) {
        natural::dword carry = 0;
        const auto w = a.words[i];
        out.words[i] = 0;
        for (natural::size_type j = 0; j < b.words.size(); ++j) {
            carry += (natural::dword)w * b.words[j];
            carry += out.words[i + j];
            out.words[i + j] = carry;
            carry >>= 64;
        }
        for (auto j = b.words.size(); carry; j++) {
            carry += out.words[i + j];
            out.words[i + j] = carry;
            carry >>= 64;
        }
    }
    out.words.normalize();
}

constexpr void square(natural& a) {
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
            natural::dword carry = static_cast<natural::dword>(ai) * aj;
            for (auto p = k; carry; p++) {
                carry += a.words[p];
                a.words[p] = carry;
                carry >>= 64;
            }
        }
    }
    a.words.normalize();
}

constexpr void mul(const natural& a, const natural& b, natural& out) {
    if (a.words.size() == 0 || b.words.size() == 0) {
        out.set_zero();
        return;
    }

    if (a.words.size() == 1) {
        out = b;
        out *= a.words[0];
        return;
    }

    if (b.words.size() == 1) {
        out = a;
        out *= b.words[0];
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
        auto w = a.words[0];
        a = b;
        a *= w;
        return;
    }

    if (&a == &b)
        square(a);
    else
        __mul(a, b, a);
}

constexpr natural operator*(const natural& a, const natural& b) { natural c; mul(a, b, /*out*/c); return c; }
constexpr natural operator*(natural a, std::integral auto b) {
    if (b < 0)
        throw std::runtime_error("multiplication of natural with negative number");
    a *= static_cast<natural::word>(b);
    return a;
}
constexpr natural operator*(std::integral auto a, natural b) { return std::move(b) * a; }

constexpr natural& operator*=(natural& a, const natural& b) { mul(a, b); return a; }

// acc += a * b
constexpr void add_product(natural& acc, const natural& a, const natural& b) {
    acc.words.resize(1 + std::max(acc.words.size(), a.words.size() * b.words.size()));
    // TODO in this algo it doesn't matter if a*b or b*a
    // TODO figure out how to choose A B order based on sizes of A and B
    for (size_t i = 0; i < a.words.size(); i++) {
        natural::dword carry = 0;
        for (natural::size_type j = 0; j < b.words.size(); j++) {
            carry += (natural::dword)a.words[i] * b.words[j];
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
    for (size_t i = 0; i < a.words.size(); i++) {
        natural::dword carry = 0;
        for (natural::size_type j = 0; j < b.words.size(); j++) {
            carry += (natural::dword)a.words[i] * b.words[j];
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

#if 0
ulong borrow = 0;
for (size_t i = 0; i < b.words.size(); ++i) {
    cent diff = (cent)words[i] - b.words[i] - borrow;
    if (diff < 0) {
        words[i] = static_cast<ulong>(diff + ((natural::dword)1 << 64));
        borrow = 1;
    } else {
        words[i] = static_cast<ulong>(diff);
        borrow = 0;
    }
}
for (size_t i = b.words.size(); borrow && i < words.size(); ++i) {
    cent diff = (cent)words[i] - borrow;
    if (diff < 0) {
        words[i] = static_cast<ulong>(diff + ((natural::dword)1 << 64));
        borrow = 1;
    } else {
        words[i] = static_cast<ulong>(diff);
        borrow = 0;
    }
}
words.normalize();
#endif

constexpr uint64_t div(const natural& dividend, uint64_t divisor, natural& quotient) {
    if (divisor == 0)
        throw std::runtime_error("division by zero");
    if (&dividend != &quotient)
        quotient.words.reset(dividend.words.size());
    natural::dword acc = 0;
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
constexpr unsigned __int128 extract_128bits(const natural& a, int e) {
    const auto bits_per_word = sizeof(natural::word) * 8;
    const auto word_shift = e / bits_per_word;
    const auto bit_shift = e % bits_per_word;

    if (word_shift >= a.words.size())
        return 0;

    unsigned __int128 res = a.words[word_shift] >> bit_shift;
    if (word_shift + 1 >= a.words.size())
        return res;
    res |= static_cast<unsigned __int128>(a.words[word_shift + 1]) << (64 - bit_shift);
    if (word_shift + 2 >= a.words.size())
        return res;
    res |= static_cast<unsigned __int128>(a.words[word_shift + 2]) << (128 - bit_shift);
    return res;
}

// returns static_cast<ulong>(a >> e) - without memory allocation
constexpr uint64_t extract_64bits(const natural& a, int e) {
    const auto bits_per_word = sizeof(natural::word) * 8;
    const auto word_shift = e / bits_per_word;
    const auto bit_shift = e % bits_per_word;

    if (word_shift >= a.words.size())
        return 0;

    uint64_t res = a.words[word_shift] >> bit_shift;
    if (word_shift + 1 >= a.words.size())
        return res;
    res |= static_cast<uint64_t>(a.words[word_shift + 1]) << (64 - bit_shift);
    return res;
}

// returns largest q such that a * q <= b (assuming a != 0)
constexpr natural::word __word_div(const natural& a, const natural& b) {
    if (a > b)
        return 0;
    if (a == b)
        return 1;
    if (b.words.size() == 1)
        return b.words[0] / a.words[0];
    if (b.words.size() == 2) {
        const natural::dword q = static_cast<natural::dword>(b) / static_cast<natural::dword>(a);
        return std::min<natural::dword>(q, std::numeric_limits<natural::word>::max());
    }

    const int e = b.num_bits() - 128;
    const natural::dword q = extract_128bits(b, e) / std::max<natural::dword>(1, extract_128bits(a, e));
    if (q > std::numeric_limits<natural::word>::max())
        return std::numeric_limits<natural::word>::max();

    natural::word g = q;
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
        const natural::word q = __word_div(divisor, remainder);
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
        const natural::word q = __word_div(divisor, remainder);
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
constexpr natural operator/(const natural& a, std::integral auto b) {
    if (b < 0)
        throw std::runtime_error("division of natural with negative number");
    natural q;
    div(a, static_cast<natural::word>(b), q);
    return q;
}

constexpr natural& operator/=(natural& a, const natural &b) { natural rem; div(a, b, /*out*/a, /*out*/rem); return a; }
constexpr natural& operator/=(natural& a, std::integral auto b) {
    if (b < 0)
        throw std::runtime_error("division of natural with negative number");
    div(a, static_cast<natural::word>(b), a);
    return a;
}

constexpr natural operator%(const natural& a, const natural& b) { natural rem; mod(a, b, /*out*/rem); return rem; }
constexpr natural& operator%=(natural& a, const natural& b) { natural rem; mod(a, b, /*out*/rem); a = rem; return a; }

constexpr natural& operator>>=(natural& a, size_t i) {
    size_t bits_per_word = sizeof(natural::word) * 8;
    size_t word_shift = i / bits_per_word;
    size_t bit_shift = i % bits_per_word;

    if (word_shift >= a.words.size()) {
        a.words.reset(0);
        return a;
    }

    a.words.erase_first_n_words(word_shift);
    if (bit_shift != 0) {
        natural::word carry = 0;
        for (auto idx = a.words.size(); idx-- > 0;) {
            auto current = a.words[idx];
            a.words[idx] = (current >> bit_shift) | carry;
            carry = (current << (bits_per_word - bit_shift));
        }
    }
    a.words.normalize();
    return a;
}

constexpr natural& operator<<=(natural& a, size_t i) {
    size_t bits_per_word = sizeof(natural::word) * 8;
    size_t word_shift = i / bits_per_word;
    size_t bit_shift = i % bits_per_word;

    if (bit_shift) {
        natural::word carry = 0;
        for (natural::size_type i = 0; i < a.words.size(); ++i) {
            auto current = a.words[i];
            a.words[i] = (current << bit_shift) | carry;
            carry = current >> (bits_per_word - bit_shift);
        }
        if (carry)
            a.words += carry;
    }
    a.words.insert_first_n_words(word_shift);
    return a;
}

constexpr natural operator>>(natural a, std::integral auto i) { a >>= i; return a; }
constexpr natural operator<<(natural a, std::integral auto i) { a <<= i; return a; }

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

constexpr natural& operator|=(natural& a, natural::word b) {
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

constexpr natural& operator&=(natural& a, natural::word b) {
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

constexpr natural& operator^=(natural& a, natural::word b) {
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

constexpr bool natural::is_uint8() const { return is_uint32() && *this <= 255u; }
constexpr bool natural::is_uint16() const { return is_uint32() && *this <= 65535u; }

constexpr natural pow(natural base, std::integral auto exp) {
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(natural, ...)");
    if (base == 2) {
        natural out;
        out.words.reset((exp + 64) / 64);
        out.words.back() = natural::word(1) << (exp % 64);
        return out;
    }

    if (exp == 0)
        return 1;

    natural result = 1;
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

constexpr natural pow(natural base, const natural& _exp) {
    if (_exp.is_uint64())
        return pow(base, static_cast<uint64_t>(_exp));

    natural result = 1;
    natural exp = _exp;
    while (exp) {
        if (exp.is_odd())
            result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

// uniformly sample from [0, (2**n)-1]
constexpr void uniform_sample_bits(const size_t n, auto& rng, natural& out) {
    static_assert(sizeof(rng()) == 8);
    auto w = (n + 63) / 64;
    out.words.reset(w, /*initialize*/false);
    while (w--)
        out.words[w] = rng();
    if ((n % 64) != 0)
        out.words.back() &= (natural::word(1) << (n % 64)) - 1;
    out.words.normalize();
}

constexpr natural uniform_sample_bits(const size_t n, auto& rng) {
    natural out;
    uniform_sample_bits(n, rng, /*out*/out);
    return out;
}

// uniformly sample from [0, count-1]
constexpr void uniform_sample(const natural& count, auto& rng, natural& out) {
    if (count.is_uint64()) {
        out = std::uniform_int_distribution<uint64_t>(0, static_cast<uint64_t>(count) - 1)(rng);
        return;
    }
    if (count.words.size() == 2 && count.words[0] == 0 && count.words[1] == 1) {
        static_assert(sizeof(rng()) == 8);
        out = rng();
        return;
    }

    const auto z = count.num_trailing_zeros();
    const auto b = count.num_bits();
    if (b == z + 1) { // if count is power of 2
        uniform_sample_bits(z, rng, /*out*/out);
        return;
    }

    const natural::size_type n = (b + 63) / 64;
    if (b == 64 * n) {
        // Note: power of two case is handled above!
        // mq = pow(2_n, 64 * n) / count
        // assert mq == 1
        while (true) {
            uniform_sample_bits(n * 64, rng, /*out*/out);
            if (out < count)
                return;
        }
    }
    // TODO avoid this division in mq == 2 case
    natural temp;
    natural mq = pow(natural(2), n * 64);
    mq /= count;
    if (mq == 2) {
        temp = count;
        temp <<= 1;
    }
    while (true) {
        if (mq == 2) {
            // optimization: avoid expensive div() for small mq
            uniform_sample_bits(n * 64, rng, /*out*/out);
            if (out < temp) {
                if (out >= count)
                    out -= count;
                return;
            }
        } else {
            uniform_sample_bits(n * 64, rng, /*out*/temp);
            div(temp, count, /*out*/temp, /*out*/out);
            if (temp < mq)
                return;
        }
    }
}

constexpr natural uniform_sample(const natural& count, auto& rng) {
    natural out;
    uniform_sample(count, rng, /*out*/out);
    return out;
}

constexpr natural uniform_sample(const natural& min, const natural& max, auto& rng) {
    natural count = max;
    count -= min;
    count += 1;
    natural out;
    uniform_sample(count, rng, /*out*/out);
    out += min;
    return out;
}

constexpr auto num_trailing_zeros(std::unsigned_integral auto a) { return a ? __builtin_ctzl(a) : 0; }

template<std::unsigned_integral T>
constexpr T __gcd_inner(T a, T b) {
    while (b) {
        b >>= __builtin_ctzl(b); // since b is non-zero
        if (a > b)
            std::swap(a, b);
        b -= a;
    }
    return a;
}

template<std::unsigned_integral T>
constexpr T gcd(T a, T b) {
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    auto az = num_trailing_zeros(a);
    if (az == 0)
        return __gcd_inner(a, b);
    auto common = std::min(az, num_trailing_zeros(b));
    return __gcd_inner(a >> az, b >> common) << common;
}

constexpr uint64_t gcd(uint64_t a, uint64_t b) {
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    auto common = __builtin_ctzl(a | b);
    a >>= __builtin_ctzl(a);
    do {
        b >>= __builtin_ctzl(b);
        if (a > b)
            std::swap(a, b);
        b -= a;
    } while (b);
    return a << common;
}

constexpr natural gcd(natural a, natural b) {
    if (a.words.size() == 1 && b.words.size() == 1)
        return gcd(a.words[0], b.words[0]);

    if (a.words.size() == 0)
        return b;
    if (b.words.size() == 0)
        return a;

    auto az = a.num_trailing_zeros();
    auto common = std::min(az, b.num_trailing_zeros());
    a >>= az;
    b >>= common;

    while (b) {
        if (a.words.size() == 1 && b.words.size() == 1) {
            a = __gcd_inner(a.words[0], b.words[0]);
            break;
        }
        b >>= b.num_trailing_zeros();
        if (a > b)
            a.swap(b);
        b -= a;
    }
    a <<= common;
    return a;
}

constexpr natural gcd(std::integral auto a, natural b) { return gcd(natural(abs_ulong(a)), std::move(b)); }
constexpr natural gcd(natural a, std::integral auto b) { return gcd(std::move(a), natural(abs_ulong(b))); }

// least common multiple
constexpr natural lcm(const natural& a, const natural& b) {
    natural m = a * b;
    m /= gcd(a, b);
    return m;
}

constexpr natural isqrt(const natural& x) {
    if (x == 0 || x == 1)
        return x;

    const auto b = x.num_bits();
    auto e = b + (b & 1);
    natural q, a;
    while (true) {
        q <<= 1;
        a = q;
        a |= 1;
        if (a * a <= (x >> e))
            q = a;
        if (e < 2)
            break;
        e -= 2;
    }
    return q;
}

constexpr bool is_prime(const uint64_t a) {
    if (a <= 3)
        return a >= 2;
    if (a % 2 == 0 || a % 3 == 0)
        return false;
    uint32_t i = 5;
    while (true) {
        // overflow here is not possible
        if (uint64_t(i) * i > a)
            return true;
        if (a % i == 0 || a % (i + 2) == 0)
            return false;
        i += 6;
        // uint32_t oveflow check
        if (i <= 5)
            return true;
    }
}

constexpr bool is_prime(const natural& a) {
    if (a.is_uint64())
        return is_prime(static_cast<uint64_t>(a));
    if (a <= 3)
        return a >= 2;
    if (a.is_even() == 0 || a.mod3() == 0)
        return false;
    uint64_t i = 5;
    while (true) {
        if (i * i > a)
            return true;
        if (a % i == 0 || a % (i + 2) == 0)
            return false;
        i += 6;
        // uint64_t oveflow check
        if (i <= 5)
            break;
    }
    natural ni = i - 6; // rollback to i before overflow
    ni += 6; // overflow safe
    natural q;
    while (true) {
        mul(ni, ni, /*out*/q); // q = ni * ni
        if (q > a)
            return true;

        mod(a, ni, /*out*/q); // q = a % ni
        if (!q)
            return false;

        ni += 2;
        mod(a, ni, /*out*/q); // q = a % ni
        if (!q)
            return false;

        ni += 4;
    }
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

// assumes are a and b are in [0, m-1] range
// a = (a * b) % m
constexpr void mul_mod(natural& a, const natural& b, const natural& m) {
    // TODO optimize
    a *= b;
    if (a >= m)
        a %= m;
}

// returns (a**n) mod p
constexpr void pow_mod(natural a, natural n, const natural& p, natural& out) {
    out = 1;
    if (a >= p)
        a %= p;
    for (size_t i = 0; i < n.num_bits(); i++) {
        if (n.bit(i))
            mul_mod(out, a, p);
        mul_mod(a, a, p);
    }
}

constexpr natural pow_mod(natural a, natural n, natural p) {
    natural out;
    pow_mod(a, n, p, /*out*/out);
    return out;
}

// Miller-Rabin algorithm
// It returns false if n is composite and returns true if n
// is probably prime.  k is an input parameter that determines
// accuracy level. Higher value of `rounds` indicates more accuracy.
constexpr bool is_likely_prime(const uint64_t n, int rounds, auto& rng) {
    if (n <= 4)
        return n == 2 || n == 3;
    if (n % 2 == 0)
        return false;

    const auto s = num_trailing_zeros(n - 1);
    const uint64_t d = (n - 1) >> s;

    for (int i = 0; i < rounds; i++) {
        const uint64_t a = std::uniform_int_distribution<uint64_t>(2, n-2)(rng);

        uint64_t x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1)
            continue;

        int r = 1;
        for (; r < s; r++) {
            x = (static_cast<__uint128_t>(x) * x) % n;
            if (x == n - 1)
                break;
        }
        if (r == s)
            return false;
    }
    return true;
}

constexpr bool is_likely_prime(const natural& n, int rounds, auto& rng) {
    if (n.is_uint64())
        return is_likely_prime(static_cast<uint64_t>(n), rounds, rng);

    if (n <= 4)
        return n == 2 || n == 3;
    if (n.is_even())
        return false;

    const natural n_minus_1 = n - 1;
    const natural n_minus_3 = n - 3;
    auto s = n_minus_1.num_trailing_zeros();
    natural d = n_minus_1;
    d >>= s;
    natural x, a;

    for (int i = 0; i < rounds; i++) {
        // sample uniformly from [2, n-2] and store in a
        uniform_sample(n_minus_3, rng, /*out*/a);
        a += 2;

        pow_mod(a, d, n, /*out*/x);
        if (x == 1 || x == n_minus_1)
            continue;

        int r = 1;
        for (; r < s; r++) {
            mul_mod(x, x, n); // x = (x * x) % n
            if (x == n_minus_1)
                break;
        }
        if (r == s)
            return false;
    }
    return true;
}

constexpr std::vector<std::pair<uint64_t, int>> factorize(uint64_t a) {
    if (a <= 1)
        return {};

    std::vector<std::pair<uint64_t, int>> out;
    auto z = num_trailing_zeros(a);
    if (z) {
        out.emplace_back(2, z);
        a >>= z;
    }

    uint64_t p = 3;
    while (a > 1) {
        int count = 0;
        while (a % p == 0) {
            a /= p;
            count += 1;
        }
        if (count)
            out.emplace_back(p, count);
        p += 2;
    }
    return out;
}

constexpr bool is_power_of_two(const natural& a) {
    return a.num_bits() == 1 + a.num_trailing_zeros();
}

constexpr auto log_lower(natural a, uint64_t base) {
    uint64_t count = 0;
    if (!a)
        return count;
    while (true) {
        a /= base;
        if (!a)
            break;
        count += 1;
    }
    return count;
}
constexpr auto log_upper(natural a, uint64_t base) {
    uint64_t count = 0;
    while (a) {
        a /= base;
        count += 1;
    }
    return count;
}

constexpr int natural::str_size_upper_bound(unsigned base) const {
    if (words.size() == 0)
        return 1;
    int m;
    switch (base) {
    case 2: m = 64; break;
    case 4: m = 32; break;
    case 8: m = 22; break;
    case 10: m = 20; break;
    case 16: m = 16; break;
    default:
        m = log_upper(std::numeric_limits<word>::max(), base);
    }
    return words.size() * m;
}

constexpr natural::operator float() const {
    if (words.size() == 0)
        return 0.0f;

    const int exponent = static_cast<int>(num_bits()) - 24;
    if (exponent > 127) return std::numeric_limits<float>::infinity();
    if (exponent <= 0) return words[0];

    const auto m = extract_64bits(*this, exponent);
    return std::ldexp(static_cast<float>(m), exponent);
}

// TODO operator long double
constexpr natural::operator double() const {
    if (words.size() == 0)
        return 0.0;

    const int exponent = static_cast<int>(num_bits()) - 53;
    if (exponent > 1023) return std::numeric_limits<double>::infinity();
    if (exponent <= 0) return words[0];

    const auto m = extract_64bits(*this, exponent);
    return std::ldexp(static_cast<double>(m), exponent);
}

static_assert(sizeof(natural) == 16);

constexpr natural::natural(std::string_view s, unsigned base) {
    const char* p = s.data();
    const char* end = s.data() + s.size();
    if (p >= end)
        throw std::runtime_error("expecting digit instead of end of string");

    natural::word acc = 0;
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
                const natural::word m = 10'000'000'000'000'000'000ull;
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
            words += acc;
        else
            words[0] |= acc;
    }
    words.normalize();
}

}
