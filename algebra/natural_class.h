#pragma once
#include "algebra/integer_backend.h"
#include "algebra/util.h"
#include "algebra/kernels.h"
#include <string_view>
#include <algorithm>
#include <vector>

namespace algebra {

struct natural;
template<> struct IsNumberClass<natural> : std::true_type {};

struct natural {
    integer_backend words;

    constexpr natural(std::initializer_list<uint64_t> a) : words(a) { }
    constexpr natural() {}
    constexpr natural(std_int auto a) : words(a) {
        Check(a >= 0, "assigning negative number to natural");
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
        Check(is_uint8(), "integer is too large to fit in uint8");
        return words[0];
    }
    constexpr operator uint16_t() const {
        Check(is_uint16(), "integer is too large to fit in uint16");
        return words[0];
    }
    constexpr operator uint32_t() const {
        Check(is_uint32(), "integer is too large to fit in uint32");
        return words[0];
    }
    constexpr operator unsigned long() const {
        static_assert(sizeof(unsigned long) == 8);
        Check(is_uint64(), "integer is too large to fit in uint64");
        return words[0];
    }
    constexpr operator unsigned long long() const {
        static_assert(sizeof(unsigned long long) == 8);
        Check(is_uint64(), "integer is too large to fit in uint64");
        return words[0];
    }

    uint128_t unsafe_u128() const {
        uint128_t a = words[0];
        if (words.size() > 1)
            a |= static_cast<uint128_t>(words[1]) << 64;
        return a;
    }

    constexpr operator uint128_t() const {
        Check(is_uint128(), "integer is too large to fit in uint128");
        return unsafe_u128();
    }

    constexpr operator cnatural() const { return {words.data(), words.size()}; }
    constexpr operator vnatural() { return {{words.data(), words.size()}, words.capacity()}; }
    constexpr operator inatural() { return {words.data(), words.size()}; }

    constexpr size_t num_trailing_zeros() const { return algebra::num_trailing_zeros(*this); }
    constexpr natural& operator+=(const uint64_t b) {
        if (__add_and_return_carry(*this, b))
            words.push_back(1);
        return *this;
    }
    constexpr natural& operator+=(const uint128_t b) {
        uint128_t carry = __add_and_return_carry(*this, b);
        if (carry) {
            if (carry >> 64)
                words.push_back(carry, carry >> 64);
            else
                words.push_back(carry);
        }
        return *this;
    }

    constexpr natural& operator+=(const natural& _b) {
        const uint64_t* b = _b.words.data();
        const int B = _b.words.size();

        if (B > words.size()) {
            // Since we need to reallocate A anyway, add one more word if there is possibility of last carry
            words.reserve((b[B - 1] == UINT64_MAX) ? B + 1 : B);
            while (words.size() < B)
                words.push_back(0);
        }
        if (__add_and_return_carry(*this, static_cast<cnatural>(_b)))
            words.push_back(1);
        return *this;
    }

    constexpr natural& operator-=(uint64_t b) {
        inatural a = *this;
        __sub(a, b);
        words.downsize(a.size);
        return *this;
    }

    constexpr natural& operator-=(uint128_t b) {
        inatural a = *this;
        __sub(a, b);
        words.downsize(a.size);
        return *this;
    }

    constexpr natural& operator-=(const natural& b) {
        inatural a = *this;
        __sub(a, static_cast<cnatural>(b));
        words.downsize(a.size);
        return *this;
    }

    constexpr void mul_add(uint64_t a, uint64_t b) {
        if (a == 0) {
            words = b;
            return;
        }
        uint64_t carry = __mul_add_return_carry(*this, a, b);
        if (carry)
            words.push_back(carry);
    }

    constexpr uint64_t operator%(std_int auto b) const {
        static_assert(sizeof(b) <= 8);
        Check(b > 0, "division of natural by zero or negative number");
        return __mod(*this, static_cast<uint64_t>(b));
    }

    constexpr uint128_t operator%(const uint128_t b) const {
        Check(b > 0, "division of natural by zero or negative number");
        if (b <= UINT64_MAX)
            return __mod(*this, static_cast<uint64_t>(b));
        return __mod(*this, b);
    }

    constexpr int mod2() const { return words[0] % 2; }

    // TODO move to kernels
    constexpr int mod3() const {
        uint64_t acc = 0;
        for (int i = words.size(); i-- > 0;)
            acc += words[i] % 3;
        return acc % 3;
    }

    constexpr int mod4() const { return words[0] % 4; }

    // TODO move to kernels
    constexpr int mod5() const {
        uint64_t acc = 0;
        // acc overflow is not possible since 4 * UINT32_MAX < UINT64_MAX
        static_assert(sizeof(words.size()) == 4);
        // (2**64) mod 5 == 1
        for (int i = words.size(); i-- > 0;)
            acc += words[i] % 5;
        return acc % 5;
    }

    // TODO move to kernels
    constexpr int mod6() const {
        uint64_t m = 0;
        int i = words.size();
        while (i > 0) {
            m = m * 4 + words[--i] % 6;
            m %= 6;
        }
        return m;
    }

    // TODO move to kernels
    constexpr int mod7() const {
        uint64_t m = 0;
        int i = 0;
        while (i + 2 < words.size()) {
            m += words[i++] % 7;
            m += words[i++] % 7 * 2;
            m += words[i++] % 7 * 4;
        }
        if (i < words.size())
            m += words[i++] % 7;
        if (i < words.size())
            m += words[i] % 7 * 2;
        return m % 7;
    }

    constexpr int mod8() const { return words[0] % 8; }

    // TODO move to kernels
    constexpr int mod9() const {
        uint64_t m = 0;
        int i = 0;
        while (i + 2 < words.size()) {
            m += words[i++] % 9;
            m += (words[i++] % 9) * 7;
            m += (words[i++] % 9) * 4;
        }
        if (i < words.size())
            m += words[i++] % 9;
        if (i < words.size())
            m += words[i] % 9 * 7;
        return m % 9;
    }

    constexpr int mod10() const { return algebra::mod10(*this); }

    constexpr natural& operator%=(std_int auto b) { *this = operator%(b); return *this; }

    constexpr natural(std::string_view s, unsigned base = 10);
    constexpr natural(const char* s, uint32_t base = 10) : natural(std::string_view(s), base) {}

    constexpr void swap(natural& o) {
        words.swap(o.words);
    }

    constexpr int str_size_upper_bound(uint32_t base = 10) const;
    constexpr int str(char* buffer, int buffer_size, uint32_t base = 10, bool upper = true) const;
    constexpr std::string str(uint32_t base = 10, bool upper = true) const {
        std::string s;
        s.resize(str_size_upper_bound(base));
        s.resize(str(s.data(), s.size(), base, upper));
        return s;
    }
    constexpr std::string hex() const { return str(16); }

    constexpr int64_t num_bits() const { return algebra::num_bits(*this); }

    constexpr bool bit(int64_t i) const {
        size_t w = i / 64;
        size_t b = i % 64;
        return w < words.size() && (words[w] & (uint64_t(1) << b));
    }

    constexpr int64_t popcount() const {
        int64_t c = 0;
        for (int i = 0; i < words.size(); i++)
            c += std::popcount(words[i]);
        return c;
    }

    constexpr int64_t size_of() const { return words.size() * 8; }

    constexpr operator bool() const { return words.size(); }

    constexpr natural& operator++() {
        if (__increment_and_return_carry(*this))
            words.push_back(1);
        return *this;
    }
    constexpr natural& operator--() {
        Check(!words.empty(), "decrementing zero natural");
        inatural a = *this;
        __decrement(a);
        words.downsize(a.size);
        return *this;
    }

    constexpr natural operator++(int) { natural a = *this; operator++(); return a; }
    constexpr natural operator--(int) { natural a = *this; operator--(); return a; }

    template<std::floating_point T> constexpr operator T() const;
};

constexpr std::string str(const natural& a) { return str(static_cast<cnatural>(a)); }

// for debugging
constexpr std::string stre(const natural& a) {
    std::string s = "[";
    for (int i = 0; i < a.words.size(); i++)
        s += std::format(" {}", a.words[i]);
    s += " ]";
    return s;
}

constexpr natural operator-(const natural& a) {
    Check(a.words.size() == 0, "natural can't be negative");
    return 0;
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
    Check(a >= b, "natural can't be negative");
    static_assert(sizeof(b) <= 16);
    return a -= static_cast<uint128_t>(b);
}

constexpr natural operator-(const std_unsigned_int auto a, natural b) {
    if (a <= UINT64_MAX) {
        Check(b.words.size() <= 1 && uint64_t(a) >= b.words[0], "natural can't be negative");
        return uint64_t(a) - b.words[0];
    }
    Check(a <= b, "natural can't be negative");
    return a - static_cast<uint128_t>(b);
}

constexpr natural operator-(natural a, const std_signed_int auto b) {
    if (b < 0)
        return a + (~make_unsigned(b) + 1);
    return a - make_unsigned(b);
}

constexpr natural operator-(const std_signed_int auto a, natural b) {
    Check(a >= 0, "natural can't be negative");
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

constexpr bool operator<(const natural& a, const natural& b) { return __less(a, b); }

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

// TODO move to kernels
constexpr void __add(natural& a, const uint64_t* b, const int B, int shift = 0) {
    // TODO optimize this
    while (a.words.size() < B + shift)
        a.words.push_back(0);

    uint128_t acc = 0;
    for (int i = 0; i < B; ++i) {
        acc += a.words[shift + i];
        acc += b[i];
        a.words[shift + i] = acc;
        acc >>= 64;
    }
    for (int i = shift + B; i < a.words.size(); ++i) {
        acc += a.words[i];
        a.words[i] = acc;
        acc >>= 64;
    }
    if (acc)
        a.words.push_back(acc);
}

constexpr natural& operator<<=(natural& a, int64_t b);

constexpr bool is_power_of_two(const natural& a) { return is_power_of_two(static_cast<cnatural>(a)); }

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

    int Q = mul_max_size(a, b);
    q.words.reset(Q);
    vnatural vq = q;
    if (std::min(A, B) <= 2 || std::max(A, B) < KARATSUBA_LIMIT) {
        __add_product(vq, static_cast<cnatural>(a), static_cast<cnatural>(b));
    } else {
        const int W = 4 * std::max(A, B);
        if (W <= 1024) {
            uint64_t w[1024];
            __mul_karatsuba_rec(a, b, vq, w, w + W);
        } else {
            auto w = new uint64_t[W];
            __mul_karatsuba_rec(a, b, vq, w, w + W);
            delete[] w;
        }
    }
    q.words.resize(vq.size);
}

constexpr natural mul_karatsuba(const natural& a, const natural& b) {
    natural q;
    mul_karatsuba(a, b, q);
    return q;
}

// supports &a == &q
constexpr void __mul(const natural& a, const natural& b, natural& q) {
    if (&a != &q)
        q.set_zero();
    auto A = a.words.size();
    auto B = b.words.size();
    q.words.resize(A + B);
    vnatural vq = q;
    __mul({a.words.data(), A}, {b.words.data(), B}, vq, /*init*/false);
    q.words.downsize(vq.size);
}

// TODO move to kernels
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

// TODO move to kernels
// supports a == &out
constexpr void __mul(const uint64_t* a, const int A, uint64_t b, uint64_t carry, natural& out) {
    if (b == 0) {
        out.words.set_zero();
        if (carry)
            out.words.push_back(carry);
        return;
    }
    if (a != out.words.data())
        out.words.reset(A, /*initialize*/false);
    for (int i = 0; i < A; ++i) {
        auto acc = __mulq(a[i], b) + carry;
        out.words[i] = acc;
        carry = acc >> 64;
    }
    if (carry)
        out.words.push_back(carry);
}

constexpr void __mul(const natural& a, uint64_t b, uint64_t carry, natural& out) {
    __mul(a.words.data(), a.words.size(), b, carry, out);
}

// TODO move to kernels
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
        auto i_min = std::max<int>(0, k - n + 1);
        auto i_max = std::min<int>(n - 1, k);
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
            uint128_t p = __mulq(a.words[0], b.words[0]);
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
        a.mul_add(b.words[0], /*carry*/0);
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
    Check(b >= 0, "multiplication of natural with negative number");
    return std::move(a) * make_unsigned(b);
}
constexpr natural operator*(std_int auto a, natural b) { return std::move(b) * a; }

constexpr natural& operator*=(natural& a, const natural& b) { mul(a, b); return a; }
constexpr natural& operator*=(natural& a, std_unsigned_int auto b) {
    if (b <= UINT64_MAX) {
        a.mul_add(static_cast<uint64_t>(b), 0);
    } else {
        static_assert(sizeof(b) == 16 || sizeof(b) <= 8);
        __mul(a, static_cast<uint128_t>(b), a);
    }
    return a;
}
constexpr natural& operator*=(natural& a, std_signed_int auto b) {
    Check(b >= 0, "multiplication of natural with negative number");
    return a *= make_unsigned(b);
}

// A += B * C (without memory allocation)
constexpr void add_product(natural& a, const natural& b, const natural& c) {
    const int B = b.words.size();
    const int C = c.words.size();
    if (B == 0 || C == 0)
        return;
    // TODO if b == 1 just call __add(a, c)
    // TODO if c == 1 just call __add(a, b)

    int A = a.words.size();
    a.words.resize(std::max(A, B + C) + 1); // TODO compute thighter bound
    vnatural va {{a.words.data(), A}, a.words.capacity()};
    if (B < C)
        __add_product(va, b, static_cast<cnatural>(c));
    else
        __add_product(va, c, static_cast<cnatural>(b));
    a.words.downsize(va.size);
}

// A += B * c (without memory allocation)
constexpr void add_product(natural& a, const natural& b, const uint64_t c) {
    const int B = b.words.size();
    if (B == 0 || c == 0)
        return;
    // TODO if b == 1 just call __add(a, c)
    // TODO if c == 1 just call __add(a, b)

    int A = a.words.size();
    a.words.resize(std::max(A, B + 1) + 1); // TODO compute thighter bound
    vnatural va {{a.words.data(), A}, a.words.capacity()};
    __add_product(va, b, c);
    a.words.downsize(va.size);
}

// Assumes A >= B * C
// A -= B * C (without memory allocation)
constexpr void sub_product(natural& a, const natural& b, const natural& c) {
    const int B = b.words.size();
    const int C = c.words.size();
    if (B == 0 || C == 0)
        return;

    inatural ia = a;
    if (B < C)
        __sub_product(ia, b, static_cast<cnatural>(c));
    else
        __sub_product(ia, c, static_cast<cnatural>(b));
    a.words.downsize(ia.size);
}

// Assumes A >= B * c
// A -= B * c (without memory allocation)
constexpr void sub_product(natural& a, const natural& b, const uint64_t c) {
    const int B = b.words.size();
    if (B == 0 || c == 0)
        return;

    inatural ia = a;
    __sub_product(ia, b, c);
    a.words.downsize(ia.size);
}

constexpr uint64_t div(const natural& a, uint64_t b, natural& q) {
    if (&a != &q)
        q.words.reset(a.words.size());
    vnatural vq = q;
    uint64_t r = __div(a, b, vq);
    q.words.downsize(vq.size);
    return r;
}

constexpr void __div(cnatural a, cnatural b, natural& q, natural& r) {
    if (b.size <= 1) {
        if (a.words != q.words.data())
            q.words.reset(a.size);
        vnatural vq = q;
        r = __div(a, b[0], vq);
        q.words.resize(vq.size);
        return;
    }
    if (__less(a, b)) {
        r.words.resize(a.size);
        std::copy(a.words, a.words + a.size, r.words.data());
        q.set_zero(); // update Q after R in case &A == &Q
        return;
    }
    Check(b.size != 0, "division by zero");

    // NOTE max word size of R is b.word.size + 1
    const int Q = std::min(a.size, a.size - b.size + 1); //div_max_size(a, A, b, B); TODO
    if (a.words != q.words.data())
        q.words.reset(Q, /*initialize*/false);

#if 0
    if (Q == 1) {
        q.words.reset(1);
        const uint64_t w = __saturated_div(a, b);
        r.words.reset(a.size - 1, /*initialize*/false);
        std::copy(a.words + 1, a.words + a.size, r.words.data());
        vnatural vr = r;
        sub_product(vr, b, w);
        r.words.downsize(vr.size);
        q.words[0] = w;
        return;
    }
#endif

    r.set_zero();
    r.words.reset(a.size - Q, /*initialize*/false);
    std::copy(a.words + Q, a.words + a.size, r.words.data());

    for (int i = Q; i-- > 0;) {
        if (r.words.size() || a[i])
            r.words.insert_first_word(a[i]);

        const uint64_t w = __saturated_div(r, b);
        vnatural vr = r;
        __sub_product(vr, b, w);
        r.words.downsize(vr.size);
        q.words[i] = w;
    }
    q.words.downsize(Q);
    q.words.normalize();
}

#if 0
    if (A == 2) {
        if (B == 2 || a[1] < b[0]) { // equivalent to a[1] < b
            r = a[1];
            q.words[1] = 0;
            r.words.insert_first_word(a[0]); // r.words.size == 2
            const uint64_t w = __saturated_div(r, _b);
            sub_product(r, _b, w); // r -= b * w
            q.words[0] = w;
            q.words.normalize();
        } else {
            r = a[1];
            const uint64_t w = a[1] / b[0];
            r.words[0] -= b[0] * w;
            q.words[1] = w;

            uint128_t rr = (uint128_t(r) << 64) | a[0];
            uint128_t m = __mulq(b[0], UINT64_MAX);
            if (rr > m) {
                q.words[0] = UINT64_MAX;
                r.words[0] = static_cast<uint64_t>(rr - m);
            } else {
                uint64_t w, rem;
                __divq(rr, b[0], w, rem);
                q.words[0] = w;
                r.words[0] = rem;
            }
            r.words.normalize();
        }
        return;
    }
#endif

constexpr void div(const natural& a, const natural& b, natural& q, natural& r) {
    __div(a, b, q, r);
}

// TODO move this kernel to util.h
constexpr void mod(const natural& a, const natural& b, natural& r) {
    if (b.words.size() <= 1) {
        r = __mod(a, b.words[0]);
        return;
    }
    Check(!b.words.empty(), "division by zero");
    r.words.set_zero();
    for (auto i = a.words.size(); i-- > 0;) {
        if (r.words.size() || a.words[i])
            r.words.insert_first_word(a.words[i]);
        const uint64_t q = __saturated_div(r, b);
        sub_product(r, b, q); // r -= b * q
    }
}

// TODO move this kernel to util.h
constexpr void mod(natural& a, const natural& b) {
    const int A = a.words.size();
    const int B = b.words.size();

    if (B <= 1) {
        a = __mod(a, b.words[0]);
        return;
    }
    Check(B != 0, "division by zero");

    uint64_t* r = a.words.data() + A;
    int R = 0;
    for (auto i = A; i-- > 0;) {
        r -= 1;
        R += 1;
        inatural ir {r, R};
        const uint64_t w = __saturated_div(ir, b);
        __sub_product(ir, b, w); // r -= b * w
        R = ir.size;
    }
    a.words.resize(R);
}

constexpr natural& operator>>=(natural& a, int64_t b);

const int BZ_BASE_CASE_SIZE = 2;

// recursively divide A by D, and accumulate quotient into Q (with shift) and return remainder R
// p is a stack of temporaries, one for each recursion depth
void __divide_2n1n(cnatural a, cnatural d, int n, natural& q, int q_shift, natural& r, natural* p) {
    if (n <= BZ_BASE_CASE_SIZE) {
        __div(a, d, *p, r);
        __add(q, p->words.data(), p->words.size(), q_shift);
        return;
    }

    const int an = std::min(n, a.size);
    __divide_2n1n({a.words, an}, d, n/2, q, q_shift + n, *p, p + 1);
    p->words.insert_first_n_words(n);
    __add(*p, a.words + an, a.size - an);
    __divide_2n1n(*p, d, n/2, q, q_shift, r, p + 1);
}

constexpr void __divide_bz(natural a, natural d, natural& q, natural& r) {
    const auto D = d.words.size();
    const auto A = a.words.size();

    natural new_r;
    q.words.reserve_and_set_zero(div_max_size(a, d));
    r.words.reserve_and_set_zero(D);
    new_r.words.reserve_and_set_zero(D);

    const int shift = std::countl_zero(d.words.back());
    a <<= shift;
    d <<= shift;

    const int target_size = ((A + 2*D - 1) / (2*D)) * 2*D;

    static_assert(sizeof(decltype(A)) == 4);
    // TODO even better: precompute total sizes of all p needed, and allocate into one contiguous work buffer instead of p_stack
    natural p_stack[32]; // is enough for natural with at most UINT32_MAX words

    for (int i = target_size; i > 0; i -= 2*D) {
        r.words.insert_first_n_words(2*D); // TODO it might be possible to optimize this
        const int start = (i >= 2*D) ? i - 2*D : 0;
        const int end = std::min(i, A);
        std::copy(a.words.data() + start, a.words.data() + end, r.words.data());

        q.words.insert_first_n_words(2*D); // TODO it might be possible to optimize this with non-zero q_shift param below (q_shoft = i-2*D)
        __divide_2n1n(r, d, D, q, 0, new_r, p_stack);
        std::swap(r, new_r);
    }

    r >>= shift;
}

constexpr void divide_bz(const natural& a, const natural& d, natural& q, natural& r) {
    Check(!d.words.empty(), "division by zero");
    if (d > a) {
        q.set_zero();
        r = a;
        return;
    }
    if (d == a) {
        q = 1;
        r.set_zero();
        return;
    }
    __divide_bz(a, d, q, r);
}

constexpr int natural::str(char* buffer, int buffer_size, unsigned base, const bool upper) const {
    char* p = buffer;
    const char* end = buffer + buffer_size;

    const char A = (upper ? 'A' : 'a') - 10;
    if (words.size() == 0) {
        Check(p < end, "buffer too small");
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
                Check(p < end, "buffer too small");
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
                Check(p < end, "buffer too small");
                *p++ = (c < 10) ? ('0' + c) : (A + c);
            }
        }
    }
    std::reverse(buffer, p);
    return p - buffer;
}

constexpr natural operator/(const natural& a, const natural& b) { natural quot, rem; div(a, b, /*out*/quot, /*out*/rem); return quot; }
constexpr natural operator/(const natural& a, std_int auto b) {
    Check(b >= 0, "division of natural with negative number");
    natural q;
    div(a, static_cast<uint64_t>(b), q);
    return q;
}

constexpr natural& operator/=(natural& a, const natural &b) { natural rem; div(a, b, /*out*/a, /*out*/rem); return a; }
constexpr natural& operator/=(natural& a, std_int auto b) {
    Check(b >= 0, "division of natural with negative number");
    div(a, static_cast<uint64_t>(b), a);
    return a;
}

constexpr natural operator%(natural a, const natural& b) { mod(a, b); return a; }
constexpr natural& operator%=(natural& a, const natural& b) { mod(a, b); return a; }

constexpr natural& operator<<=(natural& a, int64_t b) {
    if (b > 0) {
        if (a.words.size() == 0)
            return a;
        auto word_shift = b / 64;
        auto bit_shift = b % 64;

        if (bit_shift) {
            uint64_t carry = 0;
            for (int i = 0; i < a.words.size(); ++i) {
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
        b = -b; // TODO undefined behavior for INT64_MIN, but that is huge amount of shift
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

constexpr natural& operator|=(natural& a, uint64_t b);
constexpr natural operator|(natural a, uint64_t b) { return a |= b; }

constexpr natural& operator|=(natural& a, const natural& b) {
    auto bs = b.words.size();
    if (bs > a.words.size())
        a.words.resize(bs);
    for (int i = 0; i < bs; i++)
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

constexpr natural& operator&=(natural& a, uint64_t b);
constexpr natural operator&(natural a, uint64_t b) { return a &= b; }

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

constexpr natural& operator^=(natural& a, uint64_t b);
constexpr natural operator^(natural a, uint64_t b) { return a ^= b; }

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
        algebra::Check(it != ctx.end() && *it == '}', "Invalid format specifier for natural.");
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
        return std::numeric_limits<T>::infinity();
    if (exponent <= 0)
        return words[0];
    const auto m = extract_u64(*this, exponent);
    return std::ldexp(static_cast<T>(m), exponent);
}

static_assert(sizeof(natural) == 16);

constexpr natural::natural(std::string_view s, unsigned base) {
    const char* p = s.data();
    const char* end = s.data() + s.size();
    Check(p < end, "expecting digit instead of end of string");

    uint64_t acc = 0;
    unsigned count = 0;
    if (base == 10) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            Check('0' <= c && c <= '9', "expecting 0-9 for base 10");
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
            Check(c == '0' || c == '1', "expecting 0-1 for base 2");
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
            Check('0' <= c && c <= '7', "expecting 0-7 for base 8");
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
                Fail("expecting 0-9 or A-F for base 16");
            acc = acc * 16 + d;
            count += 4;
            if (count == 64) {
                words.insert_first_word(acc);
                acc = 0;
                count = 0;
            }
        }
    } else
        Fail("unsupported base");
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
