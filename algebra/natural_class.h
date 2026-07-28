#pragma once
#include "algebra/integer_backend.h"
#include "algebra/util.h"
#include "algebra/kernels.h"
#include <string_view>
#include <algorithm>
#include <vector>

namespace algebra {

struct integer;
template<> struct IsNumberClass<integer> : std::true_type {};

struct integer {
    integer_backend words;

    constexpr integer(std::initializer_list<uint64_t> a) : words(a) { }
    constexpr integer() {}
    constexpr integer(std_int auto a) : words(a) { }
    constexpr integer(integer&& o) : words(std::move(o.words)) { }
    constexpr integer(const integer& o) : words(o.words) { }

    constexpr void set_zero() { words.set_zero(); }
    constexpr void operator=(std_int auto a) { words = a; }
    constexpr void operator=(integer&& o) { words = std::move(o.words); }
    constexpr void operator=(const integer& o) { words = o.words; }

    constexpr bool is_one() const { return words.size() == 1 && words[0] == 1 && sign() >= 0; }
    constexpr bool is_even() const { return (low_word() % 2) == 0; }
    constexpr bool is_odd() const { return low_word() % 2; }

    constexpr bool is_uint8() const;
    constexpr bool is_uint16() const;
    constexpr bool is_uint32() const { return sign() >= 0 && words.size() <= 1 && low_word() <= UINT32_MAX; }
    constexpr bool is_uint64() const { return sign() >= 0 && words.size() <= 1; }
    constexpr bool is_uint128() const { return sign() >= 0 && words.size() <= 2; }

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
    constexpr integer& __mag_add_word(const uint64_t b) {
        // note: for an empty (zero) value the carry is b itself, not 1
        const uint64_t carry = __add_and_return_carry(*this, b);
        if (carry)
            words.push_back(carry);
        return *this;
    }
    constexpr integer& __mag_add_word(const uint128_t b) {
        uint128_t carry = __add_and_return_carry(*this, b);
        if (carry) {
            if (carry >> 64)
                words.push_back(carry, carry >> 64);
            else
                words.push_back(carry);
        }
        return *this;
    }


    constexpr integer& __mag_sub_word(uint64_t b) {
        // note: spelled out instead of *this >= b, because operator<(integer, unsigned) is
        // declared further down this header and is not visible from inside the class
        Check(words.size() > 1 || words[0] >= b, "integer can't be negative");
        inatural a = *this;
        __sub(a, b);
        words.downsize(a.size);
        return *this;
    }

    constexpr integer& __mag_sub_word(uint128_t b) {
        if (b <= UINT64_MAX)
            return __mag_sub_word(static_cast<uint64_t>(b));
        Check(words.size() > 2 || (words.size() == 2 && unsafe_u128() >= b), "integer can't be negative");
        inatural a = *this;
        __sub(a, b);
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
        Check(b > 0, "division of integer by zero or negative number");
        return __mod(*this, static_cast<uint64_t>(b));
    }

    constexpr uint128_t operator%(const uint128_t b) const {
        Check(b > 0, "division of integer by zero or negative number");
        if (b <= UINT64_MAX)
            return __mod(*this, static_cast<uint64_t>(b));
        return __mod(*this, b);
    }

    constexpr int mod2() const { return words[0] % 2; }
    constexpr int mod3() const { return algebra::mod3(*this); }
    constexpr int mod4() const { return words[0] % 4; }
    constexpr int mod5() const { return algebra::mod5(*this); }
    constexpr int mod6() const { return algebra::mod6(*this); }
    constexpr int mod7() const { return algebra::mod7(*this); }
    constexpr int mod8() const { return words[0] % 8; }
    constexpr int mod9() const { return algebra::mod9(*this); }
    constexpr int mod10() const { return algebra::mod10(*this); }

    constexpr integer(std::string_view s, unsigned base = 10);
    constexpr integer(const char* s, uint32_t base = 10) : integer(std::string_view(s), base) {}

    constexpr void swap(integer& o) {
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

    // these two work on the magnitude and leave the sign alone
    constexpr void __inc_magnitude() {
        if (__increment_and_return_carry(*this))
            words.push_back(1);
    }
    constexpr void __dec_magnitude() {
        Check(!words.empty(), "decrementing zero magnitude");
        inatural a = *this;
        __decrement(a);
        words.downsize(a.size);
    }

    constexpr integer& operator++() {
        if (is_negative()) {
            const bool was_one = words.size() == 1 && words[0] == 1;
            __dec_magnitude();
            if (was_one)
                words.set_negative(false); // -1 + 1 is zero, which has no sign
        } else {
            __inc_magnitude();
        }
        return *this;
    }
    constexpr integer& operator--() {
        if (is_negative()) {
            __inc_magnitude();
        } else if (words.size() > 0) {
            __dec_magnitude();
        } else {
            *this = -1;
        }
        return *this;
    }

    constexpr integer operator++(int) { integer a = *this; operator++(); return a; }
    constexpr integer operator--(int) { integer a = *this; operator--(); return a; }

    template<std::floating_point T> constexpr operator T() const;

    // ---- the sign aware half, from what used to be a separate integer class ----
    constexpr operator char() const {
        Check(is_int8(), "integer -> int8 overflow");
        return is_negative() ? -static_cast<int>(words[0]) : words[0];
    }
    constexpr operator short() const {
        Check(is_int16(), "integer -> int16 overflow");
        return is_negative() ? -static_cast<int>(words[0]) : words[0];
    }
    constexpr operator int() const {
        Check(is_int32(), "integer -> int32 overflow");
        // negate in unsigned: for INT_MIN the magnitude is 2**31, which does not fit in an int,
        // so negating after the cast is signed overflow. Converting an out of range unsigned
        // value to a signed type is well defined (modular) since C++20, same as operator int128_t.
        if (is_negative())
            return static_cast<int>(-static_cast<uint32_t>(words[0]));
        return static_cast<int>(words[0]);
    }
    constexpr operator long() const {
        Check(is_int64(), "integer -> int64 overflow");
        // see operator int() for why the negation is done in unsigned
        if (is_negative())
            return static_cast<long>(-words[0]);
        return static_cast<long>(words[0]);
    }
    constexpr operator int128_t() const {
        Check(is_int128(), "integer -> int128 overflow");
        if (sign() == 2) return (uint128_t(words[1]) << 64) | words[0];
        if (sign() == 1) return uint128_t(words[0]);
        if (sign() == -1) return -uint128_t(words[0]);
        if (sign() == -2) return -((uint128_t(words[1]) << 64) | words[0]);
        return 0;
    }

    // words[0] is only meaningful when the value is not zero
    constexpr uint64_t low_word() const { return words.size() ? words[0] : 0; }
    // the sign lives in the sign of the size, which is how integer_backend stores it
    constexpr int sign() const { return words.sign(); }
    constexpr bool is_negative() const { return sign() < 0; }
    constexpr bool is_zero() const { return words.size() == 0; }
    constexpr void negate() { words.negate(); }

    constexpr bool is_int8() const {
        if (words.size() > 1) return false;
        if (low_word() <= 127) return true;
        return sign() < 0 && low_word() == 128;
    }
    constexpr bool is_int16() const {
        if (words.size() > 1) return false;
        if (low_word() <= INT16_MAX) return true;
        return sign() < 0 && low_word() == static_cast<uint64_t>(INT16_MAX) + 1;
    }
    constexpr bool is_int32() const {
        if (words.size() > 1) return false;
        if (low_word() <= INT32_MAX) return true;
        return sign() < 0 && low_word() == static_cast<uint64_t>(INT32_MAX) + 1;
    }
    constexpr bool is_int64() const {
        if (words.size() > 1) return false;
        if (low_word() <= INT64_MAX) return true;
        return sign() < 0 && low_word() == static_cast<uint64_t>(INT64_MAX) + 1;
    }
    constexpr bool is_int128() const {
        if (words.size() > 2) return false;
        if (words.size() < 2) return true;
        const uint64_t w = words[1];
        if ((w & (uint64_t(1) << 63)) == 0) return true;
        // only INT128_MIN reaches past the positive range, and its magnitude is exactly 2**127
        return sign() < 0 && w == uint64_t(1) << 63 && words[0] == 0;
    }
};

// Transitional: natural is another name for integer now. What the type used to guarantee --
// a non-negative value -- is checked at runtime by the functions that need it.
using natural = integer;

// The sign aware operators are defined below, after the magnitude machinery they build on, so
// they are declared here for that machinery to use.
constexpr bool operator==(const integer& a, const integer& b);
constexpr bool operator==(const integer& a, std_int auto b);
constexpr integer operator-(integer a);
constexpr integer& operator+=(integer& a, const integer& b);
constexpr integer& operator-=(integer& a, const integer& b);
constexpr integer& operator+=(integer& a, std_int auto b);
constexpr integer& operator-=(integer& a, std_int auto b);
constexpr integer operator+(integer a, std_int auto b);
constexpr integer operator+(std_int auto a, integer b);
constexpr integer operator+(integer a, const integer& b);
constexpr integer operator-(integer a, const integer& b);
constexpr integer operator-(integer a, std_int auto b);
constexpr integer operator-(std_int auto a, integer b);
constexpr integer operator*(const integer& a, const integer& b);
constexpr integer operator*(const integer& a, std_int auto b);
constexpr integer& operator*=(integer& a, const integer& b);
constexpr integer& operator*=(integer& a, std_int auto b);
constexpr integer operator/(const integer& a, const integer& b);
constexpr integer operator/(const integer& a, const std_int auto b);
constexpr integer& operator/=(integer& a, const integer& b);
constexpr integer& operator/=(integer& a, const std_int auto b);
constexpr integer operator%(const integer& a, const integer& divisor);
constexpr int64_t operator%(const integer& a, int64_t b);
constexpr int operator%(const integer& a, int b);
constexpr int64_t operator%(const integer& a, unsigned b);
constexpr integer operator%(const integer& a, uint64_t b);
constexpr integer& operator%=(integer& a, const integer& b);
constexpr integer& operator%=(integer& a, std_int auto b);
constexpr bool operator<(const integer& a, const integer& b);
constexpr bool operator<(const integer& a, std_int auto b);
constexpr bool operator<(std_int auto a, const integer& b);
constexpr integer operator~(integer a);
constexpr void operator<<=(integer& a, int64_t i);
constexpr void mul(const integer& a, const integer& b, integer& c);
constexpr void mul(integer& a, const integer& b);
constexpr std::string str(const integer& a);
constexpr std::string stre(const integer& a);
constexpr void add_product(integer& a, const integer& b, const integer& c);
constexpr void sub_product(integer& a, const integer& b, const integer& c);
constexpr void add_product(integer& a, const integer& b, const std_int auto c);
constexpr void sub_product(integer& a, const integer& b, const std_int auto c);
constexpr void div(const integer& a, const integer& b, integer& quot, integer& rem);
template<std_signed_int T> constexpr T div(const integer& a, T b, integer& quot);
template<std_unsigned_int T> constexpr T div(const integer& a, T b, integer& quot);
constexpr integer mod(const integer& a, const integer& b);
constexpr uint64_t mod(const integer& a, uint64_t b);
constexpr unsigned mod(const integer& a, uint32_t b);
//@@FORWARD_DECLS_END@@

ALGEBRA_SHIFT_OP(integer)

























// q must not alias a or b
constexpr void mul_karatsuba(const integer& a, const integer& b, integer& q) {
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

    if (is_power_of_two(static_cast<cnatural>(a))) {
        const size_t z = (A - 1) * 64 + std::countr_zero(a.words[A - 1]); // = a.num_trailing_zeros() but O(1)
        const size_t bits = b.num_bits() + z;
        const size_t words = (bits + 63) / 64;
        q.words.reset(words, /*init*/false); // preallocates memory!
        // TODO this can be done directly without moving
        q = b;
        q <<= z;
        return;
    }
    if (is_power_of_two(static_cast<cnatural>(b))) {
        const size_t z = (B - 1) * 64 + std::countr_zero(b.words[B - 1]); // = b.num_trailing_zeros() but O(1)
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

constexpr integer mul_karatsuba(const integer& a, const integer& b) {
    integer q;
    mul_karatsuba(a, b, q);
    return q;
}

// supports &a == &q
constexpr void __mul(const integer& a, const integer& b, integer& q) {
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
constexpr void __mul(const integer& a, uint128_t b, integer& out) {
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
constexpr void __mul(const uint64_t* a, const int A, uint64_t b, uint64_t carry, integer& out) {
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

constexpr void __mul(const integer& a, uint64_t b, uint64_t carry, integer& out) {
    __mul(a.words.data(), a.words.size(), b, carry, out);
}

// TODO move to kernels
// assumes a.words.size() >= 2
constexpr void __square(integer& a) {
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

    // a*a == sum(a[i]*a[i] << 2i) + 2 * sum(a[i]*a[j] << (i+j)) for i<j,
    // which needs half the multiplications of a full a*a
    const int n = a.words.size();
    integer r;
    r.words.reset(n * 2); // zero initialized

    // cross products, each counted once
    for (int i = 0; i < n; i++) {
        const uint64_t ai = a.words[i];
        if (ai == 0)
            continue;
        uint128_t carry = 0;
        for (int j = i + 1; j < n; j++) {
            carry += r.words[i + j];
            carry += __mulq(ai, a.words[j]);
            r.words[i + j] = carry;
            carry >>= 64;
        }
        for (int k = i + n; carry; k++) {
            carry += r.words[k];
            r.words[k] = carry;
            carry >>= 64;
        }
    }

    // double them (the sum of cross products is less than a*a/2, so this fits)
    uint64_t shifted = 0;
    for (int k = 0; k < n * 2; k++) {
        const uint64_t w = r.words[k];
        r.words[k] = (w << 1) | shifted;
        shifted = w >> 63;
    }

    // add the squares of the individual words
    uint128_t carry = 0;
    for (int i = 0; i < n; i++) {
        carry += r.words[i * 2];
        carry += __mulq(a.words[i], a.words[i]);
        r.words[i * 2] = carry;
        carry >>= 64;
        carry += r.words[i * 2 + 1];
        r.words[i * 2 + 1] = carry;
        carry >>= 64;
    }

    r.words.normalize();
    a = std::move(r);
}

constexpr void square(integer& a) {
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

// Inputs are views so that a caller holding the words elsewhere -- integer, for one -- does not
// have to materialise a integer for them.
constexpr void mul(cnatural a, cnatural b, integer& out) {
    if (a.size == 0 || b.size == 0) {
        out.set_zero();
        return;
    }

    if (a.size == 1) {
        if (b.size == 1) {
            uint128_t p = __mulq(a[0], b[0]);
            out.words.reset_one_without_init();
            out.words[0] = p;
            uint64_t high = p >> 64;
            if (high)
                out.words.push_back(high);
            return;
        }
        __mul(b.words, b.size, a[0], /*carry*/0, out);
        return;
    }

    if (b.size == 1) {
        __mul(a.words, a.size, b[0], /*carry*/0, out);
        return;
    }

    out.set_zero();
    out.words.resize(a.size + b.size);
    vnatural vq = out;
    __mul(a, b, vq, /*init*/false);
    out.words.downsize(vq.size);
}


constexpr integer operator*(std_int auto a, integer b) { return std::move(b) * a; }



// A += B * c (without memory allocation)


// Assumes A >= B * c
// A -= B * c (without memory allocation)


constexpr void __div(cnatural a, cnatural b, integer& q, integer& r) {
    Check(b.size != 0, "division by zero");
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

    // NOTE max word size of R is b.word.size + 1
    const int Q = a.size - b.size + 1; // b.size >= 2 here, so this is less than a.size
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
        Check(__sub_product(vr, b, w), "__saturated_div() overestimated the quotient");
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


// TODO move this kernel to util.h
// The remainder, taking views so that a caller holding the words elsewhere -- integer -- needs no
// integer for them. Named apart from mod() so that adding it does not change any overload set.
constexpr void __mod_into(cnatural a, cnatural b, integer& r) {
    Check(b.size != 0, "division by zero");
    if (b.size <= 1) {
        r = __mod(a, b[0]);
        return;
    }
    r.words.set_zero();
    r.words.reserve(b.size + 1); // the remainder never exceeds this, so the loop below never grows it
    for (auto i = a.size; i-- > 0;) {
        if (r.words.size() || a[i])
            r.words.insert_first_word(a[i]);
        const uint64_t q = __saturated_div(r, b);
        inatural ir = r;
        Check(__sub_product(ir, b, q), "mod() remainder went negative");
        r.words.downsize(ir.size);
    }
}

constexpr void mod(const integer& a, const integer& b, integer& r) { __mod_into(a, b, r); }

// TODO move this kernel to util.h
constexpr void __mod_magnitude(integer& a, const integer& b) {
    const int A = a.words.size();
    const int B = b.words.size();

    Check(B != 0, "division by zero");
    if (B <= 1) {
        a = __mod(a, b.words[0]);
        return;
    }

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

constexpr integer& operator>>=(integer& a, int64_t b);

// ---- magnitude level operations ----
// The sign aware layer below borrows a magnitude and then operates on it. Now that natural and
// integer are one type, calling the operator by its usual name would resolve back to the sign
// aware overload and recurse, so the magnitude versions live under their own names.

constexpr void __mag_mul(integer& a, const integer& b) {
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

constexpr uint64_t __mag_div(const integer& a, uint64_t b, integer& q) {
    Check(b != 0, "division by zero");
    if (&a != &q)
        q.words.reset(a.words.size());
    vnatural vq = q;
    uint64_t r = __div(a, b, vq);
    q.words.downsize(vq.size);
    return r;
}

constexpr void __mag_div(const integer& a, const integer& b, integer& q, integer& r) { __div(a, b, q, r); }

constexpr integer& __mag_shl(integer& a, int64_t b) {
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

constexpr std::string __mag_stre(const integer& a) {
    std::string s = "[";
    for (int i = 0; i < a.words.size(); i++)
        s += std::format(" {}", a.words[i]);
    s += " ]";
    return s;
}

constexpr integer& __mag_add(integer& a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a.__mag_add_word(static_cast<uint64_t>(b));
    static_assert(sizeof(b) <= 16);
    return a.__mag_add_word(static_cast<uint128_t>(b));
}

constexpr integer& __mag_sub(integer& a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a.__mag_sub_word(static_cast<uint64_t>(b));
    static_assert(sizeof(b) <= 16);
    return a.__mag_sub_word(static_cast<uint128_t>(b));
}

constexpr integer& __mag_mul(integer& a, std_unsigned_int auto b) {
    if (b <= UINT64_MAX) {
        a.mul_add(static_cast<uint64_t>(b), 0);
    } else {
        static_assert(sizeof(b) == 16 || sizeof(b) <= 8);
        __mul(a, static_cast<uint128_t>(b), a);
    }
    return a;
}

constexpr integer& __mag_div(integer& a, std_unsigned_int auto b) {
    Check(b > 0, "division of a magnitude by zero");
    __mag_div(a, static_cast<uint64_t>(b), a);
    return a;
}



constexpr int integer::str(char* buffer, int buffer_size, unsigned base, const bool upper) const {
    Check(base >= 2, "str() with base less than 2");
    Check(base <= 36, "str() with base greater than 36");
    char* p = buffer;
    const char* end = buffer + buffer_size;

    // digits go in least significant first and the whole range is reversed at the end, so the
    // minus sign is written last
    const bool negative = is_negative();
    const char A = (upper ? 'A' : 'a') - 10;
    if (words.size() == 0) {
        Check(p < end, "buffer too small");
        *p++ = '0';
    } else if (words.size() == 1) {
        uint64_t a = low_word();
        if (base == 10 && buffer_size >= 20 * words.size()) {
            while (a) {
                *p++ = '0' + int(a % 10);
                a /= 10;
            }
        } else if (base == 16 && buffer_size >= 16 * words.size()) {
            uint64_t w = low_word();
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
            // divide by 10**19 (the largest power of 10 that fits into uint64_t) and split
            // the remainder with 64-bit arithmetic: 19x fewer full precision divisions
            constexpr uint64_t CHUNK = 10'000'000'000'000'000'000ull;
            integer a = *this;
            a.words.set_negative(false);
            while (a) {
                uint64_t rem = __mag_div(a, CHUNK, /*out*/a);
                if (a) {
                    for (int i = 0; i < 19; i++) {
                        *p++ = '0' + int(rem % 10);
                        rem /= 10;
                    }
                } else {
                    while (rem) {
                        *p++ = '0' + int(rem % 10);
                        rem /= 10;
                    }
                }
            }
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
            // same chunking as above, with the largest power of base that fits into uint64_t
            uint64_t chunk = base;
            int chunk_digits = 1;
            while (chunk <= UINT64_MAX / base) {
                chunk *= base;
                chunk_digits += 1;
            }
            integer n = *this;
            n.words.set_negative(false);
            while (n) {
                uint64_t rem = __mag_div(n, chunk, /*out*/n);
                for (int i = 0; (n && i < chunk_digits) || (!n && rem); i++) {
                    const int c = rem % base;
                    Check(p < end, "buffer too small");
                    *p++ = (c < 10) ? ('0' + c) : (A + c);
                    rem /= base;
                }
            }
        }
    }
    if (negative) {
        Check(p < end, "buffer too small");
        *p++ = '-';
    }
    std::reverse(buffer, p);
    return p - buffer;
}







// a = op(a, b), for a bitwise op. GROW is false when no word of the result can be
// longer than A, which is the case for AND.
template<bool GROW, typename Op>
constexpr integer& __bitwise_op(integer& a, const integer& b, Op op) {
    const int bs = b.words.size();
    if constexpr (GROW) {
        if (bs > a.words.size())
            a.words.resize(bs);
        for (int i = 0; i < bs; i++)
            a.words[i] = op(a.words[i], b.words[i]);
    } else {
        const int n = std::min<int>(a.words.size(), bs);
        if (n == 0) {
            a.words.set_zero(); // downsize(0) alone would leave a stale word behind
            return a;
        }
        a.words.downsize(n);
        for (int i = 0; i < n; i++)
            a.words[i] = op(a.words[i], b.words[i]);
    }
    a.words.normalize();
    return a;
}

constexpr integer& operator|=(integer& a, const integer& b) {
    return __bitwise_op<true>(a, b, [](uint64_t x, uint64_t y) { return x | y; });
}
constexpr integer& operator&=(integer& a, const integer& b) {
    return __bitwise_op<false>(a, b, [](uint64_t x, uint64_t y) { return x & y; });
}
constexpr integer& operator^=(integer& a, const integer& b) {
    return __bitwise_op<true>(a, b, [](uint64_t x, uint64_t y) { return x ^ y; });
}

constexpr integer operator|(integer a, const integer& b) { return a |= b; }
constexpr integer operator&(integer a, const integer& b) { return a &= b; }
constexpr integer operator^(integer a, const integer& b) { return a ^= b; }

constexpr integer& operator|=(integer& a, uint64_t b) {
    if (a.words.size() == 0)
        a = b;
    else
        a.words[0] |= b;
    return a;
}

constexpr integer& operator&=(integer& a, uint64_t b) {
    if (a.words.size() == 0)
        return a;
    a.words.downsize(1);
    a.words[0] &= b;
    a.words.normalize();
    return a;
}

constexpr integer& operator^=(integer& a, uint64_t b) {
    if (a.words.size() == 0) {
        a = b;
    } else {
        a.words[0] ^= b;
        a.words.normalize();
    }
    return a;
}

constexpr integer operator|(integer a, uint64_t b) { return a |= b; }
constexpr integer operator&(integer a, uint64_t b) { return a &= b; }
constexpr integer operator^(integer a, uint64_t b) { return a ^= b; }

// Below this divisor size (in words) recursive division falls back to schoolbook division.
constexpr int BZ_BASE_CASE_SIZE = 8;

// returns the low N words of A, i.e. A mod 2**(64*n)
constexpr integer __low_words(const integer& a, int n) {
    const int s = std::min<int>(n, a.words.size());
    integer c;
    c.words.reset(s, /*initialize*/false);
    for (int i = 0; i < s; i++)
        c.words[i] = a.words[i];
    c.words.normalize();
    return c;
}

// A / B, assuming B has 2n words with its top bit set, A has at most 3n words and A / B < 2**(64*n)
constexpr void __divide_3n2n(const integer& a, const integer& b, int n, integer& q, integer& r);

// A / B, assuming B has n words with its top bit set, A has at most 2n words and (A >> 64*n) < B
constexpr void __divide_2n1n(const integer& a, const integer& b, int n, integer& q, integer& r) {
    if ((n & 1) || n <= BZ_BASE_CASE_SIZE) {
        div(a, b, /*out*/q, /*out*/r);
        return;
    }
    const int h = n / 2;

    // A = [a1 a2 a3 a4] and B = [b1 b2], each part of h words
    integer q1, r1;
    __divide_3n2n(a >> (64 * h), b, h, /*out*/q1, /*out*/r1); // [a1 a2 a3] / B

    integer t = r1 << (64 * h);
    t += __low_words(a, h); // [r1 a4]

    integer q2;
    __divide_3n2n(t, b, h, /*out*/q2, /*out*/r);

    q = q1 << (64 * h);
    q += q2;
}

constexpr void __divide_3n2n(const integer& a, const integer& b, int n, integer& q, integer& r) {
    // A = [a1 a2 a3] and B = [b1 b2], each part of n words
    const integer b1 = b >> (64 * n);
    const integer a12 = a >> (64 * n);

    integer r1;
    if ((a12 >> (64 * n)) < b1) { // a1 < b1
        __divide_2n1n(a12, b1, n, /*out*/q, /*out*/r1);
    } else {
        // quotient would not fit into n words, so use the largest value that does
        q = 1;
        q <<= 64 * n;
        q -= 1u;
        // r1 = [a1 a2] - b1 * q = [a1 a2] - (b1 << 64*n) + b1
        r1 = a12;
        r1 -= b1 << (64 * n);
        r1 += b1;
    }

    // r = [r1 a3] - q * b2, with at most two corrections
    r = r1 << (64 * n);
    r += __low_words(a, n);
    const integer qb2 = q * __low_words(b, n);
    while (r < qb2) {
        q -= 1u;
        r += b;
    }
    r -= qb2;
}

// A / D, using recursive (Burnikel-Ziegler) division
constexpr void divide_bz(const integer& a, const integer& d, integer& q, integer& r) {
    Check(!d.words.empty(), "division by zero");
    const int n = d.words.size();
    if (n <= BZ_BASE_CASE_SIZE || a.words.size() <= n) {
        div(a, d, /*out*/q, /*out*/r);
        return;
    }

    // normalize, so that the top bit of the divisor is set
    const int shift = std::countl_zero(d.words.back());
    const integer an = a << shift;
    const integer dn = d << shift;

    // divide chunk by chunk, most significant chunk first
    const int chunks = (an.words.size() + n - 1) / n;
    q.set_zero();
    r.set_zero();
    integer qi, t;
    for (int i = chunks; i-- > 0;) {
        t = r << (64 * n);
        t += __low_words(an >> (64 * n * i), n);
        __divide_2n1n(t, dn, n, /*out*/qi, /*out*/r);
        q <<= 64 * n;
        q += qi;
    }
    r >>= shift;
}

namespace literals {
constexpr auto operator""_n(const char* s) { return integer(s); }
}

}

template <>
struct std::formatter<algebra::integer, char> {
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
        algebra::Check(it != ctx.end() && *it == '}', "Invalid format specifier for integer.");
        return it;
    }

    // shared by std::formatter<integer>
    constexpr auto format_padded(const auto& a, auto& ctx) const {
        auto it = ctx.out();
        const int bound = a.str_size_upper_bound(base);

        char c_array[100];
        std::string str;
        char* buffer;
        if (bound <= 100) {
            buffer = c_array;
        } else {
            str.resize(bound);
            buffer = str.data();
        }

        const int n = a.str(buffer, bound, base, upper);
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

    constexpr auto format(const algebra::integer& a, auto& ctx) const { return format_padded(a, ctx); }
};

constexpr std::ostream& operator<<(std::ostream& os, const algebra::integer& a) { return os << a.str(); }

namespace algebra {

constexpr bool integer::is_uint8() const { return is_uint32() && words[0] <= 255; }
constexpr bool integer::is_uint16() const { return is_uint32() && words[0] <= 65535; }

constexpr int integer::str_size_upper_bound(unsigned base) const {
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
    default: m = 16; // enough for any base above 16, avoids a dependency on log_upper()
    }
    return words.size() * m + (is_negative() ? 1 : 0); // room for the minus sign
}

template<std::floating_point T>
constexpr integer::operator T() const {
    const int exponent = static_cast<int>(num_bits()) - std::numeric_limits<T>::digits;
    T a;
    if (exponent >= std::numeric_limits<T>::max_exponent)
        a = std::numeric_limits<T>::infinity();
    else if (exponent <= 0)
        a = words[0];
    else
        a = std::ldexp(static_cast<T>(extract_u64(*this, exponent)), exponent);
    return is_negative() ? -a : a;
}

static_assert(sizeof(integer) == 16);

constexpr integer::integer(std::string_view s, unsigned base) {
    const bool negative = s.size() && s[0] == '-';
    if (negative)
        s = s.substr(1);
    const char* p = s.data();
    const char* end = s.data() + s.size();
    Check(p < end, "expecting digit instead of end of string");
    // parse the magnitude first, then apply the sign, since the digit loops below return early
    struct __apply_sign {
        integer_backend& words;
        bool negative;
        constexpr ~__apply_sign() {
            if (negative && words.size())
                words.negate();
        }
    } __sign{words, negative};

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


constexpr integer abs(integer a) { a.words.set_negative(false); return a; }

class magnitude_ref {
    integer_backend& _words;
    bool _negative;
public:
    integer value;
    constexpr magnitude_ref(integer_backend& w) : _words(w), _negative(w.sign() < 0) {
        value.words.swap(w);
        value.words.set_negative(false);
    }
    magnitude_ref(const magnitude_ref&) = delete;
    magnitude_ref(magnitude_ref&&) = delete;
    magnitude_ref& operator=(const magnitude_ref&) = delete;
    magnitude_ref& operator=(magnitude_ref&&) = delete;
    constexpr ~magnitude_ref() {
        _words.swap(value.words);
        _words.set_negative(_negative && _words.size() != 0);
    }
    constexpr integer* operator->() { return &value; }
    constexpr integer& operator*() { return value; }
};

constexpr magnitude_ref magnitude(integer& a) { return magnitude_ref(a.words); }

//@@SIGN_AWARE_LAYER_BEGINS@@


struct integer;




constexpr bool operator==(const integer& a, const integer& b) {
    return a.is_negative() == b.is_negative() && __equal(static_cast<cnatural>(a), static_cast<cnatural>(b));
}

constexpr bool operator==(const integer& a, std_int auto b) {
    if (b < 0)
        return a.is_negative() && __equal_u(static_cast<cnatural>(a), abs_unsigned(b));
    return !a.is_negative() && __equal_u(static_cast<cnatural>(a), make_unsigned(b));
}

constexpr void negate(integer& a) { a.negate(); }
constexpr integer operator-(integer a) { a.negate(); return a; }

template<bool plus>
constexpr integer& __add(integer& a, const integer& b) {

    if (a.words.size() < b.words.size())
        a.words.resize(b.words.size());
    bool a_neg = a.is_negative();
    const uint64_t carry = __add_and_return_carry(a, a_neg, b, plus == b.is_negative());
    if (carry)
        a.words.push_back(carry);
    else
        a.words.normalize();
    a.words.set_negative(a_neg);
    return a;
}

constexpr integer& operator+=(integer& a, const integer& b) { return __add<true>(a, b); }
constexpr integer& operator-=(integer& a, const integer& b) { return __add<false>(a, b); }

template<bool plus>
constexpr integer& __add(integer& a, std_int auto b) {
    auto ub = make_unsigned(b);
    if (b < 0) // a + b == a - abs(b) and a - b == a + abs(b)
        return __add<!plus>(a, static_cast<decltype(ub)>(~ub + 1));

    if (plus == !a.is_negative()) {
        __mag_add(*magnitude(a), ub);
        return a;
    }

    if (!__less_u(static_cast<cnatural>(a), ub)) {
        __mag_sub(*magnitude(a), ub);
        return a;
    }
    // here plus == a.is_negative(), so the result of the magnitude subtraction is
    // negative only for (a >= 0) - b
    a = ub - static_cast<decltype(ub)>(abs(a));
    if constexpr (!plus)
        a.negate();
    return a;
}

constexpr integer& operator+=(integer& a, std_int auto b) { return __add<true>(a, b); }
constexpr integer& operator-=(integer& a, std_int auto b) { return __add<false>(a, b); }

constexpr integer operator+(integer a, std_int auto b) { return a += b; }
constexpr integer operator+(std_int auto a, integer b) { return b += a; }
constexpr integer operator+(integer a, const integer& b) { return a += b; }

constexpr integer operator-(integer a, const integer& b) { return a -= b; }
constexpr integer operator-(integer a, std_int auto b) { return a -= b; }
constexpr integer operator-(std_int auto a, integer b) { b -= a; return -b; }

constexpr void mul(const integer& a, const integer& b, integer& c) {
    const bool negative = a.is_negative() != b.is_negative();
    // views for the operands and a borrow for the result, so nothing is copied. c may alias a or b,
    // so the views are taken before the borrow empties c.
    const cnatural ca = a, cb = b;
    {
        auto m = magnitude(c);
        mul(ca, cb, *m);
    }
    c.words.set_negative(negative);
}

constexpr void mul(integer& a, const integer& b) {
    const bool negative = a.is_negative() != b.is_negative();
    if (&a == &b) {
        auto m = magnitude(a);
        __mag_mul(*m, *m); // squaring: __mag_mul() takes the &a == &b path
    } else {
        // b is a distinct object, so a view of it stays valid while a's magnitude is borrowed
        const cnatural cb = b;
        natural out;
        mul(static_cast<cnatural>(a), cb, out);
        a.words = std::move(out.words);
    }
    a.words.set_negative(negative);
}

constexpr integer operator*(const integer& a, const integer& b) {
    integer c;
    mul(a, b, c);
    return c;
}


// TODO avoid memory allocation here for int128!
constexpr integer operator*(const integer& a, std_int auto b) {
    integer c;
    mul(a, integer(b), c);
    return c;
}

constexpr integer& operator*=(integer& a, const integer& b) {
    mul(a, b);
    return a;
}


constexpr integer& operator*=(integer& a, std_int auto b) {
    const bool negative = a.is_negative() != (b < 0);
    __mag_mul(*magnitude(a), abs_unsigned(b));
    a.words.set_negative(negative);
    return a;
}

constexpr std::string str(const integer& a) {
    const std::string m = str(static_cast<cnatural>(a));
    return a.is_negative() ? "-" + m : m;
}

constexpr std::string stre(const integer& a) {
    const std::string m = __mag_stre(a);
    return a.is_negative() ? "-" + m : m;
}

// The magnitude products, done on integer's own backend. These mirror natural's add_product and
// sub_product, which are themselves thin wrappers over the same kernels, so nothing is copied and
// no natural temporary is built -- which matters because sub_product sits on the division path.
constexpr void __magnitude_add_product(integer& a, cnatural b, cnatural c) {
    if (b.size == 0 || c.size == 0)
        return;
    const int A = a.words.size();
    a.words.resize(std::max<int>(A, b.size + c.size) + 1);
    vnatural va{{a.words.data(), A}, a.words.capacity()};
    if (b.size < c.size)
        __add_product(va, b, c);
    else
        __add_product(va, c, b);
    a.words.downsize(va.size);
}

constexpr void __magnitude_sub_product(integer& a, cnatural b, cnatural c) {
    if (b.size == 0 || c.size == 0)
        return;
    inatural ia{a.words.data(), a.words.size()};
    const bool ok = (b.size < c.size) ? __sub_product(ia, b, c) : __sub_product(ia, c, b);
    Check(ok, "sub_product() assumes A >= B * C");
    a.words.downsize(ia.size);
}

constexpr void __magnitude_add_product(integer& a, cnatural b, const uint64_t c) {
    if (b.size == 0 || c == 0)
        return;
    const int A = a.words.size();
    a.words.resize(std::max<int>(A, b.size + 1) + 1);
    vnatural va{{a.words.data(), A}, a.words.capacity()};
    __add_product(va, b, c);
    a.words.downsize(va.size);
}

constexpr void __magnitude_sub_product(integer& a, cnatural b, const uint64_t c) {
    if (b.size == 0 || c == 0)
        return;
    inatural ia{a.words.data(), a.words.size()};
    Check(__sub_product(ia, b, c), "sub_product() assumes A >= B * c");
    a.words.downsize(ia.size);
}

constexpr void __magnitude_complement(integer& a) {
    __complement(inatural{a.words.data(), a.words.size()});
    a.words.normalize();
}

template<bool plus>
constexpr void __add_product(integer& a, const integer& b, const integer& c) {
    if (b.words.empty() || c.words.empty())
        return;
    const bool a_negative = a.is_negative();
    const bool bc_negative = b.is_negative() != c.is_negative();

    const cnatural cb = b, cc = c;
    if ((plus && a_negative == bc_negative) || (!plus && a_negative != bc_negative)) {
        __magnitude_add_product(a, cb, cc);
        a.words.set_negative(a_negative);
    } else if (a.num_bits() > b.num_bits() + c.num_bits()) {
        __magnitude_sub_product(a, cb, cc);
        a.words.set_negative(a_negative);
    } else {
        const int m = mul_max_size(cb, cc);
        a.words.resize(m + 1);
        a.words[m] = 1;
        __magnitude_sub_product(a, cb, cc);
        if (a.words.size() > m) {
            a.words[m] -= 1;
            a.words.normalize();
        } else {
            a.words.resize(m);
            __magnitude_complement(a);
            a.words.set_negative(!a_negative);
        }
    }
}

constexpr void add_product(integer& a, const integer& b, const integer& c) { __add_product<true>(a, b, c); }
constexpr void sub_product(integer& a, const integer& b, const integer& c) { __add_product<false>(a, b, c); }

constexpr int __signum(const integer& a) {
    if (a.is_negative())
        return -1;
    return a.is_zero() ? 0 : 1;
}

template<bool plus>
constexpr void __add_product(integer& a, const integer& b, const uint64_t cu, const bool c_negative) {
    if (b.words.empty() || cu == 0)
        return;
    const bool a_negative = a.is_negative();
    const bool bc_negative = b.is_negative() != c_negative;

    if ((plus && a_negative == bc_negative) || (!plus && a_negative != bc_negative)) {
        __magnitude_add_product(a, static_cast<cnatural>(b), cu);
        a.words.set_negative(a_negative);
    } else if (a.num_bits() > b.num_bits() + num_bits(cu)) {
        __magnitude_sub_product(a, static_cast<cnatural>(b), cu);
        a.words.set_negative(a_negative);
    } else {
        const int m = mul_max_size(static_cast<cnatural>(b), {&cu, 1});
        a.words.resize(m + 1);
        a.words[m] = 1;
        __magnitude_sub_product(a, static_cast<cnatural>(b), cu);
        if (a.words.size() > m) {
            a.words[m] -= 1;
            a.words.normalize();
        } else {
            a.words.resize(m);
            __magnitude_complement(a);
            a.words.set_negative(!a_negative);
        }
    }
}

constexpr void add_product(integer& a, const integer& b, const std_int auto c) { __add_product<true>(a, b, abs_unsigned(c), c < 0); static_assert(sizeof(c) <= 8); }
constexpr void sub_product(integer& a, const integer& b, const std_int auto c) { __add_product<false>(a, b, abs_unsigned(c), c < 0); static_assert(sizeof(c) <= 8); }

constexpr void div(const integer& a, const integer& b, integer& quot, integer& rem) {
    if (b == 1) {
        quot = a;
        rem = 0;
        return;
    }
    if (b == -1) {
        quot = -a;
        rem = 0;
        return;
    }
    const bool a_negative = a.is_negative();
    const bool negative = a.is_negative() != b.is_negative();
    // Take the operand magnitudes before borrowing the outputs: quot or rem may be the same object
    // as a or b, and a borrow leaves the borrowed value's words empty while it is alive.
    const natural an = abs(a);
    const natural bn = abs(b);
    {
        auto q = magnitude(quot);
        auto r = magnitude(rem);
        __mag_div(an, bn, *q, *r);
    }
    quot.words.set_negative(negative);
    rem.words.set_negative(a_negative);
}

constexpr integer operator/(const integer& a, const integer& b) {
    integer quot, rem;
    div(a, b, quot, rem);
    return quot;
}

// The remainder takes the sign of a, as for the integer overload. Signed and unsigned divisors
// are separate templates so that a divisor above INT64_MAX is never read as a negative one.
template<std_signed_int T>
constexpr T div(const integer& a, T b, integer& quot) {
    if (b == 1) {
        quot = a;
        return 0;
    }
    if (b == -1) {
        quot = -a;
        return 0;
    }
    uint64_t rem;
    {
        const natural an = abs(a); // quot may alias a, see div() above
        auto q = magnitude(quot);
        rem = __mag_div(an, abs_unsigned(b), *q);
    }
    if (!quot.is_zero())
        quot.words.set_negative(a.is_negative() != (b < 0));
    return a.is_negative() ? -static_cast<T>(rem) : static_cast<T>(rem);
}

template<std_unsigned_int T>
constexpr T div(const integer& a, T b, integer& quot) {
    if (b == 1) {
        quot = a;
        return 0;
    }
    uint64_t rem;
    {
        const natural an = abs(a);
        auto q = magnitude(quot);
        rem = __mag_div(an, b, *q);
    }
    if (!quot.is_zero())
        quot.words.set_negative(a.is_negative());
    return static_cast<T>(rem);
}

constexpr integer operator/(const integer& a, const std_int auto b) {
    integer c = a;
    c.words.set_negative(false);
    __mag_div(c, abs_unsigned(b));
    c.words.set_negative(a.is_negative() != (b < 0));
    return c;
}

constexpr integer& operator/=(integer& a, const integer& b) {
    integer rem;
    div(a, b, a, rem);
    return a;
}
constexpr integer& operator/=(integer& a, const std_int auto b) {
    const bool negative = a.is_negative();
    __mag_div(*magnitude(a), abs_unsigned(b));
    a.words.set_negative(negative != (b < 0));
    return a;
}

constexpr integer operator%(const integer& a, const integer& divisor) {
    // the remainder alone: going through div() would build a quotient only to discard it
    const bool negative = a.is_negative();
    natural r;
    __mod_into(static_cast<cnatural>(a), static_cast<cnatural>(divisor), r);
    integer c = std::move(r);
    c.words.set_negative(negative && !c.is_zero());
    return c;
}

// TODO generalize for any std_int
constexpr int64_t operator%(const integer& a, int64_t b) {
    Check(b != 0, "division of integer by zero");
    uint64_t m = __mod(static_cast<cnatural>(a), abs_unsigned(b));
    return (a.sign() >= 0) ? m : -static_cast<int64_t>(m);
}

constexpr int operator%(const integer& a, int b) { return a % (int64_t)b; }
constexpr int64_t operator%(const integer& a, unsigned b) { return a % (int64_t)b; }

// Note: return type is integer instead of uint64_t, as it can be negative (can't fit into int64_t either)
constexpr integer operator%(const integer& a, uint64_t b) {
    Check(b > 0, "division of integer by zero");
    integer c = __mod(static_cast<cnatural>(a), b);
    if (a.is_negative())
        c.negate();
    return c;
}

constexpr integer mod(const integer& a, const integer& b) {
    // Note: mod() on the magnitudes would resolve back to this function, since natural
    // converts implicitly to integer. Call the natural kernel explicitly.
    natural r;
    mod(abs(a), abs(b), /*out*/r);
    if (a.is_negative() && !r.words.empty()) {
        // result is in [0, abs(b)) range
        natural e = abs(b);
        e -= r;
        return e;
    }
    return r;
}

constexpr uint64_t mod(const integer& a, uint64_t b) {
    uint64_t m = abs(a) % b;
    return (a.is_negative() && m) ? (b - m) : m;
}

// TODO this version doesn't use uint128_t %, is it faster?
// TODO if it is, move it to natural.h instead, it doesn't belong here
constexpr unsigned mod(const integer& a, uint32_t b) {
    if (b == 0)
        throw std::runtime_error("division by zero");
    uint64_t m = (uint64_t(1) << 32) % b;
    m = (m * m) % b;

    uint64_t acc = 0;
    for (auto i = a.words.size(); i-- > 0;) {
        acc *= m;
        acc += a.words[i] % b;
        acc %= b;
    }
    if (a.is_negative()) {
        acc *= b - 1;
        acc %= b;
    }
    return acc;
}

constexpr integer& operator%=(integer& a, const integer& b) { a = a % b; return a; }
// TODO issue temporary memory allocation for cent / ucent
constexpr integer& operator%=(integer& a, std_int auto b) { a = a % integer(b); return a; }

constexpr bool operator<(const integer& a, const integer& b) {
    if (a.is_negative())
        return !b.is_negative() || __less(static_cast<cnatural>(b), static_cast<cnatural>(a));
    return !b.is_negative() && __less(static_cast<cnatural>(a), static_cast<cnatural>(b));
}
// TODO issue temporary memory allocation for cent / ucent
constexpr bool operator<(const integer& a, std_int auto b) { return a < integer(b); }
constexpr bool operator<(std_int auto a, const integer& b) { return integer(a) < b; }

constexpr integer operator~(integer a) {
    if (a.sign() >= 0) {
        a += 1;
        a.negate();
        return a;
    }
    a.negate();
    a -= 1;
    return a;
}

namespace literals {
constexpr auto operator""_i(const char* s) { return integer(s); }
}

constexpr void operator<<=(integer& a, int64_t i) {
    const bool negative = a.is_negative();
    __mag_shl(*magnitude(a), i);
    a.words.set_negative(negative);
}


static_assert(sizeof(integer) == 16);

}


#if 0
template <>
struct std::formatter<algebra::neg_integer, char> {
    constexpr auto parse(auto& ctx) {
        auto it = ctx.begin();
        if (it == ctx.end() || *it != '}')
            throw std::format_error("Invalid format specifier for integer.");
        return it;
    }

    constexpr auto format(const algebra::neg_integer& a, auto& ctx) const {
        if (a.a->sign() >= 0)
            return std::format_to(ctx.out(), "{}", a.a->abs);
        return std::format_to(ctx.out(), "-{}", a.a->abs);
    }
};
#endif


template<>
struct std::hash<algebra::integer> {
    constexpr size_t operator()(const algebra::integer& a) const {
        // the backend hash includes the sign in the size, which is what makes -x hash apart from x
        return std::hash<algebra::integer_backend>()(a.words);
    }
};
