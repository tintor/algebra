#pragma once
#include "algebra/natural.h"

namespace algebra {

struct integer;
template<> struct IsNumberClass<integer> : std::true_type {};

#if 0
struct neg_integer {
    const integer* a;
};
#endif

struct integer {
    using size_type = natural::size_type;

    natural abs;

    constexpr integer() {}
    constexpr integer(std_int auto a) : abs(abs_unsigned(a)) { if (a < 0) negate(); }
    constexpr integer(integer&& o) : abs(std::move(o.abs)) { }
    constexpr integer(natural&& o) : abs(std::move(o)) { abs.words.set_negative(false); }
    constexpr integer(const integer& o) : abs(o.abs) { }
    constexpr integer(const natural& o) : abs(o) { abs.words.set_negative(false); }

#if 0
    constexpr integer(const neg_integer& o) : integer(*o.a) { negate(); }
    constexpr void operator=(const neg_integer& o) {
        if (o.a != this)
            operator=(*o.a);
        negate();
    }
#endif

    constexpr void operator=(std_int auto a) { abs.words = a; }
    constexpr void operator=(integer&& o) { abs = std::move(o.abs); }
    constexpr void operator=(natural&& o) { abs = std::move(o); abs.words.set_negative(false); }
    constexpr void operator=(const integer& o) { abs = o.abs; }
    constexpr void operator=(const natural& o) { abs = o; abs.words.set_negative(false); }

    constexpr size_type sign() const { return abs.words.sign(); }
    constexpr bool is_negative() const { return sign() < 0; }
    constexpr bool is_even() const { return abs.is_even(); }
    constexpr bool is_odd() const { return abs.is_odd(); }
    constexpr bool is_one() const { return abs.words.size() == 1 && abs.words[0] == 1 && sign() >= 0; }
    constexpr bool is_zero() const { return abs.words.size() == 0; }

    constexpr bool is_int8() const {
        if (abs.words.size() > 1)
            return false;
        if (abs.words[0] <= 127)
            return true;
        return sign() < 0 && abs.words[0] == 128;
    }

    constexpr bool is_int16() const {
        if (abs.words.size() > 1)
            return false;
        if (abs.words[0] <= INT16_MAX)
            return true;
        return sign() < 0 && abs.words[0] == static_cast<uint64_t>(INT16_MAX) + 1;
    }

    constexpr bool is_int32() const {
        if (abs.words.size() > 1)
            return false;
        if (abs.words[0] <= INT32_MAX)
            return true;
        return sign() < 0 && abs.words[0] == static_cast<uint64_t>(INT32_MAX) + 1;
    }

    constexpr bool is_int64() const {
        if (abs.words.size() > 1)
            return false;
        if (abs.words[0] <= INT64_MAX)
            return true;
        return sign() < 0 && abs.words[0] == static_cast<uint64_t>(INT64_MAX) + 1;
    }

    constexpr bool is_int128() const {
        if (abs.words.size() > 2)
            return false;
        if (abs.words.size() < 2)
            return true;

        uint64_t w = abs.words[1];
        if ((w & (uint64_t(1) << 63)) == 0)
            return true;
        return sign() < 0 && w == uint64_t(1) << 63;
    }

    constexpr bool is_uint8() const { return sign() >= 0 && abs.is_uint8(); }
    constexpr bool is_uint16() const { return sign() >= 0 && abs.is_uint16(); }
    constexpr bool is_uint32() const { return sign() >= 0 && abs.is_uint32(); }
    constexpr bool is_uint64() const { return sign() >= 0 && abs.is_uint64(); }
    constexpr bool is_uint128() const { return sign() >= 0 && abs.is_uint128(); }

    constexpr integer(std::string_view s, unsigned base = 10) : abs((s.size() && s[0] == '-') ? s.substr(1) : s) {
        if (s.size() && s[0] == '-')
            abs.words.negate();
    }
    constexpr integer(const char* s, unsigned base = 10) : integer(std::string_view(s), base) {}

    constexpr operator char() const { return operator int(); }
    constexpr operator uint8_t() const { return operator uint32_t(); }
    constexpr operator short() const { return operator int(); }
    constexpr operator uint16_t() const { return operator uint32_t(); }
    constexpr operator int() const { return is_negative() ? -static_cast<int>(abs.words[0]) : abs.words[0]; }
    constexpr operator unsigned() const { return abs.operator uint32_t(); }
    constexpr operator long() const { return is_negative() ? -static_cast<long>(abs.words[0]) : abs.words[0]; }
    constexpr operator unsigned long() const { return abs.operator unsigned long(); }
    constexpr operator unsigned long long() const { return abs.operator unsigned long long(); }
    static_assert(sizeof(long) == 8);
    static_assert(sizeof(long long) == 8);

    constexpr operator __int128() const {
        if (sign() == 2) return (uint128_t(abs.words[1]) << 64) | abs.words[0];
        if (sign() == 1) return uint128_t(abs.words[0]);
        if (sign() == -1) return -uint128_t(abs.words[0]);
        if (sign() == -2) return -((uint128_t(abs.words[1]) << 64) | abs.words[0]);
        return 0;
    }

    constexpr operator unsigned __int128() const { return abs.operator unsigned __int128(); }

    constexpr int str_size_upper_bound(unsigned base = 10) const { return is_negative() + abs.str_size_upper_bound(base); }
    constexpr int str(char* buffer, int buffer_size, unsigned base = 10, bool upper = true) const {
        int result = 0;
        if (is_negative()) {
            if (buffer_size <= 0)
                throw std::runtime_error("buffer too small");
            buffer_size -= 1;
            *buffer++ = '-';
            result = 1;
        }
        return result + abs.str(buffer, buffer_size, base, upper);
    }

    constexpr std::string str(unsigned base = 10, bool upper = true) const {
        std::string s;
        s.resize(str_size_upper_bound(base));
        s.resize(str(s.data(), s.size(), base, upper));
        return s;
    }
    constexpr std::string hex() const { return str(16); }

    constexpr void negate() { abs.words.negate(); }

    constexpr size_t popcount() const {
        if (!is_negative())
            return abs.popcount();

        size_t c = 0;
        uint64_t carry = 1;
        for (size_type i = 0; i < abs.words.size(); i++) {
            uint128_t w = (uint128_t)abs.words[i] + carry;
            carry = w >> 64;
            c += std::popcount(~static_cast<uint64_t>(w));
        }
        return c;
    }

    constexpr int size_of() const { return abs.words.size() * 8; }

    constexpr operator bool() const { return sign(); }

    constexpr integer& operator++() {
        if (!is_negative())
            abs += uint64_t(1);
        else
            abs -= uint64_t(1);
        return *this;
    }

    constexpr integer& operator--() {
        if (is_negative())
            abs += uint64_t(1);
        else if (abs.words.size() > 0)
            abs -= uint64_t(1);
        else
            *this = -1;
        return *this;
    }

    constexpr integer operator++(int) { integer a = *this; operator++(); return a; }
    constexpr integer operator--(int) { integer a = *this; operator--(); return a; }

    constexpr auto num_bits() const { return abs.num_bits(); }
    constexpr auto num_trailing_zeros() const { return abs.num_trailing_zeros(); }

    template<std::floating_point T>
    constexpr operator T() const {
        auto a = static_cast<T>(abs);
        return (sign() < 0) ? -a : a;
    }

    constexpr void swap(integer& o) { abs.swap(o.abs); }
};

constexpr void negate(integer& a) { a.negate(); }

#if 0
constexpr neg_integer operator-(const integer& a) { return {&a}; }
#else
constexpr integer operator-(integer a) { a.negate(); return a; }
#endif

constexpr void add(const integer& a, const integer& b, integer& c) {
    if (a.is_negative() == b.is_negative()) {
        c.abs = a.abs + b.abs;
        return;
    }
    if (b.abs < a.abs) {
        c.abs = a.abs - b.abs;
        c.abs.words.set_negative(a.is_negative());
        return;
    }
    c.abs = b.abs - a.abs;
    c.abs.words.set_negative(b.is_negative());
}

constexpr void sub(const integer& a, const integer& b, integer& c) {
    if (a.is_negative() != b.is_negative()) {
        c.abs = a.abs + b.abs;
        c.abs.words.set_negative(a.is_negative());
        return;
    }
    if (b.abs < a.abs) {
        c.abs = a.abs - b.abs;
        c.abs.words.set_negative(a.is_negative());
        return;
    }
    c.abs = b.abs - a.abs;
    c.abs.words.set_negative(b.is_negative());
    c.negate();
}

constexpr integer operator+(const integer& a, const integer& b) {
    integer c;
    add(a, b, c);
    return c;
}

constexpr integer operator-(const integer& a, const integer& b) {
    integer c;
    sub(a, b, c);
    return c;
}

constexpr integer operator+(integer a, std_int auto b) { return a += b; }
constexpr integer operator+(std_int auto a, integer b) { return b += a; }

constexpr integer operator-(integer a, std_int auto b) { return a -= b; }
constexpr integer operator-(std_int auto a, integer b) { b -= a; return -b; }

constexpr integer& operator+=(integer& a, const integer& b) {
    if (a.is_negative() == b.is_negative()) {
        a.abs += b.abs;
        return a;
    }
    if (b.abs < a.abs) {
        a.abs -= b.abs;
        return a;
    }
    a.abs = b.abs - a.abs; // TODO optimize this temporary
    a.abs.words.set_negative(b.is_negative());
    return a;
}

constexpr integer& operator-=(integer& a, const integer& b) {
    if (a.is_negative() != b.is_negative()) {
        a.abs += b.abs;
        return a;
    }
    if (b.abs < a.abs) {
        a.abs -= b.abs;
        return a;
    }
    a.abs = b.abs - a.abs; // TODO optimize this temporary
    a.abs.words.set_negative(b.is_negative());
    a.negate();
    return a;
}

constexpr integer& operator+=(integer& a, std_int auto b) {
    auto ub = make_unsigned(b);
    if (b < 0)
        return a -= ~ub + 1;

    if (!a.is_negative()) {
        a.abs += ub;
        return a;
    }

    if (a.abs >= ub) {
        a.abs -= ub;
        return a;
    }
    a = ub - static_cast<decltype(ub)>(a.abs);
    a.negate();
    return a;
}

constexpr integer& operator-=(integer& a, std_int auto b) {
    auto ub = make_unsigned(b);
    if (b < 0)
        return a += ~ub + 1;

    if (a.is_negative()) {
        a.abs += ub;
        return a;
    }

    if (a.abs >= ub) {
        a.abs -= ub;
        return a;
    }
    a = ub - static_cast<decltype(ub)>(a.abs);
    a.negate();
    return a;
}

constexpr bool operator==(const integer& a, const integer& b) { return a.abs == b.abs && a.is_negative() == b.is_negative(); }

constexpr bool operator==(const integer& a, std_int auto b) {
    if (b < 0)
        return a.is_negative() && a.abs == abs_unsigned(b);
    return !a.is_negative() && a.abs == make_unsigned(b);
}

constexpr void mul(const integer& a, const integer& b, integer& c) {
    c.abs = a.abs * b.abs;
    c.abs.words.set_negative(a.is_negative() != b.is_negative());
}

constexpr void mul(integer& a, const integer& b) {
    const bool negative = a.is_negative() != b.is_negative();
    a.abs *= b.abs;
    a.abs.words.set_negative(negative);
}

constexpr integer operator*(const integer& a, const integer& b) {
    integer c;
    mul(a, b, c);
    return c;
}

constexpr integer operator*(const integer& a, const natural& b) {
    integer c;
    c.abs = a.abs * b;
    c.abs.words.set_negative(a.is_negative());
    return c;
}
constexpr integer operator*(const natural& a, const integer& b) { return b * a; }

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

constexpr integer& operator*=(integer& a, const natural& b) {
    const bool negative = a.is_negative();
    a.abs *= b;
    a.abs.words.set_negative(negative);
    return a;
}

constexpr integer& operator*=(integer& a, std_int auto b) {
    const bool negative = a.is_negative() != (b < 0);
    a.abs *= abs_unsigned(b);
    a.abs.words.set_negative(negative);
    return a;
}

constexpr void add_product(integer& acc, const integer& a, const integer& b) {
    // TODO optimize temporary allocation
    acc += a * b;
}

constexpr void sub_product(integer& acc, const integer& a, const integer& b) {
    // TODO optimize temporary allocation
    acc -= a * b;
}

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
    div(a.abs, b.abs, quot.abs, rem.abs);
    quot.abs.words.set_negative(negative);
    rem.abs.words.set_negative(a_negative);
}

constexpr integer operator/(const integer& a, const integer& b) {
    integer quot, rem;
    div(a, b, quot, rem);
    return quot;
}

constexpr int64_t div(const integer& a, int64_t b, integer& quot) {
    if (b == 1) {
        quot = a;
        return 0;
    }
    if (b == -1) {
        quot = -a;
        return 0;
    }
    int64_t rem = div(a.abs, abs_unsigned(b), quot.abs);
    if (quot.abs)
        quot.abs.words.set_negative(a.is_negative() != (b < 0));
    return a.is_negative() ? -rem : rem;
}

constexpr integer operator/(const integer& a, int64_t b) {
    integer quot;
    div(a, b, quot);
    return quot;
}

constexpr integer operator/(const integer& a, int b) { return a / (int64_t)b; }
constexpr integer operator/(const integer& a, unsigned b) { return a / (int64_t)b; }
// TODO ulong?

constexpr void operator/=(integer& a, const integer& b) {
    integer rem;
    div(a, b, a, rem);
}

constexpr void operator/=(integer& a, int64_t b) { div(a, b, a); }
constexpr void operator/=(integer& a, int b) { div(a, (long)b, a); }
constexpr void operator/=(integer& a, unsigned b) { div(a, (long)b, a); }
// TODO ulong?

constexpr integer operator%(const integer& a, const integer& divisor) {
    integer quotient, remainder;
    div(a, divisor, quotient, remainder);
    return remainder;
}

constexpr int64_t operator%(const integer& a, int64_t b) {
    if (b == 0)
        throw std::runtime_error("division by zero");
    uint128_t acc = 0;
    for (integer::size_type i = a.abs.words.size(); i-- > 0;) {
        acc <<= 64;
        acc |= a.abs.words[i];
        acc %= abs_unsigned(b);
    }
    return (a.sign() >= 0) ? acc : -static_cast<int64_t>(acc);
}

constexpr int operator%(const integer& a, int b) { return a % (int64_t)b; }
constexpr int64_t operator%(const integer& a, unsigned b) { return a % (int64_t)b; }

constexpr integer operator%(const integer& a, uint64_t b) {
    if (b == 0)
        throw std::runtime_error("division by zero");
    uint128_t acc = 0;
    for (integer::size_type i = a.abs.words.size(); i-- > 0;) {
        acc <<= 64;
        acc |= a.abs.words[i];
        acc %= b;
    }
    return acc;
}

constexpr uint64_t mod(const integer& a, uint64_t b) {
    if (b == 0)
        throw std::runtime_error("division by zero");
    uint128_t acc = 0;
    for (integer::size_type i = a.abs.words.size(); i-- > 0;) {
        acc <<= 64;
        acc |= a.abs.words[i];
        acc %= b;
    }
    if (a.is_negative()) {
        acc *= b - 1;
        acc %= b;
    }
    return acc;
}

constexpr unsigned mod(const integer& a, unsigned b) {
    if (b == 0)
        throw std::runtime_error("division by zero");
    uint64_t m = (uint64_t(1) << 32) % b;
    m = (m * m) % b;

    uint64_t acc = 0;
    for (integer::size_type i = a.abs.words.size(); i-- > 0;) {
        acc *= m;
        acc += a.abs.words[i] % b;
        acc %= b;
    }
    if (a.is_negative()) {
        acc *= b - 1;
        acc %= b;
    }
    return acc;
}

constexpr integer& operator%=(integer& a, const integer& b) { a = a % b; return a; }

constexpr integer& operator%=(integer& a, std_int auto b) { a = a % integer(b); return a; }

constexpr bool operator<(const integer& a, const integer& b) {
    if (a.is_negative())
        return !b.is_negative() || a.abs > b.abs;
    return !b.is_negative() && a.abs < b.abs;
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
    a.abs <<= i;
    a.abs.words.set_negative(negative);
}

ALGEBRA_SHIFT_OP(integer)

static_assert(sizeof(integer) == 16);

}

template <>
struct std::formatter<algebra::integer, char> : public std::formatter<algebra::natural, char> {
    constexpr auto format(const algebra::integer& a, auto& ctx) const {
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

constexpr std::ostream& operator<<(std::ostream& os, const algebra::integer& a) { return os << a.str(); }

template<>
struct std::hash<algebra::integer> {
    constexpr size_t operator()(const algebra::integer& a) const {
        return std::hash<algebra::natural>()(a.abs);
    }
};
