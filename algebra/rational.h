#pragma once
#include "algebra/integer.h"
#include "algebra/integer_func.h"
#include <regex>
#include <charconv>

namespace algebra {
using namespace algebra;

struct rational {
    integer num, den;

    constexpr rational() : den(1) { }
    constexpr rational(integer a) : num(std::move(a)), den(1) { }
    constexpr rational(integer a, integer b) : num(std::move(a)), den(std::move(b)) { simplify(); }
    // no simplification
    constexpr rational(integer a, integer b, int) : num(std::move(a)), den(std::move(b)) { }

    constexpr rational(std::integral auto a) : num(a), den(1) { }
    // TODO simplify in small integers first, before converting to integer class!
    constexpr rational(std::integral auto a, std::integral auto b) : num(a), den(b) { simplify(); }

    constexpr rational(float x);
    constexpr rational(double x);

    constexpr rational(std::string_view s) {
        static std::regex pattern(R"(([-]?\d+)([./]\d+)?(e[-+]?(\d+))?)");
        std::match_results<std::string_view::const_iterator> match;
        if (!std::regex_match(s.begin(), s.end(), match, pattern))
            throw std::runtime_error("syntax error parsing rational");

        std::string_view group1(match[1].first, match[1].second - match[1].first);
        num = integer(group1);

        std::string_view group2(match[2].first, match[2].second - match[2].first);
        if (group2.empty()) {
            den = 1;
        } else if (group2[0] == '/') {
            den = integer(group2.substr(1));
        } else if (group2[0] == '.') {
            den = pow(integer(10), group2.size() - 1);
            num *= den;
            integer frac(group2.substr(1));
            if (num >= 0) num += frac; else num -= frac;
        } else
            throw std::runtime_error("unreachable");

        std::string_view group3(match[3].first, match[3].second - match[3].first);
        if (!group3.empty()) {
            std::string_view group4(match[4].first, match[4].second - match[4].first);
            unsigned long e;
            std::from_chars(match[4].first, match[4].second, e);
            integer exp = pow(integer(10), e);
            if (group3[1] == '-') den *= exp; else num *= exp;
        }

        simplify();
    }

    constexpr rational(const std::string& s) : rational(std::string_view(s)) { }
    constexpr rational(const char* s) : rational(std::string_view(s)) { }

    constexpr void simplify();

    constexpr operator float() const;
    constexpr operator double() const;

    constexpr void invert() {
        std::swap(num, den);
        if (den.is_negative()) {
            den.negate();
            num.negate();
        }
    }

    constexpr void negate() { num.negate(); }

    constexpr std::string str() const {
        // TODO precalculate number of chars needed and do one allocation only
        auto s = num.str();
        if (!den.is_one()) {
            s += '/';
            s += den.str();
        }
        return s;
    }

    constexpr auto sign() const { return num.sign(); }
    constexpr bool is_integer() const { return den.is_one(); }
    constexpr bool is_even() const { return is_integer() && num.is_even(); }
    constexpr bool is_odd() const { return is_integer() && num.is_odd(); }
};

constexpr void negate(rational& a) { a.negate(); }
constexpr rational fract(const rational&);

constexpr rational operator-(const rational& a) { return {-a.num, a.den}; }

constexpr rational& operator+=(rational& a, const rational& b);
constexpr rational operator+(const rational& a, const rational& b);
constexpr rational operator+(const rational& a, const integer& b);
constexpr rational operator+(const integer& a, const rational& b);

constexpr rational operator+(std::integral auto a, const rational& b) { return integer(a) + b; }
constexpr rational operator+(const rational& a, std::integral auto b) { return a + integer(b); }
constexpr void operator+=(rational& a, std::integral auto b) { a += integer(b); }

constexpr rational& operator-=(rational& a, const rational& b);
constexpr rational operator-(const rational& a, const rational& b);
constexpr rational operator-(const rational& a, const integer& b);
constexpr rational operator-(const integer& a, const rational& b);

constexpr rational operator-(std::integral auto a, const rational& b) { return integer(a) - b; }
constexpr rational operator-(const rational& a, std::integral auto b) { return a - integer(b); }
constexpr void operator-=(rational& a, std::integral auto b) { a -= integer(b); }

constexpr rational operator*(const rational& a, const rational& b);
constexpr rational operator*(const rational& a, const integer& b);
constexpr rational operator*(const integer& a, const rational& b);

constexpr rational operator*(std::integral auto a, const rational& b) { return integer(a) * b; }
constexpr rational operator*(const rational& a, std::integral auto b) { return a * integer(b); }

constexpr rational& operator*=(rational& a, const rational& b);
constexpr rational& operator*=(rational& a, const integer& b);
constexpr void operator*=(rational& a, std::integral auto b) { a *= integer(b); }

constexpr rational& operator/=(rational&, const integer&);
constexpr rational operator/(const rational&, const rational&);
constexpr rational operator/(const rational&, const integer&);
constexpr rational operator/(const integer&, const rational&);

constexpr rational operator/(std::integral auto a, const rational& b) { return (a == 1) ? rational{b.den, b.num} : (integer(a) / b); }
constexpr rational operator/(const rational& a, std::integral auto b) { return a / integer(b); }
constexpr void operator/=(rational& a, std::integral auto b) { a /= integer(b); }

constexpr bool operator<(const rational& a, const rational& b) { return (a.den == b.den) ? (a.num <  b.num) : (a.num * b.den <  b.num * a.den); }
constexpr bool operator>(const rational& a, const rational& b) { return (a.den == b.den) ? (a.num >  b.num) : (a.num * b.den >  b.num * a.den); }

template<typename T>
concept Integral = std::integral<T> || std::same_as<T, integer> || std::same_as<T, natural>;

constexpr bool operator<(Integral auto a, const rational& b) { return b.den.is_one() ? (a <  b.num) : (a * b.den <  b.num); }
constexpr bool operator<(const rational& a, Integral auto b) { return a.den.is_one() ? (a.num <  b) : (a.num <  b * a.den); }

constexpr bool operator==(const rational& a, const rational& b) { return a.num == b.num && a.den == b.den; }
constexpr bool operator==(const rational& a, const Integral auto b) { return a.num == b && a.den.is_one(); }
constexpr bool operator==(const Integral auto a, const rational& b) { return a == b.num && b.den.is_one(); }

namespace literals {
constexpr auto operator""_q(const char* s) { return rational(s); }
}

constexpr rational::rational(float x) {
    if (std::isnan(x))
        throw std::runtime_error("can't convert nan to rational");
    if (std::isinf(x))
        throw std::runtime_error("can't convert infinite to rational");
    if (x == 0.0f) {
        num = 0;
        den = 1;
        return;
    }

    int exponent;
    float mantissa = std::frexp(x, &exponent);

    // Convert mantissa to an exact integer representation
    const int mantissa_bits = 24;
    float scaled_mantissa = std::ldexp(mantissa, mantissa_bits); // mantissa * 2^24
    num = static_cast<long>(scaled_mantissa);
    exponent -= mantissa_bits;

    den = 1;
    if (exponent > 0) num <<= exponent;
    if (exponent < 0) den <<= -exponent;
    if (x < 0.0) num = -num;

    simplify();
}

constexpr rational::rational(double x) {
    if (std::isnan(x))
        throw std::runtime_error("can't convert nan to rational");
    if (std::isinf(x))
        throw std::runtime_error("can't convert infinite to rational");
    if (x == 0.0) {
        num = 0;
        den = 1;
        return;
    }

    int exponent;
    double mantissa = std::frexp(x, &exponent);

    // Convert mantissa to an exact integer representation
    const int mantissa_bits = 53;
    double scaled_mantissa = std::ldexp(mantissa, mantissa_bits); // mantissa * 2^53
    num = static_cast<long>(scaled_mantissa);
    exponent -= mantissa_bits;

    den = 1;
    if (exponent > 0) num <<= exponent;
    if (exponent < 0) den <<= -exponent;
    if (x < 0.0) num = -num;

    simplify();
}

constexpr void rational::simplify() {
    if (den.is_zero())
        throw std::runtime_error("rational with zero denominator");
    if (den.is_one())
        return;
    if (num.is_zero()) {
        den = 1;
        return;
    }
    if (den.is_negative()) {
        den.negate();
        num.negate();
    }

    // TODO optimize for small integers
    auto az = num.num_trailing_zeros();
    auto bz = den.num_trailing_zeros();

    auto z = std::min(az, bz);
    num >>= z;
    den >>= z;

    if (num.is_one() || den.is_one())
        return;

    integer a = abs(num);
    integer b = den;
    a >>= az - z;
    while (!b.is_zero()) {
        b >>= b.num_trailing_zeros();
        if (a > b)
            std::swap(a, b);
        b -= a;
    }
    if (!a.is_one()) {
        num /= a;
        den /= a;
    }
}
constexpr rational::operator float() const {
    const auto num_bits = den.num_bits();
    if (num_bits <= 50) return static_cast<float>(num) / static_cast<float>(den);

    auto e = num_bits - 50;
    return static_cast<float>(num >> e) / static_cast<float>(den >> e);
}

constexpr rational::operator double() const {
    const auto num_bits = den.num_bits();
    if (num_bits <= 900) return static_cast<double>(num) / static_cast<double>(den);

    auto e = num_bits - 900;
    return static_cast<double>(num >> e) / static_cast<double>(den >> e);
}

constexpr rational& operator+=(rational& a, const rational& b) {
    if (a.den == b.den) {
        a.num += b.num;
    } else {
        a.num *= b.den;
        add_product(a.num, b.num, a.den);
        a.den *= b.den;
    }
    a.simplify();
    return a;
}

constexpr rational operator+(const rational& a, const rational& b) {
    if (a.den == b.den)
        return rational{a.num + b.num, a.den};

    integer p = a.num * b.den;
    add_product(p, b.num, a.den);
    return rational{std::move(p), a.den * b.den};
}

constexpr rational operator+(const rational& a, const integer& b) {
    return {a.num + b * a.den, a.den};
}

constexpr rational operator+(const integer& a, const rational& b) {
    return {a * b.den + b.num, b.den};
}

constexpr rational& operator-=(rational& a, const rational& b) {
    if (a.den == b.den) {
        a.num -= b.num;
    } else {
        a.num *= b.den;
        sub_product(a.num, b.num, a.den);
        a.den *= b.den;
    }
    a.simplify();
    return a;
}

constexpr rational operator-(const rational& a, const rational& b) {
    if (a.den == b.den)
        return rational{a.num - b.num, a.den};

    integer p = a.num * b.den;
    sub_product(p, b.num, a.den);
    return rational{std::move(p), a.den * b.den};
}

constexpr rational operator-(const rational& a, const integer& b) {
    return {a.num - b * a.den, a.den};
}

constexpr rational operator-(const integer& a, const rational& b) {
    return {a * b.den - b.num, b.den};
}

constexpr rational& operator*=(rational& a, const rational& b) {
    a.num *= b.num;
    a.den *= b.den;
    if (&a != &b)
        a.simplify();
    return a;
}

constexpr rational& operator*=(rational& a, const integer& b) {
    a.num *= b;
    a.simplify();
    return a;
}

constexpr rational operator*(const rational& a, const rational& b) {
    return {a.num * b.num, a.den * b.den};
}

constexpr rational operator*(const rational& a, const integer& b) {
    return {a.num * b, a.den};
}

constexpr rational operator*(const integer& a, const rational& b) {
    return {a * b.num, b.den};
}

constexpr rational& operator/=(rational& a, const integer& b) {
    a.den *= b;
    a.simplify();
    return a;
}

constexpr rational operator/(const rational& a, const rational& b) {
    return {a.num * b.den, a.den * b.num};
}

constexpr rational operator/(const rational& a, const integer& b) {
    return {a.num, a.den * b};
}

constexpr rational operator/(const integer& a, const rational& b) {
    return {a * b.den, b.num};
}

static_assert(sizeof(rational) == 32);

constexpr rational operator%(const rational& a, const rational& b) {
    return a - (a.num * b.den) / (a.den * b.num) * b;
}

constexpr rational operator%(const rational& a, const integer& b) {
    return a - a.num / (a.den * b) * b;
}

constexpr rational& operator%=(rational& a, const rational& b) {
    a -= (a.num * b.den) / (a.den * b.num) * b;
    return a;
}

constexpr rational& operator%=(rational& a, const integer& b) {
    a -= a.num / (a.den * b) * b;
    return a;
}

}

template <>
struct std::formatter<algebra::rational, char> {
    std::optional<int> frac_digits;

    constexpr auto parse(auto& ctx) {
        auto it = ctx.begin();
        auto end = ctx.end();

        if (it != end && *it == '.') {
            ++it;

            if (it != end && '0' <= *it && *it <= '9') {
                int digits = 0;
                while (it != end && '0' <= *it && *it <= '9')
                    digits = digits * 10 + (*it++ - '0');
                frac_digits = digits;
            } else {
                throw std::format_error("Expected digits after '.' in format specifier.");
            }
        }

        if (it == end || *it != '}')
            throw std::format_error("Invalid format specifier for rational.");
        return it;
    }

    constexpr auto format(const algebra::rational& a, auto& ctx) const {
        using namespace algebra;
        auto it = ctx.out();
        if (frac_digits == nullopt) {
            int capacity = a.num.str_size_upper_bound();
            if (!a.den.is_one())
                capacity = std::max(capacity, a.den.str_size_upper_bound());

            std::string str;
            str.resize(capacity);

            int n = a.num.str(str.data(), capacity);
            for (int i = 0; i < n; i++)
                *it++ = str[i];

            if (!a.den.is_one()) {
                *it++ = '/';
                n = a.den.str(str.data(), capacity);
                for (int i = 0; i < n; i++)
                    *it++ = str[i];
            }
            return it;
        }

        if (a.num < 0)
            *it++ = '-';
        const int r = *frac_digits;
        integer n = abs(a.num) * pow(integer(10), r);
        integer w; // w = n/a.den + 1/2
        if (a.den.is_even())
            w = (n + (a.den >> 1)) / a.den;
        else
            w = ((n << 1) + a.den) / (a.den << 1);
        // TODO allocate on stack if small enough
        std::string s = w.str();
        for (int i = r; i < s.size(); i++)
            *it++ = s[i - r];
        if (r >= s.size())
            *it++ = '0';
        *it++ = '.';
        for (int i = 0; i < r; i++)
            *it++ = s[i + s.size() - r];
        return it;
    }
};

constexpr std::ostream& operator<<(std::ostream& os, const algebra::rational& a) { return os << a.str(); }

template<>
struct std::hash<algebra::rational> {
    constexpr size_t operator()(const algebra::rational& a) const {
        uint64_t seed = 0;
        seed = algebra::integer_backend::hash_fn_64bit(seed ^ std::hash<algebra::integer>()(a.num));
        seed = algebra::integer_backend::hash_fn_64bit(seed ^ std::hash<algebra::integer>()(a.den));
        return seed;
    }
};
