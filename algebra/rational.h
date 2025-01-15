#pragma once
#include "algebra/integer.h"
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
concept Integral = std::integral<T> || std::same_as<T, integer>;

constexpr bool operator<(Integral auto a, const rational& b) { return b.den.is_one() ? (a <  b.num) : (a * b.den <  b.num); }
constexpr bool operator<(const rational& a, Integral auto b) { return a.den.is_one() ? (a.num <  b) : (a.num <  b * a.den); }

constexpr bool operator==(const rational& a, const rational& b) { return a.num == b.num && a.den == b.den; }
constexpr bool operator==(const rational& a, const Integral auto b) { return a.num == b && a.den.is_one(); }
constexpr bool operator==(const Integral auto a, const rational& b) { return a == b.num && b.den.is_one(); }

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

namespace algebra {

constexpr auto operator""_q(const char* s) { return rational(s); }

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

constexpr int approx_log2(const rational& a) {
    return a.num.num_bits() - a.den.num_bits();
}

constexpr rational sqrt(const integer& x, unsigned iterations) {
    rational s = x;
    s.num += 1;
    s /= 2;
    integer temp;
    while (iterations--) {
        mul(s.den, s.den, /*out*/temp);
        mul(temp, x);
        add_product(temp, s.num, s.num);

        s.den *= s.num;
        s.den <<= 1;
        s.num.swap(temp);
        s.simplify();
    }
    return s;
}

constexpr rational sqrt(const rational& x, unsigned iterations) {
    if (x.is_integer())
        return sqrt(x.num, iterations);
    rational s = x;
    s.num += s.den;
    s /= 2;
    while (iterations--) {
        s += x / s;
        s /= 2;
    }
    return s;
}

constexpr rational pow(const rational& base, long exp) {
    if (exp == 0)
        return rational{1};
    if (exp == 1)
        return base;
    if (exp == 2)
        return base * base;
    if (exp == -1)
        return rational{base.den, base.num};
    if (exp == -2)
        return rational{base.den * base.den, base.num * base.num};

    const bool invert = exp < 0;
    rational result = 1;
    rational _base = base;
    while (exp) {
        if (exp % 2)
            result *= _base;
        _base *= _base;
        exp >>= 1;
    }
    if (invert)
        result.invert();
    return result;
}

constexpr void pow(const rational& base, const integer& exp, rational& out) {
    if (exp.is_zero()) {
        out = 1;
        return;
    }
    if (exp == 1) {
        out = base;
        return;
    }
    if (exp == 2) {
        out.num = base.num;
        out.num *= base.num;
        out.den = base.den;
        out.den *= base.den;
        return;
    }
    if (exp == -1) {
        out.num = base.den;
        out.den = base.num;
        return;
    }
    if (exp == -2) {
        out.num = base.den;
        out.num *= base.den;
        out.den = base.num;
        out.den *= base.num;
        return;
    }

    out = 1;
    rational _base = base;
    if (exp.abs.is_odd())
        out *= _base;
    for (size_t i = 1; i < exp.abs.num_bits(); i++) {
        _base *= _base;
        if (exp.abs.bit(i))
            out *= _base;
    }
    if (exp.sign() < 0)
        out.invert();
}

constexpr rational pow(const rational& base, const integer& exp) {
    rational out;
    pow(base, exp, out);
    return out;
}

constexpr rational nth_root(const rational& base, const integer& exp, unsigned iterations) {
    if (exp == 2)
        return sqrt(base, iterations);
    if (exp == -2)
        return 1 / sqrt(base, iterations);

    rational result = 1;
    rational q;
    integer e = 1 - exp;
    while (iterations--) {
        pow(result, e, q);
        q *= base;
        q -= result;
        q /= exp;
        result += q;
    }
    return result;
}

constexpr rational pow(const rational& base, const rational& exp, unsigned iterations) {
    if (exp.is_integer())
        return pow(base, exp.num);

    integer quot, rem;
    div(exp.num, exp.den, quot, rem);
    rational root = nth_root(base, exp.den, iterations);
    return pow(base, quot) * pow(root, rem);
}

constexpr rational fract(const rational& a) {
    return {abs(a.num) % a.den, a.den};
}

constexpr rational abs(rational a) {
    if (a.num.is_negative())
        a.num.negate();
    return a;
}

constexpr rational round(const rational& a, unsigned digits, unsigned base = 10) {
    const natural b = pow(natural(base), digits);
    return {(a.num * b) / a.den, b};
}

static_assert(sizeof(rational) == 32);

// round towards 0 to integer
constexpr integer trunc(const rational& a) { return a.is_integer() ? a.num : (a.num / a.den); }

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

struct BinarySplit {
    integer p, q, r;
};

constexpr BinarySplit __PI(unsigned a, unsigned b) {
    if (b == a + 1) {
        integer e = a;
        integer p = -(6*e - 5)*(2*e - 1)*(6*e - 1);
        integer q = 10939058860032000 * (e * e * e);
        integer r = p * (545140134*a + 13591409);
        return {p, q, r};
    }
    unsigned m = (a + b) / 2;
    BinarySplit low = __PI(a, m);
    BinarySplit high = __PI(m, b);
    return {low.p * high.p, low.q * high.q, high.q * low.r + low.p * high.r};
}

// Chudnovsky algorithm
constexpr rational PI(unsigned n) {
    BinarySplit e = __PI(1, std::max(2u, n));
    return sqrt(10005_q, n) * rational{426880 * e.q, 13591409 * e.q + e.r};
}

constexpr rational sin(rational x, unsigned n) {
    const bool negate = x.num.sign() < 0;
    x.num.abs.words.set_negative(false);
    x %= 2 * PI(10);

    rational out = x;
    rational a = x;
    rational b = x * x;
    b.negate();
    unsigned i = 1;
    while (n--) {
        a *= b;
        a /= ++i;
        a /= ++i;
        out += a;
    }
    if (negate)
        out.negate();
    return out;
}

constexpr rational cos(rational x, unsigned n) {
    x.num.abs.words.set_negative(false);
    x %= 2 * PI(10);

    rational out = 1;
    rational a = x;
    rational b = x * x;
    b.negate();
    unsigned i = 0;
    while (n--) {
        a *= b;
        a /= ++i;
        a /= ++i;
        out += a;
    }
    return out;
}

constexpr rational exp(rational x, unsigned n) {
    rational out = x;
    rational a = x;
    unsigned i = 0;
    while (++i < n) {
        a *= x;
        a /= i;
        out += a;
    }
    return out;
}

}
