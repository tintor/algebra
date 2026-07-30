#pragma once
#include "algebra/integer.h"
#include <random>

namespace algebra {

struct rational;
template<> struct IsNumberClass<rational> : std::true_type {};

template<typename T> concept integral = std::same_as<T, integer> || std_int<T>;
template<typename T> concept rational_like = integral<T> || std::same_as<T, rational>;

struct rational {
    integer num, den;

    constexpr rational() : den(1) { }
    constexpr rational(integer a) : num(std::move(a)), den(1) { }
    constexpr rational(integer a, integer b) : num(std::move(a)), den(std::move(b)) { simplify(); }
private:
    constexpr rational(integer a, integer b, int) : num(std::move(a)), den(std::move(b)) { }
public:
    constexpr static rational normalized(integer num, integer den) { return {std::move(num), std::move(den), /*dummy*/0}; }

    constexpr rational(std_int auto a) : num(a), den(1) { }

    constexpr rational(std_int auto n, std_int auto d) {
        Check(d != 0, "rational with zero denominator");
        if (n == 0) {
            num = uint64_t(0);
            den = uint64_t(1);
            return;
        }

        using T = make_unsigned_t<larger_type<decltype(n), decltype(d)>>;
        T un = abs_unsigned(n);
        T ud = abs_unsigned(d);

        auto z = countr_zero(un | ud);
        un >>= z;
        ud >>= z;
        auto e = gcd(un, ud);
        num = un / e;
        den = ud / e;
        if ((n < 0) != (d < 0))
            num.negate();
    }

    template<std::floating_point T>
    constexpr rational(T x);

    // grammar: [-] digits [ ('.' | '/') digits ] [ 'e' [+-] digits ]
    constexpr rational(std::string_view s) {
        auto is_digit = [](char c) { return '0' <= c && c <= '9'; };
        size_t i = 0;
        const bool negative = i < s.size() && s[i] == '-';
        if (negative)
            i += 1;

        const size_t int_begin = i;
        while (i < s.size() && is_digit(s[i]))
            i += 1;
        if (i == int_begin)
            throw std::runtime_error("syntax error parsing rational");
        num = integer(s.substr(0, i));

        den = 1;
        if (i < s.size() && (s[i] == '.' || s[i] == '/')) {
            const char kind = s[i++];
            const size_t frac_begin = i;
            while (i < s.size() && is_digit(s[i]))
                i += 1;
            if (i == frac_begin)
                throw std::runtime_error("syntax error parsing rational");
            const std::string_view digits = s.substr(frac_begin, i - frac_begin);
            if (kind == '/') {
                den = integer(digits);
            } else {
                den = pow(integer(10), digits.size());
                num *= den;
                const integer frac(digits);
                // integer("-0") is +0, so the sign has to come from the text
                if (negative) num -= frac; else num += frac;
            }
        }

        if (i < s.size() && s[i] == 'e') {
            i += 1;
            bool negative_exponent = false;
            if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
                negative_exponent = s[i] == '-';
                i += 1;
            }
            const size_t exp_begin = i;
            uint64_t e = 0;
            while (i < s.size() && is_digit(s[i])) {
                e = e * 10 + (s[i] - '0');
                Check(e <= 100'000'000, "exponent too large parsing rational");
                i += 1;
            }
            if (i == exp_begin)
                throw std::runtime_error("syntax error parsing rational");
            const integer p = pow(integer(10), e);
            if (negative_exponent) den *= p; else num *= p;
        }

        if (i != s.size())
            throw std::runtime_error("syntax error parsing rational");

        simplify();
    }

    constexpr rational(const std::string& s) : rational(std::string_view(s)) { }
    constexpr rational(const char* s) : rational(std::string_view(s)) { }

    constexpr void simplify();

    constexpr operator float() const;
    constexpr operator double() const;

    constexpr void invert() {
        Check(!num.is_zero(), "rational with zero denominator");
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
    constexpr bool is_negative() const { return num.is_negative(); }
    constexpr bool is_zero() const { return num.is_zero(); }
};

constexpr auto signum(const rational& a) { return signum(a.num); }
constexpr void negate(rational& a) { a.negate(); }
constexpr rational fract(const rational&);

constexpr rational operator-(const rational& a) { return {-a.num, a.den}; }

constexpr rational& operator+=(rational& a, const rational& b);
constexpr rational operator+(const rational& a, const rational& b);
// Adding a whole number keeps the denominator, and gcd(num + b*den, den) == gcd(num, den) == 1, so
// the result is already in lowest terms: normalized() skips the gcd the constructor would compute.
constexpr rational operator+(const rational& a, const integral auto& b) {
    return rational::normalized(a.num + b * a.den, a.den);
}
constexpr rational operator+(const integral auto& a, const rational& b) { return b + a; }

constexpr rational& operator-=(rational& a, const rational& b);
constexpr rational operator-(const rational& a, const rational& b);
constexpr rational operator-(const rational& a, const integral auto& b);
constexpr rational operator-(const integral auto& a, const rational& b);

constexpr rational operator*(const rational& a, const rational& b);
constexpr rational operator*(const rational& a, const integral auto& b);
constexpr rational operator*(const integral auto& a, const rational& b);

constexpr rational& operator*=(rational& a, const rational& b);
constexpr rational& operator*=(rational& a, const integral auto& b);

constexpr rational operator/(const rational&, const rational&);
constexpr rational operator/(const rational&, const integral auto&);
constexpr rational operator/(const integral auto&, const rational&);

constexpr bool operator<(const rational& a, const rational& b) { return (a.den == b.den) ? (a.num < b.num) : less_ab_cd(a.num, b.den, b.num, a.den); }
constexpr bool operator<(integral auto a, const rational& b) { return b.den.is_one() ? (a < b.num) : less_ab_c(a, b.den, b.num); }
constexpr bool operator<(const rational& a, integral auto b) { return a.den.is_one() ? (a.num < b) : less_a_bc(a.num, b, a.den); }

constexpr bool operator==(const rational& a, const rational& b) { return a.num == b.num && a.den == b.den; }
constexpr bool operator==(const rational& a, const integral auto b) { return a.num == b && a.den.is_one(); }

namespace literals {
constexpr auto operator""_q(const char* s) { return rational(s); }
}

template<std::floating_point T>
constexpr rational::rational(T x) {
    if (std::isnan(x))
        throw std::runtime_error("can't convert nan to rational");
    if (std::isinf(x))
        throw std::runtime_error("can't convert infinite to rational");
    if (x == 0) {
        num = 0;
        den = 1;
        return;
    }

    int exponent;
    T mantissa = std::frexp(x, &exponent);

    // Convert mantissa to an exact integer representation
    const int mantissa_bits = std::numeric_limits<T>::digits;
    T scaled_mantissa = std::ldexp(mantissa, mantissa_bits);
    num = static_cast<long>(scaled_mantissa); // already signed, same sign as x
    exponent -= mantissa_bits;

    den = 1;
    if (exponent > 0)
        num <<= exponent;
    if (exponent < 0) {
        den <<= -exponent;
        simplify();
    }
}

constexpr void rational::simplify() {
    Check(!den.is_zero(), "rational with zero denominator");
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

    auto az = num.num_trailing_zeros();
    auto bz = den.num_trailing_zeros();

    auto z = std::min(az, bz);
    num >>= z;
    den >>= z;

    if (den.words.size() <= 1) {
        if (is_power_of_two(den.words[0]))
            return;
        if (num.words.size() <= 1) {
            // note that num can be negative here, but we only need its absolute value
            if (is_power_of_two(num.words[0]))
                return;
            // both fit in a word, and the numerator is odd after the shift, which is what
            // __gcd_inner() wants. Any factor of two left in the denominator cannot be common.
            const uint64_t g = __gcd_inner<uint64_t>(num.words[0] >> (az - z), den.words[0]);
            if (g != 1) {
                num /= g;
                den /= g;
            }
            return;
        }
    }

    // TODO allocate the operands on stack if they are small enough
    const integer g = gcd(abs(num) >> (az - z), abs(den));
    if (!g.is_one()) {
        num /= g;
        den /= g;
    }
}

// Converts num/den to T. Each term is scaled down to at most `keep` significant bits (which is
// still far more than T's precision) so that neither of them overflows T on the way, and the
// discarded powers of two are put back with ldexp. Both terms have to be scaled independently:
// shifting them by the same amount loses the whole numerator when den is much larger than num.
template<std::floating_point T>
constexpr T __rational_to_float(const integer& num, const integer& den, int64_t keep) {
    const int64_t nb = num.num_bits();
    const int64_t db = den.num_bits();
    if (nb <= keep && db <= keep)
        return static_cast<T>(num) / static_cast<T>(den);

    const int64_t ne = (nb > keep) ? (nb - keep) : 0;
    const int64_t de = (db > keep) ? (db - keep) : 0;
    const T q = static_cast<T>(num >> ne) / static_cast<T>(den >> de);

    int64_t e = ne - de;
    // the result over/underflows T long before this, and ldexp() takes an int
    if (e > 1'000'000) e = 1'000'000;
    if (e < -1'000'000) e = -1'000'000;
    return std::ldexp(q, static_cast<int>(e));
}

constexpr rational::operator float() const { return __rational_to_float<float>(num, den, 50); }

constexpr rational::operator double() const { return __rational_to_float<double>(num, den, 900); }

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

// as with operator+ above, subtracting a whole number cannot introduce a common factor
constexpr rational operator-(const rational& a, const integral auto& b) {
    return rational::normalized(a.num - b * a.den, a.den);
}
constexpr rational operator-(const integral auto& a, const rational& b) {
    return rational::normalized(a * b.den - b.num, b.den);
}

constexpr rational& operator*=(rational& a, const rational& b) {
    a.num *= b.num;
    a.den *= b.den;
    if (&a != &b)
        a.simplify();
    return a;
}

constexpr rational& operator*=(rational& a, const integral auto& b) {
    a.num *= b;
    a.simplify();
    return a;
}

constexpr rational operator*(const rational& a, const rational& b) { return {a.num * b.num, a.den * b.den}; }
constexpr rational operator*(const rational& a, const integral auto& b) { return {a.num * b, a.den}; }
constexpr rational operator*(const integral auto& a, const rational& b) { return {a * b.num, b.den}; }

constexpr rational& operator/=(rational& a, const rational& b) {
    if (&a == &b) {
        // a.den *= b.num below would read the already updated numerator
        Check(!a.num.is_zero(), "rational with zero denominator");
        a.num = 1;
        a.den = 1;
        return a;
    }
    a.num *= b.den;
    a.den *= b.num;
    a.simplify();
    return a;
}

constexpr rational& operator/=(rational& a, const integral auto& b) {
    a.den *= b;
    a.simplify();
    return a;
}

constexpr rational operator/(const rational& a, const rational& b) { return {a.num * b.den, a.den * b.num}; }
constexpr rational operator/(const rational& a, const integral auto& b) { return {a.num, a.den * b}; }
constexpr rational operator/(const integral auto& a, const rational& b) { return {b.den * a, b.num}; }

static_assert(sizeof(rational) == 32);

constexpr rational operator%(const rational& a, const rational& b) { return a - (a.num * b.den) / (a.den * b.num) * b; }
constexpr rational operator%(const rational& a, const integral auto& b) { return a - a.num / (a.den * b) * b; }

constexpr rational& operator%=(rational& a, const rational& b) { a -= (a.num * b.den) / (a.den * b.num) * b; return a; }
constexpr rational& operator%=(rational& a, const integral auto& b) { a -= a.num / (a.den * b) * b; return a; }

constexpr rational& operator<<=(rational& a, int64_t b) {
    if (b > 0) {
        auto z = a.den.num_trailing_zeros();
        if (b > z) {
            a.num <<= b - z;
            a.den >>= z;
        } else {
            a.den >>= b;
        }
        return a;
    }
    if (b < 0) {
        b = -b;
        auto z = a.num.num_trailing_zeros();
        if (b > z) {
            a.num >>= z;
            a.den <<= b - z;
        } else {
            a.num >>= b;
        }
        if (a.num == 0)
            a.den = 1;
    }
    return a;
}

ALGEBRA_SHIFT_OP(rational)

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

            if (it != end && *it == 'f')
                ++it; // fixed notation, same as without it
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
        const std::string s = w.str();
        const int len = s.size();
        if (len > r) {
            for (int i = 0; i < len - r; i++)
                *it++ = s[i];
        } else {
            *it++ = '0';
        }
        if (r == 0)
            return it; // no fraction digits requested, and no trailing '.' either
        *it++ = '.';
        for (int i = len; i < r; i++) // fraction shorter than r digits
            *it++ = '0';
        for (int i = (len > r) ? (len - r) : 0; i < len; i++)
            *it++ = s[i];
        return it;
    }
};

inline std::ostream& operator<<(std::ostream& os, const algebra::rational& a) { return os << a.str(); }

template<>
struct std::hash<algebra::rational> {
    constexpr size_t operator()(const algebra::rational& a) const {
        uint64_t seed = 0;
        seed = algebra::hash_fn_64bit(seed ^ std::hash<algebra::integer>()(a.num));
        seed = algebra::hash_fn_64bit(seed ^ std::hash<algebra::integer>()(a.den));
        return seed;
    }
};
