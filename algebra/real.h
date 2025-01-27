#pragma once
#include "algebra/rational.h"
#include "algebra/natural_func.h"
#include "algebra/integer_func.h"

namespace algebra {

template<int B> struct real;
template<int B> struct IsNumberClass<real<B>> : std::true_type {};

// assuming exp >= 0
template<int B>
constexpr integer shift(const integer& a, std::integral auto exp) {
    if (exp == 0)
        return a;
    if constexpr (B == 2)
        return a << exp;
    return pow(integer(B), exp, a);
}

template<int Base = 2>
struct real {
    integer num;
    int exp;

    constexpr real(std::integral auto a, int exp = 0) : num(a), exp(exp) { normalize(); }
    constexpr real(integer a, int exp = 0) : num(std::move(a)), exp(exp) { normalize(); }
    constexpr real(natural a, int exp = 0) : num(std::move(a)), exp(exp) { normalize(); }
    constexpr real(integer a, int exp, int dummy) : num(std::move(a)), exp(exp) { }

    constexpr real(float a) : real(rational(a)) { }
    constexpr real(double a) : real(rational(a)) { }
    constexpr real(const rational& a); // exact conversion
    static constexpr real round(const rational& a, int digits);

    constexpr void normalize() {
        if constexpr (Base == 2) {
            auto e = num.num_trailing_zeros();
            if (e) {
                num >>= e;
                exp += e;
            }
        } else {
            integer quot;
            while (num) {
                if (div(num, Base, /*out*/quot))
                    break;
                num.swap(quot);
                exp += 1;
            }
        }
        if (num.is_zero())
            exp = 0;
    }

    std::string str() const;
};
using decimal = real<10>;

template<int Base>
constexpr real<Base>::real(const rational& s) {
    // s.den must be a non-negative power of B (otherwise exception will be thrown!
    num = s.num;
    if constexpr (Base == 2) {
        if (!is_power_of_two(s.den.abs))
            throw std::runtime_error("inexact conversion");
        exp = -s.den.num_trailing_zeros();
        auto z = num.num_trailing_zeros();
        exp += z;
    } else if constexpr (Base == 10) {
        // Check if s.den has any prime factors not in base
        integer a = s.den;
        int twos = a.num_trailing_zeros();
        a >>= twos;
        int fives = 0;
        while (a > 1 && a.abs.mod5() == 0) {
            a /= 5;
            fives += 1;
        }
        if (a > 1)
            throw std::runtime_error("inexact conversion");

        if (twos > fives)
            num *= pow(natural(5), twos - fives);
        if (fives > twos)
            num <<= fives - twos;
        exp = -std::max(twos, fives);
    } else {
        const auto base_factors = factorize(Base);
        // Check if s.den has any prime factors not in base
        natural a = s.den.abs, q;
        for (auto [factor, count] : base_factors)
            while (a > 1) {
                if (div(a, factor, q))
                    break;
                std::swap(a, q);
            }
        if (a > 1)
            throw std::runtime_error("inexact conversion");

        throw std::runtime_error("unimplemented");
    }
    normalize();
}

template<int B>
constexpr real<B> real<B>::round(const rational& a, int digits) {
    return {shift<B>(a.num, digits) / a.den, -digits};
}

template<int B>
constexpr real<B> operator-(const real<B>& a) { return {-a.num, a.exp}; }

template<int B>
constexpr real<B> operator+(const real<B>& a, const real<B>& b) {
    if (a.exp > b.exp) return {shift<B>(a.num, a.exp - b.exp) + b.num, b.exp};
    if (b.exp > a.exp) return {a.num + shift<B>(b.num, b.exp - a.exp), a.exp};
    return {a.num + b.num, a.exp};
}

template<int B>
constexpr real<B> operator+(const real<B>& a, const integral auto& b) {
    if (a.exp > 0) return {shift<B>(a.num, a.exp) + b};
    if (0 > a.exp) return {a.num + shift<B>(b, -a.exp), a.exp};
    return {a.num + b};
}

template<int B> constexpr real<B> operator+(const integral auto& a, const real<B>& b) { return b + a; }

template<int B>
constexpr real<B> operator-(const real<B>& a, const real<B>& b) {
    if (a.exp > b.exp) return {shift<B>(a.num, a.exp - b.exp) - b.num, b.exp};
    if (b.exp > a.exp) return {a.num - shift<B>(b.num, b.exp - a.exp), a.exp};
    return {a.num - b.num, a.exp};
}

template<int B>
constexpr real<B> operator-(const real<B>& a, const integral auto& b) {
    if (a.exp > 0) return {shift<B>(a.num, a.exp) - b};
    if (0 > a.exp) return {a.num - shift<B>(b, -a.exp), a.exp};
    return {a.num - b};
}

template<int B>
constexpr real<B> operator-(const integral auto& a, const real<B>& b) {
   auto c = b - a;
   c.num.negate();
   return c;
};

template<int B> constexpr real<B> operator*(const real<B>& a, const real<B>& b) { return {a.num * b.num, a.exp + b.exp}; }
template<int B> constexpr real<B> operator*(const real<B>& a, const integral auto& b) { return {a.num * b, a.exp}; }
template<int B> constexpr real<B> operator*(const integral auto& a, const real<B>& b) { return {a * b.num, b.exp}; }

template<int B, typename T> constexpr real<B>& operator+=(real<B>& a, const T& b) { a = a + b; return a; }
template<int B, typename T> constexpr real<B>& operator-=(real<B>& a, const T& b) { a = a - b; return a; }
template<int B, typename T> constexpr real<B>& operator*=(real<B>& a, const T& b) { a = a * b; return a; }

template<int B>
constexpr real<B> operator/(const real<B>& a, const real<B>& b) {
    int scale = 100;
    return {shift<B>(a.num, scale) / b.num, a.exp - b.exp - scale};
}

template<int B>
constexpr real<B>& operator/=(real<B>& a, const integral auto& b) {
    int scale = 100;
    if constexpr (B == 2) {
        a.num <<= scale;
    } else {
        a.num = pow(integer(B), scale, a.num);
    }
    a.num /= b;
    a.exp -= scale;
    a.normalize();
    return a;
}

template<int B> constexpr real<B> operator/(real<B> a, const integral auto& b) { a /= b; return a; }
template<int B> constexpr real<B> operator/(const integral auto& a, real<B> b) { return real<B>(a) / b; }

template<int B>
constexpr real<B>& operator<<=(real<B>& a, int64_t b) {
    if constexpr (B == 2) {
        a.exp += b;
        return a;
    }
    if (b > 0) {
        a.num <<= b;
        if constexpr (B % 2 == 0)
            a.normalize();
    }
    if (b < 0)
        a /= exp2(-b);
    return a;
}

ALGEBRA_SHIFT_OP(real<2>)
ALGEBRA_SHIFT_OP(real<10>)

template<int B>
constexpr bool operator==(const real<B>& a, const real<B>& b) {
    if (a.exp > b.exp) return shift<B>(a.num, a.exp - b.exp) == b.num;
    if (a.exp < b.exp) return a.num == shift<B>(b.num, b.exp - a.exp);
    return a.num == b.num;
}

template<int B>
constexpr bool operator==(const real<B>& a, const integral auto& b) {
    if (a.exp > 0) return shift<B>(a.num, a.exp) == b;
    if (a.exp < 0) return a.num == shift<B>(b, -a.exp);
    return a.num == b;
}

template<int B>
constexpr bool operator<(const real<B>& a, const real<B>& b) {
    if (a.exp > b.exp) return shift<B>(a.num, a.exp - b.exp) < b.num;
    if (a.exp < b.exp) return a.num < shift<B>(b.num, b.exp - a.exp);
    return a.num < b.num;
}

template<int B>
constexpr bool operator<(const real<B>& a, const integral auto& b) {
    if (a.exp > 0) return shift<B>(a.num, a.exp) < b;
    if (a.exp < 0) return a.num < shift<B>(b, -a.exp);
    return a.num < b;
}

template<int B>
constexpr bool operator<(const integral auto& a, const real<B>& b) {
    if (0 > b.exp) return shift<B>(a, -b.exp) < b.num;
    if (0 < b.exp) return a < shift<B>(b.num, b.exp);
    return a < b.num;
}

template <int B>
constexpr rational to_rational(const real<B>& a) {
    if (a.exp > 0)
        return {shift<B>(a.num, a.exp)};
    if (a.exp < 0)
        return {a.num, shift<B>(integer(1), -a.exp)};
    return {a.num};
}

namespace literals {
constexpr real<2> operator""_f(const char* s) { return real<2>::round(rational(s), 53); }
// TODO it would be more efficient to parse into decimal directly!
constexpr decimal operator""_d(const char* s) { return rational(s); }
}

}

std::optional<int> REAL_FRACT_DIGITS;

template <int B>
struct std::formatter<algebra::real<B>, char> : std::formatter<algebra::rational, char> {
    constexpr auto format(const algebra::real<B>& a, auto& ctx) const {
        // TODO handle frac_digits specifier from formatter::parse()
        if constexpr (B == 10) {
            auto it = ctx.out();
            if (a.exp >= 0) {
                it = std::format_to(it, "{}", a.num);
                for (int i = 0; i < a.exp; i++)
                    *it++ = '0';
            } else {
                auto s = a.num.abs.str();
                if (a.num.sign() < 0)
                    *it++ = '-';
                if (s.size() <= -a.exp) {
                    *it++ = '0';
                    *it++ = '.';
                    for (int i = s.size(); i < -a.exp; i++)
                        *it++ = '0';
                    for (char c: s)
                        *it++ = c;
                } else {
                    for (int i = 0; i < s.size(); i++) {
                        if (i == s.size() + a.exp)
                            *it++ = '.';
                        *it++ = s[i];
                    }
                }
            }
            return it;
        }

        if (frac_digits != std::nullopt || REAL_FRACT_DIGITS == std::nullopt)
            return std::formatter<algebra::rational, char>::format(to_rational(a), ctx);

        std::formatter<algebra::rational, char> f;
        f.frac_digits = REAL_FRACT_DIGITS.value();
        f.format(to_rational(a), ctx);
        return ctx.out();
    }
};

template <int B>
constexpr std::ostream& operator<<(std::ostream& os, const algebra::real<B>& a) { return os << a.str(); }

template <int B>
struct std::hash<algebra::real<B>> {
    constexpr size_t operator()(const algebra::real<B>& a) const {
        uint64_t seed = 0;
        seed = algebra::integer_backend::hash_fn_64bit(seed ^ std::hash<algebra::integer>()(a.num));
        seed = algebra::integer_backend::hash_fn_64bit(seed ^ a.exp);
        return seed;
    }
};

namespace algebra {

template<int B>
std::string real<B>::str() const {
    return std::format("{}", *this);
}

}
