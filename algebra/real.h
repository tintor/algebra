#pragma once
#include "algebra/rational.h"

namespace algebra {

template<int Base = 2>
struct real {
    integer num;
    int exp;

    constexpr real(std::integral auto a, int exp = 0) : num(a), exp(exp) { normalize(); }
    constexpr real(integer a, int exp = 0) : num(std::move(a)), exp(exp) { normalize(); }

    constexpr real(float a) : real(rational(a)) { }
    constexpr real(double a) : real(rational(a)) { }
    constexpr real(const rational& a); // exact conversion
    static constexpr real round(const rational& a, unsigned digits);

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
        num >= z;
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
        integer a = s.den;
        for (auto [factor, count] : base_factors)
            while (a > 1 && a % factor == 0)
                a /= static_cast<long>(factor);
        if (a > 1)
            throw std::runtime_error("inexact conversion");

        throw std::runtime_error("unimplemented");
    }
    normalize();
}

template<int Base>
constexpr real<Base> real<Base>::round(const rational& a, unsigned digits) {
    const natural b = pow(natural(Base), digits);
    return real<Base>((a.num * b) / a.den, -digits);
}

template<int B>
constexpr real<B> operator-(const real<B>& a) { return {-a.num, a.exp}; }

template<int B>
constexpr real<B> operator+(const real<B>& a, const real<B>& b) {
    if (a.exp > b.exp) return {a.num * pow(integer(B), a.exp - b.exp) + b.num, b.exp};
    if (b.exp > a.exp) return {a.num + b.num * pow(integer(B), b.exp - a.exp), a.exp};
    return {a.num + b.num, a.exp};
}

template<int B>
constexpr real<B> operator-(const real<B>& a, const real<B>& b) {
    if (a.exp > b.exp) return {a.num * pow(integer(B), a.exp - b.exp) + b.num, b.exp};
    if (b.exp > a.exp) return {a.num + b.num * pow(integer(B), b.exp - a.exp), a.exp};
    return {a.num + b.num, a.exp};
}

template<int B>
constexpr real<B> operator*(const real<B>& a, const real<B>& b) {
    return {a.num * b.num, a.exp + b.exp};
}

template<int B>
constexpr real<B> operator/(const real<B>& a, const real<B>& b) {
    if (b.num == 0)
        throw std::runtime_error("division by zero");
    int scale = 100;
    real<B> result(a.num * pow(integer(B), scale) / b.num, a.exp - b.exp - scale);
    result.normalize();
    return result;
}

template<int B>
constexpr bool operator==(const real<B>& a, const real<B>& b) {
    if constexpr (B == 2) {
        if (a.exp > b.exp) return a.num == (b.num << (a.exp - b.exp));
        if (a.exp < b.exp) return (a.num << (b.exp - a.exp)) == b.num;
    } else {
        if (a.exp > b.exp) return a.num == b.num * pow(integer(B), a.exp - b.exp);
        if (a.exp < b.exp) return a.num * pow(integer(B), b.exp - a.exp) == b.num;
    }
    return a.num == b.num;
}

template<int B>
constexpr bool operator<(const real<B>& a, const real<B>& b) {
    if constexpr (B == 2) {
        if (a.exp > b.exp) return a.num < (b.num << (a.exp - b.exp));
        if (a.exp < b.exp) return (a.num << (b.exp - a.exp)) < b.num;
    } else {
        if (a.exp > b.exp) return a.num < b.num * pow(integer(B), a.exp - b.exp);
        if (a.exp < b.exp) return a.num * pow(integer(B), b.exp - a.exp) < b.num;
    }
    return a.num < b.num;
}

template <int B>
constexpr rational to_rational(const real<B>& a) {
    if (a.exp > 0)
        return {a.num * pow(integer(B), a.exp)};
    if (a.exp < 0)
        return {a.num, pow(integer(B), -a.exp)};
    return {a.num};
}

namespace literals {
constexpr real<2> operator""_f(const char* s) { return real<2>::round(rational(s), 53); }
// TODO it would be more efficient to parse into decimal directly!
constexpr decimal operator""_d(const char* s) { return rational(s); }
}

}

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
        return std::formatter<algebra::rational, char>::format(to_rational(a), ctx);
    }
};

template <int B>
constexpr std::ostream& operator<<(std::ostream& os, const algebra::real<B>& a) { return os << a.str(); }
