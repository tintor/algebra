#pragma once
#include "algebra/rational.h"

namespace algebra {

template<int Base = 2>
struct real {
    integer num;
    int exp;

    template<std::integral I>
    constexpr real(I a, int exp = 0) : num(a), exp(exp) { }
    constexpr real(integer a, int exp = 0) : num(std::move(a)), exp(exp) { }

    constexpr real(float a);
    constexpr real(double a);
    constexpr real(const rational& a, int prec = -53);

    constexpr void simplify() {
        if constexpr (Base == 2) {
            if (num.sign() > 0) {
                auto e = num.num_trailing_zeros();
                num >>= e;
                exp += e;
            }
        } else {
            while (num.sign() > 0 && num % Base == 0) {
                num /= Base;
                exp += 1;
            }
        }
        if (num.is_zero()) exp = 0;
    }
};
using decimal = real<10>;

template<int Base>
constexpr real<Base>::real(const rational& s, int prec) {
    num = s.num / s.den;
    integer e = 1;
    exp = 0;
    while (approx_log2(s - rational(num, e)) > prec) {
        exp -= 1;
        e *= Base;
        num = s.num * e / s.den;
    }
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
    result.simplify();
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
constexpr bool operator!=(const real<B>& a, const real<B>& b) {
    if constexpr (B == 2) {
        if (a.exp > b.exp) return a.num != (b.num << (a.exp - b.exp));
        if (a.exp < b.exp) return (a.num << (b.exp - a.exp)) != b.num;
    } else {
        if (a.exp > b.exp) return a.num < b.num != pow(integer(B), a.exp - b.exp);
        if (a.exp < b.exp) return a.num * pow(integer(B), b.exp - a.exp) != b.num;
    }
    return a.num != b.num;
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

template<int B>
constexpr bool operator<=(const real<B>& a, const real<B>& b) {
    if constexpr (B == 2) {
        if (a.exp > b.exp) return a.num <= (b.num << (a.exp - b.exp));
        if (a.exp < b.exp) return (a.num << (b.exp - a.exp)) <= b.num;
    } else {
        if (a.exp > b.exp) return a.num <= b.num * pow(integer(B), a.exp - b.exp);
        if (a.exp < b.exp) return a.num * pow(integer(B), b.exp - a.exp) <= b.num;
    }
    return a.num <= b.num;
}

template<int B>
constexpr bool operator>(const real<B>& a, const real<B>& b) {
    if constexpr (B == 2) {
        if (a.exp > b.exp) return a.num > (b.num << (a.exp - b.exp));
        if (a.exp < b.exp) return (a.num << (b.exp - a.exp)) > b.num;
    } else {
        if (a.exp > b.exp) return a.num > b.num * pow(integer(B), a.exp - b.exp);
        if (a.exp < b.exp) return a.num * pow(integer(B), b.exp - a.exp) > b.num;
    }
    return a.num > b.num;
}

template<int B>
constexpr bool operator>=(const real<B>& a, const real<B>& b) {
    if constexpr (B == 2) {
        if (a.exp > b.exp) return a.num >= (b.num << (a.exp - b.exp));
        if (a.exp < b.exp) return (a.num << (b.exp - a.exp)) >= b.num;
    } else {
        if (a.exp > b.exp) return a.num >= b.num * pow(integer(B), a.exp - b.exp);
        if (a.exp < b.exp) return a.num * pow(integer(B), b.exp - a.exp) >= b.num;
    }
    return a.num >= b.num;
}

template <int B>
constexpr rational to_rational(const real<B>& a) {
    if (a.exp > 0)
        return {a.num * pow(integer(B), a.exp)};
    if (a.exp < 0)
        return {a.num, pow(integer(B), -a.exp)};
    return {a.num};
}

}

template <int B>
struct std::formatter<algebra::real<B>, char> : std::formatter<algebra::rational, char> {
    constexpr auto format(const algebra::real<B>& a, auto& ctx) const {
        return std::formatter<algebra::rational, char>::format(to_rational(a), ctx);
    }
};

template <int B>
constexpr std::ostream& operator<<(std::ostream& os, const algebra::real<B>& a) { return os << a.str(); }
