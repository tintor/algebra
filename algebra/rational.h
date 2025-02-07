#pragma once
#include "algebra/rational_class.h"

namespace algebra {

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
    if (exp.is_odd())
        out *= _base;
    for (size_t i = 1; i < exp.num_bits(); i++) {
        _base *= _base;
        if (exp.bit(i))
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

// returns abs(a) > abs(b), minimizing memory allocation
constexpr bool abs_greater(const rational& a, const rational& b) {
    if (a.den == b.den)
        return abs_greater(a.num, b.num);
    if (a.den == 1u)
        return abs_greater(((b.den == 1u) ? a.num : (a.num * b.den)), b.num);
    return abs_greater(((b.den == 1u) ? a.num : (a.num * b.den)), b.num * a.den);
}

constexpr rational round(const rational& a, unsigned digits, unsigned base = 10) {
    const integer b = pow(integer(base), digits);
    return {(a.num * b) / a.den, b};
}

// round towards 0 to integer
constexpr integer trunc(const rational& a) { return a.is_integer() ? a.num : (a.num / a.den); }

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
    return sqrt(rational(10005), n) * rational{426880 * e.q, 13591409 * e.q + e.r};
}

constexpr rational sin(rational x, unsigned n) {
    const bool negate = x.num.is_negative();
    if (negate)
        x.num.negate();
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
    if (x.num.is_negative())
        x.num.negate();
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

constexpr void simplify(rational& x, rational& y) {
    simplify(x.num, y.num);
    simplify(x.den, y.den);
}

constexpr void simplify(rational& x, rational& y, rational& z) {
    simplify(x.num, y.num, z.den);
    simplify(x.den, y.den, z.den);
}

}
