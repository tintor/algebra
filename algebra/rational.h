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

// returns sqrt(x) with an absolute error smaller than 2**-bits
// (much cheaper than sqrt(x, iterations), whose operands double in size every iteration)
constexpr rational sqrt_bits(const rational& x, int bits) {
    Check(!x.num.is_negative(), "sqrt() of a negative number");
    Check(bits >= 0, "sqrt_bits() with negative precision");
    if (x.num.is_zero())
        return 0;
    // floor(sqrt(num * 2**(2*bits) / den)) / 2**bits
    natural n = abs(x.num);
    n.words.set_negative(false);
    n <<= 2 * bits;
    n /= abs(x.den);
    return {integer(isqrt(n)), integer(power_of_two(bits))};
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
    // exp >>= 1 on a negative value never reaches 0
    unsigned long e = invert ? -static_cast<unsigned long>(exp) : static_cast<unsigned long>(exp);
    rational result = 1;
    rational _base = base;
    while (e) {
        if (e & 1)
            result *= _base;
        e >>= 1;
        if (e)
            _base *= _base;
    }
    if (invert)
        result.invert();
    return result;
}

constexpr void pow(const rational& base, const integer& exp, rational& out) {
    if (&out == &base) { // the branches below overwrite out before they are done reading base
        const rational b = base;
        pow(b, exp, out);
        return;
    }
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
        if (out.den.is_negative()) { // rational keeps the sign in the numerator
            out.den.negate();
            out.num.negate();
        }
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
    Check(exp != 0, "nth_root() with a zero exponent");
    if (exp == 2)
        return sqrt(base, iterations);
    if (exp == -2)
        return 1 / sqrt(base, iterations);

    // Seed with the power of two closest to the root, from log2(root) == log2(base) / exp. Starting
    // at 1 is far from the root for anything but a base near 1, and Newton diverges from there for
    // a negative exp. Every iteration also multiplies the size of the operands by abs(exp), so the
    // iterations that a poor start needs get expensive quickly.
    // abs(shift) is at most abs(approx_log2(base)), so it fits
    const int64_t shift = static_cast<int64_t>(integer(approx_log2(base)) / exp);
    rational result = (shift >= 0) ? rational(exp2(shift)) : rational(integer(1), exp2(-shift));

    rational q;
    integer e = 1 - exp;
    while (iterations--) {
        pow(result, e, q);
        q *= base;
        q -= result;
        q /= exp;
        result += q;
        result.simplify(); // the operands grow by a factor of abs(exp) per iteration otherwise
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
    return abs_greater((b.den == 1u) ? a.num : (a.num * b.den),
                       (a.den == 1u) ? b.num : (b.num * a.den));
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
        integer r = p * (545140134ull * a + 13591409); // 545140134 * a overflows 32 bits for a >= 8
        return {p, q, r};
    }
    unsigned m = (a + b) / 2;
    BinarySplit low = __PI(a, m);
    BinarySplit high = __PI(m, b);
    return {low.p * high.p, low.q * high.q, high.q * low.r + low.p * high.r};
}

// Chudnovsky algorithm
// every term of the series is worth about 14.18 decimal digits (47.1 bits)
constexpr rational PI(unsigned n) {
    BinarySplit e = __PI(1, std::max(2u, n));
    const int bits = static_cast<int>(n) * 48 + 64;
    return sqrt_bits(rational(10005), bits) * rational{426880 * e.q, 13591409 * e.q + e.r};
}

// reduce x into [0, 2*PI). Note: expensive, since the denominator of PI(n) ends up in x,
// and unnecessary for arguments that are already small - 6 < 2*PI, so the cheap test first.
constexpr void __reduce_mod_two_pi(rational& x, unsigned n) {
    if (x < 6)
        return;
    const rational two_pi = 2 * PI(n);
    if (x < two_pi)
        return;
    x %= two_pi;
    // The remainder carries the denominator of PI(n), which would make every term of the
    // series that follows enormous. Keep the digits the series can actually use.
    x = round(x, 10 * n);
}

constexpr rational sin(rational x, unsigned n) {
    const bool negate = x.num.is_negative();
    if (negate)
        x.num.negate();
    __reduce_mod_two_pi(x, 10);

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
    __reduce_mod_two_pi(x, 10);

    rational out = 1;
    rational a = 1; // cos(x) = 1 - x^2/2! + x^4/4! - ...
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

// exp(x) = 1 + x + x^2/2! + ... + x^n/n!
constexpr rational exp(rational x, unsigned n) {
    rational out = 1;
    rational a = 1;
    for (unsigned i = 1; i <= n; i++) {
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
    simplify(x.num, y.num, z.num);
    simplify(x.den, y.den, z.den);
}

}
