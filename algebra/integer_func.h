#pragma once
#include "algebra/integer.h"
#include "algebra/natural_func.h"

namespace algebra {

constexpr integer abs(integer a) {
    if (a.is_negative())
        a.negate();
    return a;
}

// returns abs(a) > abs(b), minimizing memory allocation
constexpr bool abs_greater(const integer& a, const integer& b) {
    return a.abs > b.abs;
}

constexpr integer uniform_sample(const integer& min, const integer& max, auto& rng) {
    integer max_min = max - min;
    if (max_min.sign() < 0)
        throw std::runtime_error("max smaller than min in uniform_sample()");
    ++max_min;
    return integer(uniform_sample(max_min, rng)) + min;
}

constexpr integer exp2(std_int auto exp) {
    if (exp < 0)
        throw std::runtime_error("negative exponent in exp2(...)");
    integer out = 1;
    out <<= exp;
    return out;
}

constexpr integer pow(integer base, std_int auto exp) {
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(integer, ...)");
    if (base == 2)
        return exp2(exp);
    if (base == 4)
        return exp2(static_cast<uint64_t>(exp) * 2);
    if (base == 8)
        return exp2(static_cast<uint64_t>(exp) * 3);
    if (is_power_of_two(base))
        return exp2(static_cast<uint64_t>(exp) * base.num_trailing_zeros());
    if (exp == 0)
        return 1;
    if (exp == 1)
        return base;
    if (exp == 2)
        return base * base;

    integer result = 1;
    if (exp & 1)
        result = base;
    exp >>= 1;
    while (exp) {
        base *= base;
        if (exp & 1)
            result *= base;
        exp >>= 1;
    }
    return result;
}

constexpr integer pow(integer base, std_int auto exp, integer result) {
    if (base == 2)
        return result << exp;
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(integer, ...)");
    if (exp == 0)
        return 1;
    if (exp == 1)
        return result * base;
    if (is_power_of_two(base))
        return result << (static_cast<uint64_t>(exp) * base.num_trailing_zeros());

    if (exp & 1)
        result *= base;
    exp >>= 1;
    while (exp) {
        base *= base;
        if (exp & 1)
            result *= base;
        exp >>= 1;
    }
    return result;
}

constexpr integer pow(integer base, const natural& exp) {
    if (exp.is_uint64())
        return pow(base, static_cast<uint64_t>(exp));
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(integer, ...)");

    integer result = 1;
    if (exp.is_odd())
        result = base;
    for (int i = 1; i < exp.num_bits(); i++) {
        base *= base;
        if (exp.bit(i))
            result *= base;
    }
    return result;
}

// returns x such that (a * x) mod n == 1, (or false if such number doesn't exist)
constexpr bool mod_inverse(const natural& a, const natural& n, natural& out) {
    integer t = 0;
    integer r = n;
    integer new_t = 1;
    integer new_r = a;
    integer e, q;

    while (new_r) {
        div(r, new_r, /*out*/q, /*out*/e); // remainder is discarded

        // (t, new_t) = (new_t, t − q * new_t)
        e = t;
        sub_product(e, q, new_t);
        t.swap(new_t);
        new_t.swap(e);

        // (r, new_r) = (new_r, r − quotient * new_r)
        e = r;
        sub_product(r, q, new_r);
        r.swap(new_r);
        new_r.swap(e);
    }
    if (r > 1)
        return false;
    if (t.sign() < 0)
        t += n;
    if (t.sign() < 0)
        t.negate();
    out = t;
    return true;
}

// returns (n k) mod p
constexpr void binominal_mod(const natural& n, uint64_t k, const natural& p, natural& out) {
    out = 1;
    natural e, inv;
    for (uint64_t i = 0; i < k; i++) {
        e = n;
        e -= i;
        out *= e;

        e = i;
        e += 1;
        mod_inverse(e, p, inv);

        out *= inv;
        out %= p;
    }
}

constexpr int signum(const integer& a) {
    return (a.sign() > 0) - (a.sign() < 0);
}

constexpr bool is_power_of_two(const integer& a) {
    return a.sign() > 0 && a.num_bits() == 1 + a.num_trailing_zeros();
}

// reduce vector's length, without changing vector's direction
constexpr void simplify(integer& x, integer& y) {
    integer a = gcd(x, y);
    if (a != 1) {
        x /= a;
        y /= a;
    }
}

constexpr void simplify(integer& x, integer& y, integer& z) {
    integer a = gcd(gcd(x, y), z);
    if (a != 1) {
        x /= a;
        y /= a;
        z /= a;
    }
}

}
