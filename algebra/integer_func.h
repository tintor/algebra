#pragma once
#include "algebra/integer.h"

namespace algebra {

constexpr integer abs(integer a) {
    if (a.is_negative())
        a.negate();
    return a;
}

constexpr integer uniform_sample(const integer& min, const integer& max, auto& rng) {
    integer max_min = max - min;
    if (max_min.sign() < 0)
        throw std::runtime_error("max smaller than min in uniform_sample()");
    ++max_min;
    return integer(uniform_sample(max_min.abs, rng)) + min;
}

constexpr integer exp2(std::integral auto exp) {
    if (exp < 0)
        throw std::runtime_error("negative exponent in exp2(...)");
    integer out;
    out.abs.words.reset((exp + 64) / 64);
    out.abs.words.back() = integer::word(1) << (exp % 64);
    return out;
}

constexpr integer pow(integer base, std::integral auto exp) {
    if (base == 2)
        return exp2(exp);
    if (base == 4)
        return exp2(static_cast<uint64_t>(exp) * 2);
    if (base == 8)
        return exp2(static_cast<uint64_t>(exp) * 3);
    if (base.sign() > 0 && base.is_power_of_two())
        return exp2(static_cast<uint64_t>(exp) * base.num_trailing_zeros());
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(integer, ...)");
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

constexpr integer pow(integer base, std::integral auto exp, integer result) {
    if (base == 2)
        return result << exp;
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(integer, ...)");
    if (exp == 0)
        return 1;
    if (exp == 1)
        return result * base;
    if (base.sign() > 0 && base.is_power_of_two())
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
    if (t < 0)
        t += n;
    out = t.abs;
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

}
