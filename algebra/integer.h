#pragma once
#include "algebra/integer_class.h"
#include "algebra/natural.h"

namespace algebra {



// returns x such that (a * x) mod m == 1, (or false if such number doesn't exist)
constexpr bool inverse_mod(const natural& a, const natural& m, natural& out) {
    Check(!a.is_negative(), "inverse_mod() of a negative number");
    Check(!m.is_negative(), "inverse_mod() of a negative number");
    integer t = 0;
    integer r = m;
    integer new_t = 1;
    integer new_r = a;
    integer e, q, nt;

    while (new_r) {
        div(r, new_r, /*out*/q, /*out*/e); // e = r - q * new_r

        // (t, new_t) = (new_t, t − q * new_t)
        nt = t;
        sub_product(nt, q, new_t);
        t.swap(new_t);
        new_t.swap(nt);

        // (r, new_r) = (new_r, r − q * new_r)
        r.swap(new_r);
        new_r.swap(e);
    }
    if (r > 1)
        return false;
    if (t.is_negative())
        t += m;
    out = abs(t);
    return true;
}

// returns (n k) mod m
constexpr void mod(integer& a, const integer& b) {
    const bool negative = a.is_negative();
    a.words.set_negative(false);
    {
        const natural bn = abs(b); // b may be a itself, see mul() on aliasing
        auto m = magnitude(a);
        __mod_magnitude(*m, bn);
    }
    if (negative && !a.words.empty()) {
        // result is in [0, abs(b)) range
        natural e = abs(b);
        e -= abs(a);
        a = std::move(e);
    }
}

constexpr void binominal_mod(const natural& n, uint64_t k, const natural& m, natural& out) {
    Check(!n.is_negative(), "binominal_mod() of a negative number");
    Check(!m.is_negative(), "binominal_mod() of a negative number");
    // Fast path: multiply by the modular inverse of each i+1, which keeps out below m. It needs
    // every i+1 in [1, k] to be invertible mod m, i.e. m coprime with k!, so inverse_mod can fail.
    out = 1;
    natural e, inv;
    for (uint64_t i = 0; i < k; i++) {
        e = i;
        e += 1;
        if (!inverse_mod(e, m, inv)) {
            // i+1 shares a factor with m and has no inverse. Compute the coefficient exactly
            // instead, where every division is exact because the partial product is C(n, j+1),
            // and reduce only at the end. Slower, since the intermediate is not bounded by m.
            out = 1;
            for (uint64_t j = 0; j < k; j++) {
                e = n;
                e -= j;
                out *= e;
                e = j;
                e += 1;
                out /= e;
            }
            mod(out, m);
            return;
        }

        e = n;
        e -= i;
        out *= e;
        __mul_mod(out, inv, m);
    }
    mod(out, m); // k == 0 leaves out == 1, which still has to be reduced
}


constexpr int signum(const integer& a) {
    if (a.words.empty())
        return 0;
    return a.is_negative() ? -1 : 1;
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

// returns a * b < c (cheaper than naive multiplication)
constexpr bool less_ab_c(const integer& a, const integer& b, const integer& c) {
    int ab = signum(a) * signum(b);
    int cc = signum(c);
    if (ab != cc)
        return ab < cc;
    return (ab > 0) ? __less_ab_c(a, b, c) : __less_a_bc(c, a, b); // integer converts to cnatural
}

// returns a < b * c (cheaper than naive multiplication)
constexpr bool less_a_bc(const integer& a, const integer& b, const integer& c) {
    int aa = signum(a);
    int bc = signum(b) * signum(c);
    if (aa != bc)
        return aa < bc;
    return (aa > 0) ? __less_a_bc(a, b, c) : __less_ab_c(b, c, a);
}

// returns a * b < c * d (cheaper than naive multiplication)
// this can be useful for geometry algorithms!
constexpr bool less_ab_cd(const integer& a, const integer& b, const integer& c, const integer& d) {
    int ab = signum(a) * signum(b);
    int cd = signum(c) * signum(d);
    if (ab != cd)
        return ab < cd;
    return (ab > 0) ? __less_ab_cd(a, b, c, d) : __less_ab_cd(c, d, a, b);
}

}
