#pragma once
#include "algebra/integer_class.h"
#include <array>
#include <optional>
#include <random>
#include <span>
#include <vector>

namespace algebra {

constexpr integer power_of_two(size_t e) {
    integer x;
    x.words.reset((e + 64) / 64);
    x.words.back() = uint64_t(1) << (e % 64);
    return x;
}

// Moved here from integer.h: iroot() and the isqrt family below call pow(), and picking up
// util.h's pow(uint64_t, unsigned) instead would silently truncate a big operand.

// uniformly sample from [0, (2**n)-1]
constexpr void uniform_sample_bits(const size_t n, auto& rng, integer& out) {
    static_assert(sizeof(rng()) == 8);
    auto w = (n + 64 - 1) / 64;
    out.words.reset(w, /*initialize*/false);
    while (w--)
        out.words[w] = rng();
    if ((n % 64) != 0)
        out.words.back() &= (uint64_t(1) << (n % 64)) - 1;
    out.words.normalize();
}

constexpr integer uniform_sample_bits(const size_t n, auto& rng) {
    integer out;
    uniform_sample_bits(n, rng, /*out*/out);
    return out;
}

// uniformly sample from [0, count-1]
constexpr void uniform_sample(const integer& count, auto& rng, integer& out) {
    Check(!count.is_negative(), "uniform_sample() with a negative count");
    if (count.is_uint64()) {
        out = std::uniform_int_distribution<uint64_t>(0, static_cast<uint64_t>(count) - 1)(rng);
        return;
    }
    if (count.words.size() == 2 && count.words[0] == 0 && count.words[1] == 1) {
        static_assert(sizeof(rng()) == 8);
        out = rng();
        return;
    }

    const auto z = count.num_trailing_zeros();
    const auto b = count.num_bits();
    if (b == z + 1) { // if count is power of 2
        uniform_sample_bits(z, rng, /*out*/out);
        return;
    }

    const auto n = (b + 64 - 1) / 64;
    if (b == 64 * n) {
        // Note: power of two case is handled above!
        // mq = pow(2_i, 64 * n) / count
        // assert mq == 1
        while (true) {
            uniform_sample_bits(n * 64, rng, /*out*/out);
            if (out < count)
                return;
        }
    }
    // TODO avoid this division in mq == 2 case
    integer temp;
    integer mq = power_of_two(n * 64);
    mq /= count;
    if (mq == 2) {
        temp = count;
        temp <<= 1;
    }
    while (true) {
        if (mq == 2) {
            // optimization: avoid expensive div() for small mq
            uniform_sample_bits(n * 64, rng, /*out*/out);
            if (out < temp) {
                if (out >= count)
                    out -= count;
                return;
            }
        } else {
            uniform_sample_bits(n * 64, rng, /*out*/temp);
            div(temp, count, /*out*/temp, /*out*/out);
            if (temp < mq)
                return;
        }
    }
}

constexpr integer uniform_sample(const integer& count, auto& rng) {
    integer out;
    uniform_sample(count, rng, /*out*/out);
    return out;
}


template<std_unsigned_int T>
constexpr T __gcd_inner(T a, T b) {
    while (b) {
        b >>= countr_zero(b); // since b is non-zero (__builtin_ctzl would truncate T to 64 bits)
        if (a > b)
            std::swap(a, b);
        b -= a;
    }
    return a;
}

constexpr auto gcd(std_int auto a, std_int auto b) -> make_unsigned_t<larger_type<decltype(a), decltype(b)>> {
    using T = make_unsigned_t<larger_type<decltype(a), decltype(b)>>;
    T ua = abs_unsigned(a);
    T ub = abs_unsigned(b);
    if (ua == 0)
        return ub;
    if (ub == 0)
        return ua;
    auto az = countr_zero(ua);
    if (az == 0)
        return __gcd_inner(ua, ub);
    auto common = std::min(az, countr_zero(ub));
    return __gcd_inner(ua >> az, ub >> common) << common;
}

// the greatest common divisor of the magnitudes, so the sign of either argument does not matter
constexpr integer gcd(integer a, integer b) {
    a.words.set_negative(false);
    b.words.set_negative(false);
    if (a.words.size() == 1 && b.words.size() == 1)
        return gcd(a.words[0], b.words[0]);

    if (a.words.size() == 0)
        return b;
    if (b.words.size() == 0)
        return a;

    auto az = a.num_trailing_zeros();
    auto common = std::min(az, b.num_trailing_zeros());
    a >>= az;
    b >>= common;

    while (b) {
        if (a.words.size() == 1 && b.words.size() == 1) {
            a = __gcd_inner(a.words[0], b.words[0]);
            break;
        }
        b >>= b.num_trailing_zeros();
        if (a > b)
            a.swap(b);
        b -= a;
    }
    a <<= common;
    return a;
}

constexpr integer gcd(integer a, std_int auto b) { return gcd(std::move(a), integer(abs_unsigned(b))); }
constexpr integer gcd(std_int auto a, integer b) { return gcd(integer(abs_unsigned(a)), std::move(b)); }

// least common multiple
// likewise of the magnitudes
constexpr integer lcm(const integer& a, const integer& b) {
    if (a.words.empty() || b.words.empty())
        return 0;
    // divide first: a * b would be twice as long as the result
    integer m = a;
    m /= gcd(a, b);
    m *= b;
    return m;
}

// Constrained so that a class type never matches it. With a plain uint64_t parameter, isqrt(integer)
// was ambiguous: integer converts to uint64_t just as readily as it binds to the integer overload.
template<std_unsigned_int T>
constexpr uint64_t isqrt(const T x) {
    constexpr uint64_t MAX = UINT32_MAX; // isqrt(UINT64_MAX), also the largest q with q*q <= UINT64_MAX
    uint64_t q = static_cast<uint64_t>(std::sqrt(static_cast<double>(x)));
    if (q > MAX)
        q = MAX; // q * q below must not overflow
    if (q * q > x)
        q -= 1;
    else if (q < MAX && (q + 1) * (q + 1) <= x)
        q += 1; // double can round a large perfect square down
    return q;
}

// returns the largest Q such that Q*Q <= x
constexpr uint64_t __isqrt_u128(const uint128_t x) {
    if (x <= UINT64_MAX)
        return isqrt(static_cast<uint64_t>(x));

    const double d = std::sqrt(static_cast<double>(x));
    // double has 53 bits of mantissa, so the estimate can be off by up to ~2**11
    uint64_t q = (d >= 18446744073709551616.0) ? UINT64_MAX : static_cast<uint64_t>(d);
    if (q == 0)
        return 0;
    q = static_cast<uint64_t>(std::min<uint128_t>((static_cast<uint128_t>(q) + x / q) / 2, UINT64_MAX)); // newton step
    while (q > 0 && static_cast<uint128_t>(q) * q > x)
        q -= 1;
    while (q != UINT64_MAX && static_cast<uint128_t>(q + 1) * (q + 1) <= x)
        q += 1;
    return q;
}

/*
constexpr integer isqrt(const integer& x) {
    using u128 = unsigned __int128;
    using u64 = uint64_t;

    if (x.words.size() == 0)
        return 0;

    const int bits = x.num_bits();
    const int exponent = bits - 53;
    if (exponent > 1023)
        return __slow_isqrt(x); // too large for double

    // convert x to double
    const double xd = (exponent <= 0) ? x.words[0] : std::ldexp(static_cast<double>(extract_64bits(x, exponent)), exponent);

    // initial guess using hardware double (max error of 1 for up to 107 bits)
    u64 iq = static_cast<u64>(std::sqrt(xd));

    if (bits <= 64) {
        u128 cx = x.words[0];
        if (static_cast<u128>(iq) * iq > cx)
            iq -= 1;
        return iq;
    }

    if (bits <= 77) {
        u128 cx = (static_cast<u128>(x.words[1]) << 64) | x.words[0];
        iq += 1;
        if (static_cast<u128>(iq) * iq > cx)
            iq -= 1;
        if (static_cast<u128>(iq) * iq > cx)
            iq -= 1;
        return iq;
    }

    if (bits <= 126) {
        u128 cx = (static_cast<u128>(x.words[1]) << 64) | x.words[0];
        return (iq + static_cast<u64>(cx / iq)) / 2;
    }
}
*/

constexpr uint128_t concat(const uint64_t a, const uint64_t b) {
    return (static_cast<uint128_t>(a) << 64) | b;
}

// Note: there is bool inverse_mod() in integer_func.h!

// Deterministic Miller-Rabin algorithm for small numbers
constexpr bool is_prime(const uint64_t n) {
    if (n <= 4)
        return n == 2 || n == 3;
    if (n % 2 == 0)
        return false;

    const auto s = std::countr_zero(n - 1);
    const uint64_t d = (n - 1) >> s;

    static constexpr std::array<uint32_t, 3> bases32 = {2, 7, 61};
    static constexpr std::array<uint32_t, 7> bases64 = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    for (uint32_t a : (n <= UINT32_MAX) ? std::span<const uint32_t>(bases32) : std::span<const uint32_t>(bases64)) {
        if (a >= n)
            continue;
        uint64_t x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1)
            continue;

        int r = 1;
        for (; r < s; r++) {
            x = (static_cast<__uint128_t>(x) * x) % n;
            if (x == n - 1)
                break;
        }
        if (r == s)
            return false;
    }
    return true;
}

constexpr bool is_one_of(int a, std::initializer_list<int> b) {
    for (int p : b)
        if (a == p)
            return true;
    return false;
}

inline uint64_t try_fermat_factorize(uint64_t n) {
    if (n % 2 == 0)
        return 2;

    uint64_t a = isqrt(n);
    if (a * a < n)
        a++;
    if (a * a == n)
        return a;

    for (int i = 0; i < 100'000; i++) {
        const uint128_t a_sq = static_cast<uint128_t>(a) * a;
        const uint128_t b_sq = a_sq - n;
        const uint64_t b = __isqrt_u128(b_sq);
        if (static_cast<uint128_t>(b) * b == b_sq)
            return a - b;
        if (a == UINT64_MAX)
            return 0;
        a++;
    }
    return 0;
}

constexpr bool is_likely_prime(const integer& n, int rounds);

// Constrained for the same reason as isqrt() above: with a plain uint64_t parameter,
// factorize(integer) was ambiguous between this and the integer overload.
template<std_unsigned_int T>
constexpr std::vector<std::pair<uint64_t, int>> factorize(T a) {
    if (a <= 1)
        return {};
    std::vector<std::pair<uint64_t, int>> out;
    auto z = std::countr_zero(a);
    if (z) {
        out.emplace_back(2, z);
        a >>= z;
    }
    if (a == 1)
        return out;

    int count = 0;
    while (a % 3 == 0) {
        a /= 3;
        count += 1;
    }
    if (count) {
        out.emplace_back(3, count);
        if (a == 1)
            return out;
    }

    int f = 1;
    while (true) {
        uint64_t s = isqrt(a);
        if (s * s != a)
            break;
        a = s;
        f *= 2;
    }

    if (is_prime(a)) {
        out.emplace_back(a, f);
        return out;
    }

    uint64_t p = 5;
    while (true) {
        if (a % p == 0) {
            int count = f;
            a /= p;
            while (a % p == 0) {
                a /= p;
                count += f;
            }
            out.emplace_back(p, count);
            if (a == 1)
                break;
            if (is_prime(a)) {
                out.emplace_back(a, f);
                break;
            }
        }
        p += 2;

        if (a % p == 0) {
            int count = f;
            a /= p;
            while (a % p == 0) {
                a /= p;
                count += f;
            }
            out.emplace_back(p, count);
            if (a == 1)
                break;
            if (is_prime(a)) {
                out.emplace_back(a, f);
                break;
            }
        }
        p += 4;
    }
    return out;
}

constexpr void invert_bits(integer& a) {
    Check(!a.is_negative(), "invert_bits() of a negative number");
    for (int i = 0; i < a.words.size(); i++)
        a.words[i] = ~a.words[i];
    a.words.normalize();
}

constexpr void complement(integer& a) {
    Check(!a.is_negative(), "complement() of a negative number");
    __complement(a);
    a.words.normalize();
}

}
