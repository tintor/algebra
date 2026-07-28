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
    // [0, count) is empty for a zero count. The fast path below would ask for a distribution over
    // [0, count - 1], where the subtraction wraps to the whole 64-bit range.
    Check(!count.is_zero(), "uniform_sample() with a zero count");
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

// returns the largest Q such that Q*Q <= x
constexpr uint64_t __isqrt_u128(const uint128_t x);

// Constrained so that a class type never matches it. With a plain uint64_t parameter, isqrt(integer)
// was ambiguous: integer converts to uint64_t just as readily as it binds to the integer overload.
template<std_unsigned_int T>
constexpr uint64_t isqrt(const T x) {
    // the clamp below is a one word bound, so a wider value needs the two word implementation
    if constexpr (sizeof(T) > 8) {
        if (x > UINT64_MAX)
            return __isqrt_u128(x);
    }
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




constexpr bool is_power_of_two(const integer& a) { return !a.is_negative() && is_power_of_two(static_cast<cnatural>(a)); }

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
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(integer, ...)");
    if (base == 2)
        return result << exp;
    if (exp == 0)
        return result;
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

constexpr integer pow(integer base, const integer& exp) {
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

// b = the value truncated towards zero, so a negative value keeps its sign
template<std::floating_point T>
constexpr void round_to_zero(const T& a, integer& b) {
    // frexp() leaves the exponent unspecified for these, and the cast of the mantissa below would
    // be undefined as well
    Check(!std::isnan(a), "round_to_zero() of nan");
    Check(!std::isinf(a), "round_to_zero() of infinity");

    int exponent;
    auto mantissa = std::frexp(a, &exponent);

    const int bits = std::numeric_limits<T>::digits;
    auto m = std::ldexp(mantissa, bits);
    exponent -= bits;

    const bool negative = m < 0;
    if (negative)
        m = -m;

    // the magnitude first, since a negative shift truncates the magnitude, which is what
    // truncation towards zero means on either side of zero
    b = static_cast<uint64_t>(m);
    b <<= exponent;
    if (negative)
        b.negate();
}

// very fast, but only approximate for large A
constexpr integer isqrt_hardware(const integer& a) {
    Check(!a.is_negative(), "isqrt_hardware() of a negative number");
    if (a <= 1)
        return a;

    const int FP_DIGITS = std::numeric_limits<double>::digits;

    int a_exp = static_cast<int>(a.num_bits()) - FP_DIGITS;
    int delta = 0;
    if (a_exp >= std::numeric_limits<double>::max_exponent) {
        delta = a_exp - (std::numeric_limits<double>::max_exponent - 1);
        delta += delta % 2;
        a_exp -= delta;
    }

    Check(a_exp < std::numeric_limits<double>::max_exponent);
    Check(delta % 2 == 0);

    double a_fp;
    if (a_exp <= 0) {
        a_fp = a.words[0];
    } else {
        const uint64_t m = extract_u64(a, a_exp);
        a_fp = std::ldexp(static_cast<double>(m), a_exp);
    }

    const double x_fp = std::sqrt(a_fp);

    int x_exp;
    auto x_mantissa = std::frexp(x_fp, &x_exp);

    const uint64_t x = std::ldexp(x_mantissa, FP_DIGITS);
    return integer(x) << (x_exp - FP_DIGITS + delta / 2);
}

// exact result for values of at most two words, shared by the isqrt implementations
constexpr bool __isqrt_small(const integer& a, integer& out) {
    if (a.words.size() <= 1) {
        uint64_t q = std::sqrt(static_cast<double>(a.words[0]));
        if (__mulq(q, q) > a.words[0])
            q -= 1;
        out = q;
        return true;
    }
    if (a.words.size() == 2) {
        const uint128_t ac = concat(a.words[1], a.words[0]);
        out = __isqrt_u128(ac);
        return true;
    }
    return false;
}

constexpr integer isqrt(const integer& a) {
    Check(!a.is_negative(), "isqrt() of a negative number");
    integer small;
    if (__isqrt_small(a, small))
        return small;

    integer y = power_of_two((a.num_bits() + 1) / 2);
    integer x, r;
    do {
        x = y;
        div(a, x, y, r);
        y += x;
        // TODO if (y.is_even() && r >= x/2) y += 1    | Would this 1) speed up iteration? 2) avoid mul() at the end?
        y >>= 1;
    } while (y < x);

    x += 1;
    mul(x, x, r);
    if (r <= a)
        return x;
    x -= 1;
    return x;
}

constexpr integer isqrt2(const integer& a) {
    Check(!a.is_negative(), "isqrt2() of a negative number");
    integer small;
    if (__isqrt_small(a, small))
        return small;

    integer x = power_of_two((a.num_bits() + 1) / 2);
    integer v = x * x;
    integer r;

    //int i = 0;
    while (v > a) {
        //std::print("{}\n", x);
        v -= a;
        div(v, x, v, r); // v is much smaller than a, which makes this division cheaper!
        v >>= 1;
        if (v == 0) {
            x -= 1;
            return x;
        }
        x -= v;
        mul(x, x, v);
        //Check(++i <= 10000);
    }

    if (v == a)
        return x;
    v += x;
    x += 1;
    v += x;
    if (v <= a)
        return x;
    x -= 1;
    return x;
}

constexpr integer isqrt3(const integer& a) {
    Check(!a.is_negative(), "isqrt3() of a negative number");
    integer small;
    if (__isqrt_small(a, small))
        return small;

    integer x = power_of_two((a.num_bits() + 1) / 2);
    integer x2 = x * x;
    integer r, v, m;

    while (x2 > a) {
        // newton step: v = (x*x - a) / (2*x)
        v = x2;
        v -= a;
        div(v, x, /*out*/v, /*out*/r); // v is much smaller than a, which makes this division cheaper!
        v >>= 1;
        if (v == 0) {
            x -= 1;
            return x;
        }
        // x*x decreases by (x_old + x_new) * (x_old - x_new)
        m = x;
        x -= v;
        m += x;
        m *= v;
        Check(m <= x2);
        x2 -= m;
    }

    if (x2 == a)
        return x;
    x2 += x;
    x += 1;
    x2 += x; // x2 = (x + 1) * (x + 1)
    if (x2 <= a)
        return x;
    x -= 1;
    return x;
}

constexpr integer iroot(const integer& a, uint32_t n) {
    Check(!a.is_negative(), "iroot() of a negative number");
    if (a <= 1 || n == 1)
        return a;
    if (n == 0)
        return 1;
    if (n == 2)
        return isqrt(a);

    // exact bracket from the bit length: 2**((bits-1)/n) <= root < 2**(bits/n + 1)
    const auto bits = a.num_bits();
    integer left = power_of_two((bits - 1) / n);
    integer right = power_of_two(bits / n + 1);

    integer m, mn, t, t2;

    // narrow it with a floating point estimate, but only after verifying each bound: the
    // estimate can land on either side of the root, and for small roots the window below
    // used to be empty, which left the true root outside the bracket
    const double estimate = std::pow(static_cast<double>(a), 1.0 / n);
    if (std::isfinite(estimate)) // a value above ~1e308 has no double estimate to narrow with
        round_to_zero(estimate, m);
    if (m > 1u) {
        const integer w = (m >> 19) + 2;
        integer lo = (m > w) ? (m - w) : integer(1);
        if (lo > left && pow(lo, n) <= a)
            left = std::move(lo);
        integer hi = m + w;
        if (hi < right && pow(hi, n) > a)
            right = std::move(hi);
    }

    while (left < right) {
        m = right;
        m -= left;
        ++m;
        m >>= 1;
        m += left;

        // mn = pow(m, n)
        if (n == 3) {
            mul(m, m, mn);
            mn *= m;
        } else if (n == 4) {
            mul(m, m, t);
            mul(t, t, mn);
        } else if (n == 5) {
            mul(m, m, t);
            mul(t, t, mn);
            mn *= m;
        } else if (n == 6) {
            mul(m, m, t);
            t *= m;
            mul(t, t, mn);
        } else {
            int c = n;
            mn = 1;
            if (c & 1)
                mn = m;
            t = m;
            c >>= 1;
            while (c) {
                mul(t, t, t2);
                std::swap(t, t2); // t = t * t
                if (c & 1)
                    mn *= t;
                c >>= 1;
            }
        }

        if (mn > a) {
            --m;
            right = m;
            continue;
        }
        if (mn < a) {
            left = m;
            continue;
        }
        return m;
    }
    return left;
}

// assumes are a and b are in [0, m-1] range
// a = (a + b) % m
constexpr void add_mod(integer& a, const integer& b, const integer& m) {
    a += b;
    if (a >= m)
        a -= m;
}

// assumes are a and b are in [0, m-1] range
// a = (a - b) % m
constexpr void sub_mod(integer& a, const integer& b, const integer& m) {
    if (b > a)
        a += m;
    a -= b;
}

// assumes are a and b are in [0, m-1] range
// a = (a * b) % m
constexpr void __mul_mod(integer& a, const integer& b, const integer& m) {
    // This is simple and slow implementation, for testing.
    a *= b;
    if (a >= m)
        a %= m;
}

constexpr void mul_mod(const integer& a, const integer& b, const integer& m, integer& out) {
    if (a == 1 || b == 0) {
        out = b;
        return;
    }
    if (b == 1 || a == 0) {
        out = a;
        return;
    }
    if (a.is_uint128() && b.is_uint128() && m.is_uint128()) {
        out = mul_mod(static_cast<uint128_t>(a), static_cast<uint128_t>(b), static_cast<uint128_t>(m));
        return;
    }

    if (a.num_bits() + b.num_bits() <= m.num_bits()) {
        out = a;
        out *= b;
        if (a.num_bits() + b.num_bits() == m.num_bits() && out >= m)
            out -= m;
        return;
    }

    out = 0;
    integer aa = a, bb = b;
    while (aa && bb) {
        if (aa < bb)
            std::swap(aa, bb);
        if (bb == 1) {
            add_mod(out, aa, m); // result = (result + aa) % m
            return;
        }
        if (bb.is_odd())
            add_mod(out, aa, m); // result = (result + aa) % m
        add_mod(aa, aa, m); // aa = (aa + aa) % m
        bb >>= 1;
    }
}

// returns (a**b) mod m
constexpr void pow_mod(integer a, const integer& b, const integer& m, integer& out) {
    out = 1;
    if (a >= m)
        a %= m;
    for (size_t i = 0; i < b.num_bits(); i++) {
        if (b.bit(i))
            __mul_mod(out, a, m);
        if (i == b.num_bits() - 1)
            break;
        __mul_mod(a, a, m);
    }
}

constexpr integer pow_mod(integer a, const integer& b, const integer& m) {
    integer out;
    pow_mod(std::move(a), b, m, /*out*/out);
    return out;
}

// Miller-Rabin algorithm
// It returns false if n is composite and returns true if n
// is probably prime.  k is an input parameter that determines
// accuracy level. Higher value of `rounds` indicates more accuracy.
constexpr bool is_likely_prime(const integer& n, int rounds) {
    Check(!n.is_negative(), "is_likely_prime() of a negative number");
    if (n.is_uint64())
        return is_prime(static_cast<uint64_t>(n));

    std::array<int, 40> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
        71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173};
    if (rounds > primes.size())
        throw std::runtime_error("rounds arg is too high");
    if (n.mod2() == 0 || n.mod3() == 0 || n.mod5() == 0)
        return false;

    const integer n_minus_1 = n - 1;
    auto s = n_minus_1.num_trailing_zeros();
    integer d = n_minus_1;
    d >>= s;
    integer x, a;

    for (int i = 0; i < rounds; i++) {
        a = primes[i];

        pow_mod(a, d, n, /*out*/x);
        if (x == 1 || x == n_minus_1)
            continue;

        int r = 1;
        for (; r < s; r++) {
            __mul_mod(x, x, n); // x = (x * x) % n
            if (x == n_minus_1)
                break;
        }
        if (r == s)
            return false;
    }
    return true;
}

constexpr std::pair<int, int> mod63_65(const integer& a) {
    Check(!a.is_negative(), "mod63_65() of a negative number");
    int m63 = 0;
    int m65 = 0;
    int i = 0;
    while (i < a.num_bits()) {
        uint64_t b = extract_u64(a, i);

        m63 += b % 64;
        if (m63 >= 63)
            m63 -= 63;

        m65 += b % 64;
        if (m65 >= 65)
            m65 -= 65;

        b >>= 6;
        m63 += b % 64;
        if (m63 >= 63)
            m63 -= 63;

        m65 -= b % 64;
        if (m65 < 0)
            m65 += 65;
        i += 12;
    }
    return {m63, m65};
}

// rejects ~98% of all numbers
constexpr bool is_possible_square(const integer& a) {
    Check(!a.is_negative(), "is_possible_square() of a negative number");
    if (!is_one_of(a.words[0] % 16, {0,1,4,9}))
        return false;

    auto [m63, m65] = mod63_65(a);
    return is_one_of(m63, {0,1,4,7,9,16,18,22,25,28,36,37,43,46,49,58})
        && is_one_of(m65, {0,1,4,9,10,14,16,25,26,29,30,35,36,39,40,49,51,55,56,61,64});
}

constexpr std::vector<std::pair<integer, int>> factorize(integer a) {
    Check(!a.is_negative(), "factorize() of a negative number");
    if (a <= 1)
        return {};
    std::vector<std::pair<integer, int>> out;

    auto count = a.num_trailing_zeros();
    if (count) {
        out.emplace_back(2, count);
        a >>= count;
        if (a == 1)
            return out;
    }

    if (a.mod3() == 0) {
        count = 1;
        a /= 3;
        while (a > 1 && a.mod3() == 0) {
            a /= 3;
            count += 1;
        }
        out.emplace_back(3, count);
        if (a == 1)
            return out;
    }

    int f = 1;
    while (is_possible_square(a)) {
        integer s = isqrt(a);
        if (s * s != a)
            break;
        a = s;
        f *= 2;
    }
    if (is_likely_prime(a, 40)) {
        out.emplace_back(a, f);
        return out;
    }
    if (a.is_uint64()) {
        for (auto e : factorize(static_cast<uint64_t>(a)))
            out.push_back({e.first, e.second * f});
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

            while (is_possible_square(a)) {
                integer s = isqrt(a);
                if (s * s != a)
                    break;
                a = s;
                f *= 2;
            }
            if (is_likely_prime(a, 40)) {
                out.emplace_back(a, f);
                break;
            }
            if (a.is_uint64()) {
                for (auto e : factorize(static_cast<uint64_t>(a)))
                    out.push_back({e.first, e.second * f});
                return out;
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

            while (is_possible_square(a)) {
                integer s = isqrt(a);
                if (s * s != a)
                    break;
                a = s;
                f *= 2;
            }
            if (is_likely_prime(a, 40)) {
                out.emplace_back(a, f);
                break;
            }
            if (a.is_uint64()) {
                for (auto e : factorize(static_cast<uint64_t>(a)))
                    out.push_back({e.first, e.second * f});
                return out;
            }
        }
        p += 4;
        if (p < 5)
            throw std::runtime_error("overflow");
    }
    return out;
}

// returns (n k)
constexpr void binominal(const integer& n, uint64_t k, integer& out) {
    Check(!n.is_negative(), "binominal() of a negative number");
    out = 1;
    integer e;
    for (uint64_t i = 0; i < k; i++) {
        e = n;
        e -= i;
        out *= e;
        out /= i + 1;
    }
}

constexpr bool __exact_sqrt1(const integer& a, integer& b) {
    auto z = a.num_trailing_zeros();
    if (z & 1)
        return false;

    b = a;
    b >>= z / 2;
    return is_possible_square(b);
}

constexpr bool __exact_sqrt2(const integer& a, integer& b) {
    // option 1) try to factorize with small factors (fast, but doesn't work for all numbers)
    // option 2) isqrt_nr() -> it detects if number is perfect square

    integer s = isqrt(a);
    mul(s, s, b);
    if (b != a)
        return false;
    b = std::move(s);
    return true;
}

constexpr bool exact_sqrt(const integer& a, integer& b) { return __exact_sqrt1(a, b) && __exact_sqrt2(a, b); }

// assumes that whole and root are already initialized
constexpr void exact_sqrt(integer a, integer& whole, integer& root) {
    Check(!a.is_negative(), "exact_sqrt() of a negative number");
    if (a <= 1)
        return;

    // factorize a
    auto z = a.num_trailing_zeros();
    if (z) {
        if (z > 1)
            whole <<= z / 2;
        if (z & 1) {
            if (root.is_even()) {
                root >>= 1;
                whole <<= 1;
            } else
                root <<= 1;
        }
        a >>= z;
    }

    std::optional<integer> a_sqrt;
    if (is_possible_square(a)) {
        integer s = isqrt(a);
        if (s * s == a) {
            whole *= s;
            return;
        }
        a_sqrt = std::move(s);
    }

    uint64_t p = 3;
    while (a > 1) {
        if (p > 256) {
            if (a_sqrt == std::nullopt) {
                integer s = isqrt(a);
                if (s * s == a) {
                    whole *= s;
                    return;
                }
                a_sqrt = std::move(s);
                if (is_likely_prime(a, 40))
                    break;
            }
            if (p > *a_sqrt)
                break;
        } else
            if (p >= a)
                break;

        int count = 0;
        while (a % p == 0) {
            a /= p;
            count += 1;
        }
        if (count) {
            if (count >= 2)
                whole *= pow(integer(p), count / 2);
            if (count & 1) {
                if (root % p == 0) {
                    root /= p;
                    whole *= p;
                } else
                    root *= p;
            }
            a_sqrt = std::nullopt;
        }
        p += 2;
    }
    if (a > 1) {
        if (root % a == 0) {
            root /= a;
            whole *= a;
        } else
            root *= a;
    }
}


// returns abs(a) > abs(b), minimizing memory allocation
constexpr bool abs_greater(const integer& a, const integer& b) {
    return __less(static_cast<cnatural>(b), static_cast<cnatural>(a)); // no allocation: views only
}

constexpr integer uniform_sample(const integer& min, const integer& max, auto& rng) {
    integer max_min = max - min;
    if (max_min.sign() < 0)
        throw std::runtime_error("max smaller than min in uniform_sample()");
    ++max_min;
    return uniform_sample(max_min, rng) + min;
}



// returns result * (base ** exp)


// returns x such that (a * x) mod m == 1, (or false if such number doesn't exist)
constexpr bool inverse_mod(const integer& a, const integer& m, integer& out) {
    Check(!a.is_negative() && !m.is_negative(), "inverse_mod() of a negative number");
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
        const integer bn = abs(b); // b may be a itself, see mul() on aliasing
        auto m = magnitude(a);
        __abs_mod(*m, bn); // by name: mod() here would resolve back to this function
    }
    if (negative && !a.words.empty()) {
        // result is in [0, abs(b)) range
        integer e = abs(b);
        e -= abs(a);
        a = std::move(e);
    }
}

constexpr void binominal_mod(const integer& n, uint64_t k, const integer& m, integer& out) {
    Check(!n.is_negative() && !m.is_negative(), "binominal_mod() of a negative number");
    // Fast path: multiply by the modular inverse of each i+1, which keeps out below m. It needs
    // every i+1 in [1, k] to be invertible mod m, i.e. m coprime with k!, so inverse_mod can fail.
    out = 1;
    integer e, inv;
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



constexpr uint64_t log_lower(const integer& n, uint64_t base) {
    Check(!n.is_negative(), "log_lower() of a negative number");
    // dividing by one never reaches zero, so the loop below would not terminate
    Check(base >= 2, "log_lower() with a base below two");
    integer a = abs(n);
    uint64_t count = 0;
    if (!a)
        return count;
    while (true) {
        a /= base;
        if (!a)
            break;
        count += 1;
    }
    return count;
}

constexpr uint64_t log_upper(const integer& n, uint64_t base) {
    Check(!n.is_negative(), "log_upper() of a negative number");
    Check(base >= 2, "log_upper() with a base below two");
    integer a = abs(n);
    uint64_t count = 0;
    while (a) {
        a /= base;
        count += 1;
    }
    return count;
}

constexpr bool is_power_of_three(const integer& n) {
    Check(!n.is_negative(), "is_power_of_three() of a negative number");
    integer a = n; // isqrt() takes an integer now, so no magnitude copy is needed
    if (a.words.empty())
        return false;
    integer m;
    while (a > 1) {
        if (a.mod3())
            return false;
        // a power of three that is also a perfect square stays a power of three when
        // halved in exponent, and the square root is much cheaper to keep dividing
        if (is_possible_square(a)) {
            integer s = isqrt(a);
            mul(s, s, m);
            if (m == a) {
                a = std::move(s);
                continue;
            }
        }
        a /= 3u;
    }
    return true;
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
