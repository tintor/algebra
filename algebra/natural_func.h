#pragma once
#include "algebra/natural.h"
#include <span>

namespace algebra {

constexpr natural pow(natural base, std::integral auto exp) {
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(natural, ...)");
    if (base == 2) {
        natural out;
        out.words.reset((exp + 64) / 64);
        out.words.back() = uint64_t(1) << (exp % 64);
        return out;
    }

    if (exp == 0)
        return 1;

    natural result = 1;
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

constexpr natural pow(natural base, const natural& _exp) {
    if (_exp.is_uint64())
        return pow(base, static_cast<uint64_t>(_exp));

    natural result = 1;
    natural exp = _exp;
    while (exp) {
        if (exp.is_odd())
            result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

// uniformly sample from [0, (2**n)-1]
constexpr void uniform_sample_bits(const size_t n, auto& rng, natural& out) {
    static_assert(sizeof(rng()) == 8);
    auto w = (n + 64 - 1) / 64;
    out.words.reset(w, /*initialize*/false);
    while (w--)
        out.words[w] = rng();
    if ((n % 64) != 0)
        out.words.back() &= (uint64_t(1) << (n % 64)) - 1;
    out.words.normalize();
}

constexpr natural uniform_sample_bits(const size_t n, auto& rng) {
    natural out;
    uniform_sample_bits(n, rng, /*out*/out);
    return out;
}

// uniformly sample from [0, count-1]
constexpr void uniform_sample(const natural& count, auto& rng, natural& out) {
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

    const natural::size_type n = (b + 64 - 1) / 64;
    if (b == 64 * n) {
        // Note: power of two case is handled above!
        // mq = pow(2_n, 64 * n) / count
        // assert mq == 1
        while (true) {
            uniform_sample_bits(n * 64, rng, /*out*/out);
            if (out < count)
                return;
        }
    }
    // TODO avoid this division in mq == 2 case
    natural temp;
    natural mq = pow(natural(2), n * 64);
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

constexpr natural uniform_sample(const natural& count, auto& rng) {
    natural out;
    uniform_sample(count, rng, /*out*/out);
    return out;
}

constexpr natural uniform_sample(const natural& min, const natural& max, auto& rng) {
    natural count = max;
    count -= min;
    count += 1;
    natural out;
    uniform_sample(count, rng, /*out*/out);
    out += min;
    return out;
}

template<std::unsigned_integral T>
constexpr T __gcd_inner(T a, T b) {
    while (b) {
        b >>= __builtin_ctzl(b); // since b is non-zero
        if (a > b)
            std::swap(a, b);
        b -= a;
    }
    return a;
}

template <typename A, typename B>
using larger_type = std::conditional_t<(sizeof(A) >= sizeof(B)), A, B>;

template<std::integral A, std::integral B>
constexpr auto gcd(A _a, B _b) {
    std::make_unsigned_t<larger_type<A, B>> a = (_a < 0) ? -_a : _a;
    std::make_unsigned_t<larger_type<A, B>> b = (_b < 0) ? -_b : _b;
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    auto az = std::countr_zero(a);
    if (az == 0)
        return __gcd_inner(a, b);
    auto common = std::min(az, std::countr_zero(b));
    return __gcd_inner(a >> az, b >> common) << common;
}

constexpr natural gcd(natural a, natural b) {
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

constexpr natural gcd(natural a, std::integral auto b) { return gcd(std::move(a), natural(make_unsigned((b < 0) ? -b : b))); }
constexpr natural gcd(std::integral auto a, natural b) { return gcd(std::move(b), a); }

// least common multiple
constexpr natural lcm(const natural& a, const natural& b) {
    natural m = a * b;
    m /= gcd(a, b);
    return m;
}

constexpr natural __slow_isqrt(const natural& x) {
    const auto b = x.num_bits();
    if (b <= 1)
        return x;

    auto e = b + (b & 1);
    natural q, a, a2, xe;
    int i = 0;
    uint64_t qi = 0;
    while (true) {
        xe = x;
        xe >>= e;
        if (i < 64) {
            qi <<= 1;
            if (xe.words.size() > 2 || static_cast<unsigned __int128>(qi | 1) * (qi | 1) <= xe)
                qi |= 1;
        } else {
            if (i == 64) {
                q = qi;
            }
            q <<= 1;
            a = q;
            a |= 1;
            mul(a, a, a2);
            if (a2 <= xe)
                q |= 1;
        }
        if (e < 2)
            break;
        e -= 2;
        i++;
    }
    if (i < 64)
        q = qi;
    return q;
}

constexpr uint64_t isqrt(const uint64_t x) {
    uint64_t q = static_cast<uint64_t>(std::sqrt(x));
    if (q * q > x)
        q -= 1;
    return q;
}

constexpr natural isqrt(const natural& x) {
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

    /*if (bits <= 107) {
        q = iq + 1;
        mul(q, q, e);
        if (e > x)
            q -= 1;
        mul(q, q, e);
        if (e > x)
            q -= 1;
        return q;
    }*/

    if (bits <= 126) {
        u128 cx = (static_cast<u128>(x.words[1]) << 64) | x.words[0];
        return (iq + static_cast<u64>(cx / iq)) / 2;
    }

    natural q, e;
    // this doesn't work for 129 bits!
    /*natural rem;
    div(x, q, e, rem);
    q += e;
    q >>= 1;

    div(x, q, e, rem);
    q += e;
    q >>= 1;
    return q;*/

    return __slow_isqrt(x);
}

// returns (a ** n) mod p
constexpr uint64_t pow_mod(uint64_t a, uint64_t n, uint64_t p) {
    uint64_t b = 1;
    a %= p;
        while (n) {
        if (n & 1)
            b = (static_cast<__uint128_t>(b) * a) % p;
        n >>= 1;
        a = (static_cast<__uint128_t>(a) * a) % p;
    }
    return b;
}

// assumes are a and b are in [0, m-1] range
// a = (a * b) % m
constexpr void mul_mod(natural& a, const natural& b, const natural& m) {
    // TODO optimize
    a *= b;
    if (a >= m)
        a %= m;
}

// returns (a**n) mod p
constexpr void pow_mod(natural a, natural n, const natural& p, natural& out) {
    out = 1;
    if (a >= p)
        a %= p;
    for (size_t i = 0; i < n.num_bits(); i++) {
        if (n.bit(i))
            mul_mod(out, a, p);
        mul_mod(a, a, p);
    }
}

constexpr natural pow_mod(natural a, natural n, natural p) {
    natural out;
    pow_mod(a, n, p, /*out*/out);
    return out;
}

// Deterministic Miller-Rabin algorithm for small numbers
constexpr bool is_prime(const uint64_t n) {
    if (n <= 4)
        return n == 2 || n == 3;
    if (n % 2 == 0)
        return false;

    const auto s = std::countr_zero(n - 1);
    const uint64_t d = (n - 1) >> s;

    std::array<uint32_t, 3> bases32 = {2, 7, 61};
    std::array<uint32_t, 7> bases64 = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    for (uint32_t a : (n <= UINT32_MAX) ? std::span<uint32_t>(bases32) : std::span<uint32_t>(bases64)) {
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

// Miller-Rabin algorithm
// It returns false if n is composite and returns true if n
// is probably prime.  k is an input parameter that determines
// accuracy level. Higher value of `rounds` indicates more accuracy.
constexpr bool is_likely_prime(const natural& n, int rounds) {
    if (n.is_uint64())
        return is_prime(static_cast<uint64_t>(n));

    std::array<int, 40> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
        71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173};
    if (rounds > primes.size())
        throw std::runtime_error("rounds arg is too high");
    if (n.mod2() == 0 || n.mod3() == 0 || n.mod5() == 0)
        return false;

    const natural n_minus_1 = n - 1;
    auto s = n_minus_1.num_trailing_zeros();
    natural d = n_minus_1;
    d >>= s;
    natural x, a;

    for (int i = 0; i < rounds; i++) {
        a = primes[i];

        pow_mod(a, d, n, /*out*/x);
        if (x == 1 || x == n_minus_1)
            continue;

        int r = 1;
        for (; r < s; r++) {
            mul_mod(x, x, n); // x = (x * x) % n
            if (x == n_minus_1)
                break;
        }
        if (r == s)
            return false;
    }
    return true;
}

constexpr bool is_possible_square(const natural& a) {
    auto w = a.words[0] % 128;
    if (w > 57)
        return false;
    if (w < 25)
        return w == 0 || w == 1 || w == 4 || w == 9 || w == 16 || w == 17;
    return w == 25 || w == 33 || w == 36 || w == 41 || w == 49 || w == 57;
}

uint64_t try_fermat_factorize(uint64_t n) {
    if (n % 2 == 0)
        return 2;

    uint64_t a = isqrt(n);
    if (a * a < n)
        a++;
    if (a * a == n)
        return a;

    for (int i = 0; i < 100'000; i++) {
        uint128_t a_sq = (uint128_t)a * a;
        uint128_t b_sq = a_sq - n;
        uint64_t b = isqrt(b_sq);
        if (b * b == b_sq)
            return a - b;
        if (a == UINT64_MAX)
            return 0;
        a++;
    }
    return 0;
}

constexpr bool is_likely_prime(const natural& n, int rounds);

constexpr std::vector<std::pair<uint64_t, int>> factorize(uint64_t a) {
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

constexpr std::vector<std::pair<natural, int>> factorize(natural a) {
    if (a <= 1)
        return {};
    std::vector<std::pair<natural, int>> out;

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
        natural s = isqrt(a);
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
                natural s = isqrt(a);
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
                natural s = isqrt(a);
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

constexpr bool is_power_of_two(const natural& a) {
    return a.num_bits() == 1 + a.num_trailing_zeros();
}

constexpr uint64_t log_lower(natural a, uint64_t base) {
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
constexpr uint64_t log_upper(natural a, uint64_t base) {
    uint64_t count = 0;
    while (a) {
        a /= base;
        count += 1;
    }
    return count;
}

// returns (n k)
constexpr void binominal(const natural& n, uint64_t k, natural& out) {
    out = 1;
    natural e;
    for (uint64_t i = 0; i < k; i++) {
        e = n;
        e -= i;
        out *= e;
        out /= i + 1;
    }
}

// assumes that whole and root are already initialized
constexpr void exact_sqrt(natural a, natural& whole, natural& root) {
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

    std::optional<natural> a_sqrt;
    if (is_possible_square(a)) {
        natural s = isqrt(a);
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
                natural s = isqrt(a);
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
                whole *= pow(natural(p), count / 2);
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

}
