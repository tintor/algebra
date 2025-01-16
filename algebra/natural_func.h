#pragma once
#include "algebra/natural.h"

namespace algebra {

constexpr natural pow(natural base, std::integral auto exp) {
    if (exp < 0)
        throw std::runtime_error("negative exponent in pow(natural, ...)");
    if (base == 2) {
        natural out;
        out.words.reset((exp + 64) / 64);
        out.words.back() = natural::word(1) << (exp % 64);
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
    auto w = (n + 63) / 64;
    out.words.reset(w, /*initialize*/false);
    while (w--)
        out.words[w] = rng();
    if ((n % 64) != 0)
        out.words.back() &= (natural::word(1) << (n % 64)) - 1;
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

    const natural::size_type n = (b + 63) / 64;
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

template<std::unsigned_integral T>
constexpr T gcd(T a, T b) {
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    auto az = num_trailing_zeros(a);
    if (az == 0)
        return __gcd_inner(a, b);
    auto common = std::min(az, num_trailing_zeros(b));
    return __gcd_inner(a >> az, b >> common) << common;
}

constexpr uint64_t gcd(uint64_t a, uint64_t b) {
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    auto common = __builtin_ctzl(a | b);
    a >>= __builtin_ctzl(a);
    do {
        b >>= __builtin_ctzl(b);
        if (a > b)
            std::swap(a, b);
        b -= a;
    } while (b);
    return a << common;
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

constexpr natural gcd(std::integral auto a, natural b) { return gcd(natural(abs_ulong(a)), std::move(b)); }
constexpr natural gcd(natural a, std::integral auto b) { return gcd(std::move(a), natural(abs_ulong(b))); }

// least common multiple
constexpr natural lcm(const natural& a, const natural& b) {
    natural m = a * b;
    m /= gcd(a, b);
    return m;
}

constexpr natural isqrt(const natural& x) {
    if (x == 0 || x == 1)
        return x;

    const auto b = x.num_bits();
    auto e = b + (b & 1);
    natural q, a;
    while (true) {
        q <<= 1;
        a = q;
        a |= 1;
        if (a * a <= (x >> e))
            q = a;
        if (e < 2)
            break;
        e -= 2;
    }
    return q;
}

constexpr bool is_prime(const uint64_t a) {
    if (a <= 3)
        return a >= 2;
    if (a % 2 == 0 || a % 3 == 0)
        return false;
    uint32_t i = 5;
    while (true) {
        // overflow here is not possible
        if (uint64_t(i) * i > a)
            return true;
        if (a % i == 0 || a % (i + 2) == 0)
            return false;
        i += 6;
        // uint32_t oveflow check
        if (i <= 5)
            return true;
    }
}

constexpr bool is_prime(const natural& a) {
    if (a.is_uint64())
        return is_prime(static_cast<uint64_t>(a));
    if (a <= 3)
        return a >= 2;
    if (a.is_even() == 0 || a.mod3() == 0)
        return false;
    uint64_t i = 5;
    while (true) {
        if (i * i > a)
            return true;
        if (a % i == 0 || a % (i + 2) == 0)
            return false;
        i += 6;
        // uint64_t oveflow check
        if (i <= 5)
            break;
    }
    natural ni = i - 6; // rollback to i before overflow
    ni += 6; // overflow safe
    natural q;
    while (true) {
        mul(ni, ni, /*out*/q); // q = ni * ni
        if (q > a)
            return true;

        mod(a, ni, /*out*/q); // q = a % ni
        if (!q)
            return false;

        ni += 2;
        mod(a, ni, /*out*/q); // q = a % ni
        if (!q)
            return false;

        ni += 4;
    }
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

// Miller-Rabin algorithm
// It returns false if n is composite and returns true if n
// is probably prime.  k is an input parameter that determines
// accuracy level. Higher value of `rounds` indicates more accuracy.
constexpr bool is_likely_prime(const uint64_t n, int rounds, auto& rng) {
    if (n <= 4)
        return n == 2 || n == 3;
    if (n % 2 == 0)
        return false;

    const auto s = num_trailing_zeros(n - 1);
    const uint64_t d = (n - 1) >> s;

    for (int i = 0; i < rounds; i++) {
        const uint64_t a = std::uniform_int_distribution<uint64_t>(2, n-2)(rng);

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

constexpr bool is_likely_prime(const natural& n, int rounds, auto& rng) {
    if (n.is_uint64())
        return is_likely_prime(static_cast<uint64_t>(n), rounds, rng);

    if (n <= 4)
        return n == 2 || n == 3;
    if (n.is_even())
        return false;

    const natural n_minus_1 = n - 1;
    const natural n_minus_3 = n - 3;
    auto s = n_minus_1.num_trailing_zeros();
    natural d = n_minus_1;
    d >>= s;
    natural x, a;

    for (int i = 0; i < rounds; i++) {
        // sample uniformly from [2, n-2] and store in a
        uniform_sample(n_minus_3, rng, /*out*/a);
        a += 2;

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

constexpr std::vector<std::pair<uint64_t, int>> factorize(uint64_t a) {
    if (a <= 1)
        return {};

    std::vector<std::pair<uint64_t, int>> out;
    auto z = num_trailing_zeros(a);
    if (z) {
        out.emplace_back(2, z);
        a >>= z;
    }

    uint64_t p = 3;
    while (a > 1) {
        int count = 0;
        while (a % p == 0) {
            a /= p;
            count += 1;
        }
        if (count)
            out.emplace_back(p, count);
        p += 2;
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

}
