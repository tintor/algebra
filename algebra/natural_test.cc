#include "algebra/natural.h"
#include "algebra/__test.h"
#include <catch2/benchmark/catch_benchmark.hpp>
#include <chrono>

TEST_CASE("extract_64bits") {
    std::mt19937_64 rng(233);
    for (int i = 0; i < 1000; i++) {
        int bits = std::uniform_int_distribution<int>(32, 256)(rng);
        natural x = uniform_sample_bits(bits, rng);

        REQUIRE(extract_u64(x, x.num_bits()) == 0);
        for (int j = 0; j < x.num_bits(); j++) {
            REQUIRE(extract_u64(x, j) == (x >> j).words[0]);
        }
    }
}

TEST_CASE("mod63_65") {
    REQUIRE(mod63_65(1_n) == std::pair{1, 1});
    REQUIRE(mod63_65(63_n) == std::pair{0, 63});
    REQUIRE(mod63_65(66_n) == std::pair{3, 1});
    REQUIRE(mod63_65(130_n) == std::pair{4, 0});
    REQUIRE(mod63_65(natural(UINT64_MAX)) == std::pair{15, 15});
    REQUIRE(mod63_65(natural(UINT128_MAX)) == std::pair{3, 60});

    std::mt19937_64 rng(234);
    for (int i = 0; i < 1000; i++) {
        int bits = std::uniform_int_distribution<int>(32, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        auto [m63, m65] = mod63_65(x);
        REQUIRE(m63 == x % 63);
        REQUIRE(m65 == x % 65);
    }
}

TEST_CASE("is_possible_square") {
    std::mt19937_64 rng(233);
    for (int i = 0; i < 1000; i++) {
        int bits = std::uniform_int_distribution<int>(32, 256)(rng);
        natural x = uniform_sample_bits(bits, rng);
        if (!is_possible_square(x)) {
            natural a = isqrt(x);
            REQUIRE(a * a != x);
        }
    }
}

constexpr natural diff(const natural& a, const natural& b) { return (a > b) ? a - b : (b - a); }

void verify_isqrt(const natural& x, const natural& a) {
    if (a * a > x || (a + 1) * (a + 1) <= x) {
        natural t = isqrt(x);
        natural e = diff(a, t);
        print("wrong result!\n");
        print("x       = {} ({})\n", x, x.num_bits());
        print("a       = {} ({})\n", a, a.num_bits());
        print("t       = {} ({})\n", t, t.num_bits());
        print("a^2     = {}\n", a * a);
        print("(a+1)^2 = {}\n", (a + 1) * (a + 1));
        print("abs_err = {} ({})\n", e, e.num_bits());
        REQUIRE(false);
    }
}

void test_isqrt(const auto& fn) {
    std::mt19937_64 rng(0);
    natural x;

    print("small\n");
    for (int i = 0; i < 1'000'000; i++) verify_isqrt(natural(i), fn(i));

#if 0
    print("64 bit\n");
    for (int i = 0; i < 100'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(30, 64)(rng);
        natural x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }

    print("96 bit\n");
    for (int i = 0; i < 100'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(65, 96)(rng);
        natural x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }

    print("128 bit\n");
    for (int i = 0; i < 100'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(97, 128)(rng);
        natural x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }
#endif
    natural e("4723193678752028155961467022770253813739277744022718882812203718746");
    verify_isqrt(e, fn(e));

    print("192 bit\n");
    for (int i = 0; i < 1'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(129, 192)(rng);
        natural x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }

    print("256 bit\n");
    for (int i = 0; i < 1'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(193, 256)(rng);
        natural x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 1'000'000 == 0) print("{}\n", i / 1'000'000);
    }

    print("512 bit\n");
    for (int i = 0; i < 100'000; i++) {
        int bits = std::uniform_int_distribution<int>(257, 512)(rng);
        natural x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 100'000 == 0) print("{}\n", i / 100'000);
    }

    print("1024 bit\n");
    for (int i = 0; i < 100'000; i++) {
        int bits = std::uniform_int_distribution<int>(513, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 100'000 == 0) print("{}\n", i / 10'000);
    }
}

constexpr natural isqrt_natural(const natural& x) { return isqrt(x); }
TEST_CASE("isqrt stress") { test_isqrt(isqrt_natural); }
TEST_CASE("isqrt2") { test_isqrt(isqrt2); }
TEST_CASE("isqrt3") { test_isqrt(isqrt3); }

TEST_CASE("modX()") {
    std::mt19937_64 rng(0);
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(1, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        REQUIRE(x % 2 == x.mod2());
        REQUIRE(x % 3 == x.mod3());
        REQUIRE(x % 4 == x.mod4());
        REQUIRE(x % 5 == x.mod5());
        REQUIRE(x % 6 == x.mod6());
        REQUIRE(x % 7 == x.mod7());
        REQUIRE(x % 8 == x.mod8());
        REQUIRE(x % 9 == x.mod9());
        REQUIRE(x % 10 == x.mod10());
    }
}

TEST_CASE("round_to_zero") {
    natural a;
    round_to_zero(10.1, a);
    REQUIRE(a == 10);
    round_to_zero(10.9, a);
    REQUIRE(a == 10);
    round_to_zero(-7.1, a);
    REQUIRE(a == 7);
    round_to_zero(-7.9, a);
    REQUIRE(a == 7);
}

TEST_CASE("iroot 2") {
    const auto a = 3213123_n;
    REQUIRE(isqrt(a) == iroot(a, 2));
    const auto b = 88888888888_n;
    REQUIRE(isqrt(b) == iroot(b, 2));
}

TEST_CASE("iroot 3 small") {
    std::mt19937_64 rng(3);
    int n = 3;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(4, 63)(rng);
        natural x = uniform_sample_bits(bits, rng);
        natural r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 3") {
    std::mt19937_64 rng(3);
    int n = 3;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        natural r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 4") {
    std::mt19937_64 rng(4);
    int n = 4;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        natural r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 5") {
    std::mt19937_64 rng(5);
    int n = 5;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        natural r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 6") {
    std::mt19937_64 rng(6);
    int n = 5;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        natural r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 7") {
    std::mt19937_64 rng(7);
    int n = 5;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        natural x = uniform_sample_bits(bits, rng);
        natural r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

#if 0
TEST_CASE("factorize") {
    natural a = pow(2_n, 128);
    int ms_max = 0;
    while (a > 1) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        auto factors = factorize(a);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        int ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

        if (ms >= 0) {
            ms_max = std::max(ms, ms_max);
            std::print("{} = ", a);
            for (int i = 0; i < factors.size(); i++) {
                if (i > 0)
                    std::print(" * ");
                std::print("{}", factors[i].first);
                if (factors[i].second > 1)
                    std::print("^{}", factors[i].second);
            }
            std::print(" in {} ms (max {} ms)\n", ms, ms_max);
        }

        natural m = 1;
        for (const auto& [factor, count] : factors) {
            if (!is_likely_prime(factor, 40)) {
                std::print("returned factor {} is not prime!\n", factor);
                REQUIRE(false);
            }
            for (int i = 0; i < count; i++)
                m *= factor;
        }
        REQUIRE(m == a);

        a -= 1;
    }
}
#endif

ulong doubleToLongBits(double a) {
    return *reinterpret_cast<const ulong*>(&a);
}

// UNTESTED
natural double_to_natural(double a) {
    ulong bits = doubleToLongBits(a);
    ulong e = ulong(1) << 52;
    natural b = (bits & e - 1) | e;
    b <<= ((bits >> 52) & 0x7ff) - 1075;
    return b;
}

// UNTESTED
natural fast_isqrt(const natural& x) {
    if (x == 0)
        return 0;
    double xd = static_cast<double>(x);
    natural val, q, s;
    if (xd < 2.1267e37) { // 2.12e37 largest here since sqrt(long.max*long.max) > long.max
        ulong vInt = (ulong)sqrt(xd);
        natural nInt = vInt;
        val = natural((vInt + static_cast<ulong>(x / nInt)) >> 1);
    } else if (xd < 4.3322e127) {
        val = double_to_natural(sqrt(xd));

        div(x, val, q, s);
        q += val;
        q >>= 1;
        std::swap(q, val); // val = ((x / val) + val) >> 1;

        if (xd > 2e63) {
            /// val = ((x / val) + val) >> 1;
            div(x, val, q, s);
            q += val;
            q >>= 1;
            val = q; // val = ((x / val) + val) >> 1;
        }
    } else { // handle large numbers over 4.3322e127
        uint xLen = x.num_bits();
        uint wantedPrecision = ((xLen + 1) / 2);
        uint xLenMod = xLen + (xLen & 1) + 1;

        //////// Do the first Sqrt on Hardware ////////
        ulong valLong = doubleToLongBits(std::sqrt(static_cast<ulong>(x >> (xLenMod - 63)))) & 0x1fffffffffffffL;
        if (valLong == 0)
            valLong = 1L << 53;

        //////// Classic Newton Iterations ////////
        val = valLong;
        val <<= 52;
        q = x >> (xLenMod - (3 * 53));
        Check(valLong != 0);
        q /= valLong;
        val += q;

        uint size = 106;
        for (; size < 256; size <<= 1) {
            q = x;
            q >>= xLenMod - (3 * size);
            Check(val != 0);
            div(q, val, q, s);
            val <<= size - 1;
            val += q;
        }

        if (xd > 4e254) { // 4e254 = 1<<845.77
            uint numOfNewtonSteps = 32 - std::countl_zero(wantedPrecision / size);

            ////// Apply Starting Size ////////
            uint wantedSize = (wantedPrecision >> numOfNewtonSteps) + 2;
            uint needToShiftBy = size - wantedSize;
            val >>= needToShiftBy;

            size = wantedSize;
            do {
                //////// Newton Plus Iteration ////////
                uint shiftX = xLenMod - (3 * size);
                mul(val, val, s);
                s <<= size - 1; // s = (val * val) << (size - 1);

                q = x;
                q >>= shiftX;
                q -= s; // q = (x >> shiftX) - s
                Check(val != 0);
                div(q, val, q, s);

                val <<= size;
                val += q;
                size *= 2;
            } while (size < wantedPrecision);
        }
        val >>= size - wantedPrecision;
    }

    // Detect a round ups. This function can be further optimized - see article.
    // For a ~7% speed bump the following line can be removed but round-ups will occur.
    mul(val, val, q);
    if (q > x)
        val -= 1;
    return val;
}

#if 0
TEST_CASE("fast_isqrt stress") {
    natural a = 1;
    while (true) {
        std::print("{} -> {}\n", a.num_bits(), fast_isqrt(a).num_bits());
        a <<= 1;
    }
}
#endif

TEST_CASE("uniform_sample") {
    natural a = (1_n << 128) - 1;
    std::mt19937_64 rng(0);
    for (int i = 0; i < 20; i++) {
        natural m = uniform_sample(0, a, rng);
        REQUIRE(m.words.size() <= a.words.size());
        REQUIRE(0 <= m);
        REQUIRE(m <= a);
        REQUIRE(0 <= m.words.size());
        REQUIRE(m.words.size() <= 2);
        while (m.words.size()) {
            div(m, 10ull, /*out*/m);
        }
    }
}

TEST_CASE("uniform_sample2") {
    natural a = 999'999'999'999'999'999'999'999'999'999'999'999_n;
    natural b = 500'000'000'000'000'000'000'000'000'000'000'000_n;
    std::mt19937_64 rng(0);
    natural sum;
    int n = 100000;
    for (int i = 0; i < n; i++) {
        natural m = uniform_sample(0, a, rng);
        REQUIRE(m.words.size() <= a.words.size());
        REQUIRE(0 <= m);
        REQUIRE(m <= a);
        REQUIRE(0 <= m.words.size());
        REQUIRE(m.words.size() <= 2);
        if (i < 20)
            print("{}\n", m);
        sum += m;
    }
    sum /= n;
    print("avg {}\n", sum);
    REQUIRE(diff(sum, b) <= b / 100);
}

TEST_CASE("pow") {
    REQUIRE(pow(natural(2), 3) == 8);
    REQUIRE(pow(natural(10), 30) == natural("1000000000000000000000000000000"));
}

TEST_CASE("gcd") {
    REQUIRE(std::countr_zero(0u) == 32);
    REQUIRE(std::countr_zero(1u) == 0);
    REQUIRE(std::countr_zero(2u) == 1);
    REQUIRE(std::countr_zero(8u) == 3);
    REQUIRE(std::countr_zero(24u) == 3);

    REQUIRE(gcd(5u, 5u) == 5u);
    REQUIRE(gcd(6u, 15u) == 3u);
    REQUIRE(gcd(7u, 3u) == 1u);

    REQUIRE(gcd(natural(5), natural(5)) == natural(5));
    REQUIRE(gcd(natural(6), natural(15)) == natural(3));
    REQUIRE(gcd(natural(7), natural(3)) == natural(1));

    REQUIRE(gcd(natural(5), 5) == 5);
    REQUIRE(gcd(1, natural(5)) == 1);
}

TEST_CASE("isqrt") {
    REQUIRE(isqrt(natural(0)) == 0);
    REQUIRE(isqrt(natural(1)) == 1);
    REQUIRE(isqrt(natural(2)) == 1);
    REQUIRE(isqrt(natural(3)) == 1);
    REQUIRE(isqrt(natural(4)) == 2);
    REQUIRE(isqrt(natural(5)) == 2);
    REQUIRE(isqrt(natural(9)) == 3);
    REQUIRE(isqrt(natural(9999)) == 99);
    REQUIRE(isqrt(natural(10000)) == 100);
}

TEST_CASE("is_prime") {
    REQUIRE(!is_prime(natural(0)));
    REQUIRE(!is_prime(natural(1)));
    REQUIRE(is_prime(natural(2)));
    REQUIRE(is_prime(natural(3)));
    REQUIRE(!is_prime(natural(4)));
    REQUIRE(is_prime(natural(5)));
    REQUIRE(!is_prime(natural(6)));
}

TEST_CASE("is_power_of_two") {
    REQUIRE(is_power_of_two(1_n));
    REQUIRE(is_power_of_two(2));
    REQUIRE(!is_power_of_two(3));
    REQUIRE(is_power_of_two(4));
    REQUIRE(!is_power_of_two(5));

    REQUIRE(natural(64).num_bits() == 7);
    REQUIRE(natural(64).num_trailing_zeros() == 6);
    natural a = natural(1) << 100;
    REQUIRE(a.words.size() == 2);
    REQUIRE(a.words[0] == 0);
    REQUIRE(a.words[1] == uint64_t(1) << 36);
    REQUIRE(a.num_bits() == 101);
    REQUIRE(a.num_trailing_zeros() == 100);
    REQUIRE(is_power_of_two(a));

    REQUIRE(is_power_of_two(1_n << 280));
}

TEST_CASE("is_power_of_three") {
    REQUIRE(!is_power_of_three(0));
    REQUIRE(is_power_of_three(1));
    REQUIRE(!is_power_of_three(2));
    REQUIRE(is_power_of_three(3));
    REQUIRE(!is_power_of_three(4));
    REQUIRE(is_power_of_three(9));
    REQUIRE(is_power_of_three(pow(3_n, 30)));
    REQUIRE(!is_power_of_three(pow(3_n, 30) - 1));
}

TEST_CASE("power_of_two") {
    REQUIRE(pow(2_n, 0) == 1_n);
    REQUIRE(pow(2_n, 3) == 8_n);
    REQUIRE(pow(2_n, 63) == 1_n << 63);
    REQUIRE(pow(2_n, 64) == 1_n << 64);
}

TEST_CASE("factorize(uint64_t)") {
    using f = std::vector<std::pair<uint64_t, int>>;
    REQUIRE(factorize(uint64_t(0)) == f{});
    REQUIRE(factorize(uint64_t(1)) == f{});
    REQUIRE(factorize(uint64_t(12)) == f{{2, 2}, {3, 1}});
    REQUIRE(factorize(uint64_t(16)) == f{{2, 4}});
    REQUIRE(factorize(uint64_t(13)) == f{{13, 1}});
    REQUIRE(factorize(uint64_t(30)) == f{{2, 1}, {3, 1}, {5, 1}});
    REQUIRE(factorize(uint64_t(49)) == f{{7, 2}});
}

TEST_CASE("factorize(natural)") {
    using f = std::vector<std::pair<natural, int>>;
    REQUIRE(factorize(0_n) == f{});
    REQUIRE(factorize(1_n) == f{});
    REQUIRE(factorize(12_n) == f{{2, 2}, {3, 1}});
    REQUIRE(factorize(16_n) == f{{2, 4}});
    REQUIRE(factorize(13_n) == f{{13, 1}});
    REQUIRE(factorize(30_n) == f{{2, 1}, {3, 1}, {5, 1}});
    REQUIRE(factorize(49_n) == f{{7, 2}});
    REQUIRE(factorize(340282366920938463463374607431768211453_n) == f{{11, 1}, {6949, 1}, {4451685225093714772084598273548427_n, 1}});
}

TEST_CASE("is_prime vs is_likely_prime") {
    for (uint64_t p = 50'000'000; p <= 50'100'000; p++) {
        bool m = is_likely_prime(natural(p), 10);
        bool e = is_prime(p);
        REQUIRE(m == e);
    }
}

TEST_CASE("merseinne primes vs is_likely_prime") {
    std::vector<int> mp = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279};
    for (int p = 2; p <= mp.back() + 1; p++) {
        natural a = pow(2_n, p) - 1;
        const bool actual = is_likely_prime(a, 40);
        const bool expected = (std::find(mp.begin(), mp.end(), p) != mp.end());
        if (actual != expected)
            print("p={} a.num_bits={} actual={} expected={}\n", p, a.num_bits(), actual, expected);
        CHECK(actual == expected);
    }
}

TEST_CASE("isqrt uint64") {
    REQUIRE(isqrt(static_cast<uint64_t>(0)) == 0);
    REQUIRE(isqrt(static_cast<uint64_t>(1)) == 1);
    REQUIRE(isqrt(static_cast<uint64_t>(3)) == 1);
    REQUIRE(isqrt(static_cast<uint64_t>(4)) == 2);
    REQUIRE(isqrt(static_cast<uint64_t>(8)) == 2);
    REQUIRE(isqrt(static_cast<uint64_t>(9)) == 3);
    REQUIRE(isqrt(UINT64_MAX) == 4294967295ull);
    REQUIRE(isqrt(UINT64_MAX - 1) == 4294967295ull);

    // large perfect squares and their neighbours
    for (uint64_t k : {3037000499ull, 4000000000ull, 4294967290ull, 4294967291ull, 4294967295ull}) {
        REQUIRE(isqrt(k * k) == k);
        REQUIRE(isqrt(k * k - 1) == k - 1);
        if (k * k + 1 != 0)
            REQUIRE(isqrt(k * k + 1) == k);
    }

    std::mt19937_64 rng(11);
    for (int i = 0; i < 20000; i++) {
        const uint64_t x = rng();
        const uint64_t q = isqrt(x);
        REQUIRE(static_cast<uint128_t>(q) * q <= x);
        REQUIRE(static_cast<uint128_t>(q + 1) * (q + 1) > x);
    }
}

TEST_CASE("__isqrt_u128") {
    REQUIRE(__isqrt_u128(0) == 0);
    REQUIRE(__isqrt_u128(1) == 1);
    REQUIRE(__isqrt_u128(3) == 1);
    REQUIRE(__isqrt_u128(4) == 2);
    REQUIRE(__isqrt_u128(8) == 2);
    REQUIRE(__isqrt_u128(UINT64_MAX) == 4294967295ull);
    REQUIRE(__isqrt_u128(static_cast<uint128_t>(UINT64_MAX) + 1) == 4294967296ull);

    const uint128_t one = 1;
    REQUIRE(__isqrt_u128((one << 64) + (one << 33)) == (one << 32));
    REQUIRE(__isqrt_u128((one << 64) + (one << 33) + 1) == (one << 32) + 1);
    REQUIRE(__isqrt_u128(static_cast<uint128_t>(UINT64_MAX) * UINT64_MAX) == UINT64_MAX);
    REQUIRE(__isqrt_u128(UINT128_MAX) == UINT64_MAX);

    std::mt19937_64 rng(2);
    for (int i = 0; i < 2000; i++) {
        const uint128_t x = (static_cast<uint128_t>(rng()) << 64) | rng();
        const uint64_t q = __isqrt_u128(x);
        REQUIRE(static_cast<uint128_t>(q) * q <= x);
        if (q != UINT64_MAX)
            REQUIRE(static_cast<uint128_t>(q + 1) * (q + 1) > x);
    }
    for (int i = 0; i < 2000; i++) {
        const uint128_t x = static_cast<uint128_t>(rng());
        const uint64_t q = __isqrt_u128(x);
        REQUIRE(static_cast<uint128_t>(q) * q <= x);
        REQUIRE(static_cast<uint128_t>(q + 1) * (q + 1) > x);
    }
}

TEST_CASE("try_fermat_factorize") {
    REQUIRE(try_fermat_factorize(4) == 2);
    REQUIRE(try_fermat_factorize(9) == 3);
    REQUIRE(try_fermat_factorize(25) == 5);
    REQUIRE(try_fermat_factorize(15) == 3);
    REQUIRE(try_fermat_factorize(5959) == 59); // 59 * 101
    REQUIRE(try_fermat_factorize(4294967291ull * 4294967279ull) == 4294967279ull);

    // whatever it returns has to be a proper divisor (or 0 for "not found")
    std::mt19937_64 rng(1);
    for (int i = 0; i < 300; i++) {
        const uint64_t p = std::uniform_int_distribution<uint64_t>(3, 4294967291ull)(rng) | 1;
        const uint64_t q = std::uniform_int_distribution<uint64_t>(3, 4294967291ull / p * 2 + 3)(rng) | 1;
        const uint64_t n = p * q;
        const uint64_t f = try_fermat_factorize(n);
        if (f) {
            REQUIRE(f > 1);
            REQUIRE(n % f == 0);
        }
    }
}

TEST_CASE("gcd 128 bit") {
    const uint128_t one = 1;
    REQUIRE(gcd(one << 100, one << 70) == (one << 70));
    REQUIRE(gcd(3 * (one << 70), one << 100) == (one << 70));
    REQUIRE(gcd(one << 64, one << 64) == (one << 64));
    REQUIRE(gcd(5 * (one << 64), 15 * (one << 64)) == 5 * (one << 64));
    REQUIRE(gcd(one << 127, 3 * (one << 65)) == (one << 65));
    REQUIRE(gcd((one << 100) + 1, one << 100) == 1);

    // odd 128 bit values (low word non-zero)
    const uint128_t p = (one << 89) - 1; // 618970019642690137449562111 = 618970019642690137449562111
    REQUIRE(gcd(p * 3, p * 5) == p);
    REQUIRE(gcd(p, p) == p);
}

TEST_CASE("isqrt of powers of two") {
    // squaring these needs a carry that travels through several all-ones words
    for (int e = 0; e < 900; e++) {
        natural a = 1;
        a <<= e;
        const natural s = isqrt(a);
        REQUIRE(s * s <= a);
        REQUIRE((s + 1u) * (s + 1u) > a);
    }
}

TEST_CASE("modular arithmetic") {
    std::mt19937_64 rng(3);
    for (int i = 0; i < 200; i++) {
        const int bits = std::uniform_int_distribution<int>(2, 300)(rng);
        natural m = uniform_sample_bits(bits, rng);
        if (m < 2u)
            continue;
        natural a = uniform_sample_bits(bits, rng) % m;
        const natural b = uniform_sample_bits(bits, rng) % m;

        natural s = a;
        add_mod(s, b, m);
        REQUIRE(s == (a + b) % m);

        natural d = a;
        sub_mod(d, b, m);
        REQUIRE(d == (a + m - b) % m);

        natural p;
        mul_mod(a, b, m, p);
        REQUIRE(p == (a * b) % m);

        natural q = a;
        __mul_mod(q, b, m);
        REQUIRE(q == (a * b) % m);
    }

    // pow_mod against repeated multiplication
    for (int i = 0; i < 40; i++) {
        const natural m = uniform_sample_bits(std::uniform_int_distribution<int>(2, 90)(rng), rng) + 2u;
        const natural a = uniform_sample_bits(60, rng) % m;
        const uint64_t e = std::uniform_int_distribution<uint64_t>(0, 40)(rng);
        natural expected = 1;
        for (uint64_t k = 0; k < e; k++)
            expected = (expected * a) % m;
        REQUIRE(pow_mod(a, natural(e), m) == expected);
        natural out;
        pow_mod(a, natural(e), m, out);
        REQUIRE(out == expected);
    }

    // the 128 bit helpers in util.h
    for (int i = 0; i < 2000; i++) {
        const uint128_t m = (static_cast<uint128_t>(rng()) << 64 | rng()) | 2;
        const uint128_t a = (static_cast<uint128_t>(rng()) << 64 | rng()) % m;
        const uint128_t b = (static_cast<uint128_t>(rng()) << 64 | rng()) % m;
        const natural bignum = natural(a) + natural(b);
        REQUIRE(natural(add_mod(a, b, m)) == bignum % natural(m));
        REQUIRE(natural(mul_mod(a, b, m)) == (natural(a) * natural(b)) % natural(m));
    }
    for (int i = 0; i < 500; i++) {
        const uint64_t m = rng() | 2;
        const uint64_t a = rng() % m;
        const uint64_t e = rng() % 64;
        natural expected = 1;
        for (uint64_t k = 0; k < e; k++)
            expected = (expected * a) % m;
        REQUIRE(pow_mod(a, e, m) == expected);
    }
}

TEST_CASE("number theory helpers") {
    // power_of_two
    for (int e = 0; e < 300; e++) {
        const natural p = power_of_two(e);
        REQUIRE(p == pow(natural(2), e));
        REQUIRE(p.num_bits() == e + 1);
    }

    // log_lower / log_upper
    REQUIRE(log_lower(natural(0), 10) == 0);
    REQUIRE(log_upper(natural(0), 10) == 0);
    for (uint64_t base : {2ull, 3ull, 10ull}) {
        for (int e = 0; e < 20; e++) {
            const natural p = pow(natural(base), e);
            REQUIRE(log_lower(p, base) == e);
            REQUIRE(log_upper(p, base) == e + 1);
            if (e > 0) {
                REQUIRE(log_lower(p - 1u, base) == e - 1);
                REQUIRE(log_upper(p - 1u, base) == e);
            }
        }
    }

    // binominal
    natural out;
    binominal(natural(5), 0, out);
    REQUIRE(out == 1u);
    binominal(natural(5), 2, out);
    REQUIRE(out == 10u);
    binominal(natural(10), 5, out);
    REQUIRE(out == 252u);
    binominal(natural(52), 5, out);
    REQUIRE(out == 2598960u);
    // symmetry: C(n, k) == C(n, n-k)
    natural x, y;
    binominal(natural(30), 12, x);
    binominal(natural(30), 18, y);
    REQUIRE(x == y);

    // iroot
    REQUIRE(iroot(natural(0), 3) == 0);
    REQUIRE(iroot(natural(1), 5) == 1);
    REQUIRE(iroot(natural(8), 3) == 2);
    REQUIRE(iroot(natural(9), 3) == 2);
    REQUIRE(iroot(natural(27), 3) == 3);
    for (uint32_t n : {3u, 4u, 5u, 7u}) {
        for (uint64_t v : {2ull, 5ull, 10ull, 1000ull}) {
            const natural p = pow(natural(v), n);
            REQUIRE(iroot(p, n) == v);
            REQUIRE(iroot(p - 1u, n) == v - 1);
            REQUIRE(iroot(p + 1u, n) == v);
        }
    }

    // exact_sqrt, both forms
    natural b;
    REQUIRE(exact_sqrt(natural(144), b));
    REQUIRE(b == 12u);
    REQUIRE(!exact_sqrt(natural(145), b));
    natural whole = 1, root = 1;
    exact_sqrt(natural(72), whole, root); // sqrt(72) == 6 * sqrt(2)
    REQUIRE(whole == 6u);
    REQUIRE(root == 2u);
    whole = 1;
    root = 1;
    exact_sqrt(natural(196), whole, root);
    REQUIRE(whole == 14u);
    REQUIRE(root == 1u);

    // is_possible_square never rejects a real square
    for (uint64_t v = 0; v < 3000; v++)
        REQUIRE(is_possible_square(natural(v) * natural(v)));

    // concat
    REQUIRE(concat(0, 0) == 0);
    REQUIRE(concat(1, 2) == ((static_cast<uint128_t>(1) << 64) | 2));

    // complement
    natural c = 1;
    c <<= 64;
    complement(c); // two's complement of 2**64 in 2 words
    REQUIRE(c == (natural(1) << 64));
}
