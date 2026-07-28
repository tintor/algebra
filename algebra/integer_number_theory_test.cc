#include "algebra/integer.h"
#include "algebra/__test.h"
#include <catch2/benchmark/catch_benchmark.hpp>
#include <chrono>

TEST_CASE("extract_64bits") {
    std::mt19937_64 rng(233);
    for (int i = 0; i < 1000; i++) {
        int bits = std::uniform_int_distribution<int>(32, 256)(rng);
        integer x = uniform_sample_bits(bits, rng);

        REQUIRE(extract_u64(x, x.num_bits()) == 0);
        for (int j = 0; j < x.num_bits(); j++) {
            REQUIRE(extract_u64(x, j) == (x >> j).words[0]);
        }
    }
}

TEST_CASE("mod63_65") {
    REQUIRE(mod63_65(1_i) == std::pair{1, 1});
    REQUIRE(mod63_65(63_i) == std::pair{0, 63});
    REQUIRE(mod63_65(66_i) == std::pair{3, 1});
    REQUIRE(mod63_65(130_i) == std::pair{4, 0});
    REQUIRE(mod63_65(integer(UINT64_MAX)) == std::pair{15, 15});
    REQUIRE(mod63_65(integer(UINT128_MAX)) == std::pair{3, 60});

    std::mt19937_64 rng(234);
    for (int i = 0; i < 1000; i++) {
        int bits = std::uniform_int_distribution<int>(32, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
        auto [m63, m65] = mod63_65(x);
        REQUIRE(m63 == x % 63);
        REQUIRE(m65 == x % 65);
    }
}

TEST_CASE("is_possible_square") {
    std::mt19937_64 rng(233);
    for (int i = 0; i < 1000; i++) {
        int bits = std::uniform_int_distribution<int>(32, 256)(rng);
        integer x = uniform_sample_bits(bits, rng);
        if (!is_possible_square(x)) {
            integer a = isqrt(x);
            REQUIRE(a * a != x);
        }
    }
}

constexpr integer diff(const integer& a, const integer& b) { return (a > b) ? a - b : (b - a); }

void verify_isqrt(const integer& x, const integer& a) {
    if (a * a > x || (a + 1) * (a + 1) <= x) {
        integer t = isqrt(x);
        integer e = diff(a, t);
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
    integer x;

    print("small\n");
    for (int i = 0; i < 1'000'000; i++) verify_isqrt(integer(i), fn(i));

#if 0
    print("64 bit\n");
    for (int i = 0; i < 100'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(30, 64)(rng);
        integer x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }

    print("96 bit\n");
    for (int i = 0; i < 100'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(65, 96)(rng);
        integer x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }

    print("128 bit\n");
    for (int i = 0; i < 100'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(97, 128)(rng);
        integer x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }
#endif
    integer e("4723193678752028155961467022770253813739277744022718882812203718746");
    verify_isqrt(e, fn(e));

    print("192 bit\n");
    for (int i = 0; i < 1'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(129, 192)(rng);
        integer x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 10'000'000 == 0) print("{}\n", i / 10'000'000);
    }

    print("256 bit\n");
    for (int i = 0; i < 1'000'000; i++) {
        int bits = std::uniform_int_distribution<int>(193, 256)(rng);
        integer x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 1'000'000 == 0) print("{}\n", i / 1'000'000);
    }

    print("512 bit\n");
    for (int i = 0; i < 100'000; i++) {
        int bits = std::uniform_int_distribution<int>(257, 512)(rng);
        integer x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 100'000 == 0) print("{}\n", i / 100'000);
    }

    print("1024 bit\n");
    for (int i = 0; i < 100'000; i++) {
        int bits = std::uniform_int_distribution<int>(513, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
        verify_isqrt(x, fn(x));
        if (i % 100'000 == 0) print("{}\n", i / 10'000);
    }
}

constexpr integer isqrt_integer(const integer& x) { return isqrt(x); }
TEST_CASE("isqrt stress") { test_isqrt(isqrt_integer); }
TEST_CASE("isqrt2") { test_isqrt(isqrt2); }
TEST_CASE("isqrt3") { test_isqrt(isqrt3); }

TEST_CASE("modX()") {
    std::mt19937_64 rng(0);
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(1, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
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
    integer a;
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
    const auto a = 3213123_i;
    REQUIRE(isqrt(a) == iroot(a, 2));
    const auto b = 88888888888_i;
    REQUIRE(isqrt(b) == iroot(b, 2));
}

TEST_CASE("iroot 3 small") {
    std::mt19937_64 rng(3);
    int n = 3;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(4, 63)(rng);
        integer x = uniform_sample_bits(bits, rng);
        integer r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 3") {
    std::mt19937_64 rng(3);
    int n = 3;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
        integer r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 4") {
    std::mt19937_64 rng(4);
    int n = 4;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
        integer r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 5") {
    std::mt19937_64 rng(5);
    int n = 5;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
        integer r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 6") {
    std::mt19937_64 rng(6);
    int n = 5;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
        integer r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

TEST_CASE("iroot 7") {
    std::mt19937_64 rng(7);
    int n = 5;
    for (int i = 0; i < 10'000; i++) {
        int bits = std::uniform_int_distribution<int>(64, 1024)(rng);
        integer x = uniform_sample_bits(bits, rng);
        integer r = iroot(x, n);
        REQUIRE(pow(r, n) <= x);
        REQUIRE(pow(r + 1, n) > x);
    }
}

#if 0
TEST_CASE("factorize") {
    integer a = pow(2_i, 128);
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

        integer m = 1;
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
integer double_to_integer(double a) {
    ulong bits = doubleToLongBits(a);
    ulong e = ulong(1) << 52;
    integer b = (bits & e - 1) | e;
    b <<= ((bits >> 52) & 0x7ff) - 1075;
    return b;
}

// UNTESTED
integer fast_isqrt(const integer& x) {
    if (x == 0)
        return 0;
    double xd = static_cast<double>(x);
    integer val, q, s;
    if (xd < 2.1267e37) { // 2.12e37 largest here since sqrt(long.max*long.max) > long.max
        ulong vInt = (ulong)sqrt(xd);
        integer nInt = vInt;
        val = integer((vInt + static_cast<ulong>(x / nInt)) >> 1);
    } else if (xd < 4.3322e127) {
        val = double_to_integer(sqrt(xd));

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
    integer a = 1;
    while (true) {
        std::print("{} -> {}\n", a.num_bits(), fast_isqrt(a).num_bits());
        a <<= 1;
    }
}
#endif

TEST_CASE("uniform_sample") {
    integer a = (1_i << 128) - 1;
    std::mt19937_64 rng(0);
    for (int i = 0; i < 20; i++) {
        integer m = uniform_sample(0, a, rng);
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
    integer a = 999'999'999'999'999'999'999'999'999'999'999'999_i;
    integer b = 500'000'000'000'000'000'000'000'000'000'000'000_i;
    std::mt19937_64 rng(0);
    integer sum;
    int n = 100000;
    for (int i = 0; i < n; i++) {
        integer m = uniform_sample(0, a, rng);
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
    REQUIRE(pow(integer(2), 3) == 8);
    REQUIRE(pow(integer(10), 30) == integer("1000000000000000000000000000000"));
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

    REQUIRE(gcd(integer(5), integer(5)) == integer(5));
    REQUIRE(gcd(integer(6), integer(15)) == integer(3));
    REQUIRE(gcd(integer(7), integer(3)) == integer(1));

    REQUIRE(gcd(integer(5), 5) == 5);
    REQUIRE(gcd(1, integer(5)) == 1);
}

TEST_CASE("isqrt") {
    REQUIRE(isqrt(integer(0)) == 0);
    REQUIRE(isqrt(integer(1)) == 1);
    REQUIRE(isqrt(integer(2)) == 1);
    REQUIRE(isqrt(integer(3)) == 1);
    REQUIRE(isqrt(integer(4)) == 2);
    REQUIRE(isqrt(integer(5)) == 2);
    REQUIRE(isqrt(integer(9)) == 3);
    REQUIRE(isqrt(integer(9999)) == 99);
    REQUIRE(isqrt(integer(10000)) == 100);
}

TEST_CASE("is_prime") {
    REQUIRE(!is_prime(integer(0)));
    REQUIRE(!is_prime(integer(1)));
    REQUIRE(is_prime(integer(2)));
    REQUIRE(is_prime(integer(3)));
    REQUIRE(!is_prime(integer(4)));
    REQUIRE(is_prime(integer(5)));
    REQUIRE(!is_prime(integer(6)));
}

TEST_CASE("is_power_of_two") {
    REQUIRE(is_power_of_two(1_i));
    REQUIRE(is_power_of_two(2));
    REQUIRE(!is_power_of_two(3));
    REQUIRE(is_power_of_two(4));
    REQUIRE(!is_power_of_two(5));

    REQUIRE(integer(64).num_bits() == 7);
    REQUIRE(integer(64).num_trailing_zeros() == 6);
    integer a = integer(1) << 100;
    REQUIRE(a.words.size() == 2);
    REQUIRE(a.words[0] == 0);
    REQUIRE(a.words[1] == uint64_t(1) << 36);
    REQUIRE(a.num_bits() == 101);
    REQUIRE(a.num_trailing_zeros() == 100);
    REQUIRE(is_power_of_two(a));

    REQUIRE(is_power_of_two(1_i << 280));
}


TEST_CASE("power_of_two") {
    REQUIRE(pow(2_i, 0) == 1_i);
    REQUIRE(pow(2_i, 3) == 8_i);
    REQUIRE(pow(2_i, 63) == 1_i << 63);
    REQUIRE(pow(2_i, 64) == 1_i << 64);
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

TEST_CASE("factorize(integer)") {
    using f = std::vector<std::pair<integer, int>>; // factorize() returns integers now
    REQUIRE(factorize(0_i) == f{});
    REQUIRE(factorize(1_i) == f{});
    REQUIRE(factorize(12_i) == f{{2, 2}, {3, 1}});
    REQUIRE(factorize(16_i) == f{{2, 4}});
    REQUIRE(factorize(13_i) == f{{13, 1}});
    REQUIRE(factorize(30_i) == f{{2, 1}, {3, 1}, {5, 1}});
    REQUIRE(factorize(49_i) == f{{7, 2}});
    REQUIRE(factorize(340282366920938463463374607431768211453_i) == f{{11, 1}, {6949, 1}, {4451685225093714772084598273548427_i, 1}});
}

TEST_CASE("is_prime vs is_likely_prime") {
    for (uint64_t p = 50'000'000; p <= 50'100'000; p++) {
        bool m = is_likely_prime(integer(p), 10);
        bool e = is_prime(p);
        REQUIRE(m == e);
    }
}

TEST_CASE("merseinne primes vs is_likely_prime") {
    std::vector<int> mp = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279};
    for (int p = 2; p <= mp.back() + 1; p++) {
        integer a = pow(integer(2), p) - 1; // pow() returns an integer now
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
        integer a = 1;
        a <<= e;
        const integer s = isqrt(a);
        REQUIRE(s * s <= a);
        REQUIRE((s + 1u) * (s + 1u) > a);
    }
}

TEST_CASE("modular arithmetic") {
    std::mt19937_64 rng(3);
    for (int i = 0; i < 200; i++) {
        const int bits = std::uniform_int_distribution<int>(2, 300)(rng);
        integer m = uniform_sample_bits(bits, rng);
        if (m < 2u)
            continue;
        integer a = integer(uniform_sample_bits(bits, rng)) % m;
        const integer b = integer(uniform_sample_bits(bits, rng)) % m;

        integer s = a;
        add_mod(s, b, m);
        REQUIRE(s == (a + b) % m);

        integer d = a;
        sub_mod(d, b, m);
        REQUIRE(d == (a + m - b) % m);

        integer p;
        mul_mod(a, b, m, p);
        REQUIRE(p == (a * b) % m);

        integer q = a;
        __mul_mod(q, b, m);
        REQUIRE(q == (a * b) % m);
    }

    // pow_mod against repeated multiplication
    for (int i = 0; i < 40; i++) {
        const integer m = uniform_sample_bits(std::uniform_int_distribution<int>(2, 90)(rng), rng) + 2u;
        const integer a = integer(uniform_sample_bits(60, rng)) % m;
        const uint64_t e = std::uniform_int_distribution<uint64_t>(0, 40)(rng);
        integer expected = 1;
        for (uint64_t k = 0; k < e; k++)
            expected = (expected * a) % m;
        REQUIRE(pow_mod(a, integer(e), m) == expected);
        integer out;
        pow_mod(a, integer(e), m, out);
        REQUIRE(out == expected);
    }

    // the 128 bit helpers in util.h
    for (int i = 0; i < 2000; i++) {
        const uint128_t m = (static_cast<uint128_t>(rng()) << 64 | rng()) | 2;
        const uint128_t a = (static_cast<uint128_t>(rng()) << 64 | rng()) % m;
        const uint128_t b = (static_cast<uint128_t>(rng()) << 64 | rng()) % m;
        const integer bignum = integer(a) + integer(b);
        REQUIRE(integer(add_mod(a, b, m)) == bignum % integer(m));
        REQUIRE(integer(mul_mod(a, b, m)) == (integer(a) * integer(b)) % integer(m));
    }
    for (int i = 0; i < 500; i++) {
        const uint64_t m = rng() | 2;
        const uint64_t a = rng() % m;
        const uint64_t e = rng() % 64;
        integer expected = 1;
        for (uint64_t k = 0; k < e; k++)
            expected = (expected * a) % m;
        REQUIRE(pow_mod(a, e, m) == expected);
    }
}

TEST_CASE("number theory helpers") {
    // power_of_two
    for (int e = 0; e < 300; e++) {
        const integer p = power_of_two(e);
        REQUIRE(p == pow(integer(2), e));
        REQUIRE(p.num_bits() == e + 1);
    }

    // binominal
    integer out;
    binominal(integer(5), 0, out);
    REQUIRE(out == 1u);
    binominal(integer(5), 2, out);
    REQUIRE(out == 10u);
    binominal(integer(10), 5, out);
    REQUIRE(out == 252u);
    binominal(integer(52), 5, out);
    REQUIRE(out == 2598960u);
    // symmetry: C(n, k) == C(n, n-k)
    integer x, y;
    binominal(integer(30), 12, x);
    binominal(integer(30), 18, y);
    REQUIRE(x == y);

    // iroot
    REQUIRE(iroot(integer(0), 3) == 0);
    REQUIRE(iroot(integer(1), 5) == 1);
    REQUIRE(iroot(integer(8), 3) == 2);
    REQUIRE(iroot(integer(9), 3) == 2);
    REQUIRE(iroot(integer(27), 3) == 3);
    for (uint32_t n : {3u, 4u, 5u, 7u}) {
        for (uint64_t v : {2ull, 5ull, 10ull, 1000ull}) {
            const integer p = pow(integer(v), n);
            REQUIRE(iroot(p, n) == v);
            REQUIRE(iroot(p - 1u, n) == v - 1);
            REQUIRE(iroot(p + 1u, n) == v);
        }
    }

    // exact_sqrt, both forms
    integer b;
    REQUIRE(exact_sqrt(integer(144), b));
    REQUIRE(b == 12u);
    REQUIRE(!exact_sqrt(integer(145), b));
    integer whole = 1, root = 1;
    exact_sqrt(integer(72), whole, root); // sqrt(72) == 6 * sqrt(2)
    REQUIRE(whole == 6u);
    REQUIRE(root == 2u);
    whole = 1;
    root = 1;
    exact_sqrt(integer(196), whole, root);
    REQUIRE(whole == 14u);
    REQUIRE(root == 1u);

    // is_possible_square never rejects a real square
    for (uint64_t v = 0; v < 3000; v++)
        REQUIRE(is_possible_square(integer(v) * integer(v)));

    // concat
    REQUIRE(concat(0, 0) == 0);
    REQUIRE(concat(1, 2) == ((static_cast<uint128_t>(1) << 64) | 2));

    // complement: two's complement over the words the value occupies
    integer c = 1;
    c <<= 64;
    complement(c);
    REQUIRE(c == (integer(1) << 128) - (integer(1) << 64));
    integer c2 = 5;
    complement(c2);
    REQUIRE(c2 == integer(UINT64_MAX) - 4u);
}



TEST_CASE("lcm") {
    REQUIRE(lcm(integer(0), integer(0)) == 0);
    REQUIRE(lcm(integer(0), integer(5)) == 0);
    REQUIRE(lcm(integer(5), integer(0)) == 0);
    REQUIRE(lcm(integer(4), integer(6)) == 12);
    REQUIRE(lcm(integer(21), integer(6)) == 42);
    REQUIRE(lcm(integer(7), integer(7)) == 7);
    REQUIRE(lcm(integer(1), integer(13)) == 13);

    std::mt19937_64 rng(0);
    for (int i = 0; i < 30; i++) {
        const integer a = uniform_sample_bits(std::uniform_int_distribution<int>(1, 300)(rng), rng);
        const integer b = uniform_sample_bits(std::uniform_int_distribution<int>(1, 300)(rng), rng);
        if (!a || !b)
            continue;
        const integer l = lcm(a, b);
        REQUIRE(l % a == 0);
        REQUIRE(l % b == 0);
        REQUIRE(l * gcd(a, b) == a * b);
    }
}

TEST_CASE("invert_bits") {
    integer a = UINT64_MAX;
    invert_bits(a);
    REQUIRE(a.words.size() == 0);
    REQUIRE(a == 0u);

    integer b = 0xF0F0u;
    integer c = b;
    invert_bits(c);
    REQUIRE(c == (UINT64_MAX ^ 0xF0F0u));
    // operator~ on an integer is the arithmetic complement, which is a different operation
    REQUIRE(~b == -integer(0xF0F1u));

    // a negative input has no bit pattern to flip
    integer d = -1;
    REQUIRE_THROWS(invert_bits(d));
}

TEST_CASE("isqrt2 and isqrt3 agree with isqrt") {
    // three independent implementations, so they are each other's reference
    for (int bits = 1; bits <= 200; bits++) {
        const integer a = power_of_two(bits);
        for (const integer& x : {a - 1u, a, a + 1u}) {
            const integer s = isqrt(x);
            REQUIRE(s * s <= x);
            REQUIRE((s + 1u) * (s + 1u) > x);
            REQUIRE(isqrt2(x) == s);
            REQUIRE(isqrt3(x) == s);
        }
    }
    std::mt19937_64 rng(7);
    for (int i = 0; i < 200; i++) {
        const integer x = uniform_sample_bits(std::uniform_int_distribution<int>(1, 300)(rng), rng);
        const integer s = isqrt(x);
        REQUIRE(isqrt2(x) == s);
        REQUIRE(isqrt3(x) == s);
    }
    for (unsigned i = 0; i <= 20; i++) {
        REQUIRE(isqrt2(integer(i)) == isqrt(integer(i)));
        REQUIRE(isqrt3(integer(i)) == isqrt(integer(i)));
    }
}

TEST_CASE("isqrt_hardware") {
    // documented as very fast but only approximate for large arguments, so exact only for small
    for (unsigned i = 0; i <= 100; i++)
        REQUIRE(isqrt_hardware(integer(i)) == isqrt(integer(i)));
    // for large arguments it must still land close to the true root
    for (int bits : {80, 130, 260, 500}) {
        const integer a = power_of_two(bits);
        const integer s = isqrt(a);
        const integer h = isqrt_hardware(a);
        const integer lo = s - s / integer(1000000u) - 1u;
        const integer hi = s + s / integer(1000000u) + 1u;
        REQUIRE(h >= lo);
        REQUIRE(h <= hi);
    }
}

TEST_CASE("is_one_of") {
    REQUIRE(is_one_of(4, {0, 1, 4, 9}));
    REQUIRE(is_one_of(0, {0, 1, 4, 9}));
    REQUIRE(is_one_of(9, {0, 1, 4, 9}));
    REQUIRE(!is_one_of(5, {0, 1, 4, 9}));
    REQUIRE(!is_one_of(4, {}));
}

TEST_CASE("is_possible_square - every square passes the filter") {
    // a filter: every real square must pass, non-squares may or may not
    for (unsigned i = 0; i <= 300; i++)
        REQUIRE(is_possible_square(integer(i) * integer(i)));
    for (int bits : {40, 100, 250}) {
        const integer a = power_of_two(bits);
        REQUIRE(is_possible_square(a * a));
    }
    // zero is a square, including a zero that came out of a subtraction
    REQUIRE(is_possible_square(0_i));
    integer z = 5;
    z -= 5u;
    REQUIRE(is_possible_square(z));
    // and it does reject: 2, 3, 5, 6 are not squares
    REQUIRE(!is_possible_square(2_i));
    REQUIRE(!is_possible_square(3_i));
    REQUIRE(!is_possible_square(5_i));
    REQUIRE(!is_possible_square(6_i));
}
