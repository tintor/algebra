#include "algebra/natural_func.h"
#include "algebra/__test.h"
#include <catch2/benchmark/catch_benchmark.hpp>

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
        q /= valLong;
        val += q;

        uint size = 106;
        for (; size < 256; size <<= 1) {
            q = x;
            q >>= xLenMod - (3 * size);
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

natural diff(const natural& a, const natural& b) {
    return (a > b) ? a - b : (b - a);
}

TEST_CASE("__slow_isqrt stress") {
    long i = 0;
    std::mt19937_64 rng(0);
    int bits_max = 1024;
    while (i < 100000) {
        int bits = std::uniform_int_distribution<int>(1, bits_max)(rng);
        natural x = uniform_sample_bits(bits, rng);
        bits = x.num_bits();
        natural q = __slow_isqrt(x);
        if (q * q > x || (q + 1) * (q + 1) <= x) {
            natural qe = q;
            std::print("__slow_isqrt failed for bits={} x={} q={} qe={} diff={}\n", bits, x, q, qe, diff(q, qe));
            bits_max = bits - 1;
        }
        i += 1;
        if (i % 1'000'000 == 0)
            print("{} mil, bits_max {}\n", i / 1'000'000, bits_max);
    }
}

TEST_CASE("isqrt benchmark") {
    natural a;

#define BENCH_DIGITS(N) \
    a = pow(10_n, N); \
    print("1e" #N " bits {}\n", a.num_bits()); \
    BENCHMARK("isqrt() " #N) { return isqrt(a); }; \
    BENCHMARK("__slow_isqrt() " #N) { return __slow_isqrt(a); }; \
    BENCHMARK("fast_isqrt() " #N) { return fast_isqrt(a); }; \
    BENCHMARK("std::sqrt() " #N) { return sqrt(static_cast<double>(a)); };

    BENCH_DIGITS(10);
    BENCH_DIGITS(20);
    BENCH_DIGITS(30);
    BENCH_DIGITS(38);
    BENCH_DIGITS(40);
    BENCH_DIGITS(50);
    BENCH_DIGITS(75);
    BENCH_DIGITS(100);
    BENCH_DIGITS(125);
    BENCH_DIGITS(150);
    BENCH_DIGITS(300);
}

TEST_CASE("isqrt stress") {
    long i = 0;
    std::mt19937_64 rng(0);
    int bits_max = 126;
    while (i < 100000) {
        int bits = std::uniform_int_distribution<int>(1, bits_max)(rng);
        natural x = uniform_sample_bits(bits, rng);
        bits = x.num_bits();
        natural q = isqrt(x);
        if (q * q > x || (q + 1) * (q + 1) <= x) {
            natural qe = __slow_isqrt(x);
            std::print("isqrt failed for bits={} x={} q={} qe={} diff={}\n", bits, x, q, qe, diff(q, qe));
            bits_max = bits - 1;
        }
        i += 1;
        if (i % 1'000'000 == 0)
            print("{} mil, bits_max {}\n", i / 1'000'000, bits_max);
    }
}

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
    REQUIRE(is_power_of_two(1));
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

TEST_CASE("power_of_two") {
    REQUIRE(pow(2_n, 0) == 1_n);
    REQUIRE(pow(2_n, 3) == 8_n);
    REQUIRE(pow(2_n, 63) == 1_n << 63);
    REQUIRE(pow(2_n, 64) == 1_n << 64);
}

TEST_CASE("factorize") {
    using f = std::vector<std::pair<uint64_t, int>>;
    REQUIRE(factorize(0_n) == f{});
    REQUIRE(factorize(1_n) == f{});
    REQUIRE(factorize(12_n) == f{{2, 2}, {3, 1}});
    REQUIRE(factorize(16_n) == f{{2, 4}});
    REQUIRE(factorize(13_n) == f{{13, 1}});
    REQUIRE(factorize(30_n) == f{{2, 1}, {3, 1}, {5, 1}});
}

TEST_CASE("is_prime vs is_likely_prime") {
    std::mt19937_64 rng(0);
    for (uint64_t p = 50'000'000; p <= 50'100'000; p++) {
        bool m = is_likely_prime(natural(p), 10, rng);
        bool e = is_prime(p);
        REQUIRE(m == e);
    }
}

TEST_CASE("merseinne primes vs is_likely_prime") {
    std::mt19937_64 rng(0);
    std::vector<int> mp = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279};
    for (int p = 2; p <= mp.back() + 1; p++) {
        natural a = pow(2_n, p) - 1;
        bool actual = is_likely_prime(a, 10, rng);
        bool expected = (std::find(mp.begin(), mp.end(), p) != mp.end());
        REQUIRE(actual == expected);
    }
}
