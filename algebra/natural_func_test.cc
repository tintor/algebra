#include "algebra/natural_func.h"
#include "algebra/__test.h"

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

natural diff(const natural& a, const natural& b) {
    return (a >= b) ? (a - b) : (b - a);
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
    REQUIRE(num_trailing_zeros(0u) == 0);
    REQUIRE(num_trailing_zeros(1u) == 0);
    REQUIRE(num_trailing_zeros(2u) == 1);
    REQUIRE(num_trailing_zeros(8u) == 3);
    REQUIRE(num_trailing_zeros(24u) == 3);

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
