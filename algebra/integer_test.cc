#include "algebra/integer.h"
#include "algebra/__test.h"

TEST_CASE("uniform_sample") {
    integer a;
    a.abs = pow(2_n, 128) - 1;
    std::mt19937_64 rng(0);
    for (int i = 0; i < 20; i++) {
        integer m = uniform_sample(0, a, rng);
        REQUIRE(m.sign() <= a.sign());
        REQUIRE(0 <= m);
        REQUIRE(m <= a);
        REQUIRE(0 <= m.sign());
        REQUIRE(m.sign() <= 2);
        while (m.sign())
            div(m, static_cast<long>(10), /*out*/m);
    }
}

integer random_integer(const int bits_max, std::mt19937_64& rng) {
    int bits = std::uniform_int_distribution<int>(0, bits_max)(rng);
    integer a;
    uniform_sample_bits(bits, rng, a.abs);
    if (std::uniform_int_distribution<int>(0, 1)(rng) == 0)
        a.negate();
    return a;
}

TEST_CASE("abs") {
    REQUIRE(abs(0_i) == 0);
    REQUIRE(abs(10_i) == 10);
    REQUIRE(abs(-3_i) == 3);
}

TEST_CASE("pow") {
    REQUIRE(pow(2_i, 3) == 8);
    REQUIRE(pow(10_i, 30) == 1000000000000000000000000000000_i);
}

TEST_CASE("add/sub_product") {
    int m = 64 * 2;
    int n = 64 * 2;
    std::mt19937_64 rng(904);
    for (int i = 0; i < 1000'000; i++) {
        integer a = random_integer(m, rng);
        integer b = random_integer(n, rng);
        integer c = random_integer(n, rng);
        int64_t d = std::uniform_int_distribution<int64_t>(INT64_MIN, INT64_MAX)(rng);

        integer e = a;
        add_product(e, b, c);
        REQUIRE(e == a + b * c);

        e = a;
        add_product(e, b, d);
        REQUIRE(e == a + b * d);

        e = a;
        sub_product(e, b, c);
        REQUIRE(e == a - b * c);

        e = a;
        sub_product(e, b, d);
        REQUIRE(e == a - b * d);
    }
}

TEST_CASE("inverse_mod") {
    natural out;

    REQUIRE(inverse_mod(3_n, 11_n, out));
    REQUIRE(out == 4u); // 3 * 4 == 12 == 1 (mod 11)

    REQUIRE(inverse_mod(1_n, 7_n, out));
    REQUIRE(out == 1u);

    REQUIRE(inverse_mod(5_n, 7_n, out));
    REQUIRE(out == 3u); // 5 * 3 == 15 == 1 (mod 7)

    // no inverse when gcd(a, m) != 1
    REQUIRE(!inverse_mod(4_n, 8_n, out));
    REQUIRE(!inverse_mod(6_n, 9_n, out));

    // exhaustive small check
    for (uint64_t m = 2; m < 30; m++)
        for (uint64_t a = 1; a < m; a++) {
            const bool ok = inverse_mod(natural(a), natural(m), out);
            REQUIRE(ok == (gcd(a, m) == 1));
            if (ok)
                REQUIRE((natural(a) * out) % m == 1);
        }

    // multi-word
    natural big = 1;
    big <<= 130;
    big += 7u;
    REQUIRE(inverse_mod(3_n, big, out));
    REQUIRE((3_n * out) % big == 1);
}

TEST_CASE("in place mod is euclidean") {
    auto m = [](long a, long b) { integer x = a; mod(x, integer(b)); return x; };
    REQUIRE(m(-10, 5) == 0);
    REQUIRE(m(-10, 3) == 2);
    REQUIRE(m(10, 3) == 1);
    REQUIRE(m(-9, 3) == 0);
    REQUIRE(m(9, 3) == 0);
    REQUIRE(m(0, 3) == 0);
    REQUIRE(m(-1, 3) == 2);
    REQUIRE(m(-10, -3) == 2);
    REQUIRE(m(10, -3) == 1);
}
