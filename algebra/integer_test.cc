#include "algebra/integer.h"
#include "algebra/__test.h"

TEST_CASE("uniform_sample") {
    integer a;
    a = pow(2_i, 128) - 1;
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
    uniform_sample_bits(bits, rng, a); // writes an integer now, and always a non negative one
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
    integer out;

    REQUIRE(inverse_mod(3_i, 11_i, out));
    REQUIRE(out == 4u); // 3 * 4 == 12 == 1 (mod 11)

    REQUIRE(inverse_mod(1_i, 7_i, out));
    REQUIRE(out == 1u);

    REQUIRE(inverse_mod(5_i, 7_i, out));
    REQUIRE(out == 3u); // 5 * 3 == 15 == 1 (mod 7)

    // no inverse when gcd(a, m) != 1
    REQUIRE(!inverse_mod(4_i, 8_i, out));
    REQUIRE(!inverse_mod(6_i, 9_i, out));

    // exhaustive small check
    for (uint64_t m = 2; m < 30; m++)
        for (uint64_t a = 1; a < m; a++) {
            const bool ok = inverse_mod(integer(a), integer(m), out);
            REQUIRE(ok == (gcd(a, m) == 1));
            if (ok)
                REQUIRE((integer(a) * out) % m == 1);
        }

    // multi-word
    integer big = 1;
    big <<= 130;
    big += 7u;
    REQUIRE(inverse_mod(3_i, big, out));
    REQUIRE((integer(3) * out) % big == 1);
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

TEST_CASE("pow with accumulator") {
    // pow(base, exp, result) == result * base**exp
    REQUIRE(pow(3_i, 0, 5_i) == 5);
    REQUIRE(pow(3_i, 1, 5_i) == 15);
    REQUIRE(pow(3_i, 2, 5_i) == 45);
    REQUIRE(pow(3_i, 3, 5_i) == 135);
    REQUIRE(pow(2_i, 0, 5_i) == 5);
    REQUIRE(pow(2_i, 3, 5_i) == 40);
    REQUIRE(pow(4_i, 2, 5_i) == 80);
    REQUIRE(pow(10_i, 3, 2_i) == 2000);
    REQUIRE(pow(1_i, 5, 7_i) == 7);
    REQUIRE_THROWS(pow(2_i, -1, 5_i));
    REQUIRE_THROWS(pow(3_i, -1, 5_i));
}

TEST_CASE("less_ab_c / less_a_bc / less_ab_cd") {
    // these are used by rational comparison and by the geometry predicates
    std::mt19937_64 rng(0);
    auto rnd = [&](int max_bits) {
        integer a;
        uniform_sample_bits(std::uniform_int_distribution<int>(0, max_bits)(rng), rng, a);
        if (std::uniform_int_distribution<int>(0, 1)(rng))
            a.negate();
        return a;
    };

    // exhaustive over small values, including zeros and both signs
    for (int a = -3; a <= 3; a++)
        for (int b = -3; b <= 3; b++)
            for (int c = -9; c <= 9; c++) {
                REQUIRE(less_ab_c(integer(a), integer(b), integer(c)) == (a * b < c));
                REQUIRE(less_a_bc(integer(c), integer(a), integer(b)) == (c < a * b));
            }

    for (int a = -3; a <= 3; a++)
        for (int b = -3; b <= 3; b++)
            for (int c = -3; c <= 3; c++)
                for (int d = -3; d <= 3; d++)
                    REQUIRE(less_ab_cd(integer(a), integer(b), integer(c), integer(d)) == (a * b < c * d));

    // random, multi word
    for (int i = 0; i < 400; i++) {
        const integer a = rnd(200), b = rnd(200), c = rnd(200), d = rnd(200);
        REQUIRE(less_ab_c(a, b, c) == (a * b < c));
        REQUIRE(less_a_bc(a, b, c) == (a < b * c));
        REQUIRE(less_ab_cd(a, b, c, d) == (a * b < c * d));
    }

    // products that are equal or one apart, where the bit length shortcuts are tightest
    for (int i = 0; i < 200; i++) {
        const integer a = rnd(100), b = rnd(100);
        const integer ab = a * b;
        REQUIRE(!less_ab_c(a, b, ab));
        REQUIRE(less_ab_c(a, b, ab + 1));
        REQUIRE(!less_ab_c(a, b, ab - 1));
        REQUIRE(!less_ab_cd(a, b, a, b));
        REQUIRE(less_a_bc(ab - 1, a, b) == (ab - 1 < ab));
    }
}

TEST_CASE("exp2") {
    REQUIRE(exp2(0) == 1);
    REQUIRE(exp2(1) == 2);
    REQUIRE(exp2(10) == 1024);
    REQUIRE(exp2(64) == pow(integer(2), 64));
    REQUIRE(exp2(200) == pow(integer(2), 200));
    REQUIRE_THROWS(exp2(-1));
}

TEST_CASE("pow(integer, integer)") {
    REQUIRE(pow(integer(3), 4_i) == 81);
    REQUIRE(pow(integer(-3), 3_i) == -27);
    REQUIRE(pow(integer(-3), 4_i) == 81);
    REQUIRE(pow(integer(2), 100_i) == exp2(100));
    REQUIRE(pow(integer(-2), 101_i) == -exp2(101));
    REQUIRE(pow(integer(7), 0_i) == 1);
    REQUIRE(pow(integer(0), 5_i) == 0);
}

TEST_CASE("binominal_mod") {
    // reference values: C(10,3) = 120, C(20,5) = 15504, C(6,4) = 15, C(10,5) = 252, C(8,3) = 56
    integer out;
    binominal_mod(10_i, 3, 7_i, out);
    REQUIRE(out == 120u % 7u);
    binominal_mod(20_i, 5, 101_i, out);
    REQUIRE(out == 15504u % 101u);
    binominal_mod(10_i, 0, 7_i, out);
    REQUIRE(out == 1u);

    // a modulus sharing a factor with some i+1 in [1, k] has no modular inverse for it, which
    // must not silently produce a wrong answer
    binominal_mod(6_i, 4, 4_i, out);
    REQUIRE(out == 15u % 4u);
    binominal_mod(6_i, 4, 6_i, out);
    REQUIRE(out == 15u % 6u);
    binominal_mod(10_i, 5, 10_i, out);
    REQUIRE(out == 252u % 10u);
    binominal_mod(8_i, 3, 9_i, out);
    REQUIRE(out == 56u % 9u);

    // m == 1 reduces everything to 0
    binominal_mod(10_i, 3, 1_i, out);
    REQUIRE(out == 0u);
    binominal_mod(10_i, 0, 1_i, out);
    REQUIRE(out == 0u);
}

// log_lower, log_upper and is_power_of_three moved to integer.h when they were converted to take
// an integer, so their tests move here with them.
TEST_CASE("log_lower / log_upper") {
    // log_lower / log_upper
    REQUIRE(log_lower(integer(0), 10) == 0);
    REQUIRE(log_upper(integer(0), 10) == 0);
    for (uint64_t base : {2ull, 3ull, 10ull}) {
        for (int e = 0; e < 20; e++) {
            const integer p = pow(integer(base), e);
            REQUIRE(log_lower(p, base) == e);
            REQUIRE(log_upper(p, base) == e + 1);
            if (e > 0) {
                REQUIRE(log_lower(p - 1u, base) == e - 1);
                REQUIRE(log_upper(p - 1u, base) == e);
            }
        }
    }
    // a negative input is rejected rather than silently treated as its magnitude
    REQUIRE_THROWS(log_lower(integer(-8), 2));
    REQUIRE_THROWS(log_upper(integer(-8), 2));
    // and the magnitude of a negative, taken explicitly, still works
    REQUIRE(log_lower(abs(integer(-8)), 2) == 3);

    // a base below two has no logarithm: dividing by one never reaches zero, which used to spin
    REQUIRE_THROWS(log_lower(100_i, 1));
    REQUIRE_THROWS(log_upper(100_i, 1));
    REQUIRE_THROWS(log_lower(100_i, 0));
    REQUIRE_THROWS(log_upper(100_i, 0));
    REQUIRE_THROWS(log_lower(0_i, 1));
    REQUIRE_THROWS(log_upper(1_i, 1));
}

TEST_CASE("is_power_of_three") {
    REQUIRE(!is_power_of_three(0));
    REQUIRE(is_power_of_three(1));
    REQUIRE(!is_power_of_three(2));
    REQUIRE(is_power_of_three(3));
    REQUIRE(!is_power_of_three(4));
    REQUIRE(is_power_of_three(9));
    REQUIRE(is_power_of_three(pow(3_i, 30)));
    REQUIRE(!is_power_of_three(pow(3_i, 30) - 1));
}

TEST_CASE("is_power_of_three more") {
    REQUIRE(!is_power_of_three(integer(0)));
    REQUIRE(is_power_of_three(integer(1)));
    REQUIRE(is_power_of_three(integer(3)));
    REQUIRE(is_power_of_three(integer(9)));
    REQUIRE(is_power_of_three(integer(81)));
    REQUIRE(is_power_of_three(pow(integer(3), 40)));
    REQUIRE(is_power_of_three(pow(integer(3), 64)));  // a perfect square
    REQUIRE(!is_power_of_three(integer(2)));
    REQUIRE(!is_power_of_three(integer(6)));
    REQUIRE(!is_power_of_three(integer(12)));
    REQUIRE(!is_power_of_three(pow(integer(3), 20) + 1u));
    REQUIRE(!is_power_of_three(pow(integer(3), 20) * 2u));
}

TEST_CASE("is_power_of_three rejects a negative") {
    REQUIRE_THROWS(is_power_of_three(integer(-3)));
    REQUIRE_THROWS(is_power_of_three(integer(-9)));
    REQUIRE(is_power_of_three(abs(integer(-9))));
}

// A uint64_t divisor above INT64_MAX used to bind to the int64_t overload and be read as a
// negative one, which left the quotient with the right magnitude and the wrong sign.
TEST_CASE("div by an unsigned divisor above INT64_MAX") {
    const uint64_t big = uint64_t(1) << 63;
    const integer a("100000000000000000000000");
    integer q;
    const uint64_t r = div(a, big, q);

    REQUIRE(q * integer(big) + integer(r) == a);
    REQUIRE(q > 0);
    REQUIRE(r < big);
    REQUIRE(q == a / big);

    // and with a negative dividend the quotient is negative, the remainder follows the dividend
    integer qn;
    const uint64_t rn = div(-a, big, qn);
    REQUIRE(qn == -q);
    REQUIRE(rn == r);
    REQUIRE(qn * integer(big) - integer(rn) == -a);

    // UINT64_MAX is the extreme case
    integer q2;
    const uint64_t r2 = div(a, UINT64_MAX, q2);
    REQUIRE(q2 > 0);
    REQUIRE(q2 * integer(UINT64_MAX) + integer(r2) == a);

    // a signed divisor still behaves as before
    integer q3;
    REQUIRE(div(integer(100), int64_t(7), q3) == 2);
    REQUIRE(q3 == 14);
    REQUIRE(div(integer(-100), int64_t(7), q3) == -2);
    REQUIRE(q3 == -14);
    REQUIRE(div(integer(100), int64_t(-7), q3) == 2);
    REQUIRE(q3 == -14);
}
