#include "algebra/rational.h"
#include "algebra/__test.h"

TEST_CASE("fract") {
    REQUIRE(fract(1/7_q) == 1/7_q);
    REQUIRE(fract(rational(16, 7)) == rational(2, 7));
    REQUIRE(fract(rational(0)) == rational(0));
    REQUIRE(fract(rational(-1, 7)) == rational(1, 7));
    REQUIRE(fract(rational(-16, 7)) == rational(2, 7));
}

TEST_CASE("sqrt") {
    REQUIRE(double(sqrt(2_q, 4)) == sqrt(2));
}

TEST_CASE("pow") {
    REQUIRE(pow(1_q, 1) == 1);
    REQUIRE(pow(1_q, 0) == 1);
    REQUIRE(pow(1_q, -1) == 1);

    REQUIRE(pow(2_q, 0) == 1);
    REQUIRE(pow(2_q, 1) == 2);
    REQUIRE(pow(2_q, -1) == 1/2_q);
    REQUIRE(pow(2_q, 2) == 4);
    REQUIRE(pow(2_q, -2) == 1/4_q);

    auto a = pow(2_q, 1/2_q, 8); // sqrt(2)
    auto aa = pow(a, 2);
    REQUIRE(aa.num > 0);
    REQUIRE(aa.den > 0);
    REQUIRE(round(aa, 100) == 2);

    auto b = pow(2_q, 1/3_q, 9); // cbrt(2)
    auto bb = pow(b, 3);
    REQUIRE(round(bb, 100) == 2);

    auto c = pow(2_q, 2/3_q, 9);
    REQUIRE(c == b * b);

    auto d = pow(2_q, 4/3_q, 9);
    REQUIRE(d == 2 * b);

    auto e = pow(2_q, -1/2_q, 8); // 1/sqrt(2)
    REQUIRE(e * a == 1);
}

TEST_CASE("PI") {
    REQUIRE(approx_log2(PI(11) - rational(M_PI)) <= -52);
}

TEST_CASE("exp") {
    REQUIRE(exp(0_q, 20) == 1);
    REQUIRE(approx_log2(exp(1_q, 20) - rational(exp(1))) <= -52);
    REQUIRE(approx_log2(exp(1/10_q, 20) - rational(exp(0.1))) <= -52);
    REQUIRE(approx_log2(exp(2_q, 30) - rational(exp(2))) <= -52);
    REQUIRE(approx_log2(exp(-1_q, 30) - rational(exp(-1))) <= -52);
}

TEST_CASE("cos") {
    REQUIRE(cos(0_q, 5) == 1);
    REQUIRE(approx_log2(cos(1_q, 12) - rational(cos(1.0))) <= -50);
    REQUIRE(approx_log2(cos(2_q, 14) - rational(cos(2.0))) <= -50);
    REQUIRE(approx_log2(cos(1/2_q, 12) - rational(cos(0.5))) <= -50);
    REQUIRE(approx_log2(cos(-2_q, 14) - rational(cos(2.0))) <= -50);
}

TEST_CASE("sin") {
    REQUIRE(sin(0_q, 5) == 0);
    REQUIRE(approx_log2(sin(1_q, 12) - rational(sin(1.0))) <= -50);
    REQUIRE(approx_log2(sin(2_q, 14) - rational(sin(2.0))) <= -50);
    REQUIRE(approx_log2(sin(1/2_q, 12) - rational(sin(0.5))) <= -50);
}

TEST_CASE("__PI leaf term") {
    // p = -(6a-5)(2a-1)(6a-1), q = 10939058860032000 * a^3, r = p * (545140134*a + 13591409)
    for (unsigned a : {1u, 2u, 7u, 8u, 9u, 20u, 100u}) {
        const BinarySplit e = __PI(a, a + 1);
        const integer ai = a;
        REQUIRE(e.p == -(6 * ai - 5) * (2 * ai - 1) * (6 * ai - 1));
        REQUIRE(e.q == integer(10939058860032000ll) * (ai * ai * ai));
        REQUIRE(e.r == e.p * (integer(545140134) * ai + 13591409));
    }
}
