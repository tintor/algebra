#include "algebra/rational_func.h"
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
    REQUIRE(round(pow(a, 2), 100) == 2);

    auto b = pow(2_q, 1/3_q, 9); // cbrt(2)
    REQUIRE(round(pow(b, 3), 100) == 2);

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
    REQUIRE(approx_log2(exp(1_q, 20) - rational(exp(1))) <= -52);
    //REQUIRE(approx_log2(exp(1/10_q, 20) - rational(exp(0.1))) <= -52);
    //REQUIRE(approx_log2(exp(10_q, 20) - rational(exp(10))) <= -52);
}
