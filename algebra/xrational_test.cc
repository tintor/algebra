#include "algebra/xrational.h"
#include "algebra/__test.h"

TEST_CASE("format") {
    REQUIRE(format("{}", xrational(0)) == "0");
    REQUIRE(format("{}", xrational(1)) == "1");
    REQUIRE(format("{}", xrational(1, 2)) == "sqrt(2)");
    REQUIRE(format("{}", xrational(1, 4)) == "2");
    REQUIRE(format("{}", xrational(2)) == "2");
    REQUIRE(format("{}", xrational(2, 2)) == "2*sqrt(2)");
    REQUIRE(format("{}", xrational(1/2_q, 2)) == "sqrt(2)/2");
    REQUIRE(format("{}", xrational(-1/2_q, 2)) == "-sqrt(2)/2");
    REQUIRE(format("{}", xrational(3/2_q, 2)) == "3/2*sqrt(2)");
}

TEST_CASE("sqrt") {
    REQUIRE(sqrt(1_x) == 1);
    REQUIRE(sqrt(2_x) == xrational(1, 2));
    REQUIRE(sqrt(4_x) == 2);
    REQUIRE(sqrt(8_x) == xrational(2, 2));
}

TEST_CASE("cmp") {
    REQUIRE(1 < sqrt(2_x));
    REQUIRE(sqrt(2_x) < 2);
    REQUIRE(-3/416_x > -1_x);
    REQUIRE(-3/416_x > -1);
}

TEST_CASE("add sub") {
    REQUIRE(sqrt(5_x) + 0_q == sqrt(5_x));
    REQUIRE(sqrt(5_x) - 0_q == sqrt(5_x));
    REQUIRE(0_q + sqrt(5_x) == sqrt(5_x));
    REQUIRE(0_q - sqrt(5_x) == -sqrt(5_x));

    REQUIRE(sqrt(2_x) + sqrt(8_x) == 3 * sqrt(2_x));
    REQUIRE(sqrt(6_x) * sqrt(10_x) == 2 * sqrt(15_x));

    REQUIRE(sqrt(4_x * 5) + sqrt(9_x * 5) == 5 * sqrt(5_x));
    REQUIRE(sqrt(4_x * 5) - sqrt(9_x * 5) == -sqrt(5_x));
}

TEST_CASE("div") {
    REQUIRE(2_x / 3_x == 2/3_q);
    REQUIRE(sqrt(4/9_x) == 2/3_q);
    REQUIRE(sqrt(8_x) / 3 == xrational(2/3_q, 2));
    REQUIRE(1 / sqrt(2_x) == xrational(1/2_q, 2));
    REQUIRE(sqr(1 / sqrt(2_x)) == 1/2_q);
}
