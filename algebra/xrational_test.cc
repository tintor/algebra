#include "algebra/xrational.h"
#include "algebra/__test.h"

TEST_CASE("xrational") {
    REQUIRE(format("{}", 2_x / 3_x) == "2/3");
    REQUIRE(format("{}", sqrt(1_x)) == "1");
    REQUIRE(format("{}", sqrt(2_x)) == "sqrt(2)");
    REQUIRE(format("{}", sqrt(4_x)) == "2");
    REQUIRE(format("{}", sqrt(8_x)) == "2*sqrt(2)");
    REQUIRE(format("{}", sqrt(4/9_x)) == "2/3");
    REQUIRE(format("{}", sqrt(8_x) / 3) == "2/3*sqrt(2)");
    REQUIRE(format("{}", 1 / sqrt(2_x)) == "sqrt(2)/2");
    REQUIRE(format("{}", sqr(1 / sqrt(2_x))) == "1/2");
    REQUIRE(format("{}", sqrt(2_x) + sqrt(8_x)) == "3*sqrt(2)");
    REQUIRE(format("{}", sqrt(6_x) * sqrt(10_x)) == "2*sqrt(15)");
}

TEST_CASE("cmp") {
    REQUIRE(1 < sqrt(2_x));
    REQUIRE(sqrt(2_x) < 2);
    REQUIRE(-3/416_x > -1_x);
    REQUIRE(-3/416_x > -1);
}
