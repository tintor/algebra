#include "algebra/real.h"
#include "algebra/__test.h"

TEST_CASE("format") {
    REQUIRE(format("{}", real<2>(3, 0)) == "3");
    REQUIRE(format("{}", real<2>(-3, 0)) == "-3");
    REQUIRE(format("{}", real<2>(3, 2)) == "12");
    REQUIRE(format("{}", real<2>(-3, 2)) == "-12");
    REQUIRE(format("{}", real<2>(3, -2)) == "3/4");
    REQUIRE(format("{}", real<2>(-3, -2)) == "-3/4");
    REQUIRE(format("{:.2}", real<2>(3, -2)) == "0.75");
    REQUIRE(format("{:.3}", real<2>(3, -2)) == "0.750");
    REQUIRE(format("{:.2}", real<2>(-3, -2)) == "-0.75");
}

TEST_CASE("misc") {
    REQUIRE(real<2>(2, 0).num == 1);
    REQUIRE(real<2>(2, 0).exp == 1);
    REQUIRE(real<2>(1, 1).num == 1);
    REQUIRE(real<2>(1, 1).exp == 1);

    //REQUIRE(real<2>(2, 0) == real<2>(1, 1));

    REQUIRE(format("{:.1}", real<2>::round(rational(1, 10), 53)) == "0.1");
}

TEST_CASE("real<2> literal") {
    REQUIRE(std::format("{}", 1.5_f) == "3/2");
}

TEST_CASE("decimal literal") {
    REQUIRE(std::format("{}", decimal(13/10_q)) == "1.3");
    REQUIRE(std::format("{}", 1.3_d) == "1.3");
    REQUIRE(std::format("{}", -1.3_d) == "-1.3");
    REQUIRE(std::format("{}", -1200_d) == "-1200");
    REQUIRE((-1200_d).num == -12);
    REQUIRE((-1200_d).exp == 2);
    REQUIRE(std::format("{}", -0.0012_q) == "-3/2500");
    REQUIRE((-0.0012_d).num == -12);
    REQUIRE((-0.0012_d).exp == -4);
    REQUIRE(std::format("{}", -0.0012_d) == "-0.0012");
    REQUIRE(std::format("{}", -1.234_d) == "-1.234");
}
