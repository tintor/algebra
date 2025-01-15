#include "algebra/real.h"
#include <catch2/catch_test_macros.hpp>
using std::format;

namespace algebra {

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
    REQUIRE(real<2>(2, 0).num == 2);
    REQUIRE(real<2>(2, 0).exp == 0);
    REQUIRE(real<2>(1, 1).num == 1);
    REQUIRE(real<2>(1, 1).exp == 1);

    //REQUIRE(real<2>(2, 0) == real<2>(1, 1));

    REQUIRE(format("{:.1}", real<2>(rational(1, 10))) == "0.1");
}

}
