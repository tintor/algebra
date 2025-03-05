#include "algebra/dual.h"
#include "algebra/__test.h"

TEST_CASE("basic") {
    REQUIRE(dual<float>{1, 2}.real == 1);
    REQUIRE(dual<float>{1, 2}.dual == 2);

    REQUIRE(format("{}", dual<float>{1, 2}) == "1+2*eps");
    REQUIRE(format("{}", dual<float>{1, 0}) == "1");
    REQUIRE(format("{}", dual<float>{1, -2}) == "1-2*eps");
    REQUIRE(format("{}", dual<float>{0, 2}) == "2*eps");
    REQUIRE(format("{}", dual<float>{0, 1}) == "eps");
    REQUIRE(format("{}", dual<float>{0, 0}) == "0");
}
