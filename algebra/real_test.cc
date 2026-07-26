#include "algebra/real.h"
#include "algebra/__test.h"

TEST_CASE("empty") {

}

TEST_CASE("pow with accumulator") {
    using R = real<2>;
    // pow(base, exp, result) == result * base**exp
    REQUIRE(to_rational(pow(R(3), 0, R(5))) == rational(5));
    REQUIRE(to_rational(pow(R(3), 1, R(5))) == rational(15));
    REQUIRE(to_rational(pow(R(3), 2, R(5))) == rational(45));
    REQUIRE(to_rational(pow(R(3), 3, R(5))) == rational(135));
    REQUIRE(to_rational(pow(R(2), 3, R(5))) == rational(40));
    REQUIRE(to_rational(pow(R(2), 0, R(5))) == rational(5));

    // without an accumulator
    REQUIRE(to_rational(pow(R(3), 0)) == rational(1));
    REQUIRE(to_rational(pow(R(3), 1)) == rational(3));
    REQUIRE(to_rational(pow(R(3), 3)) == rational(27));
    REQUIRE(to_rational(pow(R(2), 10)) == rational(1024));

    // negative exponent
    REQUIRE(to_rational(pow(R(2), -2)) == rational(1, 4));
    REQUIRE(to_rational(pow(R(2), -2, R(8))) == rational(2));
    REQUIRE(to_rational(pow(R(4), -1, R(8))) == rational(2));
}
