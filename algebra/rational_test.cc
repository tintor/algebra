#include "algebra/rational.h"
#include "algebra/__test.h"

TEST_CASE("xrational") {
    REQUIRE(format("{}", xrational(2) / xrational(3)) == "2/3");
    REQUIRE(format("{}", sqrt(xrational(1))) == "1");
    REQUIRE(format("{}", sqrt(xrational(2))) == "sqrt(2)");
    REQUIRE(format("{}", sqrt(xrational(4))) == "2");
    REQUIRE(format("{}", sqrt(xrational(8))) == "2*sqrt(2)");
    REQUIRE(format("{}", sqrt(xrational(4/9_q))) == "2/3");
    REQUIRE(format("{}", sqrt(xrational(8)) / xrational(3)) == "2/3*sqrt(2)");
    REQUIRE(format("{}", xrational(1) / sqrt(xrational(2))) == "sqrt(2)/2");
    REQUIRE(format("{}", sqr(xrational(1) / sqrt(xrational(2)))) == "1/2");
    REQUIRE(xrational(1) < sqrt(xrational(2)));
    REQUIRE(sqrt(xrational(2)) < xrational(2));
}

TEST_CASE("simplify") {
    REQUIRE(rational(3, -4).num == -3);
    REQUIRE(rational(3, -4).den == 4);

    REQUIRE(rational(-3, -4).num == 3);
    REQUIRE(rational(-3, -4).den == 4);

    REQUIRE(rational(-3, 4).num == -3);
    REQUIRE(rational(-3, 4).den == 4);

    REQUIRE(rational(8, 4).num == 2);
    REQUIRE(rational(8, 4).den == 1);

    REQUIRE(rational(0, 4).num == 0);
    REQUIRE(rational(0, 4).den == 1);
}

TEST_CASE("parse") {
    REQUIRE(rational("0") == rational(0));
    REQUIRE(rational("123") == rational(123));
    REQUIRE(rational("-123") == rational(-123));
    REQUIRE(rational("-4/8") == rational(-1, 2));
    REQUIRE(rational("-1.25") == rational(-5, 4));
    REQUIRE(rational("-1e3") == rational(-1000));
    REQUIRE(rational("-1e-3") == rational(-1, 1000));
}

TEST_CASE("format") {
    REQUIRE(format("{}", rational(1, 3)) == "1/3");
    REQUIRE(format("{}", rational(5)) == "5");
    REQUIRE(format("{}", rational(50, 7)) == "50/7");
    REQUIRE(format("{}", rational(-2, 3)) == "-2/3");
}

TEST_CASE("literal") {
    rational a = -1/2_q;
    REQUIRE(format("{}", a) == "-1/2");
}

TEST_CASE("format frac") {
    REQUIRE(format("{:.3}", 1/3_q) == "0.333");
    REQUIRE(format("{:.3}", 2/3_q) == "0.667");
    REQUIRE(format("{:.3}", 1 - 1/1000000_q) == "1.000");
    REQUIRE(format("{:.3}", -(1_i - 1/1000000_q)) == "-1.000");
}

TEST_CASE("<") {
    REQUIRE(rational(1, 7) < rational(2, 7));
    REQUIRE(rational(1, 4) < rational(1, 3));
    REQUIRE(rational(1, 2) < integer(1));
    REQUIRE(integer(1) < rational(3, 2));
    REQUIRE(rational(1, 2) < 1);
}

constexpr std::string os() {
    std::ostringstream s;
    s << -2/3_q;
    return s.str();
}

TEST_CASE("constexpr ostream") {
    REQUIRE(os() == "-2/3");
}
