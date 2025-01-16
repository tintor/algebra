#include "algebra/rational.h"
#include <catch2/catch_test_macros.hpp>
using std::format;
using std::print;
using namespace algebra;
using namespace algebra::literals;

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

TEST_CASE("fract") {
    REQUIRE(fract(1/7_q) == 1/7_q);
    REQUIRE(fract(rational(16, 7)) == rational(2, 7));
    REQUIRE(fract(rational(0)) == rational(0));
    REQUIRE(fract(rational(-1, 7)) == rational(1, 7));
    REQUIRE(fract(rational(-16, 7)) == rational(2, 7));
}

TEST_CASE("<") {
    REQUIRE(rational(1, 7) < rational(2, 7));
    REQUIRE(rational(1, 4) < rational(1, 3));
    REQUIRE(rational(1, 2) < integer(1));
    REQUIRE(integer(1) < rational(3, 2));
    REQUIRE(rational(1, 2) < 1);
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

constexpr std::string os() {
    std::ostringstream s;
    s << -2/3_q;
    return s.str();
}

TEST_CASE("constexpr ostream") {
    REQUIRE(os() == "-2/3");
}
