#include "algebra/rational_class.h"
#include "algebra/__test.h"

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

TEST_CASE("from floating point") {
    REQUIRE(rational(1.5) == rational(3, 2));
    REQUIRE(rational(-1.5) == rational(-3, 2));
    REQUIRE(rational(-1.5f) == rational(-3, 2));
    REQUIRE(rational(-6.0) == rational(-6));
    REQUIRE(rational(-0.25) == rational(-1, 4));
    REQUIRE(rational(0.0) == rational(0));
    REQUIRE(rational(-0.0) == rational(0));
    REQUIRE(rational(6.0) == rational(6));
    REQUIRE(rational(-1.0) == rational(-1));
    REQUIRE(rational(-0.1) == -rational(0.1));
}

TEST_CASE("format with fixed fraction digits") {
    REQUIRE(format("{:.3}", rational(1, 1000)) == "0.001");
    REQUIRE(format("{:.3}", rational(1, 100)) == "0.010");
    REQUIRE(format("{:.3}", rational(1, 2)) == "0.500");
    REQUIRE(format("{:.3}", rational(1, 8)) == "0.125");
    REQUIRE(format("{:.3}", rational(0)) == "0.000");
    REQUIRE(format("{:.5}", rational(1, 3)) == "0.33333");
    REQUIRE(format("{:.2}", rational(3, 4)) == "0.75");
    REQUIRE(format("{:.4}", rational(12345, 1000)) == "12.3450");
    REQUIRE(format("{:.2}", rational(12345, 1000)) == "12.35");
    REQUIRE(format("{:.3}", rational(-1, 1000)) == "-0.001");
    REQUIRE(format("{:.1}", rational(7)) == "7.0");
}

TEST_CASE("zero denominator throws") {
    auto q = [](int n, int d) { return rational(n, d); };
    REQUIRE_THROWS(q(1, 0));
    REQUIRE_THROWS(q(0, 0));
    REQUIRE_THROWS(q(-1, 0));
    REQUIRE_THROWS(rational(integer(1), integer(0)));
    REQUIRE_THROWS(rational("1/0"));

    REQUIRE(q(0, 5) == 0);
    REQUIRE(q(1, 2) == rational(1, 2));
    REQUIRE(q(-4, 8) == rational(-1, 2));
}

TEST_CASE("parse negative fraction below one") {
    REQUIRE(rational("-0.5") == rational(-1, 2));
    REQUIRE(rational("-0.001") == rational(-1, 1000));
    REQUIRE(rational("-0.25") == rational(-1, 4));
    REQUIRE(rational("-0.5e1") == rational(-5));
    REQUIRE(rational("-0.5e-1") == rational(-1, 20));
    REQUIRE(rational("0.5") == rational(1, 2));
    REQUIRE(rational("-0.0") == 0);
    REQUIRE(rational("-1.5") == rational(-3, 2));
    REQUIRE(rational("-0/5") == 0);
}

TEST_CASE("self assigning operators") {
    {
        rational a(1, 2);
        a /= a;
        REQUIRE(a == 1);
    }
    {
        rational a(-3, 7);
        a /= a;
        REQUIRE(a == 1);
    }
    {
        rational a(5);
        a /= a;
        REQUIRE(a == 1);
    }
    {
        rational z(0);
        REQUIRE_THROWS(z /= z);
    }
    {
        rational a(1, 2);
        a *= a;
        REQUIRE(a == rational(1, 4));
    }
    {
        rational a(1, 2);
        a += a;
        REQUIRE(a == 1);
    }
    {
        rational a(1, 2);
        a -= a;
        REQUIRE(a == 0);
    }
    {
        rational a(1, 2);
        a %= a;
        REQUIRE(a == 0);
    }
}

TEST_CASE("invert") {
    rational a(2, 3);
    a.invert();
    REQUIRE(a == rational(3, 2));

    rational b(-2, 3);
    b.invert();
    REQUIRE(b == rational(-3, 2));

    rational c(5);
    c.invert();
    REQUIRE(c == rational(1, 5));

    rational z(0);
    REQUIRE_THROWS(z.invert());
}
