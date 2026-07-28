#include "algebra/xrational.h"
#include "algebra/__test.h"
#include <sstream>
#include <unordered_map>
#include <vector>

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

TEST_CASE("compound assignment with rational") {
    {
        xrational a = sqrt(3_x);
        a /= 2;
        REQUIRE(a == xrational(1/2_q, 3));
    }
    {
        xrational a = sqrt(5_x);
        a /= 3_q;
        REQUIRE(a == xrational(1/3_q, 5));
        a *= 3_q;
        REQUIRE(a == sqrt(5_x));
    }
    {
        xrational a = sqrt(5_x);
        a /= integer(2);
        REQUIRE(a == xrational(1/2_q, 5));
    }
    {
        xrational a = 6_x;
        a /= 3;
        REQUIRE(a == 2_x);
    }
}

TEST_CASE("equality with different roots") {
    // (2/3)*sqrt(245) == (14/3)*sqrt(5), since sqrt(245) == 7*sqrt(5)
    const xrational a(2/3_q, natural(245));
    const xrational b(14/3_q, natural(5));
    REQUIRE(sqr(a) == sqr(b));
    REQUIRE(a == b);
    REQUIRE(b == a);

    const xrational c(15/3_q, natural(5));
    REQUIRE(!(a == c));
    REQUIRE(!(c == a));

    const xrational d(-2/3_q, natural(245));
    const xrational e(-14/3_q, natural(5));
    REQUIRE(d == e);
    REQUIRE(!(a == e));
}

TEST_CASE("simplify") {
    // xrational(1, r) is sqrt(r): square factors may move into base, but the value must not change
    for (int r : {1, 2, 3, 4, 5, 8, 9, 12, 16, 18, 24, 27, 32, 36, 45, 48, 50, 72, 81, 90, 98, 128, 200, 245}) {
        const xrational a(1, natural(r));
        REQUIRE(a.root > 0);
        REQUIRE(sqr(a) == r);
        REQUIRE(a > 0);
        REQUIRE(a * a == r);
    }

    REQUIRE(xrational(1, natural(4)) == 2);
    REQUIRE(xrational(1, natural(8)) == 2 * sqrt(2_x));
    REQUIRE(xrational(1, natural(16)) == 4);
    REQUIRE(xrational(1, natural(32)) == 4 * sqrt(2_x));
    REQUIRE(xrational(1, natural(128)) == 8 * sqrt(2_x));
    REQUIRE(xrational(1, natural(9)) == 3);
    REQUIRE(xrational(1, natural(9)).is_rational());
    REQUIRE(xrational(1, natural(18)) == 3 * sqrt(2_x));
    REQUIRE(xrational(1, natural(27)) == 3 * sqrt(3_x));
    REQUIRE(xrational(1, natural(36)) == 6);
    REQUIRE(xrational(1, natural(45)) == 3 * sqrt(5_x));
    REQUIRE(xrational(1, natural(72)) == 6 * sqrt(2_x));
    REQUIRE(xrational(1, natural(81)) == 9);
    REQUIRE(xrational(1, natural(90)) == 3 * sqrt(10_x));
    REQUIRE(xrational(1, natural(200)) == 10 * sqrt(2_x));
    REQUIRE(xrational(1, natural(243)) == 9 * sqrt(3_x));

    // the root is not fully simplified, but it still compares equal to the simplified form
    REQUIRE(xrational(1, natural(50)) == 5 * sqrt(2_x));
    REQUIRE(xrational(1, natural(98)) == 7 * sqrt(2_x));

    // square factors of a multi word root
    const xrational big(1, natural(1) << 101); // sqrt(2**101) == 2**50 * sqrt(2)
    REQUIRE(big == xrational(rational(natural(1) << 50), natural(2)));
    REQUIRE(sqr(big) == rational(natural(1) << 101));
    const xrational big3(1, abs(pow(integer(3), 100)));
    REQUIRE(big3 == rational(pow(integer(3), 50)));
    REQUIRE(big3.is_rational());

    // base and root of a simplified value
    REQUIRE(xrational(2, natural(9)) == 6);
    REQUIRE(xrational(1/3_q, natural(9)) == 1);
    REQUIRE(xrational(1/6_q, natural(18)) == xrational(1/2_q, natural(2)));
    REQUIRE(xrational(1/6_q, natural(9)).base == 1/2_q); // base stays in lowest terms
    REQUIRE(xrational(2/4_q, natural(3)).base == 1/2_q);
    REQUIRE(xrational(-1, natural(8)) == -2 * sqrt(2_x));
    REQUIRE(xrational(-1/6_q, natural(9)).base == -1/2_q);

    // zero always has root 1
    REQUIRE(xrational(0, natural(5)).root == 1);
    REQUIRE(xrational(0, natural(5)).is_zero());
    REQUIRE(xrational(0, natural(9)).root == 1);

    REQUIRE_THROWS_AS(xrational(1, natural(0)), std::runtime_error);
}

TEST_CASE("is_rational") {
    const xrational a(3/2_q, natural(5));
    REQUIRE(a.base == 3/2_q);
    REQUIRE(a.root == 5);
    REQUIRE(!a.is_rational());
    REQUIRE(!a.is_zero());
    REQUIRE(!a.is_negative());

    const xrational b(-3/2_q, natural(1));
    REQUIRE(b.is_rational());
    REQUIRE(b.is_negative());
    REQUIRE(b.base == -3/2_q);
    REQUIRE(b.root == 1);

    REQUIRE(xrational().is_rational());
    REQUIRE(xrational().is_zero());
    REQUIRE(xrational() == 0);
    REQUIRE(xrational(7).is_rational());
    REQUIRE(xrational(7_q).is_rational());
    REQUIRE(xrational(0).is_rational());

    REQUIRE(!sqrt(2_x).is_rational());
    REQUIRE(!(sqrt(2_x) * sqrt(3_x)).is_rational());
    REQUIRE((sqrt(2_x) * sqrt(2_x)).is_rational());
    REQUIRE((sqrt(2_x) * sqrt(8_x)).is_rational());
    REQUIRE((sqrt(2_x) / sqrt(2_x)).is_rational());
    REQUIRE((sqrt(2_x) / sqrt(8_x)).is_rational());
    REQUIRE(sqr(sqrt(2_x)).is_rational());
    REQUIRE(sqrt(4_x).is_rational());
    REQUIRE(sqrt(0_x).is_rational());
    REQUIRE(xrational(1, natural(256)).is_rational());
}

TEST_CASE("negate signum abs") {
    xrational a = sqrt(2_x);
    a.negate();
    REQUIRE(a == -sqrt(2_x));
    REQUIRE(a.is_negative());
    negate(a);
    REQUIRE(a == sqrt(2_x));
    REQUIRE(!a.is_negative());

    REQUIRE(signum(sqrt(2_x)) == 1);
    REQUIRE(signum(-sqrt(2_x)) == -1);
    REQUIRE(signum(xrational(0)) == 0);
    REQUIRE(signum(xrational(-7/3_q)) == -1);

    REQUIRE(abs(sqrt(2_x)) == sqrt(2_x));
    REQUIRE(abs(-sqrt(2_x)) == sqrt(2_x));
    REQUIRE(abs(xrational(-3/2_q)) == 3/2_q);
    REQUIRE(abs(xrational(0)) == 0);
    REQUIRE(abs(-xrational(3/2_q, natural(5))) == xrational(3/2_q, natural(5)));

    // negation and abs leave their argument alone
    const xrational b = -sqrt(3_x);
    REQUIRE(abs(b) == -b);
    REQUIRE(b.is_negative());
    REQUIRE(-(-b) == b);
    REQUIRE(-xrational(0) == 0);
}

TEST_CASE("invert") {
    xrational a = sqrt(2_x);
    a.invert();
    REQUIRE(a == xrational(1/2_q, natural(2))); // 1/sqrt(2) == sqrt(2)/2
    REQUIRE(a * sqrt(2_x) == 1);

    const xrational b(3/2_q, natural(5));
    xrational c = b;
    c.invert();
    REQUIRE(b * c == 1);
    REQUIRE(c == xrational(2/15_q, natural(5)));

    xrational d(-7/3_q);
    d.invert();
    REQUIRE(d == -3/7_q);

    xrational e = sqrt(5_x);
    e.invert();
    e.invert();
    REQUIRE(e == sqrt(5_x));

    xrational f = -sqrt(3_x);
    f.invert();
    REQUIRE(f * (-sqrt(3_x)) == 1);
    REQUIRE(f.is_negative());
}

TEST_CASE("__addition_compatible") {
    REQUIRE_NOTHROW(__addition_compatible(sqrt(2_x), sqrt(2_x)));
    REQUIRE_NOTHROW(__addition_compatible(xrational(1), xrational(2)));
    REQUIRE_NOTHROW(__addition_compatible(xrational(0), xrational(0)));
    REQUIRE_NOTHROW(__addition_compatible(xrational(3/2_q, natural(5)), xrational(-1/7_q, natural(5))));

    REQUIRE_THROWS_AS(__addition_compatible(sqrt(2_x), sqrt(3_x)), std::runtime_error);
    REQUIRE_THROWS_AS(__addition_compatible(sqrt(2_x), xrational(1)), std::runtime_error);
    REQUIRE_THROWS_AS(__addition_compatible(xrational(1), sqrt(2_x)), std::runtime_error);
    // roots are compared as stored, not fully simplified
    REQUIRE_THROWS_AS(__addition_compatible(xrational(1, natural(50)), 5 * sqrt(2_x)), std::runtime_error);
}

TEST_CASE("__less_abs") {
    // rational only: compares magnitudes
    REQUIRE(__less_abs(rational(1, 2), rational(2, 3)));
    REQUIRE(!__less_abs(rational(2, 3), rational(1, 2)));
    REQUIRE(__less_abs(rational(-1, 2), rational(2, 3)));
    REQUIRE(__less_abs(rational(1, 2), rational(-2, 3)));
    REQUIRE(!__less_abs(rational(-2, 3), rational(1, 2)));
    REQUIRE(!__less_abs(rational(1, 2), rational(1, 2)));
    REQUIRE(!__less_abs(rational(-1, 2), rational(1, 2)));
    REQUIRE(__less_abs(rational(0), rational(1, 100)));
    REQUIRE(!__less_abs(rational(0), rational(0)));
    REQUIRE(__less_abs(rational(3, 7), rational(4, 7))); // equal denominators
    REQUIRE(!__less_abs(rational(-4, 7), rational(3, 7)));

    // with roots: compares abs(base) * sqrt(root)
    // (1/2)*sqrt(3) = 0.866..  <  (1/3)*sqrt(8) = 0.942..
    REQUIRE(__less_abs(rational(1, 2), natural(3), rational(1, 3), natural(8)));
    REQUIRE(!__less_abs(rational(1, 3), natural(8), rational(1, 2), natural(3)));
    REQUIRE(__less_abs(rational(-1, 2), natural(3), rational(-1, 3), natural(8))); // signs ignored
    REQUIRE(__less_abs(rational(1, 2), natural(3), rational(-1, 3), natural(8)));

    // equal values are not less than each other, even with different roots
    REQUIRE(!__less_abs(rational(2), natural(2), rational(1), natural(8))); // 2*sqrt(2) == sqrt(8)
    REQUIRE(!__less_abs(rational(1), natural(8), rational(2), natural(2)));

    // equal bases: the larger root wins
    REQUIRE(__less_abs(rational(3, 5), natural(2), rational(3, 5), natural(3)));
    REQUIRE(!__less_abs(rational(3, 5), natural(3), rational(3, 5), natural(2)));
    REQUIRE(__less_abs(rational(-3, 5), natural(2), rational(-3, 5), natural(3)));

    // equal roots: the larger base wins
    REQUIRE(__less_abs(rational(1, 3), natural(7), rational(1, 2), natural(7)));
    REQUIRE(!__less_abs(rational(1, 2), natural(7), rational(1, 3), natural(7)));

    // zero is smaller than anything else
    REQUIRE(__less_abs(rational(0), natural(1), rational(1, 1000000), natural(2)));
    REQUIRE(!__less_abs(rational(1, 1000000), natural(2), rational(0), natural(1)));

    // multi word operands: (2**200/3)*sqrt(2) = 0.471*2**200 > (2**199/7)*sqrt(3) = 0.124*2**200
    const rational a(integer(natural(1) << 200), integer(3));
    const rational b(integer(natural(1) << 199), integer(7));
    REQUIRE(!__less_abs(a, natural(2), b, natural(3)));
    REQUIRE(__less_abs(b, natural(3), a, natural(2)));
}

TEST_CASE("__less") {
    REQUIRE(__less(rational(-1), natural(2), rational(1), natural(2)));
    REQUIRE(!__less(rational(1), natural(2), rational(-1), natural(2)));
    // -sqrt(3) < -sqrt(2)
    REQUIRE(__less(rational(-1), natural(3), rational(-1), natural(2)));
    REQUIRE(!__less(rational(-1), natural(2), rational(-1), natural(3)));
    // -2*sqrt(2) < -sqrt(3)
    REQUIRE(__less(rational(-2), natural(2), rational(-1), natural(3)));
    REQUIRE(!__less(rational(-1), natural(3), rational(-2), natural(2)));

    REQUIRE(__less(rational(0), natural(1), rational(1, 1000), natural(1)));
    REQUIRE(__less(rational(-1, 1000), natural(1), rational(0), natural(1)));
    REQUIRE(!__less(rational(0), natural(1), rational(0), natural(1)));
    REQUIRE(!__less(rational(3), natural(5), rational(3), natural(5)));
    REQUIRE(!__less(rational(2), natural(2), rational(1), natural(8))); // equal values
    REQUIRE(!__less(rational(-2), natural(2), rational(-1), natural(8)));
}

TEST_CASE("ordering") {
    const std::vector<xrational> v = {
        -3 * sqrt(2_x),                     // -4.2426
        -2_x,                               // -2
        -sqrt(3_x),                         // -1.7320
        -sqrt(2_x),                         // -1.4142
        xrational(-4/3_q),                  // -1.3333
        -1_x,                               // -1
        xrational(-1/2_q),                  // -0.5
        -xrational(1/3_q, natural(2)),      // -0.4714
        0_x,                                //  0
        xrational(1/3_q, natural(2)),       //  0.4714
        xrational(1/2_q),                   //  0.5
        xrational(1/2_q, natural(2)),       //  0.7071
        1_x,                                //  1
        sqrt(2_x),                          //  1.4142
        sqrt(3_x),                          //  1.7320
        2_x,                                //  2
        sqrt(5_x),                          //  2.2360
        xrational(3/2_q, natural(5)),       //  3.3541
        5_x,                                //  5
        xrational(1, natural(50)),          //  7.0710
    };
    for (size_t i = 0; i < v.size(); i++)
        for (size_t j = 0; j < v.size(); j++) {
            REQUIRE((v[i] < v[j]) == (i < j));
            REQUIRE((v[i] > v[j]) == (i > j));
            REQUIRE((v[i] <= v[j]) == (i <= j));
            REQUIRE((v[i] >= v[j]) == (i >= j));
            REQUIRE((v[i] == v[j]) == (i == j));
            REQUIRE((v[i] != v[j]) == (i != j));
        }
}

TEST_CASE("compare with rational") {
    REQUIRE(sqrt(2_x) > 1);
    REQUIRE(1 < sqrt(2_x));
    REQUIRE(sqrt(2_x) < 2);
    REQUIRE(2 > sqrt(2_x));
    REQUIRE(sqrt(2_x) > 7/5_q);   // 1.41421 > 1.4
    REQUIRE(sqrt(2_x) < 3/2_q);
    REQUIRE(3/2_q > sqrt(2_x));
    REQUIRE(-sqrt(2_x) < -7/5_q);
    REQUIRE(-7/5_q > -sqrt(2_x));
    REQUIRE(sqrt(2_x) > 41/29_q);      // 1.413793 < sqrt(2)
    REQUIRE(sqrt(2_x) < 99/70_q);      // 1.414285 > sqrt(2)
    REQUIRE(sqrt(2_x) < 665857/470832_q);
    REQUIRE(sqrt(2_x) > 470832/332929_q);

    REQUIRE(2_x >= 2);
    REQUIRE(2_x <= 2);
    REQUIRE(2_x != 3);
    REQUIRE(2 <= 2_x);
    REQUIRE(2 >= 2_x);
    REQUIRE(sqrt(4_x) == 2);
    REQUIRE(2 == sqrt(4_x));
    REQUIRE(sqrt(2_x) != 2);
    REQUIRE(2 != sqrt(2_x));
    REQUIRE(0_x == 0);
    REQUIRE(sqrt(2_x) > 0);
    REQUIRE(0 < sqrt(2_x));
    REQUIRE(-sqrt(2_x) < 0);
    REQUIRE(!(sqrt(2_x) < sqrt(2_x)));
    REQUIRE(xrational(1, natural(50)) == 5 * sqrt(2_x));
    REQUIRE(!(xrational(1, natural(50)) < 5 * sqrt(2_x)));
    REQUIRE(!(5 * sqrt(2_x) < xrational(1, natural(50))));
    REQUIRE(xrational(-7/3_q) < integer(-2));
    REQUIRE(xrational(-7/3_q) > integer(-3));
}

TEST_CASE("add sub with rationals") {
    REQUIRE(sqrt(2_x) + sqrt(2_x) == 2 * sqrt(2_x));
    REQUIRE((sqrt(2_x) - sqrt(2_x)).is_zero());
    REQUIRE(xrational(3/2_q, natural(5)) + xrational(1/2_q, natural(5)) == 2 * sqrt(5_x));
    REQUIRE(xrational(3/2_q, natural(5)) - xrational(1/2_q, natural(5)) == sqrt(5_x));
    REQUIRE(1_x + 2_x == 3);
    REQUIRE(1_x - 2_x == -1);
    REQUIRE(0_x + 0_x == 0);

    REQUIRE(2_x + 3 == 5);
    REQUIRE(3 + 2_x == 5);
    REQUIRE(2_x - 3 == -1);
    REQUIRE(3 - 2_x == 1);
    REQUIRE(2_x + 1/2_q == 5/2_q);
    REQUIRE(2_x - 1/2_q == 3/2_q);
    REQUIRE(1/2_q - 2_x == -3/2_q);
    REQUIRE(2_x + integer(3) == 5);
    REQUIRE(2_x - integer(3) == -1);

    // roots that share a square factor can still be added
    REQUIRE(xrational(1, natural(50)) + sqrt(2_x) == 6 * sqrt(2_x));
    REQUIRE(xrational(1, natural(50)) - sqrt(2_x) == 4 * sqrt(2_x));
    REQUIRE(sqrt(2_x) - xrational(1, natural(50)) == -4 * sqrt(2_x));
    REQUIRE(xrational(1, natural(98)) + xrational(1, natural(50)) == 12 * sqrt(2_x));

    // otherwise adding is not possible within this type
    REQUIRE_THROWS_AS(sqrt(2_x) + 1, std::runtime_error);
    REQUIRE_THROWS_AS(1 + sqrt(2_x), std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(2_x) - 1, std::runtime_error);
    REQUIRE_THROWS_AS(1 - sqrt(2_x), std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(2_x) + sqrt(3_x), std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(2_x) - sqrt(3_x), std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(6_x) + sqrt(10_x), std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(2_x) + 1/2_q, std::runtime_error);

    // identities
    const xrational a(3/2_q, natural(5));
    const xrational b(-1/7_q, natural(5));
    REQUIRE(a + b == b + a);
    REQUIRE(a - b == -(b - a));
    REQUIRE((a + b) - b == a);
    REQUIRE(a + 0_x == a);
    REQUIRE(0_x + a == a);
    REQUIRE(0_x - a == -a);
    REQUIRE(a - 0_x == a);
}

TEST_CASE("compound assignment with xrational") {
    xrational a = sqrt(5_x);
    a += sqrt(5_x);
    REQUIRE(a == 2 * sqrt(5_x));
    a -= sqrt(5_x);
    REQUIRE(a == sqrt(5_x));
    a *= sqrt(5_x);
    REQUIRE(a == 5);
    a /= sqrt(5_x);
    REQUIRE(a == sqrt(5_x));

    xrational b = 3_x;
    b += 4_x;
    REQUIRE(b == 7);
    b -= 10_x;
    REQUIRE(b == -3);
    b *= 2_x;
    REQUIRE(b == -6);
    b /= 3_x;
    REQUIRE(b == -2);

    xrational c = 5_x;
    c += 1/2_q;
    REQUIRE(c == 11/2_q);
    c -= 1/2_q;
    REQUIRE(c == 5);
    c *= 2;
    REQUIRE(c == 10);

    xrational d = sqrt(2_x);
    d += 0;
    REQUIRE(d == sqrt(2_x));
    d -= 0_q;
    REQUIRE(d == sqrt(2_x));
    REQUIRE_THROWS_AS(d += 1, std::runtime_error);
    REQUIRE_THROWS_AS(d -= 1, std::runtime_error);
    REQUIRE_THROWS_AS(d += sqrt(3_x), std::runtime_error);
    REQUIRE_THROWS_AS(d -= sqrt(3_x), std::runtime_error);

    // += and -= with a matching root, through the xrational overload
    xrational e = sqrt(3_x);
    e += xrational(2, natural(3));
    REQUIRE(e == 3 * sqrt(3_x));
    e -= xrational(2, natural(3));
    REQUIRE(e == sqrt(3_x));
}

TEST_CASE("assignment") {
    xrational a = sqrt(2_x);
    a = sqrt(3_x);
    REQUIRE(a == sqrt(3_x));
    a = 5;
    REQUIRE(a == 5);
    REQUIRE(a.root == 1);
    REQUIRE(a.is_rational());
    a = sqrt(2_x);
    a = 3/4_q;
    REQUIRE(a == 3/4_q);
    REQUIRE(a.root == 1);
    a = integer(-6);
    REQUIRE(a == -6);

    const xrational b(3/2_q, natural(5));
    xrational c(b);
    REQUIRE(c == b);
    c = b;
    REQUIRE(c == b);
}

TEST_CASE("multiply") {
    REQUIRE(sqrt(2_x) * sqrt(2_x) == 2);
    REQUIRE(sqrt(2_x) * sqrt(3_x) == sqrt(6_x));
    REQUIRE(sqrt(3_x) * sqrt(2_x) == sqrt(6_x));
    REQUIRE(sqrt(2_x) * sqrt(8_x) == 4);
    REQUIRE(sqrt(8_x) * sqrt(18_x) == 12);
    REQUIRE((-sqrt(2_x)) * sqrt(3_x) == -sqrt(6_x));
    REQUIRE(sqrt(2_x) * (-sqrt(3_x)) == -sqrt(6_x));
    REQUIRE((-sqrt(2_x)) * (-sqrt(3_x)) == sqrt(6_x));
    REQUIRE(sqrt(2_x) * 0_x == 0);
    REQUIRE(sqrt(2_x) * 3 == 3 * sqrt(2_x));
    REQUIRE(sqrt(2_x) * 1/2_q == xrational(1/2_q, natural(2)));
    REQUIRE(sqrt(2_x) * integer(3) == 3 * sqrt(2_x));
    REQUIRE(2_x * 3_x == 6);
    REQUIRE(xrational(3/2_q, natural(5)) * xrational(2/3_q, natural(5)) == 5);
    REQUIRE(xrational(1/2_q, natural(2)) * xrational(1/3_q, natural(3)) == xrational(1/6_q, natural(6)));

    // a * a and sqr(a) agree
    for (const xrational& a : {sqrt(2_x), -sqrt(3_x), xrational(3/2_q, natural(5)), xrational(-2/7_q), 0_x, 1_x}) {
        REQUIRE(a * a == sqr(a));
        REQUIRE(sqr(a) >= 0);
    }
    // associativity and commutativity
    const xrational p = sqrt(2_x), q = xrational(3/2_q, natural(5)), r = xrational(-2/7_q);
    REQUIRE(p * q == q * p);
    REQUIRE((p * q) * r == p * (q * r));
    REQUIRE(p * 1_x == p);
}

TEST_CASE("divide") {
    REQUIRE(sqrt(6_x) / sqrt(2_x) == sqrt(3_x));
    REQUIRE(sqrt(3_x) / sqrt(6_x) == xrational(1/2_q, natural(2)));
    REQUIRE(sqrt(2_x) / sqrt(2_x) == 1);
    REQUIRE(sqrt(8_x) / sqrt(2_x) == 2);
    REQUIRE(sqrt(2_x) / 2 == xrational(1/2_q, natural(2)));
    REQUIRE(2 / sqrt(2_x) == sqrt(2_x));
    REQUIRE(1 / sqrt(2_x) * sqrt(2_x) == 1);
    REQUIRE(6_x / 3 == 2);
    REQUIRE(6_x / 4_q == 3/2_q);
    REQUIRE(0_x / sqrt(2_x) == 0);
    REQUIRE(xrational(3/2_q, natural(5)) / xrational(1/2_q, natural(5)) == 3);
    REQUIRE(sqrt(15_x) / sqrt(3_x) == sqrt(5_x));
    REQUIRE(-sqrt(6_x) / sqrt(2_x) == -sqrt(3_x));
    REQUIRE(sqrt(6_x) / (-sqrt(2_x)) == -sqrt(3_x));

    // (a / b) * b == a
    const xrational a(3/2_q, natural(5));
    for (const xrational& b : {sqrt(2_x), -sqrt(3_x), xrational(2/7_q), 1_x, xrational(1, natural(50))}) {
        REQUIRE((a / b) * b == a);
        REQUIRE(a / b == a * (1 / b));
    }

    REQUIRE_THROWS_AS(sqrt(2_x) / 0_x, std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(2_x) / 0, std::runtime_error);
}

TEST_CASE("sqrt and sqr") {
    for (int n : {0, 1, 2, 3, 4, 5, 8, 9, 12, 16, 18, 20, 27, 32, 45, 48, 49, 50, 64, 72, 100, 128, 200}) {
        const xrational r = sqrt(xrational(n));
        REQUIRE(sqr(r) == n);
        REQUIRE(r * r == n);
        REQUIRE(r >= 0);
    }
    REQUIRE(sqrt(0_x) == 0);
    REQUIRE(sqrt(0_x).is_zero());
    REQUIRE(sqrt(1_x) == 1);
    REQUIRE(sqrt(9_x) == 3);
    REQUIRE(sqrt(50_x) == 5 * sqrt(2_x));
    REQUIRE(sqrt(1/2_x) == xrational(1/2_q, natural(2)));
    REQUIRE(sqrt(2/9_x) == xrational(1/3_q, natural(2)));
    REQUIRE(sqrt(1/4_x) == 1/2_q);
    REQUIRE(sqrt(9/16_x) == 3/4_q);
    REQUIRE_THROWS_AS(sqrt(-1_x), std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(-1/2_x), std::runtime_error);
    REQUIRE_THROWS_AS(sqrt(sqrt(2_x)), std::runtime_error);

    REQUIRE(sqr(xrational(0)) == 0);
    REQUIRE(sqr(xrational(-3/2_q)) == 9/4_q);
    REQUIRE(sqr(xrational(3/2_q, natural(5))) == 45/4_q);
    REQUIRE(sqr(xrational(-3/2_q, natural(5))) == 45/4_q);
    REQUIRE(sqr(sqrt(7_x)) == 7);
    REQUIRE(sqr(xrational(1, natural(50))) == 50);
    REQUIRE(sqr(sqrt(2_x)).is_rational());
}

TEST_CASE("pow") {
    REQUIRE(pow(sqrt(2_x), 0) == 1);
    REQUIRE(pow(sqrt(2_x), 1) == sqrt(2_x));
    REQUIRE(pow(sqrt(2_x), 2) == 2);
    REQUIRE(pow(sqrt(2_x), 3) == 2 * sqrt(2_x));
    REQUIRE(pow(sqrt(2_x), 4) == 4);
    REQUIRE(pow(sqrt(2_x), 6) == 8);
    REQUIRE(pow(2_x, 10) == 1024);
    REQUIRE(pow(-2_x, 3) == -8);
    REQUIRE(pow(-2_x, 4) == 16);
    REQUIRE(pow(3/2_x, 3) == 27/8_q);
    REQUIRE(pow(0_x, 3) == 0);
    REQUIRE(pow(0_x, 0) == 1);
    REQUIRE(pow(1_x, 1000) == 1);

    // negative exponents
    REQUIRE(pow(2_x, -1) == 1/2_q);
    REQUIRE(pow(2_x, -2) == 1/4_q);
    REQUIRE(pow(sqrt(2_x), -2) == 1/2_q);
    REQUIRE(pow(sqrt(2_x), -3) == xrational(1/4_q, natural(2)));
    REQUIRE(pow(-2_x, -3) == -1/8_q);
    REQUIRE(pow(3/2_x, -2) == 4/9_q);

    // pow(a, e + 1) == pow(a, e) * a
    for (int e = 0; e < 8; e++) {
        REQUIRE(pow(sqrt(2_x), e + 1) == pow(sqrt(2_x), e) * sqrt(2_x));
        REQUIRE(pow(xrational(-3/2_q, natural(5)), e + 1) == pow(xrational(-3/2_q, natural(5)), e) * xrational(-3/2_q, natural(5)));
    }
}

TEST_CASE("shift") {
    REQUIRE((sqrt(2_x) << 3) == 8 * sqrt(2_x));
    REQUIRE((sqrt(2_x) >> 1) == xrational(1/2_q, natural(2)));
    REQUIRE((3_x << 2) == 12);
    REQUIRE((3_x >> 2) == 3/4_q);
    REQUIRE((-sqrt(2_x) << 1) == -2 * sqrt(2_x));
    REQUIRE((0_x << 5) == 0);

    xrational a = sqrt(5_x);
    a <<= 2;
    REQUIRE(a == 4 * sqrt(5_x));
    a >>= 3;
    REQUIRE(a == xrational(1/2_q, natural(5)));

    // shifting is multiplication/division by a power of two
    for (int i = 0; i < 4; i++) {
        REQUIRE((sqrt(3_x) << i) == sqrt(3_x) * pow(integer(2), i));
        REQUIRE((sqrt(3_x) >> i) == sqrt(3_x) / pow(integer(2), i));
    }
}

TEST_CASE("format of sqrt forms") {
    REQUIRE(format("{}", xrational(-1)) == "-1");
    REQUIRE(format("{}", xrational(-3/2_q)) == "-3/2");
    REQUIRE(format("{}", sqrt(3_x)) == "sqrt(3)");
    REQUIRE(format("{}", sqrt(50_x)) == "5*sqrt(2)");
    REQUIRE(format("{}", xrational(1, natural(50))) == "sqrt(50)");
    REQUIRE(format("{}", xrational(-3/2_q, natural(2))) == "-3/2*sqrt(2)");
    REQUIRE(format("{}", xrational(1/3_q, natural(5))) == "sqrt(5)/3");
    REQUIRE(format("{}", xrational(-1/3_q, natural(5))) == "-sqrt(5)/3");
    REQUIRE(format("{}", -sqrt(2_x)) == "-1*sqrt(2)"); // generic form, base is an integer
    REQUIRE(format("{}", xrational(0)) == "0");
    REQUIRE(format("{}", 2_x * sqrt(2_x)) == "2*sqrt(2)");
    REQUIRE(format("{}", xrational(1, natural(9))) == "3");

    // the format specifier must be empty
    xrational one = 1_x;
    REQUIRE_THROWS_AS(std::vformat("{:.2}", std::make_format_args(one)), std::format_error);
    REQUIRE(std::vformat("{}", std::make_format_args(one)) == "1");

    std::ostringstream os;
    os << sqrt(2_x) << ' ' << xrational(-3/2_q) << ' ' << xrational(1/3_q, natural(5));
    REQUIRE(os.str() == "sqrt(2) -3/2 sqrt(5)/3");
}

TEST_CASE("hash") {
    std::hash<xrational> h;
    REQUIRE(h(sqrt(2_x)) == h(sqrt(2_x)));
    REQUIRE(h(xrational(3/2_q, natural(5))) == h(xrational(3/2_q, natural(5))));
    REQUIRE(h(sqrt(2_x)) != h(sqrt(3_x)));
    REQUIRE(h(sqrt(2_x)) != h(-sqrt(2_x)));
    REQUIRE(h(xrational(1, natural(4))) == h(2_x)); // both simplify to the same value

    std::unordered_map<xrational, int> m;
    m[sqrt(2_x)] = 1;
    m[sqrt(2_x)] += 1;
    m[sqrt(3_x)] = 3;
    m[xrational(1, natural(4))] = 4;
    REQUIRE(m.size() == 3);
    REQUIRE(m[sqrt(2_x)] == 2);
    REQUIRE(m[2_x] == 4);
}
