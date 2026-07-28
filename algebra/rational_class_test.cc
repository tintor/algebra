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

    // trailing 'f' is accepted, like for floating point
    REQUIRE(format("{:.2f}", rational(3, 4)) == "0.75");
    REQUIRE(format("{:.3f}", rational(1, 1000)) == "0.001");
    REQUIRE(format("{:.2f}", rational(15, 2)) == "7.50");
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

// The whole library is marked constexpr; these check that it can actually be evaluated at
// compile time. Nothing below compiles if a kernel accessor or the division helper is not
// usable in a constant expression.
static_assert([] { integer a = 12345; a *= 6789u; return a == 83810205u; }());
static_assert([] { integer a = 1; a <<= 200; a -= 1u; return a.num_bits() == 200; }());
static_assert([] { integer a("123456789012345678901234567890"); return a % 97u == 52u; }());
static_assert([] { integer a = 1; a <<= 128; integer q, r; div(a, integer(1000), q, r); return r == 456u; }());
static_assert([] { integer a = -5; a *= 7; return a == -35; }());
static_assert([] { integer a("-123456789012345678901234567890"); return a.is_negative() && a.num_bits() == 97; }());
static_assert([] { rational a(1, 3); a += rational(1, 6); return a == rational(1, 2); }());
static_assert([] { rational a(2, 3); a /= rational(4, 9); return a == rational(3, 2); }());

TEST_CASE("constexpr evaluation") {
    // the static_asserts above do the work; this keeps them visible in the test report
    SUCCEED("number types are usable in constant expressions");
}


TEST_CASE("from non finite floating point") {
    // note: these guards are compiled away by -ffast-math, which is why the build must not use it
    volatile double zero = 0.0;
    REQUIRE_THROWS(rational(0.0 / zero));   // nan
    REQUIRE_THROWS(rational(1.0 / zero));   // +inf
    REQUIRE_THROWS(rational(-1.0 / zero));  // -inf
    volatile float fzero = 0.0f;
    REQUIRE_THROWS(rational(0.0f / fzero));
    REQUIRE_THROWS(rational(1.0f / fzero));
}

// rational is expected to keep the sign in the numerator, the denominator strictly positive,
// the fraction in lowest terms, and zero as exactly 0/1.
static void check_reduced(const rational& a) {
    REQUIRE(!a.den.is_zero());
    REQUIRE(!a.den.is_negative());
    REQUIRE(gcd(a.num, a.den) == 1);
    if (a.num.is_zero())
        REQUIRE(a.den == 1);
}

TEST_CASE("normalized") {
    // normalized() is the fast path for callers that already know num and den are coprime
    // and den > 0: it stores them as given, without running simplify()
    const rational a = rational::normalized(integer(3), integer(4));
    REQUIRE(a.num == 3);
    REQUIRE(a.den == 4);
    REQUIRE(a == rational(3, 4));

    const rational b = rational::normalized(integer(-3), integer(4));
    REQUIRE(b.num == -3);
    REQUIRE(b.den == 4);
    REQUIRE(b == -a);

    REQUIRE(rational::normalized(integer(5), integer(1)).is_integer());
    REQUIRE(rational::normalized(integer(0), integer(1)).is_zero());

    // no reduction happens, so a reducible pair keeps its terms: it is equal in value but
    // not equal by operator==, which compares num and den
    const rational c = rational::normalized(integer(2), integer(4));
    REQUIRE(c.num == 2);
    REQUIRE(c.den == 4);
    REQUIRE(!(c == rational(1, 2)));
    REQUIRE(!(c < rational(1, 2)));
    REQUIRE(!(rational(1, 2) < c));

}

TEST_CASE("is_integer") {
    REQUIRE(rational().is_integer());
    REQUIRE(rational(0).is_integer());
    REQUIRE(rational(7).is_integer());
    REQUIRE(rational(-7).is_integer());
    REQUIRE(rational(8, 4).is_integer());
    REQUIRE(rational(0, 5).is_integer());
    REQUIRE(rational(integer(9)).is_integer());
    REQUIRE(!rational(1, 2).is_integer());
    REQUIRE(!rational(-1, 2).is_integer());
    REQUIRE(!rational(7, 2).is_integer());

    // integral only once the value actually reduces to a whole number
    REQUIRE((rational(1, 3) + rational(2, 3)).is_integer());
    REQUIRE(!(rational(1, 3) + rational(1, 3)).is_integer());
    REQUIRE((rational(1, 2) * rational(2)).is_integer());
    REQUIRE((rational(2, 3) / rational(2, 3)).is_integer());
    REQUIRE((rational(3, 4) - rational(1, 4) - rational(1, 2)).is_integer());
    REQUIRE(rational("1.0").is_integer());
    REQUIRE(rational("10/5").is_integer());
    REQUIRE(!rational("1.5").is_integer());
    REQUIRE(rational(4.0).is_integer());
    REQUIRE(!rational(0.5).is_integer());

    // is_even / is_odd are false for everything that is not an integer
    REQUIRE(!rational(1, 2).is_even());
    REQUIRE(!rational(1, 2).is_odd());
    REQUIRE(!rational(3, 2).is_even());
    REQUIRE(!rational(3, 2).is_odd());
    REQUIRE(rational(0).is_even());
    REQUIRE(rational(4, 2).is_even());
    REQUIRE(rational(9, 3).is_odd());
    REQUIRE(rational(-9, 3).is_odd());
}

TEST_CASE("sign is kept in the numerator") {
    REQUIRE(rational(1, -2).num == -1);
    REQUIRE(rational(1, -2).den == 2);
    REQUIRE(rational(integer(1), integer(-2)).num == -1);
    REQUIRE(rational(integer(1), integer(-2)).den == 2);
    REQUIRE(rational(-1, -2).num == 1);
    REQUIRE(rational(-1, -2).den == 2);
    REQUIRE(rational("-1/2").num == -1);
    REQUIRE(rational("-1/2").den == 2);
    REQUIRE(rational(-0.5).num == -1);
    REQUIRE(rational(-0.5).den == 2);

    rational a(1, 2);
    a.negate();
    REQUIRE(a.num == -1);
    REQUIRE(a.den == 2);
    negate(a);
    REQUIRE(a.num == 1);
    REQUIRE((-rational(1, 2)).num == -1);
    REQUIRE((-rational(1, 2)).den == 2);

    // invert moves a negative sign back into the numerator
    rational b(-2, 3);
    b.invert();
    REQUIRE(b.num == -3);
    REQUIRE(b.den == 2);

    // dividing by a negative value leaves den positive
    REQUIRE((rational(1, 2) / rational(-1, 3)).num == -3);
    REQUIRE((rational(1, 2) / rational(-1, 3)).den == 2);
    REQUIRE((rational(1, 2) / integer(-3)).num == -1);
    REQUIRE((rational(1, 2) / integer(-3)).den == 6);
    REQUIRE((integer(1) / rational(-1, 2)).num == -2);
    REQUIRE((integer(1) / rational(-1, 2)).den == 1);
    REQUIRE((rational(1, 2) * integer(-6)).num == -3);
    REQUIRE((rational(1, 2) * integer(-6)).den == 1);
    {
        rational c(2, 3);
        c /= rational(-4, 9);
        REQUIRE(c.num == -3);
        REQUIRE(c.den == 2);
    }

    REQUIRE(signum(rational(-1, 2)) == -1);
    REQUIRE(signum(rational(1, 2)) == 1);
    REQUIRE(signum(rational(0)) == 0);
    REQUIRE(rational(1, -2).is_negative());
    REQUIRE(!rational(-1, -2).is_negative());

    // zero is canonical: 0/1
    REQUIRE(rational(0, -5).num == 0);
    REQUIRE(rational(0, -5).den == 1);
    REQUIRE(rational(integer(0), integer(-5)).den == 1);
    REQUIRE((rational(1, 2) - rational(1, 2)).den == 1);
    REQUIRE((rational(0) * rational(3, 4)).den == 1);
    REQUIRE((rational(0) / rational(3, 4)).den == 1);
    REQUIRE((rational(3, 4) >> 10).num == 3);
    REQUIRE((rational(0) >> 10).den == 1);
    REQUIRE(rational(-0.0).den == 1);
}

TEST_CASE("lowest terms invariant") {
    REQUIRE(rational(6, 4).num == 3);
    REQUIRE(rational(6, 4).den == 2);
    REQUIRE(rational(integer(6), integer(-4)).num == -3);
    REQUIRE(rational(integer(6), integer(-4)).den == 2);
    REQUIRE(rational(12, 18) == rational::normalized(integer(2), integer(3)));

    // 2^64*3 / 2^64*5 reduces to 3/5 through the multi word path
    const integer p = integer(1) << 64;
    REQUIRE(rational(p * 3, p * 5) == rational::normalized(integer(3), integer(5)));
    REQUIRE(rational(p * 3, -(p * 5)) == rational::normalized(integer(-3), integer(5)));

    // 1/6 + 1/3 has 18 as the naive denominator, but must come out as 1/2
    REQUIRE((rational(1, 6) + rational(1, 3)).den == 2);
    REQUIRE((rational(1, 6) + rational(1, 3)).num == 1);
    // 3/10 * 5/9 = 1/6
    REQUIRE((rational(3, 10) * rational(5, 9)) == rational::normalized(integer(1), integer(6)));
    // (5/6) / (10/3) = 1/4
    REQUIRE((rational(5, 6) / rational(10, 3)) == rational::normalized(integer(1), integer(4)));

    check_reduced(rational());
    check_reduced(rational(0));
    check_reduced(rational(-7));
    check_reduced(rational(p * 3, p * 5));

    for (const char* s : {"0", "-0", "0.0", "-0.000", "0/7", "-0/7", "12/8", "-12/8", "1.5",
                          "-2.50", "1e-3", "6/4e1", "0.000000000000000000004", "-1/1024"})
        check_reduced(rational(s));

    for (double x : {0.0, -0.0, 1.0, -1.0, 0.1, -0.1, 1e-20, 1e20, 12345.6789, 0.5})
        check_reduced(rational(x));

    Random rng(42);
    for (int i = 0; i < 200; i++) {
        auto r = [&] { return rational(rng.Uniform<int>(-60, 60), rng.Uniform<int>(1, 60)); };
        const rational a = r(), b = r();
        check_reduced(a);
        check_reduced(b);
        check_reduced(-a);
        check_reduced(a + b);
        check_reduced(a - b);
        check_reduced(a * b);
        check_reduced(a + 3);
        check_reduced(a - 3);
        check_reduced(a * 3);
        check_reduced(a / 3);
        check_reduced(a << 2);
        check_reduced(a >> 3);
        {
            rational c = a;
            c += b;
            check_reduced(c);
            REQUIRE(c == a + b);
        }
        {
            rational c = a;
            c -= b;
            check_reduced(c);
            REQUIRE(c == a - b);
        }
        {
            rational c = a;
            c *= b;
            check_reduced(c);
            REQUIRE(c == a * b);
        }
        if (!b.is_zero()) {
            check_reduced(a / b);
            check_reduced(a % b);
            rational c = a;
            c /= b;
            check_reduced(c);
            REQUIRE(c == a / b);
        }
        if (!a.is_zero()) {
            rational c = a;
            c.invert();
            check_reduced(c);
            REQUIRE(c * a == 1);
        }
    }
}

TEST_CASE("str") {
    REQUIRE(rational(0).str() == "0");
    REQUIRE(rational(6).str() == "6");
    REQUIRE(rational(-5).str() == "-5");
    REQUIRE(rational(1, 3).str() == "1/3");
    REQUIRE(rational(-22, 8).str() == "-11/4");
    REQUIRE(rational("123456789012345678901234567890").str() == "123456789012345678901234567890");
    REQUIRE(rational("-1/123456789012345678901234567890").str() == "-1/123456789012345678901234567890");
}

TEST_CASE("ostream") {
    std::ostringstream s;
    s << rational(0) << ' ' << rational(-5) << ' ' << rational(3, 4) << ' ' << rational(-22, 8);
    REQUIRE(s.str() == "0 -5 3/4 -11/4");
}

TEST_CASE("format into other output iterators") {
    // the formatter has to write through the iterator it is given and return the new position
    std::string s;
    std::format_to(std::back_inserter(s), "{}", rational(-22, 8));
    REQUIRE(s == "-11/4");

    char buf[32] = {};
    auto r = std::format_to_n(buf, sizeof(buf), "{}", rational(-1, 3));
    REQUIRE(r.size == 4);
    REQUIRE(std::string(buf, r.out) == "-1/3");

    char small[2] = {};
    auto r2 = std::format_to_n(small, sizeof(small), "{}", rational(-1, 3));
    REQUIRE(r2.size == 4); // what would have been written
    REQUIRE(std::string(small, r2.out) == "-1");

    // the buffer capacity is computed from str_size_upper_bound of both terms
    const rational big(integer("123456789012345678901234567890"), integer("98765432109876543210987654321"));
    REQUIRE(format("{}", big) == big.str());
    REQUIRE(format("{}", rational(integer("100000000000000000000"), integer(3))) == "100000000000000000000/3");
    REQUIRE(format("{}", rational(integer(3), integer("100000000000000000000"))) == "3/100000000000000000000");

    REQUIRE(format("[{}]", rational(1, 2)) == "[1/2]");
    REQUIRE(format("{}{}", rational(1, 2), rational(3)) == "1/23");
    REQUIRE(format("{:.2}|{}", rational(1, 4), rational(-5)) == "0.25|-5");
}

TEST_CASE("format specifier is rejected when malformed") {
    // std::format checks the specifier at compile time, so these need the runtime interface
    rational a(1, 3);
    REQUIRE_THROWS_AS(std::vformat("{:.}", std::make_format_args(a)), std::format_error);
    REQUIRE_THROWS_AS(std::vformat("{:x}", std::make_format_args(a)), std::format_error);
    REQUIRE_THROWS_AS(std::vformat("{:.3g}", std::make_format_args(a)), std::format_error);
    REQUIRE_THROWS_AS(std::vformat("{:.3ff}", std::make_format_args(a)), std::format_error);
    REQUIRE_THROWS_AS(std::vformat("{:.-3}", std::make_format_args(a)), std::format_error);
    REQUIRE_THROWS_AS(std::vformat("{:3}", std::make_format_args(a)), std::format_error);

    REQUIRE(std::vformat("{}", std::make_format_args(a)) == "1/3");
    REQUIRE(std::vformat("{:.4}", std::make_format_args(a)) == "0.3333");
    REQUIRE(format("{:.10}", rational(1, 3)) == "0.3333333333");
    REQUIRE(format("{:.12}", rational(1, 7)) == "0.142857142857");
}

TEST_CASE("format with zero fraction digits") {
    // rounding is half away from zero, and there is no trailing '.'
    REQUIRE(format("{:.0}", rational(0)) == "0");
    REQUIRE(format("{:.0}", rational(7)) == "7");
    REQUIRE(format("{:.0}", rational(-7)) == "-7");
    REQUIRE(format("{:.0}", rational(1, 2)) == "1");
    REQUIRE(format("{:.0}", rational(3, 2)) == "2");
    REQUIRE(format("{:.0}", rational(-3, 2)) == "-2");
    REQUIRE(format("{:.0}", rational(1, 3)) == "0");
    REQUIRE(format("{:.0}", rational(-1, 3)) == "-0");
    REQUIRE(format("{:.0}", rational(14, 5)) == "3");
    REQUIRE(format("{:.0f}", rational(14, 5)) == "3");
    REQUIRE(format("{:.0}|", rational(5, 2)) == "3|");
}

TEST_CASE("conversion to float") {
    REQUIRE(static_cast<float>(rational(0)) == 0.0f);
    REQUIRE(static_cast<float>(rational(1)) == 1.0f);
    REQUIRE(static_cast<float>(rational(-7)) == -7.0f);
    REQUIRE(static_cast<float>(rational(1, 2)) == 0.5f);
    REQUIRE(static_cast<float>(rational(-3, 4)) == -0.75f);
    REQUIRE(static_cast<float>(rational(1, 3)) == 1.0f / 3.0f);
    REQUIRE(static_cast<float>(rational(-2, 3)) == -2.0f / 3.0f);
    REQUIRE(static_cast<float>(rational(12345, 1000)) == 12.345f);

    // 2^-e is exactly representable for every e in range, whatever the size of den
    for (int e = 0; e <= 120; e++) {
        const rational a(integer(1), integer(1) << e);
        REQUIRE(static_cast<float>(a) == std::ldexp(1.0f, -e));
        REQUIRE(static_cast<float>(-a) == -std::ldexp(1.0f, -e));
    }
    // and so is 2^e
    for (int e = 0; e <= 120; e++)
        REQUIRE(static_cast<float>(rational(integer(1) << e)) == std::ldexp(1.0f, e));

    // ratios of large numbers do not overflow on the way
    REQUIRE(static_cast<float>(rational(integer(1) << 200, integer(1) << 200)) == 1.0f);
    REQUIRE(static_cast<float>(rational(integer(3) << 200, integer(1) << 200)) == 3.0f);
    REQUIRE(static_cast<float>(rational(integer(1) << 200, integer(1) << 190)) == 1024.0f);
    REQUIRE(static_cast<float>(rational(integer(1) << 190, integer(1) << 200)) == std::ldexp(1.0f, -10));

    // out of range values saturate
    REQUIRE(std::isinf(static_cast<float>(rational(integer(1) << 200))));
    REQUIRE(static_cast<float>(rational(integer(1), integer(1) << 200)) == 0.0f);
}

TEST_CASE("conversion to double") {
    REQUIRE(static_cast<double>(rational(0)) == 0.0);
    REQUIRE(static_cast<double>(rational(1)) == 1.0);
    REQUIRE(static_cast<double>(rational(-7)) == -7.0);
    REQUIRE(static_cast<double>(rational(1, 2)) == 0.5);
    REQUIRE(static_cast<double>(rational(-3, 4)) == -0.75);
    REQUIRE(static_cast<double>(rational(1, 3)) == 1.0 / 3.0);
    REQUIRE(static_cast<double>(rational(12345, 1000)) == 12.345);

    for (int e = 0; e <= 1000; e++) {
        const rational a(integer(1), integer(1) << e);
        REQUIRE(static_cast<double>(a) == std::ldexp(1.0, -e));
        REQUIRE(static_cast<double>(rational(integer(1) << e)) == std::ldexp(1.0, e));
    }

    REQUIRE(static_cast<double>(rational(integer(1) << 2000, integer(1) << 2000)) == 1.0);
    REQUIRE(static_cast<double>(rational(integer(3) << 2000, integer(1) << 2000)) == 3.0);
    REQUIRE(static_cast<double>(rational(integer(1) << 2000, integer(1) << 1990)) == 1024.0);
    REQUIRE(static_cast<double>(rational(integer(1) << 1990, integer(1) << 2000)) == std::ldexp(1.0, -10));

    REQUIRE(std::isinf(static_cast<double>(rational(integer(1) << 2000))));
    REQUIRE(static_cast<double>(rational(integer(1), integer(1) << 2000)) == 0.0);
}

TEST_CASE("floating point round trip through rational") {
    // every finite float is a rational, so the round trip must be exact
    for (float x : {0.0f, 1.0f, -1.0f, 0.5f, 0.1f, -12345.678f, 1e-30f, 1e30f, 3.4e38f, 1.5e-38f,
                    7.0f, 1.0f / 3.0f, -0.2f})
        REQUIRE(static_cast<float>(rational(x)) == x);

    for (double x : {0.0, 1.0, -1.0, 0.1, -12345.6789, 1e-300, 1e300, 1e-100, 1e100,
                     2.2250738585072014e-308, 1.0 / 3.0, -0.2})
        REQUIRE(static_cast<double>(rational(x)) == x);

    Random rng(11);
    for (int i = 0; i < 100; i++) {
        const double x = rng.Uniform<double>(-1e6, 1e6);
        REQUIRE(static_cast<double>(rational(x)) == x);
    }
}

TEST_CASE("parse and format round trip") {
    for (const char* s : {"0", "1", "-1", "7", "-7", "1/3", "-1/3", "22/7", "-22/7",
                          "123456789012345678901234567890",
                          "-1/123456789012345678901234567890"})
        REQUIRE(format("{}", rational(s)) == s);

    Random rng(7);
    for (int i = 0; i < 200; i++) {
        const rational a(rng.Uniform<int>(-1000, 1000), rng.Uniform<int>(1, 1000));
        const std::string s = format("{}", a);
        REQUIRE(a.str() == s);
        REQUIRE(rational(s) == a);
        REQUIRE(rational(std::string_view(s)) == a);
    }

    // a fixed point rendering of an exactly representable value parses back to it
    REQUIRE(rational(format("{:.4}", rational(1, 8))) == rational(1, 8));
    REQUIRE(rational(format("{:.4}", rational(-1, 8))) == rational(-1, 8));
    REQUIRE(rational(format("{:.3}", rational(-1, 1000))) == rational(-1, 1000));
}

TEST_CASE("parse exponent and fraction forms") {
    REQUIRE(rational("1e3") == 1000);
    REQUIRE(rational("1e+3") == 1000);
    REQUIRE(rational("1e007") == 10000000);
    REQUIRE(rational("0e5") == 0);
    REQUIRE(rational("0e-5") == 0);
    REQUIRE(rational("-1e0") == -1);
    REQUIRE(rational("2/5e2") == 40);                  // (2/5) * 100
    REQUIRE(rational("2/5e-2") == rational(1, 250));   // (2/5) / 100
    REQUIRE(rational("1.5e2") == 150);
    REQUIRE(rational("-1.5e-2") == rational(-3, 200));
    REQUIRE(rational("1e100") == pow(integer(10), 100));
    REQUIRE(rational("1e-100").num == 1);
    REQUIRE(rational("1e-100").den == pow(integer(10), 100));
    REQUIRE(rational("123.456").num == 15432);         // 123456/1000
    REQUIRE(rational("123.456").den == 125);
    REQUIRE(rational("-123.456").num == -15432);
    REQUIRE(rational("-123.456").den == 125);
    REQUIRE(rational("0.1") == rational(1, 10));
    REQUIRE(rational("00.100") == rational(1, 10));
    REQUIRE(rational("10/5") == 2);
    REQUIRE(rational("-0.000") == 0);
    REQUIRE(rational("12345678901234567890.5") ==
            rational(integer("24691357802469135781"), integer(2)));
}

TEST_CASE("parse rejects malformed exponents and separators") {
    for (const char* s : {"1E3", "1e3.5", "1/2.3", "1.2/3", "1'000", "1_000", "0x10", "1e 3",
                          "1e3 ", "e3", "1.e3", "1./2", "\t1", "1\n", "1e--3", "1e+-3", "1e+ 3",
                          "1/-2", "1.-2", "-.1", "- 1", "1e1e1", "/2", "1//2", "nan", "inf"})
        REQUIRE_THROWS(rational(s));

    REQUIRE_THROWS(rational("0/0"));
    REQUIRE_THROWS(rational("0/000"));
    REQUIRE_THROWS(rational("1e100000001")); // exponent limit
}

TEST_CASE("parse rejects malformed input") {
    for (const char* s : {"", "-", "abc", "1.", "1/", "1..2", "1.2.3", "1e", "1e+", "1e-", "+1",
                          "1 ", " 1", "1x", "--1", "1/2/3", "1.2e", ".5", "1e2x"})
        REQUIRE_THROWS(rational(s));

    REQUIRE(rational("0") == 0);
    REQUIRE(rational("-0") == 0);
    REQUIRE(rational("007") == 7);
    REQUIRE(rational("1e0") == 1);
    REQUIRE(rational("12345678901234567890123456789") == integer("12345678901234567890123456789"));
    REQUIRE(rational("-12/8") == rational(-3, 2));
    REQUIRE(rational("1e+3") == 1000);
    REQUIRE(rational("2.5e-2") == rational(1, 40));
    REQUIRE_THROWS(rational("1e999999999"));
}
