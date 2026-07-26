#include "algebra/real_class.h"
#include "algebra/__test.h"
#include <sstream>
#include <unordered_map>

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

TEST_CASE("hash") {
    std::unordered_map<real<2>, int> m;
    m[1.5_f] = 0;
}

TEST_CASE("+") {
    REQUIRE(1000_d + 0 == 1000_d);
    REQUIRE(1_d + 0 == 1_d);
    REQUIRE(0.0001_d + 0 == 0.0001_d);
}

TEST_CASE("real<2> from rational") {
    REQUIRE(to_rational(real<2>(rational(4))) == rational(4));
    REQUIRE(to_rational(real<2>(rational(1))) == rational(1));
    REQUIRE(to_rational(real<2>(rational(6))) == rational(6));
    REQUIRE(to_rational(real<2>(rational(-12))) == rational(-12));
    REQUIRE(to_rational(real<2>(rational(3, 4))) == rational(3, 4));
    REQUIRE(to_rational(real<2>(rational(-3, 8))) == rational(-3, 8));
    REQUIRE(to_rational(real<2>(rational(0))) == rational(0));

    REQUIRE(to_rational(real<2>(4.0f)) == rational(4));
    REQUIRE(to_rational(real<2>(4.0)) == rational(4));
    REQUIRE(to_rational(real<2>(-6.0)) == rational(-6));
    REQUIRE(to_rational(real<2>(1.5)) == rational(3, 2));

    REQUIRE(to_rational(decimal(rational(4))) == rational(4));
    REQUIRE(to_rational(decimal(rational(-20))) == rational(-20));
    REQUIRE(to_rational(decimal(rational(1, 4))) == rational(1, 4));
}

TEST_CASE("shift") {
    // shift<B>(a, e) == a * B**e
    REQUIRE(shift<2>(integer(3), 0) == 3);
    REQUIRE(shift<2>(integer(3), 4) == 48);
    REQUIRE(shift<2>(integer(-3), 4) == -48);
    REQUIRE(shift<2>(integer(0), 7) == 0);
    REQUIRE(shift<2>(integer(1), 100) == exp2(100));

    REQUIRE(shift<10>(integer(3), 0) == 3);
    REQUIRE(shift<10>(integer(3), 2) == 300);
    REQUIRE(shift<10>(integer(-3), 3) == -3000);
    REQUIRE(shift<10>(integer(0), 4) == 0);
    REQUIRE(shift<10>(integer(7), 20) == integer("700000000000000000000"));

    REQUIRE(shift<3>(integer(2), 4) == 162); // 2 * 81
    REQUIRE(shift<3>(integer(-2), 1) == -6);
}

TEST_CASE("normalize") {
    // the representation is canonical: for base 2 the numerator is odd (or the value is zero)
    for (int n : {1, 2, 3, 4, 6, 8, 12, 48, 1024, -1, -2, -6, -48, -1024}) {
        for (int e : {-3, 0, 5}) {
            const real<2> a(n, e);
            REQUIRE(a.num % 2 != 0);
            REQUIRE(to_rational(a) == rational(n) * (rational(1) << e));
        }
    }
    // for base 10 the numerator is not divisible by 10
    for (int n : {1, 2, 5, 10, 40, 100, 1200, -5, -10, -1200}) {
        for (int e : {-3, 0, 5}) {
            const decimal a(n, e);
            REQUIRE(a.num % 10 != 0);
            REQUIRE(to_rational(a) == rational(n) * to_rational(decimal(1, e)));
        }
    }

    // zero is always num=0, exp=0
    for (int e : {-7, 0, 3}) {
        REQUIRE(real<2>(0, e).num == 0);
        REQUIRE(real<2>(0, e).exp == 0);
        REQUIRE(decimal(0, e).num == 0);
        REQUIRE(decimal(0, e).exp == 0);
    }

    // exact exponent bookkeeping
    REQUIRE(real<2>(48, -2).num == 3);
    REQUIRE(real<2>(48, -2).exp == 2);
    REQUIRE(decimal(1200, -1).num == 12);
    REQUIRE(decimal(1200, -1).exp == 1);
    // base 10 only removes factors of 10, not of 2 or 5
    REQUIRE(decimal(25, 0).num == 25);
    REQUIRE(decimal(25, 0).exp == 0);
    REQUIRE(decimal(16, 0).num == 16);
    REQUIRE(decimal(16, 0).exp == 0);

    // the 3-argument constructor skips normalize()
    REQUIRE(real<2>(integer(48), -2, 0).num == 48);
    REQUIRE(real<2>(integer(48), -2, 0).exp == -2);
    REQUIRE(to_rational(real<2>(integer(48), -2, 0)) == rational(12));
}

TEST_CASE("real<2> arithmetic") {
    using R = real<2>;
    const R a(3, -2); // 3/4
    const R b(5, 1);  // 10
    const R c(-7, 0);

    REQUIRE(to_rational(a) == rational(3, 4));
    REQUIRE(to_rational(b) == rational(10));

    REQUIRE(to_rational(-a) == rational(-3, 4));
    REQUIRE(to_rational(-b) == rational(-10));
    REQUIRE(to_rational(-R(0)) == rational(0));

    REQUIRE(to_rational(a + b) == rational(43, 4));
    REQUIRE(to_rational(b + a) == rational(43, 4));
    REQUIRE(to_rational(a - b) == rational(-37, 4));
    REQUIRE(to_rational(b - a) == rational(37, 4));
    REQUIRE(to_rational(a * b) == rational(15, 2));
    REQUIRE(to_rational(a * a) == rational(9, 16));
    REQUIRE(to_rational(a + a) == rational(3, 2));
    REQUIRE(to_rational(a - a) == rational(0));
    REQUIRE((a - a).exp == 0);

    // identities
    REQUIRE(a + b == b + a);
    REQUIRE(a * b == b * a);
    REQUIRE((a + b) + c == a + (b + c));
    REQUIRE((a * b) * c == a * (b * c));
    REQUIRE(a - b == -(b - a));
    REQUIRE((a + b) - b == a);
    REQUIRE((a + b) * c == a * c + b * c);
    REQUIRE(a + R(0) == a);
    REQUIRE(a * R(1) == a);
    REQUIRE(a * R(0) == R(0));

    // exponents line up
    REQUIRE((a + R(1, -2)).num == 1);
    REQUIRE((a + R(1, -2)).exp == 0);
    REQUIRE((R(1, 3) - R(1, 2)).num == 1);
    REQUIRE((R(1, 3) - R(1, 2)).exp == 2);

    // compound assignment
    R d = a;
    d += b;
    REQUIRE(to_rational(d) == rational(43, 4));
    d -= b;
    REQUIRE(d == a);
    d *= b;
    REQUIRE(to_rational(d) == rational(15, 2));
}

TEST_CASE("decimal arithmetic") {
    const decimal a(3, -2); // 0.03
    const decimal b(5, 1);  // 50
    const decimal c(-7, 0);

    REQUIRE(to_rational(a) == rational(3, 100));
    REQUIRE(to_rational(a + b) == rational(5003, 100));
    REQUIRE(to_rational(b + a) == rational(5003, 100));
    REQUIRE(to_rational(a - b) == rational(-4997, 100));
    REQUIRE(to_rational(b - a) == rational(4997, 100));
    REQUIRE(to_rational(a * b) == rational(3, 2));
    REQUIRE(to_rational(-a) == rational(-3, 100));

    REQUIRE(a + b == b + a);
    REQUIRE(a * b == b * a);
    REQUIRE((a + b) + c == a + (b + c));
    REQUIRE(a - b == -(b - a));
    REQUIRE((a + b) - b == a);
    REQUIRE((a + b) * c == a * c + b * c);
    REQUIRE(a * decimal(0) == decimal(0));

    decimal d = a;
    d += b;
    REQUIRE(to_rational(d) == rational(5003, 100));
    d -= b;
    REQUIRE(d == a);
    d *= 100;
    REQUIRE(to_rational(d) == rational(3));
}

TEST_CASE("real comparisons") {
    using R = real<2>;
    // values in increasing order, with a mix of exponents and signs
    const R v[] = {R(-5, 1), R(-3, 0), R(-1, -2), R(0), R(1, -3), R(3, -2), R(1), R(3, 1), R(5, 2)};
    const int n = sizeof(v) / sizeof(v[0]);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            REQUIRE((v[i] < v[j]) == (to_rational(v[i]) < to_rational(v[j])));
            REQUIRE((v[i] == v[j]) == (i == j));
            REQUIRE((v[i] != v[j]) == (i != j));
            REQUIRE((v[i] <= v[j]) == (i <= j));
            REQUIRE((v[i] > v[j]) == (i > j));
            REQUIRE((v[i] >= v[j]) == (i >= j));
        }

    // against integers, both directions
    REQUIRE(R(3, -2) < 1);
    REQUIRE(1 > R(3, -2));
    REQUIRE(R(3, 1) > 5);
    REQUIRE(5 < R(3, 1));
    REQUIRE(R(6, 0) == 6);
    REQUIRE(6 == R(6, 0));
    REQUIRE(R(3, 1) == 6);
    REQUIRE(R(3, -1) != 6);
    REQUIRE(R(-3, 1) < -5);
    REQUIRE(-5 > R(-3, 1));
    REQUIRE(R(0) == 0);
    REQUIRE(R(1, -2) > 0);
    REQUIRE(R(-1, -2) < 0);
    REQUIRE(R(6, 0) <= 6);
    REQUIRE(R(6, 0) >= 6);

    REQUIRE(decimal(3, -2) < 1);
    REQUIRE(1 > decimal(3, -2));
    REQUIRE(decimal(5, 1) == 50);
    REQUIRE(50 == decimal(5, 1));
    REQUIRE(decimal(5, 1) > 49);
    REQUIRE(49 < decimal(5, 1));
    REQUIRE(decimal(-5, 1) < -49);
    REQUIRE(-49 > decimal(-5, 1));
    REQUIRE(decimal(1, -3) > 0);
}

TEST_CASE("real division") {
    using R = real<2>;
    // exact when the quotient is representable within 100 digits
    REQUIRE(to_rational(R(3) / R(2)) == rational(3, 2));
    REQUIRE(to_rational(R(1) / R(4)) == rational(1, 4));
    REQUIRE(to_rational(R(-9) / R(4)) == rational(-9, 4));
    REQUIRE(to_rational(R(3, -2) / R(3, 1)) == rational(1, 8));
    REQUIRE(to_rational(R(6) / 3) == rational(2));
    REQUIRE(to_rational(R(1) / 8) == rational(1, 8));
    REQUIRE(to_rational(R(-1) / 8) == rational(-1, 8));
    REQUIRE(to_rational(3 / R(2)) == rational(3, 2));

    REQUIRE(to_rational(decimal(6) / decimal(3)) == rational(2));
    REQUIRE(to_rational(decimal(1) / decimal(8)) == rational(1, 8));
    REQUIRE(to_rational(decimal(1) / 4) == rational(1, 4));
    REQUIRE(to_rational(decimal(3, -1) / 3) == rational(1, 10));
    REQUIRE(to_rational(decimal(-1) / 8) == rational(-1, 8));

    // inexact quotients are truncated, with a relative error under B**-100
    const rational third = rational(1, 3);
    const rational q2 = to_rational(R(1) / R(3));
    REQUIRE(q2 < third);
    REQUIRE(third - q2 < to_rational(R(1, -100)));
    const rational q10 = to_rational(decimal(1) / decimal(3));
    REQUIRE(q10 < third);
    REQUIRE(third - q10 < to_rational(decimal(1, -100)));

    R e = R(7);
    e /= 2;
    REQUIRE(to_rational(e) == rational(7, 2));
    decimal f = decimal(7);
    f /= 2;
    REQUIRE(to_rational(f) == rational(7, 2));
}

TEST_CASE("real shift operators") {
    using R = real<2>;
    REQUIRE(to_rational(R(3) << 4) == rational(48));
    REQUIRE(to_rational(R(3) >> 2) == rational(3, 4));
    REQUIRE(to_rational(R(3, -2) << 2) == rational(3));
    REQUIRE(to_rational(R(-3) << 1) == rational(-6));
    REQUIRE((R(3) << 4) == R(48));
    R a(5);
    a <<= 3;
    REQUIRE(to_rational(a) == rational(40));
    a >>= 4;
    REQUIRE(to_rational(a) == rational(5, 2));

    // for base 10 a shift multiplies/divides by a power of two
    REQUIRE(to_rational(decimal(5) << 3) == rational(40));
    REQUIRE(to_rational(decimal(5) >> 1) == rational(5, 2));
    REQUIRE(to_rational(decimal(1) >> 3) == rational(1, 8));
    REQUIRE(to_rational(decimal(-3) << 2) == rational(-12));
    REQUIRE(to_rational(decimal(3, -2) >> 1) == rational(3, 200));
    decimal b(5);
    b <<= 4;
    REQUIRE(to_rational(b) == rational(80));
    b >>= 2;
    REQUIRE(to_rational(b) == rational(20));
}

TEST_CASE("decimal format") {
    REQUIRE(format("{}", decimal(0)) == "0");
    REQUIRE(format("{}", decimal(7)) == "7");
    REQUIRE(format("{}", decimal(-7)) == "-7");
    REQUIRE(format("{}", decimal(7, 3)) == "7000");
    REQUIRE(format("{}", decimal(-7, 3)) == "-7000");
    REQUIRE(format("{}", decimal(5, -1)) == "0.5");
    REQUIRE(format("{}", decimal(-5, -1)) == "-0.5");
    REQUIRE(format("{}", decimal(1, -5)) == "0.00001");
    REQUIRE(format("{}", decimal(-1, -5)) == "-0.00001");
    REQUIRE(format("{}", decimal(123456, -3)) == "123.456");
    REQUIRE(format("{}", decimal(-123456, -3)) == "-123.456");
    REQUIRE(format("{}", decimal(123456, -6)) == "0.123456");
    REQUIRE(format("{}", decimal(123456, -7)) == "0.0123456");
    REQUIRE(format("{}", decimal(rational(1, 8))) == "0.125");
    REQUIRE(format("{}", decimal(rational(-1, 8))) == "-0.125");
}

TEST_CASE("str and ostream") {
    REQUIRE(decimal(123456, -3).str() == "123.456");
    REQUIRE(real<2>(3, -2).str() == "3/4");

    std::ostringstream os;
    os << decimal(5, -1) << ' ' << real<2>(3, -2) << ' ' << real<2>(-6);
    REQUIRE(os.str() == "0.5 3/4 -6");
}

TEST_CASE("REAL_FRACT_DIGITS") {
    struct Guard {
        ~Guard() { REAL_FRACT_DIGITS = std::nullopt; }
    } guard;

    REQUIRE(format("{}", real<2>(3, -2)) == "3/4"); // unset: exact rational form
    REAL_FRACT_DIGITS = 3;
    REQUIRE(format("{}", real<2>(3, -2)) == "0.750");
    REQUIRE(format("{}", real<2>(-3, -2)) == "-0.750");
    REQUIRE(format("{}", real<2>(5)) == "5.000");
    REQUIRE(real<2>(3, -2).str() == "0.750");
    REQUIRE(format("{:.1}", real<2>(3, -2)) == "0.8"); // an explicit precision wins
    REAL_FRACT_DIGITS = 0;
    REQUIRE(format("{}", real<2>(3, -2)) == "1"); // no separator when there are no fraction digits
    REAL_FRACT_DIGITS = std::nullopt;
    REQUIRE(format("{}", real<2>(3, -2)) == "3/4");
}

TEST_CASE("round") {
    // exact when the value fits in the requested number of digits
    REQUIRE(to_rational(real<2>::round(rational(3, 8), 3)) == rational(3, 8));
    REQUIRE(to_rational(real<2>::round(rational(3, 8), 10)) == rational(3, 8));
    REQUIRE(to_rational(real<2>::round(rational(-3, 8), 5)) == rational(-3, 8));
    REQUIRE(to_rational(real<2>::round(rational(0), 5)) == rational(0));
    REQUIRE(to_rational(decimal::round(rational(1, 4), 2)) == rational(1, 4));
    REQUIRE(to_rational(decimal::round(rational(-1, 4), 4)) == rational(-1, 4));
    REQUIRE(to_rational(decimal::round(rational(123), 0)) == rational(123));

    // otherwise a multiple of B**-digits within one unit of the last place
    for (int digits : {1, 5, 20}) {
        const rational ulp2 = to_rational(real<2>(1, -digits));
        const rational ulp10 = to_rational(decimal(1, -digits));
        for (const rational& a : {rational(1, 3), rational(-1, 3), rational(22, 7), rational(-22, 7)}) {
            const real<2> r2 = real<2>::round(a, digits);
            REQUIRE(r2.exp >= -digits);
            REQUIRE(abs(to_rational(r2) - a) < ulp2);
            const decimal r10 = decimal::round(a, digits);
            REQUIRE(r10.exp >= -digits);
            REQUIRE(abs(to_rational(r10) - a) < ulp10);
        }
    }
}

TEST_CASE("hash of equal values") {
    std::hash<real<2>> h2;
    REQUIRE(h2(real<2>(6)) == h2(real<2>(3, 1)));
    REQUIRE(h2(real<2>(3, -2)) == h2(real<2>(12, -4)));
    std::hash<decimal> h10;
    REQUIRE(h10(decimal(50)) == h10(decimal(5, 1)));

    std::unordered_map<decimal, int> m;
    m[decimal(5, 1)] = 1;
    m[decimal(50)] += 1;
    m[decimal(1, -3)] = 3;
    REQUIRE(m.size() == 2);
    REQUIRE(m[decimal(50)] == 2);
}

TEST_CASE("to_rational") {
    REQUIRE(to_rational(real<2>(0)) == rational(0));
    REQUIRE(to_rational(real<2>(3, 5)) == rational(96));
    REQUIRE(to_rational(real<2>(3, -5)) == rational(3, 32));
    REQUIRE(to_rational(real<2>(-3, -5)) == rational(-3, 32));
    REQUIRE(to_rational(decimal(-7)) == rational(-7));
    REQUIRE(to_rational(decimal(7, 2)) == rational(700));
    REQUIRE(to_rational(decimal(7, -2)) == rational(7, 100));

    // round trip through rational
    for (int num : {-9, -1, 0, 1, 3, 40}) {
        for (int e : {-4, -1, 0, 2}) {
            const real<2> a(num, e);
            REQUIRE(real<2>(to_rational(a)) == a);
            const decimal b(num, e);
            REQUIRE(decimal(to_rational(b)) == b);
        }
    }
}

TEST_CASE("random cross check against rational") {
    Random rng(12345);
    for (int i = 0; i < 300; i++) {
        const int64_t na = rng.Uniform<int64_t>(-1000000, 1000000);
        const int64_t nb = rng.Uniform<int64_t>(-1000000, 1000000);
        const int ea = rng.Uniform<int>(-30, 30);
        const int eb = rng.Uniform<int>(-30, 30);
        // multi word numerators
        const integer ba = integer(na) << rng.Uniform<int>(0, 200);
        const integer bb = integer(nb) << rng.Uniform<int>(0, 200);

        const real<2> a(ba, ea), b(bb, eb);
        const rational ra = to_rational(a), rb = to_rational(b);
        REQUIRE(to_rational(-a) == -ra);
        REQUIRE(to_rational(a + b) == ra + rb);
        REQUIRE(to_rational(a - b) == ra - rb);
        REQUIRE(to_rational(a * b) == ra * rb);
        // note: the builtin operand cases are in #38, with the fix they need
        REQUIRE(to_rational(a * nb) == ra * rational(nb));
        REQUIRE((a < b) == (ra < rb));
        REQUIRE((a == b) == (ra == rb));
        REQUIRE((a < nb) == (ra < rational(nb)));
        REQUIRE((nb < a) == (rational(nb) < ra));
        REQUIRE((a == nb) == (ra == rational(nb)));
        REQUIRE(a.num % 2 != 0);

        const decimal c(ba, ea), d(bb, eb);
        const rational rc = to_rational(c), rd = to_rational(d);
        REQUIRE(to_rational(-c) == -rc);
        REQUIRE(to_rational(c + d) == rc + rd);
        REQUIRE(to_rational(c - d) == rc - rd);
        REQUIRE(to_rational(c * d) == rc * rd);
        REQUIRE((c < d) == (rc < rd));
        REQUIRE((c == d) == (rc == rd));
        REQUIRE((c < nb) == (rc < rational(nb)));
        REQUIRE((nb < c) == (rational(nb) < rc));
        REQUIRE(c.num % 10 != 0);
        REQUIRE(decimal(rc) == c);
        REQUIRE(real<2>(ra) == a);
    }
}

TEST_CASE("inexact conversion from rational") {
    REQUIRE_THROWS_AS(real<2>(rational(1, 3)), std::runtime_error);
    REQUIRE_THROWS_AS(real<2>(rational(1, 10)), std::runtime_error);
    REQUIRE_THROWS_AS(decimal(rational(1, 3)), std::runtime_error);
    REQUIRE_THROWS_AS(decimal(rational(1, 6)), std::runtime_error);

    // denominators built only from the base's prime factors are exact
    REQUIRE(to_rational(decimal(rational(1, 20))) == rational(1, 20));  // more twos than fives
    REQUIRE(to_rational(decimal(rational(1, 8))) == rational(1, 8));    // only twos
    REQUIRE(to_rational(decimal(rational(1, 25))) == rational(1, 25));  // more fives than twos
    REQUIRE(to_rational(decimal(rational(1, 1000))) == rational(1, 1000));
    REQUIRE(to_rational(decimal(rational(-3, 40))) == rational(-3, 40));
    REQUIRE(decimal(rational(1, 20)).num == 5);
    REQUIRE(decimal(rational(1, 20)).exp == -2);
    REQUIRE(decimal(rational(1, 25)).num == 4);
    REQUIRE(decimal(rational(1, 25)).exp == -2);
}
