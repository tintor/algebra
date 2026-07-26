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

TEST_CASE("pow of real<2>") {
    using R = real<2>;
    // negative base
    REQUIRE(to_rational(pow(R(-2), 0)) == rational(1));
    REQUIRE(to_rational(pow(R(-2), 1)) == rational(-2));
    REQUIRE(to_rational(pow(R(-2), 2)) == rational(4));
    REQUIRE(to_rational(pow(R(-2), 3)) == rational(-8));
    REQUIRE(to_rational(pow(R(-3), 5)) == rational(-243));
    REQUIRE(to_rational(pow(R(-2), -3)) == rational(-1, 8));

    // fractional base
    REQUIRE(to_rational(pow(R(3, -2), 2)) == rational(9, 16));
    REQUIRE(to_rational(pow(R(3, -2), 3)) == rational(27, 64));
    REQUIRE(to_rational(pow(R(1, -1), 10)) == rational(1, 1024));
    REQUIRE(to_rational(pow(R(1, -1), -3)) == rational(8));

    // 16/9 is not dyadic, so a negative exponent is only approximated (by a truncating divide)
    REQUIRE(abs(to_rational(pow(R(3, -2), -2)) - rational(16, 9)) < to_rational(R(1, -90)));

    // zero and one
    REQUIRE(to_rational(pow(R(0), 3)) == rational(0));
    REQUIRE(to_rational(pow(R(0), 0)) == rational(1));
    REQUIRE(to_rational(pow(R(1), 1000)) == rational(1));

    // pow(base, exp) == base * pow(base, exp - 1)
    for (int e = 0; e < 8; e++) {
        REQUIRE(pow(R(3, -1), e + 1) == pow(R(3, -1), e) * R(3, -1));
        REQUIRE(pow(R(-5), e + 1) == pow(R(-5), e) * R(-5));
    }

    // the exponent of the base itself is folded into real::exp, without overflowing it
    REQUIRE(pow(R(2), 40).exp == 40);
    REQUIRE(pow(R(2), 40, R(3)).exp == 40);
    REQUIRE(to_rational(pow(R(2), 40)) == rational(exp2(40)));
    REQUIRE_THROWS_AS(pow(R(2), int64_t(1) << 40), std::runtime_error);
    REQUIRE_THROWS_AS(pow(R(2), int64_t(1) << 32), std::runtime_error);
    REQUIRE_THROWS_AS(pow(R(2), int64_t(1) << 31), std::runtime_error);
    REQUIRE_THROWS_AS(pow(R(2), 1 << 30, R(1, 1 << 30)), std::runtime_error);
    REQUIRE_THROWS_AS(pow(decimal(10), int64_t(1) << 40), std::runtime_error);
    // just below the limit the exponent still fits
    REQUIRE(pow(R(2), int64_t(1) << 30).exp == 1 << 30);
}

TEST_CASE("pow of decimal") {
    REQUIRE(to_rational(pow(decimal(10), 3)) == rational(1000));
    REQUIRE(pow(decimal(10), 3).num == 1);
    REQUIRE(pow(decimal(10), 3).exp == 3);
    REQUIRE(to_rational(pow(decimal(10), -2)) == rational(1, 100));
    REQUIRE(to_rational(pow(decimal(10), 2, decimal(3))) == rational(300));
    REQUIRE(to_rational(pow(decimal(3), 4)) == rational(81));
    REQUIRE(to_rational(pow(decimal(-3), 3)) == rational(-27));
    REQUIRE(to_rational(pow(decimal(5, -1), 3)) == rational(1, 8));
    REQUIRE(to_rational(pow(decimal(2), -3)) == rational(1, 8));
}

TEST_CASE("abs") {
    REQUIRE(to_rational(abs(real<2>(3, -2))) == rational(3, 4));
    REQUIRE(to_rational(abs(real<2>(-3, -2))) == rational(3, 4));
    REQUIRE(to_rational(abs(real<2>(-6))) == rational(6));
    REQUIRE(to_rational(abs(real<2>(0))) == rational(0));
    REQUIRE(to_rational(abs(decimal(-123456, -3))) == rational(123456, 1000));
    REQUIRE(abs(decimal(-5, 1)) == decimal(5, 1));

    // abs() does not modify its argument, and abs(-a) == abs(a)
    const real<2> a(-7, -3);
    REQUIRE(abs(a) == abs(-a));
    REQUIRE(a.num == -7);
    REQUIRE(abs(a) > 0);
}

TEST_CASE("division by zero") {
    REQUIRE_THROWS_AS(real<2>(1) / real<2>(0), std::runtime_error);
    REQUIRE_THROWS_AS(real<2>(1) / 0, std::runtime_error);
    REQUIRE_THROWS_AS(decimal(1) / decimal(0), std::runtime_error);
}

TEST_CASE("conversion from floating point") {
    // float and double values are dyadic, so real<2> holds them exactly
    REQUIRE(to_rational(real<2>(0.375)) == rational(3, 8));
    REQUIRE(to_rational(real<2>(-0.375)) == rational(-3, 8));
    REQUIRE(to_rational(real<2>(0.375f)) == rational(3, 8));
    REQUIRE(to_rational(real<2>(-0.375f)) == rational(-3, 8));
    REQUIRE(to_rational(real<2>(0.0)) == rational(0));
    REQUIRE(to_rational(real<2>(1024.0)) == rational(1024));
    REQUIRE(real<2>(0.375).num == 3);
    REQUIRE(real<2>(0.375).exp == -3);

    // 0.1 is not exactly one tenth
    REQUIRE(to_rational(real<2>(0.1)) == rational(0.1));
    REQUIRE(to_rational(real<2>(0.1)) != rational(1, 10));
    REQUIRE(to_rational(real<2>(0.1)) > rational(1, 10));

    REQUIRE(to_rational(decimal(0.375)) == rational(3, 8));
    REQUIRE(format("{}", decimal(0.375)) == "0.375");
    REQUIRE(format("{}", decimal(-0.25)) == "-0.25");
    REQUIRE(to_rational(decimal(0.1)) == rational(0.1));
}
