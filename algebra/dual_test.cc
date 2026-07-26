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

using D = dual<double>;

static bool near(double a, double b) { return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(b)); }
static bool near(const D& a, double real, double dual) { return near(a.real, real) && near(a.dual, dual); }

// A dual number carries a value and its derivative, so evaluating f at {x, 1} must give
// {f(x), f'(x)}. Each case below states f' independently of the implementation.
TEST_CASE("dual arithmetic") {
    const D x{3, 1}; // the variable x at x == 3
    const D c{5, 0}; // a constant

    REQUIRE(near(x + c, 8, 1));
    REQUIRE(near(x - c, -2, 1));
    REQUIRE(near(c - x, 2, -1));

    // (x * x)' == 2x
    REQUIRE(near(x * x, 9, 6));
    // (5x)' == 5
    REQUIRE(near(c * x, 15, 5));
    // (x / 5)' == 1/5
    REQUIRE(near(x / c, 0.6, 0.2));
    // (5 / x)' == -5 / x^2
    REQUIRE(near(c / x, 5.0 / 3, -5.0 / 9));
    // (x / x)' == 0
    REQUIRE(near(x / x, 1, 0));
}

TEST_CASE("dual functions") {
    const D x{4, 1};

    // sqrt' == 1 / (2 sqrt(x))
    REQUIRE(near(sqrt(x), 2, 0.25));
    // (x^3)' == 3 x^2
    REQUIRE(near(pow(x, 3.0), 64, 48));
    // (x^0)' == 0
    REQUIRE(near(pow(x, 0.0), 1, 0));
    // exp' == exp
    REQUIRE(near(exp(D{0, 1}), 1, 1));
    REQUIRE(near(exp(x), std::exp(4.0), std::exp(4.0)));
    // log' == 1/x
    REQUIRE(near(log(x), std::log(4.0), 0.25));
    // abs' == signum
    REQUIRE(near(abs(x), 4, 1));
    REQUIRE(near(abs(D{-4, 1}), 4, -1));
    REQUIRE(near(abs(D{0, 1}), 0, 0)); // not differentiable at 0, reported as 0
}

TEST_CASE("dual trigonometry") {
    const D z{0, 1};
    // sin' == cos, cos' == -sin
    REQUIRE(near(sin(z), 0, 1));
    REQUIRE(near(cos(z), 1, 0));
    const D h{0.5, 1};
    REQUIRE(near(sin(h), std::sin(0.5), std::cos(0.5)));
    REQUIRE(near(cos(h), std::cos(0.5), -std::sin(0.5)));
    // tan' == 1 / cos^2
    REQUIRE(near(tan(h), std::tan(0.5), 1 / (std::cos(0.5) * std::cos(0.5))));
    // atan' == 1 / (1 + x^2)
    REQUIRE(near(atan(z), 0, 1));
    REQUIRE(near(atan(h), std::atan(0.5), 1 / 1.25));
}

TEST_CASE("dual chain rule") {
    // f(x) = sin(x^2) at x == 2, f'(x) = 2x cos(x^2)
    const D x{2, 1};
    const D f = sin(x * x);
    REQUIRE(near(f, std::sin(4.0), 4 * std::cos(4.0)));

    // f(x) = exp(2x) / x at x == 1, f'(x) = exp(2x) * (2x - 1) / x^2
    const D y{1, 1};
    const D g = exp(y * D{2, 0}) / y;
    REQUIRE(near(g, std::exp(2.0), std::exp(2.0)));
}
