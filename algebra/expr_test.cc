#include "algebra/expr.h"
#include "algebra/__test.h"

TEST_CASE("basic") {
    REQUIRE((sqrt(3_e) - sqrt(3_e))->sign() == 0);
    REQUIRE((sqrt(3_e) - sqrt(2_e))->sign() == 1);
    REQUIRE((sqrt(2_e) - sqrt(5_e))->sign() == -1);

    REQUIRE(sqrt(2_e) < sqrt(5_e));

    REQUIRE(2 * sqrt(2_e) == sqrt(8_e));
    REQUIRE(sqrt(2_e) > 1);
    REQUIRE(sqrt(2_e) < 2);

    REQUIRE(1 + sqrt(2_e) + 1 < 2 + sqrt(3_e));

    REQUIRE(sqrt(2_e) + sqrt(3_e) >= sqrt(5_e));

    REQUIRE(sqrt(125_e) == 5 * sqrt(5_e));

    REQUIRE(E_EXPR < 3);
    REQUIRE(E_EXPR > 2);

    REQUIRE(sin(sqrt(2_e)) > -2);
    //REQUIRE(pow(sin(E_EXPR), 2) >= 0);
    //REQUIRE(pow(sin(E_EXPR), 2) + pow(cos(E_EXPR), 2) == 1);
    //REQUIRE(sqrt(2_e) + sqrt(3_e) + sqrt(4_e) >= sqrt(5_e));
}

TEST_CASE("str") {
    REQUIRE(format("{}", sqrt(4_e)) == "2");
    REQUIRE(format("{}", sqrt(5_e)) == "sqrt(5)");
    REQUIRE(format("{}", sqrt(125_e)) == "sqrt(125)");
    REQUIRE(format("{}", pow(sqrt(5_e) + sqrt(7_e), 2)) == "12 + 2*sqrt(5)*sqrt(7)");
    REQUIRE(format("{}", sqrt(5_e) - 2) == "sqrt(5) - 2");
    REQUIRE(format("{}", sqrt(5_e) * -1) == "-sqrt(5)");
    REQUIRE(format("{}", sqrt(2_e) * sqrt(3_e) * sqrt(2_e)) == "2*sqrt(3)");
    REQUIRE(format("{}", sqrt(2_e) + sqrt(3_e) + sqrt(2_e)) == "2*sqrt(2) + sqrt(3)");
    REQUIRE(format("{}", 1 + E_EXPR + PI_EXPR) == "1 + e + π");
    REQUIRE(format("{}", E_EXPR - E_EXPR) == "0");
    REQUIRE(format("{}", 1 + E_EXPR + sqrt(2_e) - E_EXPR) == "1 + sqrt(2)");
}

#if 0
TEST_CASE("sqrt") {
    auto a = 2_q;
    rational lower = 1;
    rational upper = (1 + a) / 2;
    for (int i = 0; i < 10; i++) {
        REQUIRE(lower * lower < a);
        REQUIRE(a < upper * upper);
        info("{:.30} - {:.30}", lower, upper);
        while (true) {
            auto e = (lower + upper) / 2;
            if (e * e >= a)
                break;
            lower = e;
        }
    }
}
#endif

constexpr std::string os() {
    std::ostringstream s;
    s << -15_e;
    return s.str();
}

TEST_CASE("constexpr ostream") {
    REQUIRE(os() == "-15");
}

TEST_CASE("structural simplification") {
    // structurally equal but distinct nodes have to simplify too
    REQUIRE(format("{}", sqrt(2_e) + sqrt(2_e)) == "2*sqrt(2)");
    REQUIRE(format("{}", sqrt(2_e) * sqrt(2_e)) == "2");
    REQUIRE(format("{}", sin(2_e) + sin(2_e)) == "2*sin(2)");
    REQUIRE(format("{}", sin(2_e) - sin(2_e)) == "0");
    REQUIRE(format("{}", cos(2_e) - cos(2_e)) == "0");
    REQUIRE(format("{}", sin(2_e) - sin(3_e)) == "sin(2) - sin(3)");
    REQUIRE(format("{}", PI_EXPR + PI_EXPR) == "2*π");
    REQUIRE(format("{}", 2 * sqrt(2_e) + sqrt(2_e)) == "3*sqrt(2)");
}

// expr_var is an aggregate with a base class, so make_shared() cannot forward the name to it
static expr_ptr make_var(std::string name) {
    auto v = std::make_shared<expr_var>();
    v->name = std::move(name);
    return v;
}

TEST_CASE("node predicates and accessors") {
    using namespace algebra::literals;

    // make_rational collapses an integral rational into an expr_integer node
    REQUIRE(is_integer(make_rational(rational(4))));
    REQUIRE(!is_rational(sqrt(2_e)));
    REQUIRE(is_rational(make_rational(rational(1, 2))));
    REQUIRE(!is_integer(make_rational(rational(1, 2))));
    REQUIRE(integer_value(make_integer(7)) == 7);
    REQUIRE(rational_value(make_rational(rational(1, 2))) == rational(1, 2));
    // rational_value also accepts an integer node
    REQUIRE(rational_value(make_integer(7)) == 7);

    // sqrt and cbrt are powers with exponent 1/2 and 1/3
    const expr_ptr r5 = sqrt(5_e);
    REQUIRE(is_power(r5));
    REQUIRE(is_sqrt(r5));
    REQUIRE(!is_cbrt(r5));
    REQUIRE(power_exp(r5) == rational(1, 2));
    REQUIRE(is_integer(power_base(r5)));
    REQUIRE(integer_value(power_base(r5)) == 5);

    const expr_ptr c5 = pow(5_e, rational(1, 3));
    REQUIRE(is_cbrt(c5));
    REQUIRE(!is_sqrt(c5));

    // a perfect square is folded, so it is no longer a power
    REQUIRE(!is_power(sqrt(4_e)));
    REQUIRE(is_integer(sqrt(4_e)));

    // sum
    const expr_ptr s = sqrt(2_e) + sqrt(3_e);
    REQUIRE(is_sum(s));
    REQUIRE(sum_values(s).size() == 2);
    REQUIRE(!is_sum(r5));

    // product, and the rational * x shape
    const expr_ptr p = 3 * sqrt(2_e);
    REQUIRE(is_product(p));
    REQUIRE(is_product_rx(p));
    REQUIRE(product_values(p).size() == 2);
    REQUIRE(is_rational(product_values(p)[0]));
    REQUIRE(rational_value(product_values(p)[0]) == 3);
    REQUIRE(!is_product_rx(s));
    REQUIRE(!is_product(s));
}

TEST_CASE("identical") {
    using namespace algebra::literals;

    // the same node
    const expr_ptr a = sqrt(2_e);
    REQUIRE(identical(a, a));
    // distinct nodes with equal value
    REQUIRE(identical(sqrt(2_e), sqrt(2_e)));
    REQUIRE(identical(make_integer(3), make_integer(3)));
    REQUIRE(identical(make_rational(rational(1, 2)), make_rational(rational(1, 2))));
    REQUIRE(!identical(make_integer(3), make_integer(4)));
    REQUIRE(!identical(make_rational(rational(1, 2)), make_rational(rational(1, 3))));
    // different node kinds are never identical
    REQUIRE(!identical(make_integer(3), sqrt(3_e)));
    REQUIRE(!identical(E_EXPR, PI_EXPR));
    REQUIRE(identical(E_EXPR, std::make_shared<expr_e>()));
    REQUIRE(identical(PI_EXPR, std::make_shared<expr_pi>()));
    // nested
    REQUIRE(identical(sqrt(2_e) + sqrt(3_e), sqrt(2_e) + sqrt(3_e)));
    REQUIRE(!identical(sqrt(2_e) + sqrt(3_e), sqrt(2_e) + sqrt(5_e)));
    REQUIRE(identical(sin(2_e), sin(2_e)));
    REQUIRE(!identical(sin(2_e), sin(3_e)));
    REQUIRE(!identical(sin(2_e), cos(2_e)));
    REQUIRE(identical(cos(2_e), cos(2_e)));
    // a variable compares by name
    REQUIRE(identical(make_var("x"), make_var("x")));
    REQUIRE(!identical(make_var("x"), make_var("y")));
}

TEST_CASE("interval arithmetic") {
    using I = interval<rational>;
    auto eq = [](const I& a, const rational& lo, const rational& hi) { return a.min == lo && a.max == hi; };

    REQUIRE(eq(-I{2, 5}, -5, -2));
    REQUIRE(eq(I{1, 2} + I{10, 20}, 11, 22));
    REQUIRE(eq(I{-3, 2} + I{-1, 1}, -4, 3));

    // multiplication has to consider all four products, not just the corners in order
    REQUIRE(eq(I{2, 3} * I{4, 5}, 8, 15));
    REQUIRE(eq(I{-3, -2} * I{-5, -4}, 8, 15));
    REQUIRE(eq(I{-2, 3} * I{-5, 4}, -15, 12));
    REQUIRE(eq(I{-2, 3} * I{2, 4}, -8, 12));
    REQUIRE(eq(I{0, 0} * I{-5, 4}, 0, 0));

    // pow by repeated multiplication, so an even power is a sound over-approximation
    REQUIRE(eq(pow(I{2, 3}, 2), 4, 9));
    REQUIRE(eq(pow(I{2, 3}, 0), 1, 1));
    REQUIRE(eq(pow(I{2, 3}, 3), 8, 27));
    REQUIRE(eq(pow(I{-3, 2}, 2), -6, 9)); // the true range is [0, 9]; wider is still correct
    REQUIRE(eq(pow(I{2, 4}, -1), rational(1, 4), rational(1, 2)));
    REQUIRE(eq(pow(I{-4, -2}, -1), rational(-1, 2), rational(-1, 4)));
    // a negative power of an interval straddling zero is a division by zero
    REQUIRE_THROWS(pow(I{-1, 1}, -1));
}

TEST_CASE("bounds") {
    using namespace algebra::literals;

    auto b = [](expr_ptr e) { auto o = bounds(e); REQUIRE(o != std::nullopt); return *o; };

    // a rational is a point interval
    REQUIRE(b(make_rational(rational(3, 2))).min == rational(3, 2));
    REQUIRE(b(make_rational(rational(3, 2))).max == rational(3, 2));

    // the two constants must bracket their true values
    REQUIRE(b(E_EXPR).min < rational(271828, 100000));
    REQUIRE(b(E_EXPR).max > rational(271829, 100000));
    REQUIRE(b(PI_EXPR).min < rational(314159, 100000));
    REQUIRE(b(PI_EXPR).max > rational(314160, 100000));

    // sin and cos are bounded by [-1, 1] without knowing the argument
    REQUIRE(b(sin(sqrt(2_e))).min == -1);
    REQUIRE(b(sin(sqrt(2_e))).max == 1);
    REQUIRE(b(cos(E_EXPR)).min == -1);
    REQUIRE(b(cos(E_EXPR)).max == 1);

    // a sum adds the bounds, a negation flips them
    REQUIRE(b(1 + sin(2_e)).min == 0);
    REQUIRE(b(1 + sin(2_e)).max == 2);
    REQUIRE(b(-sin(2_e)).min == -1);
    // a product multiplies them
    REQUIRE(b(2 * sin(2_e)).min == -2);
    REQUIRE(b(2 * sin(2_e)).max == 2);

    // an unbounded node has no bounds
    REQUIRE(bounds(make_var("x")) == std::nullopt);
    // and neither does a sum containing one
    REQUIRE(bounds(1 + make_var("x")) == std::nullopt);
    // a non integer exponent is not handled
    REQUIRE(bounds(sqrt(2_e)) == std::nullopt);
}

TEST_CASE("safe_sign") {
    using namespace algebra::literals;
    REQUIRE(safe_sign(make_integer(5)) == 1);
    REQUIRE(safe_sign(make_integer(-5)) == -1);
    REQUIRE(safe_sign(make_integer(0)) == 0);
    REQUIRE(safe_sign(E_EXPR) == 1);
    REQUIRE(safe_sign(PI_EXPR) == 1);
    // the sign of sin and cos is not known without evaluating the argument
    REQUIRE(safe_sign(sin(2_e)) == std::nullopt);
    REQUIRE(safe_sign(cos(2_e)) == std::nullopt);
    // a sum whose terms all have a known sign
    REQUIRE(safe_sign(sqrt(2_e) + sqrt(3_e)) == 1);
    // a sum that needs bounds to decide: 2 + sin(x) is positive whatever sin is
    REQUIRE(safe_sign(2 + sin(2_e)) == 1);
    REQUIRE(safe_sign(-2 + sin(2_e)) == -1);
}

TEST_CASE("sign of a power of a negative base") {
    using namespace algebra::literals;
    // a negative base that is not a rational, so pow() keeps it as an expr_power node
    const expr_ptr neg = -sqrt(2_e);
    REQUIRE(neg->sign() == -1);

    // an odd exponent keeps the sign of the base, an even one makes the power positive
    REQUIRE(pow(neg, 3)->sign() == -1);
    REQUIRE(pow(neg, 5)->sign() == -1);
    REQUIRE(pow(neg, 2)->sign() == 1);
    REQUIRE(pow(neg, 4)->sign() == 1);
    REQUIRE(safe_sign(pow(neg, 3)) == -1);

    // and the comparisons that are built on the sign of a difference
    REQUIRE(pow(neg, 3) < 0);
    REQUIRE(pow(neg, 3) != 0);
    REQUIRE(pow(neg, 3) < pow(neg, 2));
    REQUIRE(pow(neg, 3) + pow(neg, 2) < pow(neg, 2));

    // a fractional exponent of a negative base is not a real number, and stays an error
    REQUIRE_THROWS(pow(neg, rational(1, 2))->sign());
}

TEST_CASE("expr_var has no sign") {
    REQUIRE_THROWS(make_var("x")->sign());
    REQUIRE(safe_sign(make_var("x")) == std::nullopt);
}

TEST_CASE("expressions containing a variable can be built") {
    const expr_ptr x = make_var("x");

    // the simplification code asks for signs through safe_sign, so building these only works if
    // an unknown sign is reported as unknown_sign_error rather than a plain std::runtime_error
    const expr_ptr s = 1 + x;
    REQUIRE(is_sum(s));
    REQUIRE(sum_values(s).size() == 2);

    const expr_ptr p = 2 * x;
    REQUIRE(is_product(p));
    REQUIRE(is_product_rx(p));

    REQUIRE(is_sum(x + x + 1));
    REQUIRE(is_power(pow(x, rational(1, 2))));

    // the sign still cannot be decided, and that is reported as unknown rather than escaping
    REQUIRE_THROWS_AS(x->sign(), unknown_sign_error);
    REQUIRE_THROWS_AS(s->sign(), unknown_sign_error);
    REQUIRE(safe_sign(s) == std::nullopt);
}

TEST_CASE("sign of a matrix") {
    // a 1x1 matrix is its single element, and every other shape has no sign of its own
    expr_matrix one;
    one.rows = 1;
    one.cols = 1;
    one.data.push_back(make_integer(-5));
    REQUIRE(one.sign() == -1);
    REQUIRE(safe_sign(std::make_shared<expr_matrix>(one)) == -1);

    expr_matrix empty;
    empty.rows = 0;
    empty.cols = 0;
    REQUIRE_THROWS_AS(empty.sign(), unknown_sign_error);

    expr_matrix row;
    row.rows = 1;
    row.cols = 2;
    row.data.push_back(make_integer(1));
    row.data.push_back(make_integer(2));
    REQUIRE_THROWS_AS(row.sign(), unknown_sign_error);
}

TEST_CASE("formatting into a counted buffer") {
    using namespace algebra::literals;
    // These formatters write with std::format_to(ctx.out(), ...) and drop the iterator it returns,
    // which looks like it should lose the position. It does not: the iterator a format context hands
    // out is a proxy for a sink that holds the position, and every copy shares it. This pins that,
    // so formatting into a counted buffer keeps working if the formatters are ever rewritten.
    char buf[64] = {};
    const auto r = std::format_to_n(buf, sizeof(buf) - 1, "{}", sqrt(5_e) + 2);
    REQUIRE(std::string(buf, r.out) == std::format("{}", sqrt(5_e) + 2));
    REQUIRE(r.size == static_cast<std::ptrdiff_t>(std::format("{}", sqrt(5_e) + 2).size()));

    // a nested expression, so format_expression recurses
    const expr_ptr deep = pow(sqrt(2_e) + sqrt(3_e), 2) * cos(PI_EXPR) - 4;
    char buf2[128] = {};
    const auto r2 = std::format_to_n(buf2, sizeof(buf2) - 1, "{}", deep);
    REQUIRE(std::string(buf2, r2.out) == std::format("{}", deep));
}

TEST_CASE("roots do not distribute over an operand of unknown sign") {
    using namespace algebra::literals;
    const expr_ptr x = make_var("x");
    const expr_ptr y = make_var("y");

    // sqrt(x^2) is |x|, not x, so the exponents must not be folded for a base that may be negative
    REQUIRE(!identical(sqrt(pow(x, 2)), x));
    REQUIRE(is_power(sqrt(pow(x, 2))));
    REQUIRE(is_power(power_base(sqrt(pow(x, 2)))));
    // and sqrt(x*y) is not sqrt(x)*sqrt(y)
    REQUIRE(!is_product(sqrt(x * y)));
    REQUIRE(is_power(sqrt(x * y)));

    // an integer exponent is sound for either sign, and still folds
    REQUIRE(identical(pow(pow(x, 2), rational(3)), pow(x, 6)));
    REQUIRE(is_product(pow(x * y, rational(2))));

    // a base known to be non negative folds as before
    REQUIRE(identical(sqrt(sqrt(2_e)), pow(2_e, rational(1, 4))));
    REQUIRE(identical(sqrt(pow(PI_EXPR, 2)), PI_EXPR));
    REQUIRE(format("{}", sqrt(2_e * PI_EXPR)) == "sqrt(2)*sqrt(π)");
    REQUIRE(sqrt(2_e) * sqrt(3_e) == sqrt(6_e));
}
