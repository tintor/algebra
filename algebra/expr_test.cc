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
