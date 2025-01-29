#include "algebra/natural.h"
#include "algebra/integer_func.h"
#include "algebra/__test.h"

TEST_CASE("uniform_sample") {
    integer a;
    a.abs = pow(2_n, 128) - 1;
    std::mt19937_64 rng(0);
    for (int i = 0; i < 20; i++) {
        integer m = uniform_sample(0, a, rng);
        REQUIRE(m.sign() <= a.sign());
        REQUIRE(0 <= m);
        REQUIRE(m <= a);
        REQUIRE(0 <= m.sign());
        REQUIRE(m.sign() <= 2);
        while (m.sign())
            div(m, static_cast<long>(10), /*out*/m);
    }
}

TEST_CASE("abs") {
    REQUIRE(abs(0_i) == 0);
    REQUIRE(abs(10_i) == 10);
    REQUIRE(abs(-3_i) == 3);
}

TEST_CASE("pow") {
    REQUIRE(pow(2_i, 3) == 8);
    REQUIRE(pow(10_i, 30) == 1000000000000000000000000000000_i);
}
