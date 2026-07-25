// Compiles the portable (non x86-64) implementation of __divq() and friends.
#define ALGEBRA_NO_ASM
#include "algebra/util.h"
#include "algebra/__test.h"

TEST_CASE("portable __divq") {
    static_assert(!std::is_void_v<decltype(__divq_mod(uint128_t(1), 1))>);

    uint64_t q = 0, r = 0;
    __divq(static_cast<uint128_t>(100), 7, q, r);
    REQUIRE(q == 14);
    REQUIRE(r == 2);
    REQUIRE(__divq(static_cast<uint128_t>(100), 7) == 14);
    REQUIRE(__divq_mod(static_cast<uint128_t>(100), 7) == 2);

    __divq(static_cast<uint128_t>(0), 3, q, r);
    REQUIRE(q == 0);
    REQUIRE(r == 0);

    std::mt19937_64 rng(4);
    for (int i = 0; i < 10000; i++) {
        const uint64_t b = rng() | 1;
        // keep the quotient below 2**64, which is what the callers guarantee
        const uint128_t a = (static_cast<uint128_t>(rng() % b) << 64) | rng();
        __divq(a, b, q, r);
        REQUIRE(q == static_cast<uint64_t>(a / b));
        REQUIRE(r == static_cast<uint64_t>(a % b));
        REQUIRE(__divq(a, b) == q);
        REQUIRE(__divq_mod(a, b) == r);
    }
}
