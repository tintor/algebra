#include "algebra/kernels.h"
#include "algebra/__test.h"

TEST_CASE("__diff - simple") {
    uint64_t a[] = {75};
    uint64_t b[] = {4};
    REQUIRE(!__diff({a, 1}, {b, 1}));
    REQUIRE(a[0] == 71);
}

TEST_CASE("__diff - negative 1 word") {
    uint64_t a[] = {10};
    uint64_t b[] = {23};
    REQUIRE(__diff({a, 1}, {b, 1}));
    REQUIRE(a[0] == 13);
}

TEST_CASE("__diff - negative 2 words") {
    uint64_t a[] = {10, 0};
    uint64_t b[] = {23, 2};
    REQUIRE(__diff({a, 2}, {b, 2}));
    REQUIRE(a[0] == 13);
    REQUIRE(a[1] == 2);
}

TEST_CASE("__diff - equal") {
    uint64_t a[] = {10, 10};
    uint64_t b[] = {10, 10};
    REQUIRE(!__diff({a, 2}, {b, 2}));
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
}

TEST_CASE("__diff - borrowing") {
    uint64_t a[] = {4, 0, 1};
    uint64_t b[] = {4, 2};
    REQUIRE(!__diff({a, 3}, {b, 2}));
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == UINT64_MAX - 1);
    REQUIRE(a[2] == 0);
}

TEST_CASE("__add_and_return_carry - subtract") {
    uint64_t a[] = {10, 10};
    uint64_t b[] = {10, 10};
    bool a_neg = false;
    auto carry = __add_and_return_carry({a, 2}, a_neg, {b, 2}, true);
    REQUIRE(carry == 0);
    REQUIRE(a_neg == false);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
}

TEST_CASE("__add_and_return_carry - add") {
    uint64_t a[] = {10, 10};
    uint64_t b[] = {10, 10};
    bool a_neg = false;
    auto carry = __add_and_return_carry({a, 2}, a_neg, {b, 2}, false);
    REQUIRE(carry == 0);
    REQUIRE(a_neg == false);
    REQUIRE(a[0] == 20);
    REQUIRE(a[1] == 20);
}

TEST_CASE("__add_and_return_carry - add with carry") {
    uint64_t a[] = {1, 0};
    uint64_t b[] = {UINT64_MAX, UINT64_MAX};
    bool a_neg = false;
    auto carry = __add_and_return_carry({a, 2}, a_neg, {b, 2}, false);
    REQUIRE(carry == 1);
    REQUIRE(a_neg == false);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
}

TEST_CASE("__add_and_return_carry - add negative") {
    uint64_t a[] = {8, 10};
    uint64_t b[] = {7, 11};
    bool a_neg = true;
    auto carry = __add_and_return_carry({a, 2}, a_neg, {b, 2}, true);
    REQUIRE(carry == 0);
    REQUIRE(a_neg == true);
    REQUIRE(a[0] == 15);
    REQUIRE(a[1] == 21);
}

constexpr uint128_t operator"" _u128(const char* str) {
    uint128_t result = 0;
    for ( ; *str; ++str)
        result = result * 10 + *str - '0';
    return result;
}

TEST_CASE("__add_and_return_carry - regression") {
    uint64_t a[] = {38300863014223413, 0};
    uint128_t bb = 222189650051079730985154229019380080095_u128;
    uint64_t b[] = {static_cast<uint64_t>(bb), static_cast<uint64_t>(bb >> 64)};
    bool a_neg = false;
    auto carry = __add_and_return_carry({a, 2}, a_neg, {b, 2}, true);
    REQUIRE(carry == 0);
    REQUIRE(a_neg == true);
    uint128_t c = 222189650051079730985115928156365856682_u128;
    REQUIRE(a[0] == static_cast<uint64_t>(c));
    REQUIRE(a[1] == static_cast<uint64_t>(c >> 64));
}

TEST_CASE("__add_and_return_carry - regression alt") {
    uint64_t aa = 38300863014223413;
    uint128_t bb = 222189650051079730985154229019380080095_u128;

    uint64_t a[] = {aa, 0};
    uint64_t b[] = {static_cast<uint64_t>(bb), static_cast<uint64_t>(bb >> 64)};
    bool a_neg = false;
    auto carry = __add_and_return_carry({a, 2}, a_neg, {b, 2}, false);
    REQUIRE(carry == 0);
    REQUIRE(a_neg == false);

    auto cc = bb + aa;
    REQUIRE(a[0] == static_cast<uint64_t>(cc));
    REQUIRE(a[1] == static_cast<uint64_t>(cc >> 64));
}
