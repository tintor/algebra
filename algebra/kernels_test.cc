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

constexpr uint128_t operator""_u128(const char* str) {
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

TEST_CASE("__sub scalar - low word becomes zero") {
    uint64_t a[] = {5, 1}; // 2**64 + 5
    inatural ia {a, 2};
    __sub(ia, static_cast<uint64_t>(5));
    REQUIRE(ia.size == 2);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 1);
}

TEST_CASE("__sub scalar - result normalized after borrow") {
    uint64_t a[] = {0, 1}; // 2**64
    inatural ia {a, 2};
    __sub(ia, static_cast<uint64_t>(1));
    REQUIRE(ia.size == 1);
    REQUIRE(a[0] == UINT64_MAX);
}

TEST_CASE("__sub scalar - result is zero") {
    uint64_t a[] = {7};
    inatural ia {a, 1};
    __sub(ia, static_cast<uint64_t>(7));
    REQUIRE(ia.size == 0);
}

TEST_CASE("__sub uint128 - no borrow out of low word") {
    uint64_t a[] = {5, 9}; // 9 * 2**64 + 5
    inatural ia {a, 2};
    __sub(ia, static_cast<uint128_t>(1) << 64 | 3); // 2**64 + 3
    REQUIRE(ia.size == 2);
    REQUIRE(a[0] == 2);
    REQUIRE(a[1] == 8);
}

TEST_CASE("__sub uint128 - borrow into third word") {
    uint64_t a[] = {0, 0, 1}; // 2**128
    inatural ia {a, 3};
    __sub(ia, static_cast<uint128_t>(1)); // 2**128 - 1
    REQUIRE(ia.size == 2);
    REQUIRE(a[0] == UINT64_MAX);
    REQUIRE(a[1] == UINT64_MAX);
}

TEST_CASE("pow uint64") {
    REQUIRE(algebra::pow(static_cast<uint64_t>(2), 0) == 1);
    REQUIRE(algebra::pow(static_cast<uint64_t>(2), 10) == 1024);
    REQUIRE(algebra::pow(static_cast<uint64_t>(2), 31) == (static_cast<uint64_t>(1) << 31));
    REQUIRE(algebra::pow(static_cast<uint64_t>(2), 40) == (static_cast<uint64_t>(1) << 40));
    REQUIRE(algebra::pow(static_cast<uint64_t>(2), 63) == (static_cast<uint64_t>(1) << 63));
    REQUIRE(algebra::pow(static_cast<uint64_t>(3), 5) == 243);
    REQUIRE(algebra::pow(static_cast<uint64_t>(10), 19) == 10000000000000000000ull);
}

TEST_CASE("pow of a narrower unsigned") {
    // the result type is uint64_t, so a narrower argument must not cap the arithmetic
    REQUIRE(algebra::pow(static_cast<uint32_t>(100'000), 3) == 1'000'000'000'000'000ull);
    REQUIRE(algebra::pow(static_cast<uint16_t>(10), 4) == 10'000);
    REQUIRE(algebra::pow(static_cast<uint8_t>(3), 6) == 729);
    REQUIRE(algebra::pow(static_cast<uint8_t>(2), 20) == 1'048'576); // the power of two path

    // against repeated multiplication in 64 bits, for a narrow and a wide argument type
    for (uint64_t base = 0; base < 40; base++)
        for (unsigned e = 0; e < 12; e++) {
            uint64_t r = 1;
            for (unsigned i = 0; i < e; i++)
                r *= base;
            REQUIRE(algebra::pow(base, e) == r);
            REQUIRE(algebra::pow(static_cast<uint32_t>(base), e) == r);
            REQUIRE(algebra::pow(static_cast<uint16_t>(base), e) == r);
        }
}

TEST_CASE("__add with zero carry into a full buffer") {
    uint64_t w[] = {1, 2};
    vnatural a {{w, 2}, 2}; // size == capacity
    uint64_t b[] = {1};
    __add(a, cnatural{b, 1}); // fits, no carry out
    REQUIRE(a.size == 2);
    REQUIRE(w[0] == 2);
    REQUIRE(w[1] == 2);
}

TEST_CASE("__add with zero carry into a full buffer, shifted") {
    uint64_t w[] = {1, 2, 3};
    vnatural a {{w, 3}, 3};
    uint64_t b[] = {1};
    __add(a, cnatural{b, 1}, 2);
    REQUIRE(a.size == 3);
    REQUIRE(w[2] == 4);
}

TEST_CASE("__sub_product by zero succeeds") {
    uint64_t a[] = {7};
    uint64_t b[] = {1, 2, 3};
    inatural ia {a, 1};
    // b is longer than a, but b * 0 is zero, so there is nothing to subtract
    REQUIRE(__sub_product(ia, cnatural{b, 3}, static_cast<uint64_t>(0)));
    REQUIRE(ia.size == 1);
    REQUIRE(a[0] == 7);

    // and it still reports a violated precondition for a non-zero multiplier
    REQUIRE(!__sub_product(ia, cnatural{b, 3}, static_cast<uint64_t>(1)));
}

TEST_CASE("inatural back() constness") {
    static_assert(std::is_same_v<decltype(std::declval<const inatural&>().back()), uint64_t>);
    static_assert(std::is_same_v<decltype(std::declval<inatural&>().back()), uint64_t&>);

    uint64_t w[] = {1, 2, 3};
    inatural a {w, 3};
    REQUIRE(a.back() == 3);
    a.back() = 7;
    REQUIRE(w[2] == 7);
}

TEST_CASE("extract_u128 - word aligned") {
    uint64_t w[] = {5, 6, 7};
    cnatural a {w, 3};
    REQUIRE(extract_u128(a, 0) == (static_cast<uint128_t>(6) << 64 | 5));
    REQUIRE(extract_u128(a, 64) == (static_cast<uint128_t>(7) << 64 | 6));
    REQUIRE(extract_u128(a, 128) == 7);
    REQUIRE(extract_u128(a, 192) == 0);
    REQUIRE(extract_u128(a, 64 * 100) == 0);
}

TEST_CASE("extract_u128 - inside a single word") {
    uint64_t w[] = {160}; // 0b1010'0000
    cnatural a {w, 1};
    REQUIRE(extract_u128(a, 0) == 160);
    REQUIRE(extract_u128(a, 4) == 10);
    REQUIRE(extract_u128(a, 7) == 1);
    REQUIRE(extract_u128(a, 8) == 0);
    REQUIRE(extract_u128(a, 64) == 0);
}

TEST_CASE("extract_u128 - unaligned across two words") {
    // 0x99aabbccddeeff001122334455667788 >> 8
    uint64_t w[] = {0x1122334455667788, 0x99aabbccddeeff00};
    cnatural a {w, 2};
    const uint128_t expected = static_cast<uint128_t>(0x0099aabbccddeeff) << 64 | 0x0011223344556677;
    REQUIRE(extract_u128(a, 8) == expected);
}

TEST_CASE("extract_u128 - unaligned across three words") {
    // 0x0f1e2d3c4b5a697899aabbccddeeff001122334455667788 >> 8, truncated to 128 bits
    uint64_t w[] = {0x1122334455667788, 0x99aabbccddeeff00, 0x0f1e2d3c4b5a6978};
    cnatural a {w, 3};
    const uint128_t expected = static_cast<uint128_t>(0x7899aabbccddeeff) << 64 | 0x0011223344556677;
    REQUIRE(extract_u128(a, 8) == expected);
    REQUIRE(extract_u128(a, 8 + 64) == (static_cast<uint128_t>(0x000f1e2d3c4b5a69) << 64 | 0x7899aabbccddeeff));
}

TEST_CASE("extract_u128 - low word agrees with extract_u64") {
    uint64_t w[] = {0x1122334455667788, 0x99aabbccddeeff00, 0x0f1e2d3c4b5a6978};
    for (int size = 0; size <= 3; size++) {
        cnatural a {w, size};
        for (int64_t e = 0; e <= 64 * 4; e++)
            REQUIRE(static_cast<uint64_t>(extract_u128(a, e)) == extract_u64(a, e));
    }
}

TEST_CASE("extract_u128 - negative shift is rejected") {
    uint64_t w[] = {1};
    REQUIRE_THROWS(extract_u128(cnatural{w, 1}, -1));
}

TEST_CASE("__complement - single word") {
    uint64_t a[] = {7};
    __complement({a, 1});
    REQUIRE(a[0] == UINT64_MAX - 6); // 2**64 - 7
}

TEST_CASE("__complement - zero stays zero") {
    uint64_t a[] = {0, 0, 0};
    __complement({a, 3});
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 0);
}

TEST_CASE("__complement - borrow out of the low word") {
    uint64_t a[] = {1, 0};
    __complement({a, 2}); // 2**128 - 1
    REQUIRE(a[0] == UINT64_MAX);
    REQUIRE(a[1] == UINT64_MAX);
}

TEST_CASE("__complement - low word is zero") {
    uint64_t a[] = {0, 1};
    __complement({a, 2}); // 2**128 - 2**64
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == UINT64_MAX);
}

TEST_CASE("__complement - all ones") {
    uint64_t a[] = {UINT64_MAX, UINT64_MAX, UINT64_MAX};
    __complement({a, 3}); // 2**192 - (2**192 - 1)
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 0);
}

TEST_CASE("__complement - twice is identity") {
    uint64_t a[] = {3, 0, 9};
    __complement({a, 3});
    __complement({a, 3});
    REQUIRE(a[0] == 3);
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 9);
}

TEST_CASE("__increment - no carry") {
    uint64_t w[] = {5, 9, 0};
    vnatural a {{w, 2}, 3};
    __increment(a);
    REQUIRE(a.size == 2);
    REQUIRE(w[0] == 6);
    REQUIRE(w[1] == 9);
}

TEST_CASE("__increment - carry into a new word") {
    uint64_t w[] = {UINT64_MAX, 0};
    vnatural a {{w, 1}, 2};
    __increment(a); // 2**64
    REQUIRE(a.size == 2);
    REQUIRE(w[0] == 0);
    REQUIRE(w[1] == 1);
}

TEST_CASE("__increment - carry through several words") {
    uint64_t w[] = {UINT64_MAX, UINT64_MAX, 0};
    vnatural a {{w, 2}, 3};
    __increment(a); // 2**128
    REQUIRE(a.size == 3);
    REQUIRE(w[0] == 0);
    REQUIRE(w[1] == 0);
    REQUIRE(w[2] == 1);
}

TEST_CASE("__increment - carry stops at the first non-max word") {
    uint64_t w[] = {UINT64_MAX, 4, 7};
    vnatural a {{w, 3}, 3};
    __increment(a);
    REQUIRE(a.size == 3);
    REQUIRE(w[0] == 0);
    REQUIRE(w[1] == 5);
    REQUIRE(w[2] == 7);
}

TEST_CASE("__increment - zero") {
    uint64_t w[] = {0};
    vnatural a {{w, 0}, 1};
    __increment(a);
    REQUIRE(a.size == 1);
    REQUIRE(w[0] == 1);
}

TEST_CASE("__increment - overflowing the buffer is reported") {
    uint64_t w[] = {UINT64_MAX};
    vnatural a {{w, 1}, 1};
    REQUIRE_THROWS(__increment(a));
}

TEST_CASE("__increment_and_return_carry") {
    uint64_t a[] = {5};
    REQUIRE(!__increment_and_return_carry({a, 1}));
    REQUIRE(a[0] == 6);
}

TEST_CASE("__increment_and_return_carry - single word wraps") {
    uint64_t a[] = {UINT64_MAX};
    REQUIRE(__increment_and_return_carry({a, 1}));
    REQUIRE(a[0] == 0);
}

TEST_CASE("__increment_and_return_carry - carry absorbed by a higher word") {
    uint64_t a[] = {UINT64_MAX, 5};
    REQUIRE(!__increment_and_return_carry({a, 2}));
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 6);
}

TEST_CASE("__increment_and_return_carry - all words wrap") {
    uint64_t a[] = {UINT64_MAX, UINT64_MAX};
    REQUIRE(__increment_and_return_carry({a, 2}));
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
}

TEST_CASE("__increment_and_return_carry - no words at all") {
    uint64_t a[] = {7};
    REQUIRE(__increment_and_return_carry({a, 0})); // 1 does not fit in zero words
    REQUIRE(a[0] == 7);
}

TEST_CASE("__decrement - no borrow") {
    uint64_t w[] = {5, 9};
    inatural a {w, 2};
    __decrement(a);
    REQUIRE(a.size == 2);
    REQUIRE(w[0] == 4);
    REQUIRE(w[1] == 9);
}

TEST_CASE("__decrement - result becomes zero") {
    uint64_t w[] = {1};
    inatural a {w, 1};
    __decrement(a);
    REQUIRE(a.size == 0);
    REQUIRE(w[0] == 0);
}

TEST_CASE("__decrement - top word becomes zero") {
    uint64_t w[] = {0, 1};
    inatural a {w, 2};
    __decrement(a); // 2**64 - 1
    REQUIRE(a.size == 1);
    REQUIRE(w[0] == UINT64_MAX);
    REQUIRE(w[1] == 0);
}

TEST_CASE("__decrement - borrow through several words") {
    uint64_t w[] = {0, 0, 1};
    inatural a {w, 3};
    __decrement(a); // 2**128 - 1
    REQUIRE(a.size == 2);
    REQUIRE(w[0] == UINT64_MAX);
    REQUIRE(w[1] == UINT64_MAX);
    REQUIRE(w[2] == 0);
}

TEST_CASE("__decrement - borrow keeps the size when a higher word survives") {
    uint64_t w[] = {0, 1, 1};
    inatural a {w, 3};
    __decrement(a);
    REQUIRE(a.size == 3);
    REQUIRE(w[0] == UINT64_MAX);
    REQUIRE(w[1] == 0);
    REQUIRE(w[2] == 1);
}

TEST_CASE("__decrement - zero is rejected") {
    uint64_t w[] = {0};
    inatural a {w, 0};
    REQUIRE_THROWS(__decrement(a));
}

TEST_CASE("__equal") {
    uint64_t a[] = {4, 5, 6};
    uint64_t b[] = {4, 5, 6};
    uint64_t c[] = {4, 5, 7};
    uint64_t d[] = {5, 5, 6};
    REQUIRE(__equal({a, 3}, {b, 3}));
    REQUIRE(!__equal({a, 3}, {c, 3})); // differs in the top word
    REQUIRE(!__equal({a, 3}, {d, 3})); // differs in the low word
    REQUIRE(!__equal({a, 3}, {a, 2})); // sizes have to match
    REQUIRE(__equal({a, 0}, {b, 0})); // zero equals zero
    REQUIRE(!__equal({a, 0}, {b, 1}));
    REQUIRE(__equal({a, 1}, {b, 1}));
}

TEST_CASE("__less - by size") {
    uint64_t a[] = {UINT64_MAX};
    uint64_t b[] = {1, 1};
    REQUIRE(__less({a, 1}, {b, 2}));
    REQUIRE(!__less({b, 2}, {a, 1}));
    REQUIRE(__less({a, 0}, {a, 1}));
    REQUIRE(!__less({a, 1}, {a, 0}));
    REQUIRE(!__less({a, 0}, {b, 0}));
}

TEST_CASE("__less - same size") {
    uint64_t a[] = {5, 1};
    uint64_t b[] = {6, 1};
    uint64_t c[] = {1, 2};
    REQUIRE(__less({a, 2}, {b, 2}));
    REQUIRE(!__less({b, 2}, {a, 2}));
    REQUIRE(!__less({a, 2}, {a, 2})); // not less than itself
    REQUIRE(__less({b, 2}, {c, 2})); // the top word decides
    REQUIRE(!__less({c, 2}, {b, 2}));
}

TEST_CASE("__mod uint64") {
    uint64_t a[] = {10};
    REQUIRE(__mod(cnatural{a, 1}, static_cast<uint64_t>(7)) == 3);
    REQUIRE(__mod(cnatural{a, 1}, static_cast<uint64_t>(10)) == 0);
    REQUIRE(__mod(cnatural{a, 1}, static_cast<uint64_t>(1)) == 0);
    REQUIRE(__mod(cnatural{a, 0}, static_cast<uint64_t>(7)) == 0);
}

TEST_CASE("__mod uint64 - multi word") {
    uint64_t a[] = {0, 1}; // 2**64 == 18446744073709551616
    REQUIRE(__mod(cnatural{a, 2}, static_cast<uint64_t>(10)) == 6);
    REQUIRE(__mod(cnatural{a, 2}, static_cast<uint64_t>(3)) == 1);
    REQUIRE(__mod(cnatural{a, 2}, static_cast<uint64_t>(7)) == 2);

    uint64_t b[] = {0, 0, 1}; // 2**128, and 2**128 % 7 == 2 * 2 == 4
    REQUIRE(__mod(cnatural{b, 3}, static_cast<uint64_t>(7)) == 4);
    REQUIRE(__mod(cnatural{b, 3}, static_cast<uint64_t>(2)) == 0);

    uint64_t c[] = {UINT64_MAX, UINT64_MAX};
    REQUIRE(__mod(cnatural{c, 2}, static_cast<uint64_t>(2)) == 1);
    REQUIRE(__mod(cnatural{c, 2}, UINT64_MAX) == 0); // (2**128 - 1) == (2**64 + 1) * (2**64 - 1)
}

TEST_CASE("__mod uint128") {
    const uint128_t p127 = static_cast<uint128_t>(1) << 127;
    const uint128_t p100 = static_cast<uint128_t>(1) << 100;
    const uint128_t b = (static_cast<uint128_t>(1) << 64) + 1;

    uint64_t small[] = {7};
    REQUIRE(__mod(cnatural{small, 1}, p100) == 7);
    REQUIRE(__mod(cnatural{small, 0}, p100) == 0);

    uint64_t a[] = {5, 0, 1}; // 2**128 + 5 == 2 * 2**127 + 5
    REQUIRE(__mod(cnatural{a, 3}, p127) == 5);
    a[0] = 0;
    REQUIRE(__mod(cnatural{a, 3}, p127) == 0);
    // 2**128 - 1 == (2**64 - 1) * (2**64 + 1), so 2**128 == 1 (mod 2**64 + 1)
    REQUIRE(__mod(cnatural{a, 3}, b) == 1);

    uint64_t c[] = {UINT64_MAX, UINT64_MAX}; // 2**128 - 1
    REQUIRE(__mod(cnatural{c, 2}, b) == 0);
    REQUIRE(__mod(cnatural{c, 2}, p127) == p127 - 1);
}

TEST_CASE("__mod uint64 - agrees with the mod3/5/7/9 kernels") {
    uint64_t a[] = {0x1122334455667788, 0x99aabbccddeeff00, 0x0f1e2d3c4b5a6978};
    for (int size = 0; size <= 3; size++) {
        cnatural x {a, size};
        REQUIRE(__mod(x, static_cast<uint64_t>(3)) == static_cast<uint64_t>(mod3(x)));
        REQUIRE(__mod(x, static_cast<uint64_t>(5)) == static_cast<uint64_t>(mod5(x)));
        REQUIRE(__mod(x, static_cast<uint64_t>(6)) == static_cast<uint64_t>(mod6(x)));
        REQUIRE(__mod(x, static_cast<uint64_t>(7)) == static_cast<uint64_t>(mod7(x)));
        REQUIRE(__mod(x, static_cast<uint64_t>(9)) == static_cast<uint64_t>(mod9(x)));
        REQUIRE(__mod(x, static_cast<uint64_t>(10)) == static_cast<uint64_t>(mod10(x)));
    }
}

TEST_CASE("__mul_add_return_carry - single word") {
    uint64_t a[] = {5};
    REQUIRE(__mul_add_return_carry({a, 1}, 3, 2) == 0);
    REQUIRE(a[0] == 17);
}

TEST_CASE("__mul_add_return_carry - carry out") {
    uint64_t a[] = {UINT64_MAX};
    REQUIRE(__mul_add_return_carry({a, 1}, 2, 0) == 1); // 2**65 - 2
    REQUIRE(a[0] == UINT64_MAX - 1);
}

TEST_CASE("__mul_add_return_carry - multiply by zero") {
    uint64_t a[] = {5, 7};
    REQUIRE(__mul_add_return_carry({a, 2}, 0, 9) == 0);
    REQUIRE(a[0] == 9);
    REQUIRE(a[1] == 0);
}

TEST_CASE("__mul_add_return_carry - multiply by one keeps the value") {
    uint64_t a[] = {1, 2, 3};
    REQUIRE(__mul_add_return_carry({a, 3}, 1, 0) == 0);
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    REQUIRE(a[2] == 3);
}

TEST_CASE("__mul_add_return_carry - maximal carry") {
    // (2**128 - 1) * (2**64 - 1) + (2**64 - 1) == (2**64 - 1) * 2**128
    uint64_t a[] = {UINT64_MAX, UINT64_MAX};
    REQUIRE(__mul_add_return_carry({a, 2}, UINT64_MAX, UINT64_MAX) == UINT64_MAX);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
}

TEST_CASE("__mul_add_return_carry - carry propagates between words") {
    uint64_t a[] = {UINT64_MAX, 1}; // 2**65 - 1
    REQUIRE(__mul_add_return_carry({a, 2}, 2, 1) == 0); // 2**66 - 1
    REQUIRE(a[0] == UINT64_MAX);
    REQUIRE(a[1] == 3);
}

TEST_CASE("__normalized_size") {
    uint64_t a[] = {1, 0, 0};
    REQUIRE(__normalized_size({a, 3}) == 1);
    REQUIRE(__normalized_size({a, 1}) == 1);
    REQUIRE(__normalized_size({a, 0}) == 0);

    uint64_t b[] = {0, 1, 0};
    REQUIRE(__normalized_size({b, 3}) == 2);

    uint64_t c[] = {0, 0, 0};
    REQUIRE(__normalized_size({c, 3}) == 0);

    uint64_t d[] = {1, 2, 3};
    REQUIRE(__normalized_size({d, 3}) == 3);
}

TEST_CASE("inatural normalize") {
    uint64_t w[] = {1, 0, 0};
    inatural a {w, 3};
    a.normalize();
    REQUIRE(a.size == 1);
    REQUIRE(w[0] == 1); // the words themselves are untouched

    inatural b {w, 0};
    b.normalize();
    REQUIRE(b.size == 0);

    uint64_t z[] = {0, 0};
    inatural c {z, 2};
    c.normalize();
    REQUIRE(c.size == 0);

    uint64_t d[] = {0, 5};
    inatural e {d, 2};
    e.normalize();
    REQUIRE(e.size == 2);
}

TEST_CASE("add_max_size") {
    uint64_t one[] = {1};
    uint64_t max[] = {UINT64_MAX};
    uint64_t big[] = {5, UINT64_MAX};
    uint64_t two_words[] = {1, 1};

    // no carry out of the top word
    REQUIRE(add_max_size({one, 1}, {one, 1}) == 1);
    REQUIRE(add_max_size({two_words, 2}, {one, 1}) == 2);
    REQUIRE(add_max_size({one, 1}, {two_words, 2}) == 2);

    // a carry out of the top word is possible
    REQUIRE(add_max_size({max, 1}, {one, 1}) == 2);
    REQUIRE(add_max_size({max, 1}, {max, 1}) == 2);
    REQUIRE(add_max_size({big, 2}, {one, 1}) == 3);
    REQUIRE(add_max_size({one, 1}, {big, 2}) == 3);

    // adding zero can never grow the result
    REQUIRE(add_max_size({max, 1}, {one, 0}) == 1);
    REQUIRE(add_max_size({one, 0}, {max, 1}) == 1);
    REQUIRE(add_max_size({big, 2}, {one, 0}) == 2);
    REQUIRE(add_max_size({one, 0}, {one, 0}) == 0);
}

TEST_CASE("add_max_size - is an upper bound") {
    uint64_t w[] = {0, 0, 0};
    const uint64_t values[] = {0, 1, 2, UINT64_MAX - 1, UINT64_MAX};
    for (uint64_t a1 : values)
        for (uint64_t a0 : values)
            for (uint64_t b0 : values) {
                w[0] = a0;
                w[1] = a1;
                w[2] = b0;
                cnatural a {w, a1 ? 2 : (a0 ? 1 : 0)};
                cnatural b {w + 2, b0 ? 1 : 0};
                // hand computed size of a + b
                const uint128_t low = static_cast<uint128_t>(a0) + b0;
                const uint128_t high = static_cast<uint128_t>(a1) + static_cast<uint64_t>(low >> 64);
                int size = 0;
                if (high >> 64)
                    size = 3;
                else if (static_cast<uint64_t>(high))
                    size = 2;
                else if (static_cast<uint64_t>(low))
                    size = 1;
                REQUIRE(add_max_size(a, b) >= size);
                REQUIRE(add_max_size(b, a) >= size);
            }
}

TEST_CASE("mul_max_size") {
    uint64_t one[] = {1};
    uint64_t four[] = {4};
    uint64_t max[] = {UINT64_MAX};
    uint64_t p32[] = {static_cast<uint64_t>(1) << 32};
    uint64_t three_words[] = {0, 0, 1};
    uint64_t two_words[] = {0, 1};

    REQUIRE(mul_max_size({one, 1}, {one, 1}) == 1);
    REQUIRE(mul_max_size({four, 1}, {four, 1}) == 1);
    REQUIRE(mul_max_size({max, 1}, {max, 1}) == 2);
    REQUIRE(mul_max_size({p32, 1}, {p32, 1}) == 2); // 2**64 needs a second word
    REQUIRE(mul_max_size({three_words, 3}, {two_words, 2}) == 4); // 2**128 * 2**64
    REQUIRE(mul_max_size({three_words, 3}, {one, 1}) == 3);

    // multiplying by zero
    REQUIRE(mul_max_size({max, 1}, {one, 0}) == 0);
    REQUIRE(mul_max_size({one, 0}, {max, 1}) == 0);
    REQUIRE(mul_max_size({one, 0}, {one, 0}) == 0);
}

TEST_CASE("div_max_size") {
    uint64_t one[] = {1};
    uint64_t two[] = {2};
    uint64_t seven[] = {7};
    uint64_t max[] = {UINT64_MAX};
    uint64_t two_words[] = {0, 1};

    REQUIRE(div_max_size({seven, 1}, {two, 1}) == 1);
    REQUIRE(div_max_size({two_words, 2}, {two, 1}) == 1); // 2**63 fits in one word
    REQUIRE(div_max_size({two_words, 2}, {one, 1}) == 2); // 2**64 needs two words
    REQUIRE(div_max_size({max, 1}, {one, 1}) == 1);

    // a smaller dividend gives a zero quotient
    REQUIRE(div_max_size({two, 1}, {two_words, 2}) == 0);
    REQUIRE(div_max_size({one, 0}, {two, 1}) == 0);

    // there is no quotient by zero, and no word of the divisor may be read
    REQUIRE(div_max_size({max, 1}, {two, 0}) == 0);
    REQUIRE(div_max_size({one, 0}, {two, 0}) == 0);
}

template<typename A, typename B>
concept has_min = requires (A a, B b) { algebra::min(a, b); };
template<typename A, typename B>
concept has_max = requires (A a, B b) { algebra::max(a, b); };

TEST_CASE("min and max of two builtin integers") {
    REQUIRE(algebra::min(uint32_t(3), uint64_t(5)) == 3);
    REQUIRE(algebra::max(uint32_t(3), uint64_t(5)) == 5);
    static_assert(std::same_as<decltype(algebra::min(uint32_t(3), uint64_t(5))), uint64_t>);

    REQUIRE(algebra::min(int32_t(-1), int64_t(5)) == -1);
    REQUIRE(algebra::max(int32_t(-1), int64_t(5)) == 5);
    REQUIRE(algebra::min(int32_t(-1), int32_t(-5)) == -5);
    REQUIRE(algebra::max(int64_t(-1), int64_t(-5)) == -1);

    // Mixed signedness does not compile. The comparison would go through the usual arithmetic
    // conversions, where a negative value becomes a huge unsigned one -- min(-1, 5u) returned 5 --
    // and there is no return type that holds both operands either.
    static_assert(has_min<uint32_t, uint64_t> && has_max<uint32_t, uint64_t>);
    static_assert(has_min<int32_t, int64_t> && has_max<int32_t, int64_t>);
    static_assert(!has_min<int32_t, uint32_t> && !has_max<int32_t, uint32_t>);
    static_assert(!has_min<uint32_t, int32_t> && !has_max<uint32_t, int32_t>);
    static_assert(!has_min<int64_t, uint64_t> && !has_max<int64_t, uint64_t>);
}
