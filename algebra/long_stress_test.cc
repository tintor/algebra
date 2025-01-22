#include <catch2/catch_test_macros.hpp>

#include "algebra/real_func.h"
#include "algebra/rational_func.h"
#include "algebra/integer_func.h"
#include <print>
#include <random>
#include <vector>
#include <string>
using std::print;
using namespace algebra;
using namespace algebra::literals;

constexpr integer sample_integer(auto& rng) {
    int bits = std::uniform_int_distribution<int>(0, 512)(rng);
    if (bits == 0)
        return 0;

    integer a;
    a.abs = pow(2_n, bits - 1);
    a.abs |= uniform_sample_bits(bits - 1, rng);
    REQUIRE(a.abs.num_bits() == bits);
    if (std::uniform_int_distribution<int>(0, 1)(rng) == 0)
        a.negate();
    return a;
}

constexpr natural sample_positive_natural(auto& rng) {
    int bits = std::uniform_int_distribution<int>(1, 512)(rng);
    natural a;
    a = pow(2_n, bits - 1);
    a |= uniform_sample_bits(bits - 1, rng);
    REQUIRE(a.num_bits() == bits);
    return a;
}

constexpr rational sample_rational(auto& rng) {
    return {sample_integer(rng), sample_positive_natural(rng)};
}

template<int B>
constexpr real<B> sample_real(auto& rng) {
    return {sample_integer(rng), std::uniform_int_distribution<int>(-1000000, 1000000)(rng)};
}

void integer_test(uint64_t seed) {
    std::mt19937_64 rng(seed);

    const integer a = sample_integer(rng);
    CHECK(a + 0 == a);
    CHECK(0 + a == a);
    CHECK(a - 0 == a);
    CHECK(0 - a == -a);
    CHECK(a - a == 0);

    CHECK(a * 1 == a);
    CHECK(1 * a == a);
    CHECK(a * -1 == -a);
    CHECK(a / 1 == a);
    CHECK(a / -1 == -a);
    CHECK(a / a == 1);
    CHECK(a / -a == -1);

    const integer b = sample_integer(rng);
    CHECK(a + b == b + a);
    CHECK(a - b == -(b - a));
    CHECK(a * b == b * a);
    if (a != 0) {
        CHECK(a * b / a == b);
        CHECK(a * b % a == 0);
    }
    if (b != 0) {
        integer q, r;
        div(a, b, q, r);
        //CHECK(a == b * q + r);
        CHECK(a * b / b == a);
        CHECK(a * b % b == 0);
    }
    if (a < b)
        CHECK(-a > -b);
    if (a > b)
        CHECK(-a < -b);

    uint64_t m = rng();
    while (m == 0)
        m = rng();
    //CHECK(mod(mod(a, m) + mod(b, m), m) == mod(a + b, m));
    //CHECK(mod(mod(a, m) * mod(b, m), m) == mod(a * b, m));

    const integer c = sample_integer(rng);
    CHECK(a - (b + c) == a - b - c);
    CHECK(a * b * c == a * (b * c));
    CHECK(a * (b + c) == a * b + a * c);
    CHECK(a * (b - c) == a * b - a * c);

    const natural k = sample_positive_natural(rng);
    if (a < b) {
        CHECK(a + k < b + k);
        CHECK(a * k < b * k);
    }
    if (a > b) {
        CHECK(a + k > b + k);
        CHECK(a * k > b * k);
    }
}

void rational_test(uint64_t seed) {
    std::mt19937_64 rng(seed);

    const rational a = sample_rational(rng);
    CHECK(a + 0 == a);
    CHECK(0 + a == a);
    CHECK(a - 0 == a);
    CHECK(0 - a == -a);
    CHECK(a - a == 0);

    CHECK(a * 1 == a);
    CHECK(1 * a == a);
    CHECK(a * -1 == -a);
    CHECK(a / 1 == a);
    CHECK(a / -1 == -a);
    CHECK(a / a == 1);
    CHECK(a / -a == -1);

    const rational b = sample_rational(rng);
    CHECK(a + b == b + a);
    CHECK(a - b == -(b - a));
    CHECK(a * b == b * a);
    if (a < b)
        CHECK(-a > -b);
    if (a > b)
        CHECK(-a < -b);

    const rational c = sample_rational(rng);
    CHECK(a - (b + c) == a - b - c);
    CHECK(a * b * c == a * (b * c));
    CHECK(a * (b + c) == a * b + a * c);
    CHECK(a * (b - c) == a * b - a * c);

    const natural k = sample_positive_natural(rng);
    if (a < b) {
        CHECK(a + k < b + k);
        CHECK(a * k < b * k);
    }
    if (a > b) {
        CHECK(a + k > b + k);
        CHECK(a * k > b * k);
    }
}

template<int B>
void real_test(uint64_t seed) {
    std::mt19937_64 rng(seed);

    const auto a = sample_real<B>(rng);
    CHECK(a + 0 == a);
    CHECK(0 + a == a);
    CHECK(a - 0 == a);
    CHECK(0 - a == -a);
    CHECK(a - a == 0);

    CHECK(a * 1 == a);
    CHECK(1 * a == a);
    CHECK(a * -1 == -a);
    CHECK(a / 1 == a);
    CHECK(a / -1 == -a);
    CHECK(a / a == 1);
    CHECK(a / -a == -1);

    const auto b = sample_real<B>(rng);
    CHECK(a + b == b + a);
    CHECK(a - b == -(b - a));
    CHECK(a * b == b * a);
    if (a < b)
        CHECK(-a > -b);
    if (a > b)
        CHECK(-a < -b);

    const auto c = sample_real<B>(rng);
    CHECK(a - (b + c) == a - b - c);
    CHECK(a * b * c == a * (b * c));
    CHECK(a * (b + c) == a * b + a * c);
    CHECK(a * (b - c) == a * b - a * c);

    const auto k = sample_positive_natural(rng);
    if (a < b) {
        CHECK(a + k < b + k);
        CHECK(a * k < b * k);
    }
    if (a > b) {
        CHECK(a + k > b + k);
        CHECK(a * k > b * k);
    }
}

TEST_CASE("main") {
    std::random_device rd;
    uint64_t seed = (uint64_t(rd()) << 32) + rd();

    while (true) {
        if (seed % 100000 == 0)
            print("seed={}\n", seed);
        try {
            integer_test(seed);
            rational_test(seed);
            real_test<2>(seed);
            real_test<10>(seed);
        } catch (...) {
            print("exception seed {}\n", seed);
            throw;
        }
        seed += 1;
    }
}

namespace std {
inline namespace __1 {
bool __is_posix_terminal(__sFILE*) { return true; }
}
}
