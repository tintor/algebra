#include "algebra/real_func.h"
#include "algebra/rational_func.h"
#include "algebra/integer_func.h"
#include "algebra/__stress_test.h"
using namespace algebra;
using namespace algebra::literals;

constexpr integer sample_integer(auto& rng) {
    int bits = std::uniform_int_distribution<int>(0, 256)(rng);
    if (bits == 0)
        return 0;

    integer a;
    a.abs = pow(2_n, bits - 1);
    a.abs |= uniform_sample_bits(bits - 1, rng);
    Check(a.abs.num_bits() == bits);
    if (std::uniform_int_distribution<int>(0, 1)(rng) == 0)
        a.negate();
    return a;
}

constexpr natural sample_positive_natural(auto& rng) {
    int bits = std::uniform_int_distribution<int>(1, 256)(rng);
    natural a;
    a = pow(2_n, bits - 1);
    a |= uniform_sample_bits(bits - 1, rng);
    Check(a.num_bits() == bits);
    return a;
}

constexpr rational sample_rational(auto& rng) {
    return {sample_integer(rng), sample_positive_natural(rng)};
}

template<int B>
constexpr real<B> sample_real(auto& rng) {
    int e = 0;
    if constexpr (B == 2) e = 256;
    else if constexpr (B == 10) e = 72;
    else throw std::runtime_error("unsupported base");
    return {sample_integer(rng), std::uniform_int_distribution<int>(-e, e)(rng)};
}

#define mono_zero_identities(a, z) \
    TEST(a + z == a); \
    TEST(z + a == a); \
    TEST(a - z == a); \
    TEST(z - a == -a); \
    TEST(a - a == z); \
    TEST(a * z == z); \
    TEST(z * a == z); \

#define mono_one_identities(a, o) \
    TEST(a * o == a); \
    TEST(o * a == a); \
    TEST(a * -o == -a); \
    TEST(a / o == a); \
    TEST(a / -o == -a); \
    if (a != 0) TEST(a / a == o); \
    if (a != 0) TEST(a / -a == -o); \

#define duo_identities(a, b) \
    TEST(a + b == b + a); \
    TEST(a - b == -(b - a)); \
    TEST(a - b == -(b) + a); \
    TEST(a - b == a + -(b)); \
    TEST(a * b == b * a); \
    TEST(a * b == -(a) * -(b)); \
    TEST(-a * b == a * -(b)); \
    TEST(pow(a + b, 2) == a*a + 2*a*b + b*b); \
    { auto q = a; q += b; TEST(q == a + b); } \
    { auto q = a; q -= b; TEST(q == a - b); } \
    { auto q = a; q *= b; TEST(q == a * b); } \
    TEST((a < b || a == b || a > b)); \
    if (a < b) { \
        TEST(-a > -(b)); \
        TEST(a != b); \
        TEST(a <= b); \
        TEST(!(a > b)); \
        TEST(!(a >= b)); \
    } \
    if (a > b) { \
        TEST(-a < -(b)); \
        TEST(a != b); \
        TEST(a >= b); \
        TEST(!(a < b)); \
        TEST(!(a <= b)); \
    }

#define mono_identities(a) \
    mono_zero_identities(a, 0) \
    mono_zero_identities(a, zero) \
    mono_one_identities(a, 1) \
    mono_one_identities(a, one) \
    TEST(a + a == a * 2); \
    TEST(a + a == a << 1); \
    TEST(a / 2 == a >> 1); \
    TEST(a << 8 == a * 256); \
    TEST(a >> 8 == a / 256); \
    TEST(pow(a, 0) == 1); \
    TEST(pow(a, 1) == a); \
    TEST(pow(a, 2) == a * a); \
    TEST(pow(a, 3) == a * a * a); \
    TEST(!(a < a)); \
    TEST(!(a > a)); \
    TEST(a <= a); \
    TEST(a >= a); \
    TEST(abs(a) >= 0); \
    TEST(a * a >= 0); \
    duo_identities(a, 0); \
    duo_identities(a, 1);

#define trio_identities(a, b, c) \
    TEST(a + b + c == a + c + b); \
    TEST(a - b - c == a - c - b); \
    TEST(a - (b + c) == a - b - c); \
    TEST(a * b * c == a * (b * c)); \
    TEST(a * (b + c) == a * b + a * c); \
    TEST(a * (b - c) == a * b - a * c); \
    if (a < b) \
        TEST(a + c < b + c); \
    if (a > b) \
        TEST(a + c > b + c); \
    if (c != 0) { \
        if (a < b) { \
            TEST(a * abs(c) < b * abs(c)); \
            TEST(a * -abs(c) > b * -abs(c)); \
        } \
        if (a > b) { \
            TEST(a * abs(c) > b * abs(c)); \
            TEST(a * -abs(c) < b * -abs(c)); \
        } \
    }

#define STR(A) format(#A "={}\n", A)

const bool test_add_product = false;
const bool test_sub_product = false;

void integer_test(uint64_t seed) {
    std::mt19937_64 rng(seed);
    integer q, r;
    const integer zero = 0;
    const integer one = 1;

    const integer a = sample_integer(rng);
    mono_identities(a);

    const integer b = sample_integer(rng);
    duo_identities(a, b);

    if (b != 0 && a >= 0) { integer e = a; mod(e, b); TEST(e == a % b); }

    if (test_add_product) {
        { integer e = a; add_product(e, b, one); TEST2(e == a + b, STR(a) + STR(b)); }
        { integer e = a; add_product(e, one, b); TEST(e == a + b); }
        { integer e = a; add_product(e, b, zero); TEST(e == a); }
        { integer e = a; add_product(e, zero, b); TEST(e == a); }
        { integer e = a; add_product(e, b, 1); TEST(e == a + b); }
        { integer e = a; add_product(e, 1, b); TEST(e == a + b); }
        { integer e = a; add_product(e, b, 0); TEST(e == a); }
        { integer e = a; add_product(e, 0, b); TEST(e == a); }
    }

    if (test_sub_product) {
        { integer e = a; sub_product(e, b, one); TEST(e == a - b); }
        { integer e = a; sub_product(e, one, b); TEST(e == a - b); }
        { integer e = a; sub_product(e, b, zero); TEST(e == a); }
        { integer e = a; sub_product(e, zero, b); TEST(e == a); }
        { integer e = a; sub_product(e, b, 1); TEST(e == a - b); }
        { integer e = a; sub_product(e, 1, b); TEST(e == a - b); }
        { integer e = a; sub_product(e, b, 0); TEST(e == a); }
        { integer e = a; sub_product(e, 0, b); TEST(e == a); }
    }

    if (a != 0) {
        TEST(a * b / a == b);
        TEST(a * b % a == 0);
    }
    if (b != 0) {
        TEST(a * b / b == a);
        TEST(a * b % b == 0);
        integer q, r;
        div(a, b, q, r);
        TEST(abs(r) < abs(b));
        if (a > 0) TEST(r >= 0);
        if (a == 0) TEST(r == 0);
        if (a < 0) TEST(r <= 0);
        TEST2(a == b * q + r, STR(a) + STR(b) + STR(q) + STR(r));

        q = a; q /= b; TEST(q == a / b);
        q = a; q %= b; TEST(q == a % b);
    }

    const integer c = sample_integer(rng);
    trio_identities(a, b, c);

    if (test_add_product) { integer e = a; add_product(e, b, c); TEST(e == a + b * c); }
    if (test_sub_product) { integer e = a; sub_product(e, b, c); TEST(e == a - b * c); }

    uint64_t w = rng();
    if (test_add_product) { integer e = a; add_product(e, b, w); TEST(e == a + b * w); }
    if (test_sub_product) { integer e = a; sub_product(e, b, w); TEST(e == a - b * w); }

    uint64_t m = rng();
    while (m == 0)
        m = rng();
    const uint64_t am = mod(a, m);
    const uint64_t bm = mod(b, m);
    TEST(0 <= am);
    TEST(am < m);
    using U = unsigned __int128;
    TEST(mod(U(am) + U(bm), m) == mod(a + b, m));
    TEST(mod(U(am) * U(bm), m) == mod(a * b, m));
}

void rational_test(uint64_t seed) {
    std::mt19937_64 rng(seed);
    const rational zero = 0;
    const rational one = 1;

    const rational a = sample_rational(rng);
    mono_identities(a);
    if (a != 0) {
        TEST(pow(a, -1) == 1 / a);
        TEST(pow(a, -2) == 1 / (a * a));
    }

    const rational b = sample_rational(rng);
    duo_identities(a, b);
    if (b != 0)
        { auto q = a; q /= b; TEST(q == a / b); }

    const rational c = sample_rational(rng);
    trio_identities(a, b, c);
    if (b != 0 && c != 0) {
        TEST(a / b / c == a / (b * c));
        TEST(a / b / c == a / c / b);
        TEST(a * b / c == a / c * b);
    }
    if (a != 0 && b != 0)
        TEST(a / b == 1 / (b / a));
}

template<int B>
void real_test(uint64_t seed) {
    std::mt19937_64 rng(seed);
    const real<B> zero = 0;
    const real<B> one = 1;

    const auto a = sample_real<B>(rng);
    mono_identities(a);

    const auto b = sample_real<B>(rng);
    duo_identities(a, b);

    const auto c = sample_real<B>(rng);
    trio_identities(a, b, c);
}

int main(int argc, char* argv[]) {
    stress_test([](uint64_t seed){
        integer_test(seed);
        //rational_test(seed);
        // TODO xrational
        //real_test<2>(seed);
        //real_test<10>(seed);
        // TODO expr
    });
    return 0;
}
