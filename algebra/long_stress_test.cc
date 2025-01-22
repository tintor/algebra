#include <catch2/catch_test_macros.hpp>

#include "algebra/real_func.h"
#include "algebra/rational_func.h"
#include "algebra/integer_func.h"
#include <print>
#include <random>
#include <vector>
#include <string>
#include <source_location>
using std::print;
using std::format;
using namespace algebra;
using namespace algebra::literals;

constexpr integer sample_integer(auto& rng) {
    int bits = std::uniform_int_distribution<int>(0, 256)(rng);
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
    int bits = std::uniform_int_distribution<int>(1, 256)(rng);
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
    int e = 0;
    if constexpr (B == 2) e = 256;
    else if constexpr (B == 10) e = 72;
    else throw std::runtime_error("unsupported base");
    return {sample_integer(rng), std::uniform_int_distribution<int>(-e, e)(rng)};
}

void reset_color() {
    print("\033[1;0m");
}

void set_red() {
    print("\033[1;31m");
}

void set_blue() {
    print("\033[1;34m");
}

void set_yellow() {
    print("\033[1;33m");
}

#define TEST(E) { \
    M __m{#E, seed, std::source_location::current()}; \
    try { \
        __m <=> E; \
    } catch (std::runtime_error& e) { \
        __m.print_failure(); \
        print("due to unexpected exception with message:\n"); \
        print("  {}\n", e.what()); \
        exit(0); \
    } catch (...) { \
        __m.print_failure(); \
        print("due to unexpected exception\n"); \
        exit(0); \
    }\
}

#define TEST2(E, MSG) { \
    M __m{#E, seed, std::source_location::current(), MSG}; \
    try { \
        __m <=> E; \
    } catch (std::runtime_error& e) { \
        __m.print_failure(); \
        print("due to unexpected exception with message:\n"); \
        print("  {}\n", e.what()); \
        exit(0); \
    } catch (...) { \
        __m.print_failure(); \
        print("due to unexpected exception\n"); \
        exit(0); \
    }\
}

struct M {
    const char* expr;
    uint64_t seed;
    std::source_location loc;
    std::string msg;

    void print_failure() const {
        print("{}:{}: ", loc.file_name(), loc.line());
        set_red();
        print("FAILED:\n");
        set_blue();
        print("  TEST( {} )\n", expr);
        if (msg.size()) {
            reset_color();
            print("with message:\n");
            set_yellow();
            std::stringstream ss(msg);
            std::string line;
            while (getline(ss, line, '\n'))
                print("  {}\n", line);
        }
        reset_color();
        print("with seed:\n");
        set_yellow();
        print("  {}\n\n", seed);
        reset_color();
    }
};

std::string shorten(const std::string& a, int pre, int post) {
    return (a.size() <= pre + post) ? a : (a.substr(0, pre) + "..." + a.substr(a.size() - post));
}

template<typename A>
struct M1 {
    M m;
    A a;

    template<typename B>
    void print_failure(const char* op, const B& b) const {
        m.print_failure();
        print("with expansion:\n");
        set_yellow();
        auto as = shorten(std::format("{}", a), 100, 100);
        auto bs = shorten(std::format("{}", b), 100, 100);
        if (as.size() >= 20 || bs.size() >= 20)
            print("  {}\n  {}\n  {}\n", as, op, bs);
        else
            print("  {} {} {}\n", as, op, bs);
        reset_color();
    }
};

template<typename A>
M1<A> operator<=>(const M& m, const A& a) { return {m, a}; }

#define M_OP(OP) \
template<typename A, typename B> \
void operator OP (const M1<A>& m1, const B& b) { \
    try { \
        if (!(m1.a OP b)) { \
            m1.print_failure(#OP, b); \
            exit(0); \
        } \
    } catch (std::runtime_error& e) { \
        m1.print_failure(#OP, b); \
        print("due to unexpected exception with message:\n"); \
        print("  {}\n", e.what()); \
        exit(0); \
    } catch (...) { \
        m1.print_failure(#OP, b); \
        print("due to unexpected exception\n"); \
        exit(0); \
    } \
}

M_OP(==)
M_OP(!=)
M_OP(<)
M_OP(>)
M_OP(<=)
M_OP(>=)

#define mono_identities(a) \
    TEST(a + 0 == a); \
    TEST(0 + a == a); \
    TEST(a - 0 == a); \
    TEST(0 - a == -a); \
    TEST(a - a == 0); \
    TEST(a * 0 == 0); \
    TEST(0 * a == 0); \
    TEST(a * 1 == a); \
    TEST(1 * a == a); \
    TEST(a * -1 == -a); \
    TEST(a / 1 == a); \
    TEST(a / -1 == -a); \
    if (a != 0) TEST(a / a == 1); \
    if (a != 0) TEST(a / -a == -1);

void integer_test(uint64_t seed) {
    std::mt19937_64 rng(seed);
    integer q, r;

    const integer a = sample_integer(rng);
    mono_identities(a);

    const integer b = sample_integer(rng);
    TEST(a + b == b + a);
    TEST(a - b == -(b - a));
    TEST(a * b == b * a);
    if (a != 0) {
        TEST(a * b / a == b);
        TEST(a * b % a == 0);
    }
    if (b != 0) {
        integer q, r;
        div(a, b, q, r);
        //TEST(a == b * q + r);
        TEST(a * b / b == a);
        TEST(a * b % b == 0);
    }
    if (a < b)
        TEST(-a > -b);
    if (a > b)
        TEST(-a < -b);

    q = a;
    q += b;
    TEST(q == a + b);
    q = a;
    q -= b;
    TEST(q == a - b);
    q = a;
    q *= b;
    TEST2(q == a * b, format("a {}\nb {}", a, b));
    if (b != 0) {
        q = a;
        q /= b;
        TEST2(q == a / b, format("a {}\nb {}\nq {} | {} {}", a, b, q, q.abs.words[0], q.abs.words.size()));
        q = a;
        q %= b;
        TEST(q == a % b);
    }

    uint64_t m = rng();
    while (m == 0)
        m = rng();
    //TEST(mod(mod(a, m) + mod(b, m), m) == mod(a + b, m));
    //TEST(mod(mod(a, m) * mod(b, m), m) == mod(a * b, m));

    const integer c = sample_integer(rng);
    TEST(a - (b + c) == a - b - c);
    TEST(a * b * c == a * (b * c));
    TEST(a * (b + c) == a * b + a * c);
    TEST(a * (b - c) == a * b - a * c);

    const natural k = sample_positive_natural(rng);
    if (a < b) {
        TEST(a + k < b + k);
        TEST(a * k < b * k);
    }
    if (a > b) {
        TEST(a + k > b + k);
        TEST(a * k > b * k);
    }
}

void rational_test(uint64_t seed) {
    std::mt19937_64 rng(seed);

    const rational a = sample_rational(rng);
    mono_identities(a);

    const rational b = sample_rational(rng);
    TEST(a + b == b + a);
    TEST(a - b == -(b - a));
    TEST(a * b == b * a);
    if (a < b)
        TEST(-a > -b);
    if (a > b)
        TEST(-a < -b);

    const rational c = sample_rational(rng);
    TEST(a - (b + c) == a - b - c);
    TEST(a * b * c == a * (b * c));
    TEST(a * (b + c) == a * b + a * c);
    TEST(a * (b - c) == a * b - a * c);

    const natural k = sample_positive_natural(rng);
    if (a < b) {
        TEST(a + k < b + k);
        TEST(a * k < b * k);
    }
    if (a > b) {
        TEST(a + k > b + k);
        TEST(a * k > b * k);
    }
}

template<int B>
void real_test(uint64_t seed) {
    std::mt19937_64 rng(seed);

    const auto a = sample_real<B>(rng);
    mono_identities(a);

    const auto b = sample_real<B>(rng);
    TEST(a + b == b + a);
    TEST(a - b == -(b - a));
    TEST(a * b == b * a);
    if (a < b)
        TEST(-a > -b);
    if (a > b)
        TEST(-a < -b);

    const auto c = sample_real<B>(rng);
    TEST(a - (b + c) == a - b - c);
    TEST(a * b * c == a * (b * c));
    TEST(a * (b + c) == a * b + a * c);
    TEST(a * (b - c) == a * b - a * c);

    const auto k = sample_positive_natural(rng);
    if (a < b) {
        TEST2(a + k < b + k, format("a={}\nb={}\nk={}", a, b, k));
        TEST(a * k < b * k);
    }
    if (a > b) {
        TEST2(a + k > b + k, format("a={}\nb={}\nk={}", a, b, k));
        TEST(a * k > b * k);
    }
}

TEST_CASE("main") {
    std::random_device rd;
    //uint64_t seed = (uint64_t(rd()) << 32) + rd();
    uint64_t seed = 0;

    while (true) {
        if (seed % 1000 == 0)
            print("seed {}\n", seed);
        try {
            integer_test(seed);
            rational_test(seed);
            //real_test<2>(seed);
            //real_test<10>(seed);
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
