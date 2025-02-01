#include "algebra/real_func.h"
#include "algebra/rational_func.h"
#include "algebra/integer_func.h"
#include <print>
#include <random>
#include <vector>
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
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

std::mutex g_mutex;

#define TEST(E) { M __m{#E, seed, std::source_location::current()}; catch_exceptions(__m, [&](){__m <=> E;}); }
#define TEST2(E, MSG) { M __m{#E, seed, std::source_location::current(), MSG}; catch_exceptions(__m, [&](){__m <=> E;}); }

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

void catch_exceptions(const M& m, const auto& fn) {
    try {
        fn();
    } catch (std::runtime_error& e) {
        g_mutex.lock();
        m.print_failure();
        print("due to unexpected exception with message:\n");
        print("  {}\n", e.what());
        exit(0);
    } catch (...) {
        g_mutex.lock();
        m.print_failure();
        print("due to unexpected exception\n");
        exit(0);
    }
}

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

void operator<=>(const M& m, const bool value) {
    if (!value) {
        g_mutex.lock();
        m.print_failure();
        exit(0);
    }
}

template<typename A>
void catch_exceptions(const M1<A>& m1, const char* op, const auto& fn, const auto& b) {
    try {
        if (!fn()) {
            g_mutex.lock();
            m1.print_failure(op, b);
            exit(0);
        }
    } catch (std::runtime_error& e) {
        g_mutex.lock();
        m1.print_failure(op, b);
        print("due to unexpected exception with message:\n");
        print("  {}\n", e.what());
        exit(0);
    } catch (...) {
        g_mutex.lock();
        m1.print_failure(op, b);
        print("due to unexpected exception\n");
        exit(0);
    }
}

#define M_OP(OP) template<typename A> void operator OP (const M1<A>& m1, const auto& b) { catch_exceptions(m1, #OP, [&](){ return m1.a OP b; }, b); }

M_OP(==)
M_OP(!=)
M_OP(<)
M_OP(>)
M_OP(<=)
M_OP(>=)

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

void integer_test(uint64_t seed) {
    std::mt19937_64 rng(seed);
    integer q, r;
    const integer zero = 0;
    const integer one = 1;

    const integer a = sample_integer(rng);
    mono_identities(a);

    const integer b = sample_integer(rng);
    duo_identities(a, b);

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
        TEST(a == b * q + r);

        q = a; q /= b; TEST(q == a / b);
        q = a; q %= b; TEST(q == a % b);
    }

    const integer c = sample_integer(rng);
    trio_identities(a, b, c);

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

void run_test(uint64_t seed) {
    try {
        integer_test(seed);
        //rational_test(seed);
        // TODO xrational
        //real_test<2>(seed);
        //real_test<10>(seed);
        // TODO expr
    } catch (std::runtime_error& e) {
        g_mutex.lock();
        print("exception raised for seed {}\n{}\n", seed, e.what());
        exit(0);
    } catch (...) {
        g_mutex.lock();
        print("unknown exception for seed {}\n", seed);
        exit(0);
    }
}

int main(int argc, char* argv[]) {
    std::random_device rd;
    std::atomic<uint64_t> seed = (uint64_t(rd()) << 32) + rd();

    auto func = [&seed](){
        while (true) {
            uint64_t s = seed++;
            run_test(s);
            if (s % 1000 == 0) {
                std::lock_guard g(g_mutex);
                print("seed {}\n", s);
            }
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < std::thread::hardware_concurrency(); i++)
        threads.push_back(std::thread(func));
    threads[0].join();
    return 0;
}

namespace std {
inline namespace __1 {
bool __is_posix_terminal(__sFILE*) { return true; }
}
}
