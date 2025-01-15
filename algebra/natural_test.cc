#include "algebra/natural.h"
#include <catch2/catch_test_macros.hpp>
#include <print>
#include <random>
using std::print;
using std::format;
using ucent = unsigned __int128;
using namespace algebra;

class Random {
public:
    Random() : _rng(std::random_device{}()) {}
    Random(unsigned long seed) : _rng(seed) {}
    operator std::mt19937_64&() { return _rng; }
    std::mt19937_64& get() { return _rng; }

    template<typename T>
    T Uniform(T min, T max) {
        static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>);

        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(_rng);
        }
        if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min, max);
            return dist(_rng);
        }
    }
private:
    std::mt19937_64 _rng;
};

natural rand_natural(int min_size, int max_size, Random& rng) {
    natural a;
    a.words.reset(rng.Uniform<int>(min_size, max_size));
    for (int i = 0; i < a.words.size(); i++)
        a.words[i] = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
    return a;
}

void rand_natural(natural& a, int min_size, int max_size, Random& rng) {
    a.words.reset(rng.Uniform<int>(min_size, max_size));
    for (int i = 0; i < a.words.size(); i++)
        a.words[i] = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
}

TEST_CASE("__word_div") {
    Random rng(0);
    natural a, b;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 10, rng);
        rand_natural(b, 1, 10, rng);
        const uint64_t q = __word_div(a, b);
        REQUIRE(a * q <= b);
        if (q != std::numeric_limits<uint64_t>::max())
            REQUIRE(a * (natural(q) + 1) > b);
    }
}

TEST_CASE("words") {
    natural b(1);
    REQUIRE(b.words.size() == 1);

    REQUIRE(natural(0).words.size() == 0);
    REQUIRE(natural(1).words.size() == 1);
    REQUIRE(natural(2).words.size() == 1);

    REQUIRE(natural(1).words[0] == 1);
    REQUIRE(natural(0).words[0] == 0);
}

TEST_CASE("str") {
    REQUIRE(natural(0).str() == "0");
    REQUIRE(natural(1).str() == "1");
    REQUIRE(natural(12).str() == "12");
    REQUIRE(natural(450).str() == "450");
}

constexpr std::string os() {
    std::ostringstream s;
    s << 15_n;
    return s.str();
}

TEST_CASE("constexpr ostream") {
    REQUIRE(os() == "15");
}

TEST_CASE("format") {
    REQUIRE(format("{:b}", 15_n) == "1111");
    REQUIRE(format("{:o}", 15_n) == "17");
    REQUIRE(format("{:d}", 15_n) == "15");
    REQUIRE(format("{:x}", 15_n) == "f");
    REQUIRE(format("{:X}", 15_n) == "F");
    REQUIRE(format("{:4d}", 15_n) == "  15");
    REQUIRE(format("{:1d}", 15_n) == "15");
    REQUIRE(format("{:*>4d}", 15_n) == "**15");
    REQUIRE(format("{:*<4d}", 15_n) == "15**");
    REQUIRE(format("{:*^4d}", 15_n) == "*15*");
}

TEST_CASE("hex") {
    REQUIRE(natural(0).hex() == "0");
    REQUIRE(natural(16).hex() == "10");
    REQUIRE(natural(255).hex() == "FF");
    REQUIRE(natural(256).hex() == "100");
}

TEST_CASE("parse") {
    REQUIRE(natural("0") == natural(0));
    REQUIRE(natural("1") == natural(1));
    REQUIRE(natural("12") == natural(12));
    REQUIRE(natural("450") == natural(450));
    const char* a = "18446744073709551617"; // UINT64_MAX + 2
    REQUIRE(natural(a).str() == a);

    REQUIRE(natural("1100", 2) == natural(12));
    REQUIRE(natural("111", 2) == natural(7));
    REQUIRE(natural("FF", 16) == natural(255));
    REQUIRE(natural("ff", 16) == natural(255));
}

TEST_CASE("static_cast<uint>") {
    REQUIRE(static_cast<uint>(natural(0)) == 0);
    REQUIRE(static_cast<uint>(natural(1)) == 1);
    uint a = std::numeric_limits<uint>::max();
    REQUIRE(static_cast<uint>(natural(a)) == a);
}

TEST_CASE("static_cast<uint64_t>") {
    REQUIRE(static_cast<uint64_t>(natural(0)) == 0);
    REQUIRE(static_cast<uint64_t>(natural(1)) == 1);
    uint64_t a = std::numeric_limits<uint64_t>::max();
    REQUIRE(static_cast<uint64_t>(natural(a)) == a);
}

TEST_CASE("ucent") {
    natural::dword a = 1;
    for (int i = 0; i < 128; i++) {
        natural b(a);
        REQUIRE(b == a);
        REQUIRE(a == b);
        a <<= 1;
    }
}

TEST_CASE("is_even / is_odd") {
    REQUIRE(natural(0).is_even());
    REQUIRE(!natural(0).is_odd());
    REQUIRE(!natural(7).is_even());
    REQUIRE(natural(7).is_odd());
}

TEST_CASE("cmp") {
    REQUIRE(natural(0) == natural(0));
    REQUIRE(natural(5) == natural(5));

    REQUIRE(natural(0) < natural(1));
    REQUIRE(natural(5) < natural(6));

    REQUIRE(natural(6) > natural(5));

    REQUIRE(natural(1) <= natural(1));
    REQUIRE(0_n <= 888089631791237197_n);
}

TEST_CASE("add") {
    REQUIRE(natural(0) + natural(0) == natural(0));
    REQUIRE(natural(5) + natural(0) == natural(5));
    REQUIRE(natural(5) + natural(6) == natural(11));
}

TEST_CASE("sub") {
    REQUIRE(natural(0) - natural(0) == natural(0));
    REQUIRE(natural(5) - natural(0) == natural(5));
    REQUIRE(natural(5) - natural(5) == natural(0));
    REQUIRE(natural(5) - natural(4) == natural(1));
}

TEST_CASE("+=") {
    natural a(5);
    a += natural(4);
    REQUIRE(a == natural(9));
}

TEST_CASE("-=") {
    natural a(5);
    a -= natural(4);
    REQUIRE(a == natural(1));
}

TEST_CASE("big add 1") {
    uint64_t a = std::numeric_limits<uint64_t>::max();
    natural b = a;
    REQUIRE(b == a);

    b += a;
    REQUIRE(b > a);
    REQUIRE(b.words.size() == 2);
    REQUIRE(b.words[1] == 1);
    REQUIRE(b.words[0] == a - 1);

    b += b;
    REQUIRE(b.words.size() == 2);
    REQUIRE(b.words[1] == 3);
    REQUIRE(b.words[0] == a - 3);
}

TEST_CASE("big add 2") {
    const uint64_t m = std::numeric_limits<uint64_t>::max();

    natural a;
    a.words.reset(4);
    a.words[0] = m;
    a.words[1] = m;
    a.words[2] = m;
    a.words[3] = 1;

    natural b;
    b.words.reset(4);
    b.words[0] = 1;
    b.words[1] = 0;
    b.words[2] = 0;
    b.words[3] = 1;

    natural c = a + b;
    REQUIRE(c.words.size() == 4);
    REQUIRE(c.words[0] == 0);
    REQUIRE(c.words[1] == 0);
    REQUIRE(c.words[2] == 0);
    REQUIRE(c.words[3] == 3);
}


TEST_CASE("mul") {
    REQUIRE(natural(0) * natural(0) == natural(0));
    REQUIRE(natural(5) * natural(0) == natural(0));
    REQUIRE(natural(0) * natural(2) == natural(0));
    REQUIRE(natural(5) * natural(2) == natural(10));
}

TEST_CASE("*=") {
    natural a(5);
    a *= natural(3);
    REQUIRE(a == natural(15));
    a *= 3;
    REQUIRE(a == 45);

    natural c(3);
    c *= c;
    REQUIRE(c == 9);

    natural d(3), e("1000000000000000000000000000000000000");
    REQUIRE(e.words.size() == 2);
    d *= e;
    REQUIRE(d.words.size() == 2);
    REQUIRE(d == natural("3000000000000000000000000000000000000"));
}

TEST_CASE("/") {
    REQUIRE(natural(0) / natural(7) == natural(0));
    REQUIRE(natural(7) / natural(7) == natural(1));
    REQUIRE(natural(7) / natural(3) == natural(2));
    REQUIRE(natural(7) / natural(8) == natural(0));
}

TEST_CASE("%") {
    REQUIRE(natural(0) % 7u == 0);
    REQUIRE(natural(1) % 7u == 1);
    REQUIRE(natural(7) % 7u == 0);
    REQUIRE(natural(8) % 7u == 1);
    REQUIRE(natural(14) % 7u == 0);

    REQUIRE(natural(0) % 7ul == 0);
    REQUIRE(natural(1) % 7ul == 1);
    REQUIRE(natural(7) % 7ul == 0);
    REQUIRE(natural(8) % 7ul == 1);
    REQUIRE(natural(14) % 7ul == 0);

    REQUIRE(natural(5) % 3u == 2);
    REQUIRE(natural(5) % 3ul == 2);
}

TEST_CASE("add stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        ucent b = rng.Uniform<ucent>(0, m - a);
        REQUIRE(natural(a) + natural(b) == a + b);
    }
}

TEST_CASE("sub stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        ucent b = rng.Uniform<ucent>(0, a);
        REQUIRE(natural(a).str() == format("{}", a));
        REQUIRE(natural(b).str() == format("{}", b));
        REQUIRE(natural(a) - natural(b) == a - b);
    }
}

TEST_CASE("mul stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(1, m);
        ucent b = rng.Uniform<ucent>(0, m / a);
        REQUIRE(natural(a) * natural(b) == a * b);
    }
}

TEST_CASE("div stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        long b = rng.Uniform<long>(0, std::numeric_limits<long>::max());
        natural q;
        REQUIRE(div(natural(a), b, q) == a % b);
        REQUIRE(q == a / b);
    }
}

TEST_CASE("str stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        natural b = a;
        if (a > std::numeric_limits<uint64_t>::max()) {
            REQUIRE(b.words.size() == 2);
            REQUIRE(b.words[0] == uint64_t(a));
            REQUIRE(b.words[1] == uint64_t(a >> 64));
        }
        REQUIRE(format("{}", a) == natural(a).str());
    }
}

TEST_CASE("div10 stress with ucent 2") {
    Random rng;
    const ucent m = std::numeric_limits<uint64_t>::max();
    for (int i = 0; i < 1'000'000; i++) {
        ucent a = rng.Uniform<ucent>(m + 1, m * 10);
        natural q;
        uint64_t r = div(natural(a), 10ull, q);
        REQUIRE(q * 10 + r == a);
    }
}

natural rand_natural(int size, Random& rng) {
    natural a;
    a.words.reset(size);
    for (int i = 0; i < size; i++)
        a.words[i] = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
    return a;
}

TEST_CASE("stress + and -") {
    Random rng(0);
    for (int i = 0; i < 1000'000; i++) {
        const natural a = rand_natural(rng.Uniform<int>(1, 10), rng);
        const natural b = rand_natural(rng.Uniform<int>(1, 10), rng);
        const natural c = rand_natural(rng.Uniform<int>(1, 10), rng);

        const natural sab = a + b;
        const natural sac = a + c;
        const natural sbc = b + c;
        REQUIRE(sab == b + a);
        REQUIRE(sbc == c + b);
        REQUIRE(sac == c + a);
        const natural s = sab + c;
        REQUIRE(s == c + sab);
        REQUIRE(s == sbc + a);
        REQUIRE(s == a + sbc);
        REQUIRE(s == sac + b);
        REQUIRE(s == b + sac);

        REQUIRE(s - a == sbc);
        REQUIRE(s - b == sac);
        REQUIRE(s - c == sab);
    }
}

TEST_CASE("stress a += a") {
    Random rng(0);
    for (int i = 0; i < 1000'000; i++) {
        const natural a = rand_natural(rng.Uniform<int>(1, 10), rng);
        natural m = a;
        m += m;
        REQUIRE(m == a + a);
        REQUIRE(m == a * 2);
    }
}

natural safe_mul(const natural& a, const natural& b) {
    natural out;
    for (int i = 0; i < a.words.size(); i++)
        for (int j = 0; j < b.words.size(); j++) {
            natural e = ucent(a.words[i]) * b.words[j];
            out += e << ((i + j) * 64);
        }
    return out;
}

TEST_CASE("stress mul") {
    Random rng(0);
    for (int i = 0; i < 1000'000; i++) {
        const natural a = rand_natural(rng.Uniform<int>(1, 5), rng);
        const natural b = rand_natural(rng.Uniform<int>(1, 5), rng);
        const natural ab = safe_mul(a, b);

        natural c;
        mul(a, b, c);
        REQUIRE(c == ab);
        c = 0;
        mul(b, a, c);
        REQUIRE(c == ab);
        c = a;
        mul(c, b);
        REQUIRE(c == ab);
        c = b;
        mul(c, a);
        REQUIRE(c == ab);
    }
}

TEST_CASE("square") {
    natural a;
    a.words.reset(2);
    a.words[0] = 2;
    a.words[1] = 10;
    natural m = a;
    square(m);
    natural e;
    __mul(a, a, e);
    REQUIRE(m == e);
}

TEST_CASE("stress square in-place") {
    Random rng(0);
    for (int i = 0; i < 100'000; i++) {
        natural a = rand_natural(rng.Uniform<int>(1, 5), rng);
        natural m = a;
        m *= m;
        REQUIRE(m == a * a);
    }
}

TEST_CASE("stress div with remainder") {
    Random rng;
    for (int i = 0; i < 100'000; i++) {
        const natural a = rand_natural(rng.Uniform<int>(2, 10), rng);

        const natural b = rand_natural(rng.Uniform<int>(1, 5), rng);
        natural quot, rem;
        div(a, b, quot, rem);
        REQUIRE(rem < b);
        REQUIRE(quot * b + rem == a);

        const uint64_t c = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
        uint64_t m = div(a, c, quot);
        REQUIRE(m < c);
        REQUIRE(quot * c + m == a);
    }
}

TEST_CASE("uniform_int") {
    natural a = (1_n << 128) - 1;
    Random rng;
    for (int i = 0; i < 20; i++) {
        natural m = uniform_int(0, a, rng.get());
        REQUIRE(m.words.size() <= a.words.size());
        REQUIRE(0 <= m);
        REQUIRE(m <= a);
        REQUIRE(0 <= m.words.size());
        REQUIRE(m.words.size() <= 2);
        while (m.words.size()) {
            div(m, 10ull, /*out*/m);
        }
    }
}

// TODO randomized long division test against cpp_int for big naturals!

TEST_CASE("is_x") {
    natural z = 0;
    REQUIRE(z.is_uint());
    REQUIRE(z.is_uint64());

    natural o = 1;
    REQUIRE(o.is_uint());
    REQUIRE(o.is_uint64());

    natural p = std::numeric_limits<uint>::max();
    REQUIRE(p.is_uint());
    REQUIRE(p.is_uint64());

    natural q = (uint64_t)std::numeric_limits<uint>::max() + 1;
    REQUIRE(!q.is_uint());
    REQUIRE(q.is_uint64());

    natural a = std::numeric_limits<long>::max();
    REQUIRE(!a.is_uint());
    REQUIRE(a.is_uint64());

    natural b = (uint64_t)std::numeric_limits<long>::max() + 1;
    REQUIRE(!b.is_uint());
    REQUIRE(b.is_uint64());

    natural c = std::numeric_limits<uint64_t>::max();
    REQUIRE(!c.is_uint());
    REQUIRE(c.is_uint64());

    natural d = c + 1;
    REQUIRE(!d.is_uint());
    REQUIRE(!d.is_uint64());
}

TEST_CASE("<<=") {
    for (int i = 0; i <= 10; i++) {
        natural a(i);
        a <<= 1;
        REQUIRE(a == natural(i << 1));
    }
}

TEST_CASE(">>=") {
    for (int i = 0; i <= 10; i++) {
        natural a(i);
        a >>= 1;
        REQUIRE(a == natural(i >> 1));
    }
}

TEST_CASE("pow") {
    REQUIRE(pow(natural(2), 3) == 8);
    REQUIRE(pow(natural(10), 30) == natural("1000000000000000000000000000000"));
}

TEST_CASE("gcd") {
    REQUIRE(num_trailing_zeros(0u) == 0);
    REQUIRE(num_trailing_zeros(1u) == 0);
    REQUIRE(num_trailing_zeros(2u) == 1);
    REQUIRE(num_trailing_zeros(8u) == 3);
    REQUIRE(num_trailing_zeros(24u) == 3);

    REQUIRE(gcd(5u, 5u) == 5u);
    REQUIRE(gcd(6u, 15u) == 3u);
    REQUIRE(gcd(7u, 3u) == 1u);

    REQUIRE(gcd(natural(5), natural(5)) == natural(5));
    REQUIRE(gcd(natural(6), natural(15)) == natural(3));
    REQUIRE(gcd(natural(7), natural(3)) == natural(1));

    REQUIRE(gcd(natural(5), 5) == 5);
    REQUIRE(gcd(1, natural(5)) == 1);
}

TEST_CASE("factorial") {
    natural a(1);
    for (int i = 2; i <= 50; i++) {
        natural b = a;
        for (int j = 1; j < i; j++)
            b += a;
        REQUIRE(a * i == b);
        REQUIRE(b / i == a);
        a *= natural(i);
        REQUIRE(a == b);
        if (i == 30) REQUIRE(a.str() == "265252859812191058636308480000000");
        if (i == 50) REQUIRE(a.str() == "30414093201713378043612608166064768844377641568960512000000000000");
    }
}

TEST_CASE("num_bits") {
    REQUIRE(natural(0).num_bits() == 0);
    REQUIRE(natural(1).num_bits() == 1);
    REQUIRE(natural(2).num_bits() == 2);
    REQUIRE(natural(3).num_bits() == 2);
    REQUIRE(natural(4).num_bits() == 3);
    REQUIRE(natural(15).num_bits() == 4);
    REQUIRE(natural(16).num_bits() == 5);
}

#include <bit>

TEST_CASE("popcount") {
    for (uint i: {0, 4, 31231, -3123121})
        REQUIRE(natural(i).popcount() == std::popcount(i));
}

TEST_CASE("isqrt") {
    REQUIRE(isqrt(natural(0)) == 0);
    REQUIRE(isqrt(natural(1)) == 1);
    REQUIRE(isqrt(natural(2)) == 1);
    REQUIRE(isqrt(natural(3)) == 1);
    REQUIRE(isqrt(natural(4)) == 2);
    REQUIRE(isqrt(natural(5)) == 2);
    REQUIRE(isqrt(natural(9)) == 3);
    REQUIRE(isqrt(natural(9999)) == 99);
    REQUIRE(isqrt(natural(10000)) == 100);
}

TEST_CASE("is_prime") {
    REQUIRE(!is_prime(natural(0)));
    REQUIRE(!is_prime(natural(1)));
    REQUIRE(is_prime(natural(2)));
    REQUIRE(is_prime(natural(3)));
    REQUIRE(!is_prime(natural(4)));
    REQUIRE(is_prime(natural(5)));
    REQUIRE(!is_prime(natural(6)));
}

TEST_CASE("<<") {
    natural(1) << 64; // regression test
}

TEST_CASE("is_power_of_two") {
    REQUIRE(is_power_of_two(1));
    REQUIRE(is_power_of_two(2));
    REQUIRE(!is_power_of_two(3));
    REQUIRE(is_power_of_two(4));
    REQUIRE(!is_power_of_two(5));

    REQUIRE(natural(64).num_bits() == 7);
    REQUIRE(natural(64).num_trailing_zeros() == 6);
    natural a = natural(1) << 100;
    REQUIRE(a.words.size() == 2);
    REQUIRE(a.words[0] == 0);
    REQUIRE(a.words[1] == uint64_t(1) << 36);
    REQUIRE(a.num_bits() == 101);
    REQUIRE(a.num_trailing_zeros() == 100);
    REQUIRE(is_power_of_two(a));
}


TEST_CASE("literal") {
    natural a = 1'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890_n;
    REQUIRE(a.str() == "1234567890234567890234567890234567890234567890234567890234567890234567890234567890");
    REQUIRE(is_power_of_two(1_n << 280));
}

TEST_CASE("power_of_two") {
    REQUIRE(pow(2_n, 0) == 1_n);
    REQUIRE(pow(2_n, 3) == 8_n);
    REQUIRE(pow(2_n, 63) == 1_n << 63);
    REQUIRE(pow(2_n, 64) == 1_n << 64);
}

TEST_CASE("merseinne primes") {
    std::vector<int> p = {2, 3, 5, 7, 13, 17, 19, 31, 61};
    for (int i = 0; i <= 70; i++) {
        REQUIRE(is_prime(pow(2_n, i) - 1) == (std::find(p.begin(), p.end(), i) != p.end()));
    }
}
