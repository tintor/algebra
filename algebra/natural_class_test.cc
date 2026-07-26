#include "algebra/natural_class.h"
#include "algebra/__test.h"
#include <catch2/benchmark/catch_benchmark.hpp>

TEST_CASE("div 10") {
    natural a = 10;
    int b = 10;
    REQUIRE(div(a, b, a) == 0);
}

TEST_CASE("div test") {
    const natural a = {0, 1'000'000};
    const natural b = {0, 1};
    natural q, r;
    div(a, b, q, r);
    REQUIRE(q > 0);
    REQUIRE(b * q + r == a);
    REQUIRE(r < b);
}

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

natural rand_natural(int size, Random& rng) {
    natural a;
    a.words.reset(size);
    for (int i = 0; i < size; i++)
        a.words[i] = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
    return a;
}

#if 0
TEST_CASE("mul benchmark") {
    Random rng(1);
    natural a, b;

    a = rand_natural(4, rng);
    b = rand_natural(4, rng);
    BENCHMARK("a * b 4") { return a * b; };
    BENCHMARK("karatsuba 4") { return mul_karatsuba(a, b); };

    a = rand_natural(8, rng);
    b = rand_natural(8, rng);
    BENCHMARK("a * b 8") { return a * b; };
    BENCHMARK("karatsuba 8") { return mul_karatsuba(a, b); };

    a = rand_natural(16, rng);
    b = rand_natural(16, rng);
    BENCHMARK("a * b 16") { return a * b; };
    BENCHMARK("karatsuba 16") { return mul_karatsuba(a, b); };

    a = rand_natural(32, rng);
    b = rand_natural(32, rng);
    BENCHMARK("a * b 32") { return a * b; };
    BENCHMARK("karatsuba 32") { return mul_karatsuba(a, b); };

    a = rand_natural(64, rng);
    b = rand_natural(64, rng);
    BENCHMARK("a * b 64") { return a * b; };
    BENCHMARK("karatsuba 64") { return mul_karatsuba(a, b); };

    a = rand_natural(128, rng);
    b = rand_natural(128, rng);
    BENCHMARK("a * b 128") { return a * b; };
    BENCHMARK("karatsuba 128") { return mul_karatsuba(a, b); };

    a = rand_natural(256, rng);
    b = rand_natural(256, rng);
    BENCHMARK("a * b 256") { return a * b; };
    BENCHMARK("karatsuba 256") { return mul_karatsuba(a, b); };

    a = rand_natural(512, rng);
    b = rand_natural(512, rng);
    BENCHMARK("a * b 512") { return a * b; };
    BENCHMARK("karatsuba 512") { return mul_karatsuba(a, b); };

    a = rand_natural(1024, rng);
    b = rand_natural(1024, rng);
    BENCHMARK("a * b 1024") { return a * b; };
    BENCHMARK("karatsuba 1024") { return mul_karatsuba(a, b); };

    a = rand_natural(2048, rng);
    b = rand_natural(2048, rng);
    BENCHMARK("a * b 2048") { return a * b; };
    BENCHMARK("karatsuba 2048") { return mul_karatsuba(a, b); };

    a = rand_natural(4096, rng);
    b = rand_natural(4096, rng);
    BENCHMARK("a * b 4096") { return a * b; };
    BENCHMARK("karatsuba 4096") { return mul_karatsuba(a, b); };

    a = rand_natural(8192, rng);
    b = rand_natural(8192, rng);
    BENCHMARK("a * b 8192") { return a * b; };
    BENCHMARK("karatsuba 8192") { return mul_karatsuba(a, b); };

    a = rand_natural(16384, rng);
    b = rand_natural(16384, rng);
    BENCHMARK("a * b 16384") { return a * b; };
    BENCHMARK("karatsuba 16384") { return mul_karatsuba(a, b); };
}
#endif

TEST_CASE("mul_karatsuba easy") {
    natural a;
    natural b;
    a.words.reset(4);
    a.words[0] = 1;
    a.words[1] = 2;
    a.words[2] = 3;
    a.words[3] = 4;
    b.words.reset(4);
    b.words[0] = 5;
    b.words[1] = 6;
    b.words[2] = 7;
    b.words[3] = 8;
    REQUIRE(a * b == mul_karatsuba(a, b));
}

TEST_CASE("mul_karatsuba ones") {
    std::vector<natural> p;
    natural a = 1;
    for (int i = 0; i < 256; i++) {
        p.push_back(a);
        a <<= 1;
    }
    Random rng(0);
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            const natural& a = p.at(i);
            const natural& b = p.at(j);
            REQUIRE(a * b == mul_karatsuba(a, b));
        }
    }
}

#if 0
TEST_CASE("mul_karatsuba 4") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        natural a = rand_natural(4, rng);
        natural b = rand_natural(4, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 8") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        natural a = rand_natural(8, rng);
        natural b = rand_natural(8, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 16") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        natural a = rand_natural(16, rng);
        natural b = rand_natural(16, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 32") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        natural a = rand_natural(32, rng);
        natural b = rand_natural(32, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 64") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        natural a = rand_natural(64, rng);
        natural b = rand_natural(64, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 128") {
    Random rng(0);
    for (int i = 0; i < 5'000'000; i++) {
        natural a = rand_natural(128, rng);
        natural b = rand_natural(128, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 256") {
    Random rng(0);
    for (int i = 0; i < 2'500'000; i++) {
        natural a = rand_natural(256, rng);
        natural b = rand_natural(256, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba general") {
    Random rng(0);
    for (int i = 0; i < 1'000'000; i++) {
        natural a = rand_natural(rng.Uniform<int>(0, 512), rng);
        natural b = rand_natural(rng.Uniform<int>(0, 512), rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}
#endif

TEST_CASE("__less_a_bc_scalar") {
    Random rng(555);
    natural a, b;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 2, rng);
        rand_natural(b, 1, 1, rng);
        uint64_t c = rng.Uniform<uint64_t>(0, UINT64_MAX);
        REQUIRE(__less_a_bc_scalar(a, b, c) == (a < b * c));
    }
}

TEST_CASE("__less_a_bc") {
    Random rng(666);
    natural a, b, c;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 10, rng);
        rand_natural(b, 1, 5, rng);
        rand_natural(c, 1, 5, rng);
        REQUIRE(__less_a_bc(a, b, c) == a < b * c);
    }
}

TEST_CASE("__less_ab_c") {
    Random rng(999);
    natural a, b, c;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 5, rng);
        rand_natural(b, 1, 5, rng);
        rand_natural(c, 1, 10, rng);
        REQUIRE(__less_ab_c(a, b, c) == a * b < c);
    }
}

TEST_CASE("__less_ab_cd") {
    Random rng(888);
    natural a, b, c, d;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 5, rng);
        rand_natural(b, 1, 5, rng);
        rand_natural(c, 1, 5, rng);
        rand_natural(d, 1, 5, rng);
        REQUIRE(__less_ab_cd(a, b, c, d) == a * b < c * d);
    }
}

constexpr int signum(const natural& a, const natural& b) {
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

TEST_CASE("__det_ab_cd") {
    Random rng(777);
    natural a, b, c, d;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 5, rng);
        rand_natural(b, 1, 5, rng);
        rand_natural(c, 1, 5, rng);
        rand_natural(d, 1, 5, rng);
        REQUIRE(__det_ab_cd(a, b, c, d) == signum(a * b, c * d));
    }
}

TEST_CASE("__saturated_div") {
    Random rng(0);
    natural a, b;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 10, rng);
        rand_natural(b, 1, 10, rng);
        REQUIRE(a * b == b * a);
        const uint64_t q = __saturated_div(b, a);
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
    uint128_t a = 1;
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

TEST_CASE("-= 2") {
    natural       a = {0, 4, 0, 1};
    const natural b = {1, 2, 1};
    const natural c = {UINT64_MAX, 1, UINT64_MAX};
    a -= b;
    if (a != c) {
        print("a={}\n", stre(a));
        print("c={}\n", stre(c));
    }
    REQUIRE(a == c);
}

TEST_CASE("-= 3") {
    natural       a = {0, 0, 1};
    const natural b = {2};
    const natural c = {UINT64_MAX - 1, UINT64_MAX};
    a -= b;
    if (a != c) {
        print("a={}\n", stre(a));
        print("c={}\n", stre(c));
    }
    REQUIRE(a == c);
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
    Random rng(6);
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        ucent b = rng.Uniform<ucent>(0, m - a);
        REQUIRE(natural(a) + natural(b) == a + b);
    }
}

TEST_CASE("sub stress with ucent") {
    Random rng(5);
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
    Random rng(4);
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(1, m);
        ucent b = rng.Uniform<ucent>(0, m / a);
        REQUIRE(natural(a) * natural(b) == a * b);
    }
}

TEST_CASE("div stress with ucent") {
    Random rng(3);
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
    Random rng(1);
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
    Random rng(2);
    const ucent m = std::numeric_limits<uint64_t>::max();
    for (int i = 0; i < 1'000'000; i++) {
        ucent a = rng.Uniform<ucent>(m + 1, m * 10);
        natural q;
        uint64_t r = div(natural(a), 10ull, q);
        REQUIRE(q * 10 + r == a);
    }
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
        square(m);
        REQUIRE(m == a * a);
    }
}

TEST_CASE("stress div with remainder") {
    Random rng(10);
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

// TODO randomized long division test against cpp_int for big naturals!

TEST_CASE("is_x") {
    natural z = 0;
    REQUIRE(z.is_uint32());
    REQUIRE(z.is_uint64());

    natural o = 1;
    REQUIRE(o.is_uint32());
    REQUIRE(o.is_uint64());

    natural p = std::numeric_limits<uint>::max();
    REQUIRE(p.is_uint32());
    REQUIRE(p.is_uint64());

    natural q = (uint64_t)std::numeric_limits<uint>::max() + 1;
    REQUIRE(!q.is_uint32());
    REQUIRE(q.is_uint64());

    natural a = std::numeric_limits<long>::max();
    REQUIRE(!a.is_uint32());
    REQUIRE(a.is_uint64());

    natural b = (uint64_t)std::numeric_limits<long>::max() + 1;
    REQUIRE(!b.is_uint32());
    REQUIRE(b.is_uint64());

    natural c = std::numeric_limits<uint64_t>::max();
    REQUIRE(!c.is_uint32());
    REQUIRE(c.is_uint64());

    natural d = c + 1;
    REQUIRE(!d.is_uint32());
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

TEST_CASE("factorial") {
    uint64_t _a = 1;
    uint64_t _b = 1;
    inatural aa {&_a, 1};
    cnatural bb {&_b, 1};
    auto carry = __add_and_return_carry(aa, bb);
    REQUIRE(_a == 2);
    REQUIRE(carry == 0);

    natural a(1);
    natural b = a;
    b += a;
    REQUIRE(b == 2);

    for (int i = 2; i <= 50; i++) {
        natural b = a;
        for (int j = 1; j < i; j++)
            b += a;
        if (a * i != b) {
            print("a={}\n", str(a));
            print("b={}\n", str(b));
            print("i={}\n", i);
        }
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

TEST_CASE("popcount") {
    for (uint i: {0, 4, 31231, -3123121})
        REQUIRE(natural(i).popcount() == std::popcount(i));
}

TEST_CASE("<<") {
    natural(1) << 64; // regression test
}

TEST_CASE("literal") {
    natural a = 1'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890_n;
    REQUIRE(a.str() == "1234567890234567890234567890234567890234567890234567890234567890234567890234567890");
}

TEST_CASE("add_product") {
    natural a = 1;
    add_product(a, 2_n, 3_n);
    REQUIRE(a.words.size() == 1);
    REQUIRE(a == 7);
    add_product(a, 0_n, 3_n);
    REQUIRE(a == 7);
    add_product(a, 2_n, 0_n);
    REQUIRE(a == 7);
    add_product(a, 2_n, 1_n);
    REQUIRE(a == 9);
    add_product(a, 1_n, 3_n);
    REQUIRE(a == 12);
}

TEST_CASE("sub_product") {
    natural a = 10;
    sub_product(a, 2_n, 3_n);
    REQUIRE(a == 4);

    natural b = static_cast<uint128_t>(UINT64_MAX) + 1;
    sub_product(b, 2_n, 3_n);
    REQUIRE(b == UINT64_MAX - 5);

    a = 340282366920938463463412908294782434869_n;
    b = 5850;
    uint64_t d = 2746327603956567807;
    natural e;
    e = a;
    sub_product(e, b, d);
    REQUIRE(e == a - b * d);

    a = 642925181765695293749472128589009400496477258957537017768413464458702803216754593095075074170080266581351009920531714606644854070062962_n;
    b = 3067758199959723904076027525013189935007362036353598278716_n;
    natural c = 105070284311530410717944298959725463208120064695905449166363534139015629378424_n;
    e = a;
    sub_product(e, b, c);
    REQUIRE(e == a - b * c);
}

TEST_CASE("add_product scalar") {
    natural a = 1;
    add_product(a, 2_n, static_cast<uint64_t>(3));
    REQUIRE(a == 7);
}

TEST_CASE("sub_product scalar") {
    natural a = 10;
    sub_product(a, 2_n, static_cast<uint64_t>(3));
    REQUIRE(a == 4);
}

TEST_CASE("add/sub_product scalar stress") {
    Random rng(31231);
    for (int i = 0; i < 500'000; i++) {
        natural a = rand_natural(1, 8, rng);
        natural b = rand_natural(1, 4, rng);
        uint64_t d = rng.Uniform<uint64_t>(0, INT64_MAX);
        natural e;

        e = a;
        add_product(e, b, d);
        REQUIRE(e == a + b * d);

        e = a;
        sub_product(e, b, 0);
        REQUIRE(e == a);

        if (a >= b) {
            e = a;
            sub_product(e, b, 1);
            REQUIRE(e == a - b);
        }

        if (a >= b * d) {
            e = a;
            sub_product(e, b, d);
            REQUIRE(e == a - b * d);
        }
    }
}

TEST_CASE("add/sub_product stress") {
    Random rng(31231);
    for (int i = 0; i < 1000'000; i++) {
        natural a = rand_natural(1, 8, rng);
        natural b = rand_natural(1, 4, rng);
        natural c = rand_natural(1, 4, rng);
        natural e;

        e = a;
        add_product(e, b, c);
        REQUIRE(e == a + b * c);

        if (a >= b * c) {
            e = a;
            sub_product(e, b, c);
            if (e != a - b * c) {
                print("a={}\n{}\n", a, stre(a));
                print("b={}\n{}\n", b, stre(b));
                print("c={}\n{}\n", c, stre(c));
                print("e={}\n{}\n", e, stre(e));
                print("#={}\n{}\n", a-b*c, stre(a-b*c));
            }
            REQUIRE(e == a - b * c);
        }

        e = a;
        sub_product(e, b, 0);
        REQUIRE(e == a);

        if (a >= b) {
            e = a;
            sub_product(e, b, 1);
            REQUIRE(e == a - b);
        }
    }
}

TEST_CASE("operator-= uint64 keeps high words") {
    natural a = 1;
    a <<= 64;
    a += 5u; // 2**64 + 5
    natural b = a;
    b -= static_cast<uint64_t>(5);
    REQUIRE(b == (natural(1) << 64));
    natural c = a;
    c -= static_cast<uint64_t>(6);
    REQUIRE(c == (natural(1) << 64) - 1u);
}

TEST_CASE("operator% uint128 with more than two words") {
    natural a = 1;
    a <<= 200;
    a += 12345u; // 2**200 + 12345
    natural bn = 1;
    bn <<= 100;
    bn += 7u; // 2**100 + 7
    const uint128_t b = (static_cast<uint128_t>(1) << 100) + 7;

    REQUIRE(natural(a % b) == a % bn);

    // and a case where the modulus needs the full 128 bits
    natural c = 1;
    c <<= 300;
    c -= 1u;
    natural dn = 1;
    dn <<= 127;
    dn += 12345u;
    const uint128_t d = (static_cast<uint128_t>(1) << 127) + 12345;
    REQUIRE(natural(c % d) == c % dn);
}

TEST_CASE("operator-- normalizes") {
    natural a = 1;
    --a;
    REQUIRE(a.words.size() == 0);
    REQUIRE(!a);
    REQUIRE(a == 0u);

    natural b = 1;
    b <<= 64; // 2**64
    --b;
    REQUIRE(b.words.size() == 1);
    REQUIRE(b == UINT64_MAX);

    natural c = 1;
    c <<= 128;
    --c;
    REQUIRE(c.words.size() == 2);
}

TEST_CASE("mul_karatsuba with power of two operand") {
    natural a = 1;
    a <<= 1300;
    a += 2u; // not a power of two, and a.words[10] == 0

    natural b = 1;
    b <<= 640; // power of two, 11 words

    REQUIRE(mul_karatsuba(a, b) == a * b);
    REQUIRE(mul_karatsuba(b, a) == a * b);

    natural c = 1;
    c <<= 64;
    REQUIRE(mul_karatsuba(a, c) == a * c);
    REQUIRE(mul_karatsuba(c, a) == a * c);

    natural d = 3;
    REQUIRE(mul_karatsuba(a, d) == a * d);
    REQUIRE(mul_karatsuba(d, a) == a * d);
}

TEST_CASE("mul_karatsuba with sparse operands") {
    for (int w : {33, 40, 64, 100}) {
        natural a = 1;
        a <<= 64 * w;
        a += 1u; // interior words are all zero
        natural b = 1;
        b <<= 64 * w;
        b += 3u;
        REQUIRE(mul_karatsuba(a, b) == a * b);
    }

    natural a = 1;
    a <<= 64 * 50;
    a += 1u;
    a <<= 64 * 50;
    a += 7u;
    natural b = 5;
    b <<= 64 * 70;
    b += 9u;
    REQUIRE(mul_karatsuba(a, b) == a * b);
    REQUIRE(mul_karatsuba(b, a) == a * b);
}

TEST_CASE("divide_bz matches div") {
    Random rng(7);
    natural q1, r1, q2, r2;

    // small / trivial cases
    for (auto [an, dn] : {std::pair{1, 1}, {2, 1}, {5, 3}, {9, 4}, {20, 8}, {33, 8}, {40, 16}, {64, 9}, {17, 17}, {18, 17}}) {
        for (int rep = 0; rep < 8; rep++) {
            natural a = rand_natural(an, an, rng);
            natural d = rand_natural(dn, dn, rng);
            if (!d)
                continue;
            div(a, d, q1, r1);
            divide_bz(a, d, q2, r2);
            REQUIRE(q1 == q2);
            REQUIRE(r1 == r2);
            REQUIRE(q2 * d + r2 == a);
            REQUIRE(r2 < d);
        }
    }

    // random sizes
    for (int rep = 0; rep < 60; rep++) {
        natural a = rand_natural(1, 60, rng);
        natural d = rand_natural(1, 30, rng);
        if (!d)
            continue;
        div(a, d, q1, r1);
        divide_bz(a, d, q2, r2);
        REQUIRE(q1 == q2);
        REQUIRE(r1 == r2);
    }

    // divisor with a single top bit set (shift == 0 path) and with a low top bit
    natural a = rand_natural(40, 40, rng);
    natural d = 1;
    d <<= 64 * 8;
    divide_bz(a, d, q2, r2);
    div(a, d, q1, r1);
    REQUIRE(q1 == q2);
    REQUIRE(r1 == r2);

    d = rand_natural(8, 8, rng);
    d.words.back() = 1; // maximal normalization shift
    divide_bz(a, d, q2, r2);
    div(a, d, q1, r1);
    REQUIRE(q1 == q2);
    REQUIRE(r1 == r2);

    // deeper recursion
    for (auto [an, dn] : {std::pair{200, 40}, {130, 32}, {300, 64}}) {
        natural aa = rand_natural(an, an, rng);
        natural dd = rand_natural(dn, dn, rng);
        div(aa, dd, q1, r1);
        divide_bz(aa, dd, q2, r2);
        REQUIRE(q1 == q2);
        REQUIRE(r1 == r2);
    }
}

TEST_CASE("compare with uint128") {
    natural big = 1;
    big <<= 100;

    REQUIRE(static_cast<uint128_t>(5) < big);
    REQUIRE(!(big < static_cast<uint128_t>(5)));
    REQUIRE(static_cast<uint128_t>(0) < big);

    natural small = 5;
    REQUIRE(!(static_cast<uint128_t>(5) < small));
    REQUIRE(static_cast<uint128_t>(4) < small);
    REQUIRE(!(static_cast<uint128_t>(6) < small));

    const uint128_t huge = (static_cast<uint128_t>(1) << 100) + 1;
    REQUIRE(big < huge);
    REQUIRE(!(huge < big));
    REQUIRE(small < huge);
    REQUIRE(!(huge < small));
}

// Counts outstanding new[] allocations, to detect leaks in integer_backend.
static int64_t g_array_allocs = 0;
void* operator new[](std::size_t n) {
    void* p = std::malloc(n ? n : 1);
    if (!p)
        throw std::bad_alloc();
    g_array_allocs += 1;
    return p;
}
void operator delete[](void* p) noexcept {
    if (p) {
        g_array_allocs -= 1;
        std::free(p);
    }
}
void operator delete[](void* p, std::size_t) noexcept { operator delete[](p); }

TEST_CASE("move assignment does not leak") {
    const int64_t before = g_array_allocs;
    {
        natural a = 1;
        a <<= 200; // heap allocated
        natural b = 1;
        b <<= 300; // heap allocated
        a = std::move(b);
        REQUIRE(a == (natural(1) << 300));
    }
    REQUIRE(g_array_allocs == before);

    {
        natural a = 1;
        a <<= 200;
        natural b = 7; // small buffer, no allocation
        a = std::move(b);
        REQUIRE(a == 7u);
    }
    REQUIRE(g_array_allocs == before);
}

TEST_CASE("subtraction below zero throws") {
    REQUIRE_THROWS(natural(3) - natural(5));
    REQUIRE_THROWS(natural(3) - 5u);
    REQUIRE_THROWS(natural(3) - static_cast<uint64_t>(5));
    REQUIRE_THROWS(3u - natural(5));
    REQUIRE_THROWS(static_cast<uint128_t>(3) - natural(5));
    REQUIRE_THROWS((natural(1) << 64) - (natural(1) << 65));
    REQUIRE_THROWS(natural(0) - 1u);

    natural a = 1;
    a <<= 100;
    REQUIRE_THROWS((static_cast<uint128_t>(1) << 90) - a);

    {
        natural b = 5;
        REQUIRE_THROWS(b -= natural(6));
    }
    {
        natural b = 5;
        REQUIRE_THROWS(b -= static_cast<uint64_t>(6));
    }

    // valid subtractions keep working
    REQUIRE(natural(5) - natural(3) == 2u);
    REQUIRE(natural(5) - 3u == 2u);
    REQUIRE(5u - natural(3) == 2u);
    REQUIRE((natural(1) << 65) - (natural(1) << 64) == (natural(1) << 64));
    REQUIRE(static_cast<uint128_t>((static_cast<uint128_t>(1) << 100) - natural(5)) == (static_cast<uint128_t>(1) << 100) - 5);
    REQUIRE(natural(5) - natural(5) == 0u);
}

TEST_CASE("division by zero throws") {
    natural a = 1;
    a <<= 200;
    a += 12345u;
    const natural zero = 0;
    natural q, r;

    REQUIRE_THROWS(a / zero);
    REQUIRE_THROWS(a % zero);
    REQUIRE_THROWS(div(a, zero, q, r));
    REQUIRE_THROWS(mod(a, zero, r));
    q = a;
    REQUIRE_THROWS(mod(q, zero));
    REQUIRE_THROWS(div(a, static_cast<uint64_t>(0), q));
    REQUIRE_THROWS(a / 0);
    REQUIRE_THROWS(a % 0);
    q = a;
    REQUIRE_THROWS(q /= zero);
    q = a;
    REQUIRE_THROWS(q /= 0);

    natural small = 7;
    REQUIRE_THROWS(small / zero);
    REQUIRE_THROWS(small % zero);
}

// simple digit at a time reference for natural::str()
static std::string ref_str(natural a, unsigned base) {
    if (!a)
        return "0";
    std::string s;
    while (a)
        s += "0123456789ABCDEF"[div(a, static_cast<uint64_t>(base), /*out*/a)];
    std::reverse(s.begin(), s.end());
    return s;
}

TEST_CASE("str matches digit at a time conversion") {
    REQUIRE(natural(0).str() == "0");
    REQUIRE(natural(1).str() == "1");
    REQUIRE(natural(10).str() == "10");
    REQUIRE(natural(19).str(10) == "19");

    // exactly at a chunk boundary: 10**19 and 10**19 - 1
    natural chunk = 10;
    for (int i = 1; i < 19; i++)
        chunk *= 10u;
    REQUIRE(chunk.str() == "10000000000000000000");
    REQUIRE((chunk - 1u).str() == "9999999999999999999");
    REQUIRE((chunk + 1u).str() == "10000000000000000001");
    REQUIRE((chunk * chunk).str() == ref_str(chunk * chunk, 10));

    Random rng(21);
    for (int i = 0; i < 100; i++) {
        const natural a = rand_natural(1, 12, rng);
        for (unsigned base : {2u, 3u, 7u, 8u, 10u, 15u, 16u})
            REQUIRE(a.str(base) == ref_str(a, base));
    }
    for (int i = 0; i < 5; i++) {
        const natural a = rand_natural(60, 80, rng);
        REQUIRE(a.str() == ref_str(a, 10));
        REQUIRE(a.str(7) == ref_str(a, 7));
    }
}

TEST_CASE("bitwise and with operands of different size") {
    natural b = 1;
    b <<= 200;
    b += 0xFFFFu;

    natural a = 0xF0F0F0F0ull;
    natural c = a;
    c &= b;
    REQUIRE(c == 0xF0F0u);
    REQUIRE(c == (a & b));

    natural d = b;
    d &= a;
    REQUIRE(d == 0xF0F0u);

    natural e = b;
    e &= natural(0);
    REQUIRE(e == 0u);
    REQUIRE(e.words.size() == 0);

    Random rng(31);
    for (int i = 0; i < 100; i++) {
        const natural x = rand_natural(1, 8, rng);
        const natural y = rand_natural(1, 8, rng);
        natural z = x;
        z &= y;
        REQUIRE(z == (x & y));
        REQUIRE(z <= x);
        REQUIRE(z <= y);
    }
}


TEST_CASE("sub_product rejects a violated precondition") {
    natural b = 1;
    b <<= 200;
    natural c = 3;

    natural a = 5; // much smaller than b * c
    REQUIRE_THROWS(sub_product(a, b, c));

    natural a2 = 5;
    REQUIRE_THROWS(sub_product(a2, b, static_cast<uint64_t>(7)));

    // valid uses keep working
    natural d = b * c;
    d += 11u;
    sub_product(d, b, c);
    REQUIRE(d == 11u);

    natural e = b;
    e *= 7u;
    e += 3u;
    sub_product(e, b, static_cast<uint64_t>(7));
    REQUIRE(e == 3u);
}

TEST_CASE("multiplication propagates a long carry") {
    // isqrt() produced this value; squaring it needs a carry that travels through
    // two all-ones words of the partial product
    natural s;
    s.words.reset(5);
    const uint64_t w[5] = {8296379691479107416ull, 8433445100458127462ull, 6328652237515477287ull,
                           14740434617336491689ull, 24879108095803ull};
    for (int i = 0; i < 5; i++)
        s.words[i] = w[i];

    natural expected;
    add_product(expected, s, s); // independent implementation
    REQUIRE(s * s == expected);
    REQUIRE(mul_karatsuba(s, s) == expected);
    natural sq = s;
    square(sq);
    REQUIRE(sq == expected);

    // the same for a * b with b != a
    natural t = s;
    t += 1u;
    natural expected2;
    add_product(expected2, s, t);
    REQUIRE(s * t == expected2);
    REQUIRE(mul_karatsuba(s, t) == expected2);

    // cross check the general multiplication against add_product for values that are
    // likely to produce long carry chains
    Random rng(51);
    for (int i = 0; i < 300; i++) {
        natural a = 1;
        a <<= rng.Uniform<int>(64, 400);
        a -= rng.Uniform<uint64_t>(1, 1000);
        natural b = 1;
        b <<= rng.Uniform<int>(64, 400);
        b -= rng.Uniform<uint64_t>(1, 1000);
        natural e;
        add_product(e, a, b);
        REQUIRE(a * b == e);
    }
}
