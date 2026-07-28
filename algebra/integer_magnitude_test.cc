#include "algebra/integer_class.h"
#include "algebra/__test.h"
#include <catch2/benchmark/catch_benchmark.hpp>

TEST_CASE("div 10") {
    integer a = 10;
    int b = 10;
    REQUIRE(div(a, b, a) == 0);
}

TEST_CASE("div test") {
    const integer a = {0, 1'000'000};
    const integer b = {0, 1};
    integer q, r;
    div(a, b, q, r);
    REQUIRE(q > 0);
    REQUIRE(b * q + r == a);
    REQUIRE(r < b);
}

integer rand_natural(int min_size, int max_size, Random& rng) {
    integer a;
    a.words.reset(rng.Uniform<int>(min_size, max_size));
    for (int i = 0; i < a.words.size(); i++)
        a.words[i] = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
    return a;
}

void rand_natural(integer& a, int min_size, int max_size, Random& rng) {
    a.words.reset(rng.Uniform<int>(min_size, max_size));
    for (int i = 0; i < a.words.size(); i++)
        a.words[i] = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
}

integer rand_natural(int size, Random& rng) {
    integer a;
    a.words.reset(size);
    for (int i = 0; i < size; i++)
        a.words[i] = rng.Uniform<uint64_t>(0, std::numeric_limits<uint64_t>::max());
    return a;
}

#if 0
TEST_CASE("mul benchmark") {
    Random rng(1);
    integer a, b;

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
    integer a;
    integer b;
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
    std::vector<integer> p;
    integer a = 1;
    for (int i = 0; i < 256; i++) {
        p.push_back(a);
        a <<= 1;
    }
    Random rng(0);
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            const integer& a = p.at(i);
            const integer& b = p.at(j);
            REQUIRE(a * b == mul_karatsuba(a, b));
        }
    }
}

#if 0
TEST_CASE("mul_karatsuba 4") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        integer a = rand_natural(4, rng);
        integer b = rand_natural(4, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 8") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        integer a = rand_natural(8, rng);
        integer b = rand_natural(8, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 16") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        integer a = rand_natural(16, rng);
        integer b = rand_natural(16, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 32") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        integer a = rand_natural(32, rng);
        integer b = rand_natural(32, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 64") {
    Random rng(0);
    for (int i = 0; i < 10'000'000; i++) {
        integer a = rand_natural(64, rng);
        integer b = rand_natural(64, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 128") {
    Random rng(0);
    for (int i = 0; i < 5'000'000; i++) {
        integer a = rand_natural(128, rng);
        integer b = rand_natural(128, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba 256") {
    Random rng(0);
    for (int i = 0; i < 2'500'000; i++) {
        integer a = rand_natural(256, rng);
        integer b = rand_natural(256, rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}

TEST_CASE("mul_karatsuba general") {
    Random rng(0);
    for (int i = 0; i < 1'000'000; i++) {
        integer a = rand_natural(rng.Uniform<int>(0, 512), rng);
        integer b = rand_natural(rng.Uniform<int>(0, 512), rng);
        REQUIRE(a * b == mul_karatsuba(a, b));
    }
}
#endif

TEST_CASE("__less_a_bc_scalar") {
    Random rng(555);
    integer a, b;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 2, rng);
        rand_natural(b, 1, 1, rng);
        uint64_t c = rng.Uniform<uint64_t>(0, UINT64_MAX);
        REQUIRE(__less_a_bc_scalar(a, b, c) == (a < b * c));
    }
}

TEST_CASE("__less_a_bc") {
    Random rng(666);
    integer a, b, c;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 10, rng);
        rand_natural(b, 1, 5, rng);
        rand_natural(c, 1, 5, rng);
        REQUIRE(__less_a_bc(a, b, c) == a < b * c);
    }
}

TEST_CASE("__less_ab_c") {
    Random rng(999);
    integer a, b, c;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 5, rng);
        rand_natural(b, 1, 5, rng);
        rand_natural(c, 1, 10, rng);
        REQUIRE(__less_ab_c(a, b, c) == a * b < c);
    }
}

TEST_CASE("__less_ab_cd") {
    Random rng(888);
    integer a, b, c, d;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 5, rng);
        rand_natural(b, 1, 5, rng);
        rand_natural(c, 1, 5, rng);
        rand_natural(d, 1, 5, rng);
        REQUIRE(__less_ab_cd(a, b, c, d) == a * b < c * d);
    }
}

constexpr int signum(const integer& a, const integer& b) {
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

TEST_CASE("__det_ab_cd") {
    Random rng(777);
    integer a, b, c, d;
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
    integer a, b;
    for (int i = 0; i < 1000'000; i++) {
        rand_natural(a, 1, 10, rng);
        rand_natural(b, 1, 10, rng);
        REQUIRE(a * b == b * a);
        const uint64_t q = __saturated_div(b, a);
        REQUIRE(a * q <= b);
        if (q != std::numeric_limits<uint64_t>::max())
            REQUIRE(a * (integer(q) + 1) > b);
    }
}

TEST_CASE("words") {
    integer b(1);
    REQUIRE(b.words.size() == 1);

    REQUIRE(integer(0).words.size() == 0);
    REQUIRE(integer(1).words.size() == 1);
    REQUIRE(integer(2).words.size() == 1);

    REQUIRE(integer(1).words[0] == 1);
    REQUIRE(integer(0).words[0] == 0);
}

TEST_CASE("str") {
    REQUIRE(integer(0).str() == "0");
    REQUIRE(integer(1).str() == "1");
    REQUIRE(integer(12).str() == "12");
    REQUIRE(integer(450).str() == "450");
}

constexpr std::string os() {
    std::ostringstream s;
    s << 15_i;
    return s.str();
}

TEST_CASE("constexpr ostream") {
    REQUIRE(os() == "15");
}

TEST_CASE("format") {
    REQUIRE(format("{:b}", 15_i) == "1111");
    REQUIRE(format("{:o}", 15_i) == "17");
    REQUIRE(format("{:d}", 15_i) == "15");
    REQUIRE(format("{:x}", 15_i) == "f");
    REQUIRE(format("{:X}", 15_i) == "F");
    REQUIRE(format("{:4d}", 15_i) == "  15");
    REQUIRE(format("{:1d}", 15_i) == "15");
    REQUIRE(format("{:*>4d}", 15_i) == "**15");
    REQUIRE(format("{:*<4d}", 15_i) == "15**");
    REQUIRE(format("{:*^4d}", 15_i) == "*15*");
}

TEST_CASE("hex") {
    REQUIRE(integer(0).hex() == "0");
    REQUIRE(integer(16).hex() == "10");
    REQUIRE(integer(255).hex() == "FF");
    REQUIRE(integer(256).hex() == "100");
}

TEST_CASE("parse") {
    REQUIRE(integer("0") == integer(0));
    REQUIRE(integer("1") == integer(1));
    REQUIRE(integer("12") == integer(12));
    REQUIRE(integer("450") == integer(450));
    const char* a = "18446744073709551617"; // UINT64_MAX + 2
    REQUIRE(integer(a).str() == a);

    REQUIRE(integer("1100", 2) == integer(12));
    REQUIRE(integer("111", 2) == integer(7));
    REQUIRE(integer("FF", 16) == integer(255));
    REQUIRE(integer("ff", 16) == integer(255));
}

TEST_CASE("static_cast<uint>") {
    REQUIRE(static_cast<uint>(integer(0)) == 0);
    REQUIRE(static_cast<uint>(integer(1)) == 1);
    uint a = std::numeric_limits<uint>::max();
    REQUIRE(static_cast<uint>(integer(a)) == a);
}

TEST_CASE("static_cast<uint64_t>") {
    REQUIRE(static_cast<uint64_t>(integer(0)) == 0);
    REQUIRE(static_cast<uint64_t>(integer(1)) == 1);
    uint64_t a = std::numeric_limits<uint64_t>::max();
    REQUIRE(static_cast<uint64_t>(integer(a)) == a);
}

TEST_CASE("ucent") {
    uint128_t a = 1;
    for (int i = 0; i < 128; i++) {
        integer b(a);
        REQUIRE(b == a);
        REQUIRE(a == b);
        a <<= 1;
    }
}

TEST_CASE("is_even / is_odd") {
    REQUIRE(integer(0).is_even());
    REQUIRE(!integer(0).is_odd());
    REQUIRE(!integer(7).is_even());
    REQUIRE(integer(7).is_odd());
}

TEST_CASE("cmp") {
    REQUIRE(integer(0) == integer(0));
    REQUIRE(integer(5) == integer(5));

    REQUIRE(integer(0) < integer(1));
    REQUIRE(integer(5) < integer(6));

    REQUIRE(integer(6) > integer(5));

    REQUIRE(integer(1) <= integer(1));
    REQUIRE(0_i <= 888089631791237197_i);
}

TEST_CASE("add") {
    REQUIRE(integer(0) + integer(0) == integer(0));
    REQUIRE(integer(5) + integer(0) == integer(5));
    REQUIRE(integer(5) + integer(6) == integer(11));
}

TEST_CASE("sub") {
    REQUIRE(integer(0) - integer(0) == integer(0));
    REQUIRE(integer(5) - integer(0) == integer(5));
    REQUIRE(integer(5) - integer(5) == integer(0));
    REQUIRE(integer(5) - integer(4) == integer(1));
}

TEST_CASE("+=") {
    integer a(5);
    a += integer(4);
    REQUIRE(a == integer(9));
}

TEST_CASE("-=") {
    integer a(5);
    a -= integer(4);
    REQUIRE(a == integer(1));
}

TEST_CASE("-= 2") {
    integer       a = {0, 4, 0, 1};
    const integer b = {1, 2, 1};
    const integer c = {UINT64_MAX, 1, UINT64_MAX};
    a -= b;
    if (a != c) {
        print("a={}\n", stre(a));
        print("c={}\n", stre(c));
    }
    REQUIRE(a == c);
}

TEST_CASE("-= 3") {
    integer       a = {0, 0, 1};
    const integer b = {2};
    const integer c = {UINT64_MAX - 1, UINT64_MAX};
    a -= b;
    if (a != c) {
        print("a={}\n", stre(a));
        print("c={}\n", stre(c));
    }
    REQUIRE(a == c);
}

TEST_CASE("big add 1") {
    uint64_t a = std::numeric_limits<uint64_t>::max();
    integer b = a;
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

    integer a;
    a.words.reset(4);
    a.words[0] = m;
    a.words[1] = m;
    a.words[2] = m;
    a.words[3] = 1;

    integer b;
    b.words.reset(4);
    b.words[0] = 1;
    b.words[1] = 0;
    b.words[2] = 0;
    b.words[3] = 1;

    integer c = a + b;
    REQUIRE(c.words.size() == 4);
    REQUIRE(c.words[0] == 0);
    REQUIRE(c.words[1] == 0);
    REQUIRE(c.words[2] == 0);
    REQUIRE(c.words[3] == 3);
}

TEST_CASE("mul") {
    REQUIRE(integer(0) * integer(0) == integer(0));
    REQUIRE(integer(5) * integer(0) == integer(0));
    REQUIRE(integer(0) * integer(2) == integer(0));
    REQUIRE(integer(5) * integer(2) == integer(10));
}

TEST_CASE("*=") {
    integer a(5);
    a *= integer(3);
    REQUIRE(a == integer(15));
    a *= 3;
    REQUIRE(a == 45);

    integer c(3);
    c *= c;
    REQUIRE(c == 9);

    integer d(3), e("1000000000000000000000000000000000000");
    REQUIRE(e.words.size() == 2);
    d *= e;
    REQUIRE(d.words.size() == 2);
    REQUIRE(d == integer("3000000000000000000000000000000000000"));
}

TEST_CASE("/") {
    REQUIRE(integer(0) / integer(7) == integer(0));
    REQUIRE(integer(7) / integer(7) == integer(1));
    REQUIRE(integer(7) / integer(3) == integer(2));
    REQUIRE(integer(7) / integer(8) == integer(0));
}

TEST_CASE("%") {
    REQUIRE(integer(0) % 7u == 0);
    REQUIRE(integer(1) % 7u == 1);
    REQUIRE(integer(7) % 7u == 0);
    REQUIRE(integer(8) % 7u == 1);
    REQUIRE(integer(14) % 7u == 0);

    REQUIRE(integer(0) % 7ul == 0);
    REQUIRE(integer(1) % 7ul == 1);
    REQUIRE(integer(7) % 7ul == 0);
    REQUIRE(integer(8) % 7ul == 1);
    REQUIRE(integer(14) % 7ul == 0);

    REQUIRE(integer(5) % 3u == 2);
    REQUIRE(integer(5) % 3ul == 2);
}

TEST_CASE("add stress with ucent") {
    Random rng(6);
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        ucent b = rng.Uniform<ucent>(0, m - a);
        REQUIRE(integer(a) + integer(b) == a + b);
    }
}

TEST_CASE("sub stress with ucent") {
    Random rng(5);
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        ucent b = rng.Uniform<ucent>(0, a);
        REQUIRE(integer(a).str() == format("{}", a));
        REQUIRE(integer(b).str() == format("{}", b));
        REQUIRE(integer(a) - integer(b) == a - b);
    }
}

TEST_CASE("mul stress with ucent") {
    Random rng(4);
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(1, m);
        ucent b = rng.Uniform<ucent>(0, m / a);
        REQUIRE(integer(a) * integer(b) == a * b);
    }
}

TEST_CASE("div stress with ucent") {
    Random rng(3);
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        long b = rng.Uniform<long>(0, std::numeric_limits<long>::max());
        integer q;
        REQUIRE(div(integer(a), b, q) == a % b);
        REQUIRE(q == a / b);
    }
}

TEST_CASE("str stress with ucent") {
    Random rng(1);
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        integer b = a;
        if (a > std::numeric_limits<uint64_t>::max()) {
            REQUIRE(b.words.size() == 2);
            REQUIRE(b.words[0] == uint64_t(a));
            REQUIRE(b.words[1] == uint64_t(a >> 64));
        }
        REQUIRE(format("{}", a) == integer(a).str());
    }
}

TEST_CASE("div10 stress with ucent 2") {
    Random rng(2);
    const ucent m = std::numeric_limits<uint64_t>::max();
    for (int i = 0; i < 1'000'000; i++) {
        ucent a = rng.Uniform<ucent>(m + 1, m * 10);
        integer q;
        uint64_t r = div(integer(a), 10ull, q);
        REQUIRE(q * 10 + r == a);
    }
}

TEST_CASE("stress + and -") {
    Random rng(0);
    for (int i = 0; i < 1000'000; i++) {
        const integer a = rand_natural(rng.Uniform<int>(1, 10), rng);
        const integer b = rand_natural(rng.Uniform<int>(1, 10), rng);
        const integer c = rand_natural(rng.Uniform<int>(1, 10), rng);

        const integer sab = a + b;
        const integer sac = a + c;
        const integer sbc = b + c;
        REQUIRE(sab == b + a);
        REQUIRE(sbc == c + b);
        REQUIRE(sac == c + a);
        const integer s = sab + c;
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
        const integer a = rand_natural(rng.Uniform<int>(1, 10), rng);
        integer m = a;
        m += m;
        REQUIRE(m == a + a);
        REQUIRE(m == a * 2);
    }
}

integer safe_mul(const integer& a, const integer& b) {
    integer out;
    for (int i = 0; i < a.words.size(); i++)
        for (int j = 0; j < b.words.size(); j++) {
            integer e = ucent(a.words[i]) * b.words[j];
            out += e << ((i + j) * 64);
        }
    return out;
}

TEST_CASE("stress mul") {
    Random rng(0);
    for (int i = 0; i < 1000'000; i++) {
        const integer a = rand_natural(rng.Uniform<int>(1, 5), rng);
        const integer b = rand_natural(rng.Uniform<int>(1, 5), rng);
        const integer ab = safe_mul(a, b);

        integer c;
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
    integer a;
    a.words.reset(2);
    a.words[0] = 2;
    a.words[1] = 10;
    integer m = a;
    square(m);
    integer e;
    __mul(a, a, e);
    REQUIRE(m == e);
}

TEST_CASE("stress square in-place") {
    Random rng(0);
    for (int i = 0; i < 100'000; i++) {
        integer a = rand_natural(rng.Uniform<int>(1, 5), rng);
        integer m = a;
        square(m);
        REQUIRE(m == a * a);
    }
}

TEST_CASE("stress div with remainder") {
    Random rng(10);
    for (int i = 0; i < 100'000; i++) {
        const integer a = rand_natural(rng.Uniform<int>(2, 10), rng);

        const integer b = rand_natural(rng.Uniform<int>(1, 5), rng);
        integer quot, rem;
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
    integer z = 0;
    REQUIRE(z.is_uint32());
    REQUIRE(z.is_uint64());

    integer o = 1;
    REQUIRE(o.is_uint32());
    REQUIRE(o.is_uint64());

    integer p = std::numeric_limits<uint>::max();
    REQUIRE(p.is_uint32());
    REQUIRE(p.is_uint64());

    integer q = (uint64_t)std::numeric_limits<uint>::max() + 1;
    REQUIRE(!q.is_uint32());
    REQUIRE(q.is_uint64());

    integer a = std::numeric_limits<long>::max();
    REQUIRE(!a.is_uint32());
    REQUIRE(a.is_uint64());

    integer b = (uint64_t)std::numeric_limits<long>::max() + 1;
    REQUIRE(!b.is_uint32());
    REQUIRE(b.is_uint64());

    integer c = std::numeric_limits<uint64_t>::max();
    REQUIRE(!c.is_uint32());
    REQUIRE(c.is_uint64());

    integer d = c + 1;
    REQUIRE(!d.is_uint32());
    REQUIRE(!d.is_uint64());
}

TEST_CASE("<<=") {
    for (int i = 0; i <= 10; i++) {
        integer a(i);
        a <<= 1;
        REQUIRE(a == integer(i << 1));
    }
}

TEST_CASE(">>=") {
    for (int i = 0; i <= 10; i++) {
        integer a(i);
        a >>= 1;
        REQUIRE(a == integer(i >> 1));
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

    integer a(1);
    integer b = a;
    b += a;
    REQUIRE(b == 2);

    for (int i = 2; i <= 50; i++) {
        integer b = a;
        for (int j = 1; j < i; j++)
            b += a;
        if (a * i != b) {
            print("a={}\n", str(a));
            print("b={}\n", str(b));
            print("i={}\n", i);
        }
        REQUIRE(a * i == b);
        REQUIRE(b / i == a);
        a *= integer(i);
        REQUIRE(a == b);
        if (i == 30) REQUIRE(a.str() == "265252859812191058636308480000000");
        if (i == 50) REQUIRE(a.str() == "30414093201713378043612608166064768844377641568960512000000000000");
    }
}

TEST_CASE("num_bits") {
    REQUIRE(integer(0).num_bits() == 0);
    REQUIRE(integer(1).num_bits() == 1);
    REQUIRE(integer(2).num_bits() == 2);
    REQUIRE(integer(3).num_bits() == 2);
    REQUIRE(integer(4).num_bits() == 3);
    REQUIRE(integer(15).num_bits() == 4);
    REQUIRE(integer(16).num_bits() == 5);
}

TEST_CASE("popcount") {
    for (uint i: {0, 4, 31231, -3123121})
        REQUIRE(integer(i).popcount() == std::popcount(i));
}

TEST_CASE("<<") {
    integer(1) << 64; // regression test
}

TEST_CASE("literal") {
    integer a = 1'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890_i;
    REQUIRE(a.str() == "1234567890234567890234567890234567890234567890234567890234567890234567890234567890");
}

TEST_CASE("add_product") {
    integer a = 1;
    add_product(a, 2_i, 3_i);
    REQUIRE(a.words.size() == 1);
    REQUIRE(a == 7);
    add_product(a, 0_i, 3_i);
    REQUIRE(a == 7);
    add_product(a, 2_i, 0_i);
    REQUIRE(a == 7);
    add_product(a, 2_i, 1_i);
    REQUIRE(a == 9);
    add_product(a, 1_i, 3_i);
    REQUIRE(a == 12);
}

TEST_CASE("sub_product") {
    integer a = 10;
    sub_product(a, 2_i, 3_i);
    REQUIRE(a == 4);

    integer b = static_cast<uint128_t>(UINT64_MAX) + 1;
    sub_product(b, 2_i, 3_i);
    REQUIRE(b == UINT64_MAX - 5);

    a = 340282366920938463463412908294782434869_i;
    b = 5850;
    uint64_t d = 2746327603956567807;
    integer e;
    e = a;
    sub_product(e, b, d);
    REQUIRE(e == a - b * d);

    a = 642925181765695293749472128589009400496477258957537017768413464458702803216754593095075074170080266581351009920531714606644854070062962_i;
    b = 3067758199959723904076027525013189935007362036353598278716_i;
    integer c = 105070284311530410717944298959725463208120064695905449166363534139015629378424_i;
    e = a;
    sub_product(e, b, c);
    REQUIRE(e == a - b * c);
}

TEST_CASE("add_product scalar") {
    integer a = 1;
    add_product(a, 2_i, static_cast<uint64_t>(3));
    REQUIRE(a == 7);
}

TEST_CASE("sub_product scalar") {
    integer a = 10;
    sub_product(a, 2_i, static_cast<uint64_t>(3));
    REQUIRE(a == 4);
}

TEST_CASE("add/sub_product scalar stress") {
    Random rng(31231);
    for (int i = 0; i < 500'000; i++) {
        integer a = rand_natural(1, 8, rng);
        integer b = rand_natural(1, 4, rng);
        uint64_t d = rng.Uniform<uint64_t>(0, INT64_MAX);
        integer e;

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
        integer a = rand_natural(1, 8, rng);
        integer b = rand_natural(1, 4, rng);
        integer c = rand_natural(1, 4, rng);
        integer e;

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
    integer a = 1;
    a <<= 64;
    a += 5u; // 2**64 + 5
    integer b = a;
    b -= static_cast<uint64_t>(5);
    REQUIRE(b == (integer(1) << 64));
    integer c = a;
    c -= static_cast<uint64_t>(6);
    REQUIRE(c == (integer(1) << 64) - 1u);
}

TEST_CASE("operator% uint128 with more than two words") {
    integer a = 1;
    a <<= 200;
    a += 12345u; // 2**200 + 12345
    integer bn = 1;
    bn <<= 100;
    bn += 7u; // 2**100 + 7
    const uint128_t b = (static_cast<uint128_t>(1) << 100) + 7;

    REQUIRE(integer(a.__abs_mod_word(b)) == a % bn); // % by a uint128 magnitude, by name

    // and a case where the modulus needs the full 128 bits
    integer c = 1;
    c <<= 300;
    c -= 1u;
    integer dn = 1;
    dn <<= 127;
    dn += 12345u;
    const uint128_t d = (static_cast<uint128_t>(1) << 127) + 12345;
    REQUIRE(integer(c.__abs_mod_word(d)) == c % dn);
}

TEST_CASE("operator-- normalizes") {
    integer a = 1;
    --a;
    REQUIRE(a.words.size() == 0);
    REQUIRE(!a);
    REQUIRE(a == 0u);

    integer b = 1;
    b <<= 64; // 2**64
    --b;
    REQUIRE(b.words.size() == 1);
    REQUIRE(b == UINT64_MAX);

    integer c = 1;
    c <<= 128;
    --c;
    REQUIRE(c.words.size() == 2);
}

TEST_CASE("mul_karatsuba with power of two operand") {
    integer a = 1;
    a <<= 1300;
    a += 2u; // not a power of two, and a.words[10] == 0

    integer b = 1;
    b <<= 640; // power of two, 11 words

    REQUIRE(mul_karatsuba(a, b) == a * b);
    REQUIRE(mul_karatsuba(b, a) == a * b);

    integer c = 1;
    c <<= 64;
    REQUIRE(mul_karatsuba(a, c) == a * c);
    REQUIRE(mul_karatsuba(c, a) == a * c);

    integer d = 3;
    REQUIRE(mul_karatsuba(a, d) == a * d);
    REQUIRE(mul_karatsuba(d, a) == a * d);
}

TEST_CASE("mul_karatsuba with sparse operands") {
    for (int w : {33, 40, 64, 100}) {
        integer a = 1;
        a <<= 64 * w;
        a += 1u; // interior words are all zero
        integer b = 1;
        b <<= 64 * w;
        b += 3u;
        REQUIRE(mul_karatsuba(a, b) == a * b);
    }

    integer a = 1;
    a <<= 64 * 50;
    a += 1u;
    a <<= 64 * 50;
    a += 7u;
    integer b = 5;
    b <<= 64 * 70;
    b += 9u;
    REQUIRE(mul_karatsuba(a, b) == a * b);
    REQUIRE(mul_karatsuba(b, a) == a * b);
}

TEST_CASE("divide_bz matches div") {
    Random rng(7);
    integer q1, r1, q2, r2;

    // small / trivial cases
    for (auto [an, dn] : {std::pair{1, 1}, {2, 1}, {5, 3}, {9, 4}, {20, 8}, {33, 8}, {40, 16}, {64, 9}, {17, 17}, {18, 17}}) {
        for (int rep = 0; rep < 8; rep++) {
            integer a = rand_natural(an, an, rng);
            integer d = rand_natural(dn, dn, rng);
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
        integer a = rand_natural(1, 60, rng);
        integer d = rand_natural(1, 30, rng);
        if (!d)
            continue;
        div(a, d, q1, r1);
        divide_bz(a, d, q2, r2);
        REQUIRE(q1 == q2);
        REQUIRE(r1 == r2);
    }

    // divisor with a single top bit set (shift == 0 path) and with a low top bit
    integer a = rand_natural(40, 40, rng);
    integer d = 1;
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
        integer aa = rand_natural(an, an, rng);
        integer dd = rand_natural(dn, dn, rng);
        div(aa, dd, q1, r1);
        divide_bz(aa, dd, q2, r2);
        REQUIRE(q1 == q2);
        REQUIRE(r1 == r2);
    }
}

TEST_CASE("compare with uint128") {
    integer big = 1;
    big <<= 100;

    REQUIRE(static_cast<uint128_t>(5) < big);
    REQUIRE(!(big < static_cast<uint128_t>(5)));
    REQUIRE(static_cast<uint128_t>(0) < big);

    integer small = 5;
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
        integer a = 1;
        a <<= 200; // heap allocated
        integer b = 1;
        b <<= 300; // heap allocated
        a = std::move(b);
        REQUIRE(a == (integer(1) << 300));
    }
    REQUIRE(g_array_allocs == before);

    {
        integer a = 1;
        a <<= 200;
        integer b = 7; // small buffer, no allocation
        a = std::move(b);
        REQUIRE(a == 7u);
    }
    REQUIRE(g_array_allocs == before);
}

// Subtraction is signed now, so going below zero is a negative result rather than an error. Only
// the magnitude level subtraction still refuses, since it has nowhere to put the sign.
TEST_CASE("subtraction below zero is signed") {
    REQUIRE(integer(3) - integer(5) == -2);
    REQUIRE(integer(3) - 5u == -2);
    REQUIRE(integer(3) - static_cast<uint64_t>(5) == -2);
    REQUIRE(3u - integer(5) == -2);
    REQUIRE((integer(1) << 64) - (integer(1) << 65) == -(integer(1) << 64));
    REQUIRE(integer(0) - 1u == -1);

    {
        integer b = 5;
        b -= integer(6);
        REQUIRE(b == -1);
    }
    {
        integer b = 5;
        b -= static_cast<uint64_t>(6);
        REQUIRE(b == -1);
    }
    {
        integer b = 5;
        REQUIRE_THROWS([&] { auto m = magnitude(b); __abs_sub(*m, 6u); }());
    }

    // valid subtractions keep working
    REQUIRE(integer(5) - integer(3) == 2u);
    REQUIRE(integer(5) - 3u == 2u);
    REQUIRE(5u - integer(3) == 2u);
    REQUIRE((integer(1) << 65) - (integer(1) << 64) == (integer(1) << 64));
    REQUIRE(static_cast<uint128_t>((static_cast<uint128_t>(1) << 100) - integer(5)) == (static_cast<uint128_t>(1) << 100) - 5);
    REQUIRE(integer(5) - integer(5) == 0u);
}

TEST_CASE("division by zero throws") {
    integer a = 1;
    a <<= 200;
    a += 12345u;
    const integer zero = 0;
    integer q, r;

    REQUIRE_THROWS(a / zero);
    REQUIRE_THROWS(a % zero);
    REQUIRE_THROWS(__abs_div(a, zero, q, r));
    REQUIRE_THROWS(__abs_mod(a, zero, r));
    q = a;
    REQUIRE_THROWS(__abs_mod(q, zero));
    REQUIRE_THROWS(__abs_div(a, static_cast<uint64_t>(0), q));
    REQUIRE_THROWS(a / 0);
    REQUIRE_THROWS(a % 0);
    q = a;
    REQUIRE_THROWS(q /= zero);
    q = a;
    REQUIRE_THROWS(q /= 0);

    integer small = 7;
    REQUIRE_THROWS(small / zero);
    REQUIRE_THROWS(small % zero);
}

// simple digit at a time reference for integer::str()
static std::string ref_str(integer a, unsigned base) {
    if (!a)
        return "0";
    std::string s;
    while (a)
        s += "0123456789ABCDEF"[div(a, static_cast<uint64_t>(base), /*out*/a)];
    std::reverse(s.begin(), s.end());
    return s;
}

TEST_CASE("str matches digit at a time conversion") {
    REQUIRE(integer(0).str() == "0");
    REQUIRE(integer(1).str() == "1");
    REQUIRE(integer(10).str() == "10");
    REQUIRE(integer(19).str(10) == "19");

    // exactly at a chunk boundary: 10**19 and 10**19 - 1
    integer chunk = 10;
    for (int i = 1; i < 19; i++)
        chunk *= 10u;
    REQUIRE(chunk.str() == "10000000000000000000");
    REQUIRE((chunk - 1u).str() == "9999999999999999999");
    REQUIRE((chunk + 1u).str() == "10000000000000000001");
    REQUIRE((chunk * chunk).str() == ref_str(chunk * chunk, 10));

    Random rng(21);
    for (int i = 0; i < 100; i++) {
        const integer a = rand_natural(1, 12, rng);
        for (unsigned base : {2u, 3u, 7u, 8u, 10u, 15u, 16u})
            REQUIRE(a.str(base) == ref_str(a, base));
    }
    for (int i = 0; i < 5; i++) {
        const integer a = rand_natural(60, 80, rng);
        REQUIRE(a.str() == ref_str(a, 10));
        REQUIRE(a.str(7) == ref_str(a, 7));
    }
}

// The bitwise operations are on magnitudes and are spelled __abs_and and friends: integer has no
// two's complement bitwise layer, so an operator spelling would promise something it does not do.
TEST_CASE("bitwise and with operands of different size") {
    integer b = 1;
    b <<= 200;
    b += 0xFFFFu;

    integer a = 0xF0F0F0F0ull;
    integer c = a;
    __abs_and(c, b);
    REQUIRE(c == 0xF0F0u);
    REQUIRE(c == __abs_and_copy(a, b));

    integer d = b;
    __abs_and(d, a);
    REQUIRE(d == 0xF0F0u);

    integer e = b;
    __abs_and(e, integer(0));
    REQUIRE(e == 0u);
    REQUIRE(e.words.size() == 0);

    Random rng(31);
    for (int i = 0; i < 100; i++) {
        const integer x = rand_natural(1, 8, rng);
        const integer y = rand_natural(1, 8, rng);
        integer z = x;
        __abs_and(z, y);
        REQUIRE(z == __abs_and_copy(x, y));
        REQUIRE(z <= x);
        REQUIRE(z <= y);
    }
}


TEST_CASE("sub_product rejects a violated precondition") {
    integer b = 1;
    b <<= 200;
    integer c = 3;

    // sub_product is the sign aware one, so an underflow gives a negative result
    integer a = 5; // much smaller than b * c
    sub_product(a, b, c);
    REQUIRE(a == integer(5) - b * c);

    integer a2 = 5;
    sub_product(a2, b, static_cast<uint64_t>(7));
    REQUIRE(a2 == integer(5) - b * 7u);

    // valid uses keep working
    integer d = b * c;
    d += 11u;
    sub_product(d, b, c);
    REQUIRE(d == 11u);

    integer e = b;
    e *= 7u;
    e += 3u;
    sub_product(e, b, static_cast<uint64_t>(7));
    REQUIRE(e == 3u);
}

TEST_CASE("multiplication propagates a long carry") {
    // isqrt() produced this value; squaring it needs a carry that travels through
    // two all-ones words of the partial product
    integer s;
    s.words.reset(5);
    const uint64_t w[5] = {8296379691479107416ull, 8433445100458127462ull, 6328652237515477287ull,
                           14740434617336491689ull, 24879108095803ull};
    for (int i = 0; i < 5; i++)
        s.words[i] = w[i];

    integer expected;
    add_product(expected, s, s); // independent implementation
    REQUIRE(s * s == expected);
    REQUIRE(mul_karatsuba(s, s) == expected);
    integer sq = s;
    square(sq);
    REQUIRE(sq == expected);

    // the same for a * b with b != a
    integer t = s;
    t += 1u;
    integer expected2;
    add_product(expected2, s, t);
    REQUIRE(s * t == expected2);
    REQUIRE(mul_karatsuba(s, t) == expected2);

    // cross check the general multiplication against add_product for values that are
    // likely to produce long carry chains
    Random rng(51);
    for (int i = 0; i < 300; i++) {
        integer a = 1;
        a <<= rng.Uniform<int>(64, 400);
        a -= rng.Uniform<uint64_t>(1, 1000);
        integer b = 1;
        b <<= rng.Uniform<int>(64, 400);
        b -= rng.Uniform<uint64_t>(1, 1000);
        integer e;
        add_product(e, a, b);
        REQUIRE(a * b == e);
    }
}



TEST_CASE("square matches multiplication") {
    for (uint64_t v : {uint64_t(0), uint64_t(1), uint64_t(2), uint64_t(3), UINT64_MAX, UINT64_MAX - 1}) {
        integer a = v;
        integer x = a;
        square(x);
        REQUIRE(x == a * a);
    }

    Random rng(41);
    for (int size = 1; size <= 40; size++) {
        for (int rep = 0; rep < 3; rep++) {
            const integer a = rand_natural(size, size, rng);
            integer x = a;
            square(x);
            REQUIRE(x == a * a);
        }
    }

    // sparse operands (interior zero words)
    for (int shift : {64, 128, 1000}) {
        integer a = 1;
        a <<= shift;
        a += 3u;
        integer x = a;
        square(x);
        REQUIRE(x == a * a);
    }

    // in place multiplication of a value with itself goes through square()
    integer b = rand_natural(9, 9, rng);
    integer c = b;
    mul(c, c);
    REQUIRE(c == b * b);
}

// operator~ is the arithmetic complement -(a+1). Flipping the words of a magnitude is invert_bits().
TEST_CASE("operator~ is the arithmetic complement") {
    integer a = UINT64_MAX;
    REQUIRE(~a == -(a + 1u));
    REQUIRE(~integer(0) == -1);
    REQUIRE(~~a == a);

    integer d = 0xF0F0u;
    REQUIRE(~d == -integer(0xF0F1u));
}

TEST_CASE("resize clears the inline word") {
    integer a = 12345;
    a.words.downsize(0); // the word itself is left behind
    a.words.resize(1);
    REQUIRE(a.words[0] == 0);
    REQUIRE(a == 0u);

    integer b = 999;
    b.words.downsize(0);
    b.words.resize(3); // grows onto the heap
    REQUIRE(b.words[0] == 0);
    REQUIRE(b.words[1] == 0);
    REQUIRE(b.words[2] == 0);
}

TEST_CASE("str rejects an unusable base") {
    integer a = 12345;
    REQUIRE_THROWS(a.str(0));
    REQUIRE_THROWS(a.str(1));
    REQUIRE_THROWS(a.str(37));
    REQUIRE(a.str(2) == "11000000111001");
    REQUIRE(a.str(36) == "9IX");
    REQUIRE(a.str(36, false) == "9ix");
    REQUIRE(integer(0).str(2) == "0");
}

TEST_CASE("operator+= with a builtin operand when the value is zero") {
    // __add_and_return_carry() returns the carry, which for an empty value is the operand itself
    // and not 1, so the carry has to be pushed rather than a literal 1
    integer a = 0;
    a += uint64_t(5);
    REQUIRE(a == 5u);

    integer b = 0;
    b += uint128_t(5);
    REQUIRE(b == 5u);

    // a zero that came out of a subtraction behaves the same
    integer c = 7;
    c -= 7u;
    REQUIRE(c.words.size() == 0);
    c += uint64_t(9);
    REQUIRE(c == 9u);

    // and the ordinary carry out of a full word still works
    integer d = UINT64_MAX;
    d += uint64_t(1);
    REQUIRE(d == (uint128_t(UINT64_MAX) + 1));
    integer e = UINT64_MAX;
    e += uint64_t(UINT64_MAX);
    REQUIRE(e == (uint128_t(UINT64_MAX) * 2));

    // nonzero start, no carry
    integer f = 10;
    f += uint64_t(5);
    REQUIRE(f == 15u);
    integer g = 0;
    g += uint64_t(0);
    REQUIRE(g == 0u);
}
