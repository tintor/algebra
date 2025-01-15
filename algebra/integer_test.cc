#include "algebra/integer.h"
#include <catch2/catch_test_macros.hpp>
#include <print>
#include <random>
using std::print;
using std::format;
using ucent = unsigned __int128;
using ulong = unsigned long;
using uint = unsigned int;
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

TEST_CASE("operator-") {
    integer a = 20;
    a = -a;
    REQUIRE(a == -20);

    integer b = 3;
    a = -b;
    REQUIRE(a == -3);

    REQUIRE(-(20_i) + 20_i == 0);
}

TEST_CASE("str") {
    REQUIRE(integer(0).str() == "0");
    REQUIRE(integer(1).str() == "1");
    REQUIRE(integer(-1).str() == "-1");
    REQUIRE(integer(12).str() == "12");
    REQUIRE(integer(450).str() == "450");
    REQUIRE(integer(-3692).str() == "-3692");
}

TEST_CASE("format") {
    REQUIRE(format("{}", 15_i) == "15");

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

    REQUIRE(format("{:b}", -15_i) == "-1111");
    REQUIRE(format("{:o}", -15_i) == "-17");
    REQUIRE(format("{:d}", -15_i) == "-15");
    REQUIRE(format("{:x}", -15_i) == "-f");
    REQUIRE(format("{:X}", -15_i) == "-F");
    REQUIRE(format("{:4d}", -15_i) == " -15");
    REQUIRE(format("{:1d}", -15_i) == "-15");
    REQUIRE(format("{:*>4d}", -15_i) == "*-15");
    REQUIRE(format("{:*<4d}", -15_i) == "-15*");
    REQUIRE(format("{:*^4d}", -15_i) == "-15*");
}

TEST_CASE("parse") {
    REQUIRE(integer("0") == integer(0));
    REQUIRE(integer("1") == integer(1));
    REQUIRE(integer("-1") == integer(-1));
    REQUIRE(integer("12") == integer(12));
    REQUIRE(integer("450") == integer(450));
    REQUIRE(integer("-3692") == integer(-3692));
    const char* a = "18446744073709551617"; // UINT64_MAX + 2
    REQUIRE(integer(a).str() == a);
}

TEST_CASE("static_cast<int>") {
    REQUIRE(static_cast<int>(integer(0)) == 0);
    REQUIRE(static_cast<int>(integer(1)) == 1);
    REQUIRE(static_cast<int>(integer(-1)) == -1);
    int a = std::numeric_limits<int>::max();
    int b = std::numeric_limits<int>::min();
    REQUIRE(static_cast<int>(integer(a)) == a);
    integer e(b);
    REQUIRE(e.abs.words[0] == 2147483648);
    REQUIRE(static_cast<int>(integer(b)) == b);
}

TEST_CASE("static_cast<uint>") {
    REQUIRE(static_cast<uint>(integer(0)) == 0);
    REQUIRE(static_cast<uint>(integer(1)) == 1);
    uint a = std::numeric_limits<uint>::max();
    REQUIRE(static_cast<uint>(integer(a)) == a);
}

TEST_CASE("static_cast<long>") {
    REQUIRE(static_cast<long>(integer(0)) == 0);
    REQUIRE(static_cast<long>(integer(1)) == 1);
    REQUIRE(static_cast<long>(integer(-1)) == -1);
    long a = std::numeric_limits<long>::max();
    long b = std::numeric_limits<long>::min();
    REQUIRE(static_cast<long>(integer(a)) == a);
    REQUIRE(static_cast<long>(integer(b)) == b);
}

TEST_CASE("static_cast<ulong>") {
    REQUIRE(static_cast<ulong>(integer(0)) == 0);
    REQUIRE(static_cast<ulong>(integer(1)) == 1);
    ulong a = std::numeric_limits<ulong>::max();
    REQUIRE(static_cast<ulong>(integer(a)) == a);
}

TEST_CASE("ucent") {
    ucent a = 1;
    for (int i = 0; i < 128; i++) {
        integer b(a);
        REQUIRE(b == a);
        REQUIRE(a == b);
        a <<= 1;
    }
}

TEST_CASE("sign") {
    REQUIRE(integer(0).sign() == 0);
    REQUIRE(integer(100).sign() > 0);
    REQUIRE(integer(-2).sign() < 0);

    REQUIRE(!integer(0).is_negative());
    REQUIRE(!integer(5).is_negative());
    REQUIRE(integer(-5).is_negative());

    REQUIRE(integer(0).is_even());
    REQUIRE(!integer(0).is_odd());
    REQUIRE(!integer(7).is_even());
    REQUIRE(integer(7).is_odd());
}

TEST_CASE("cmp") {
    REQUIRE(integer(0) == integer(0));
    REQUIRE(integer(5) == integer(5));
    REQUIRE(integer(-5) == integer(-5));

    REQUIRE(integer(-1) < integer(0));
    REQUIRE(integer(0) < integer(1));
    REQUIRE(integer(5) < integer(6));
    REQUIRE(integer(-5) < integer(6));
    REQUIRE(integer(-6) < integer(5));
    REQUIRE(integer(-6) < integer(-5));

    REQUIRE(integer(6) > integer(5));

    REQUIRE(integer(1) <= integer(1));
}

TEST_CASE("add") {
    REQUIRE(integer(0) + integer(0) == integer(0));
    REQUIRE(integer(5) + integer(0) == integer(5));
    REQUIRE(integer(5) + integer(6) == integer(11));
    REQUIRE(integer(5) + integer(-6) == integer(-1));
    REQUIRE(integer(-5) + integer(6) == integer(1));
    REQUIRE(integer(-5) + integer(-6) == integer(-11));
}

TEST_CASE("sub") {
    REQUIRE(integer(0) - integer(0) == integer(0));
    REQUIRE(integer(5) - integer(0) == integer(5));
    REQUIRE(integer(0) - integer(5) == integer(-5));
    REQUIRE(integer(5) - integer(5) == integer(0));

    REQUIRE(integer(5) - integer(6) == integer(-1));
    REQUIRE(integer(5) - integer(-6) == integer(11));
    REQUIRE(integer(-5) - integer(6) == integer(-11));
    REQUIRE(integer(-5) - integer(-6) == integer(1));
}

TEST_CASE("+=") {
    integer a(5);
    a += integer(4);
    REQUIRE(a == integer(9));

    integer b(-5);
    b += integer(4);
    REQUIRE(b == integer(-1));

    integer c(4);
    c += integer(-1);
    REQUIRE(c == integer(3));

    integer d(-4);
    d += integer(-3);
    REQUIRE(d == integer(-7));

    integer e(-4);
    e += integer(4);
    REQUIRE(e == integer(0));

    integer f(4);
    f += integer(-4);
    REQUIRE(f == integer(0));
}

TEST_CASE("-=") {
    integer a(5);
    a -= integer(4);
    REQUIRE(a == integer(1));

    integer b(-5);
    b -= integer(4);
    REQUIRE(b == integer(-9));

    integer c(4);
    c -= integer(-1);
    REQUIRE(c == integer(5));

    integer d(-4);
    d -= integer(-3);
    REQUIRE(d == integer(-1));

    integer e(-4);
    e -= integer(4);
    REQUIRE(e == integer(-8));

    integer f(4);
    f -= integer(-4);
    REQUIRE(f == integer(8));
}

TEST_CASE("big add 1") {
    ulong a = std::numeric_limits<ulong>::max();
    integer b = a;
    REQUIRE(b == a);

    b += a;
    REQUIRE(b > a);
    REQUIRE(b.sign() == 2);
    REQUIRE(b.abs.words[1] == 1);
    REQUIRE(b.abs.words[0] == a - 1);

    b += b;
    REQUIRE(b.sign() == 2);
    REQUIRE(b.abs.words[1] == 3);
    REQUIRE(b.abs.words[0] == a - 3);
}

TEST_CASE("big add 2") {
    const ulong m = std::numeric_limits<ulong>::max();

    integer a;
    a.abs.words.reset(4);
    a.abs.words[0] = m;
    a.abs.words[1] = m;
    a.abs.words[2] = m;
    a.abs.words[3] = 1;

    integer b;
    b.abs.words.reset(4);
    b.abs.words[0] = 1;
    b.abs.words[1] = 0;
    b.abs.words[2] = 0;
    b.abs.words[3] = 1;

    integer c = a + b;
    REQUIRE(c.sign() == 4);
    REQUIRE(c.abs.words[0] == 0);
    REQUIRE(c.abs.words[1] == 0);
    REQUIRE(c.abs.words[2] == 0);
    REQUIRE(c.abs.words[3] == 3);
}

TEST_CASE("add stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        ucent b = rng.Uniform<ucent>(0, m - a);
        REQUIRE(integer(a) + integer(b) == a + b);
    }
}

TEST_CASE("sub stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        ucent b = rng.Uniform<ucent>(0, a);
        REQUIRE(integer(a) - integer(b) == a - b);
    }
}

TEST_CASE("mul stress with ucent") {
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(1, m);
        ucent b = rng.Uniform<ucent>(0, m / a);
        REQUIRE(integer(a) * integer(b) == a * b);
    }
}

TEST_CASE("div stress with ucent") {
    Random rng;
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
    Random rng;
    const ucent m = std::numeric_limits<ucent>::max();
    for (int i = 0; i < 100000; i++) {
        ucent a = rng.Uniform<ucent>(0, m);
        integer b = a;
        if (a > std::numeric_limits<ulong>::max()) {
            REQUIRE(b.sign() == 2);
            REQUIRE(b.abs.words[0] == ulong(a));
            REQUIRE(b.abs.words[1] == ulong(a >> 64));
        }
        REQUIRE(format("{}", a) == integer(a).str());
    }
}

TEST_CASE("div10 stress with ucent 2") {
    Random rng;
    const ucent m = std::numeric_limits<ulong>::max();
    for (int i = 0; i < 1'000'000; i++) {
        ucent a = rng.Uniform<ucent>(m + 1, m * 10);
        integer q;
        long r = div(integer(a), static_cast<long>(10), q);
        REQUIRE(q * 10 + r == a);
    }
}

TEST_CASE("uniform_int") {
    integer a;
    a.abs = pow(2_n, 128) - 1;
    Random rng;
    for (int i = 0; i < 20; i++) {
        integer m = uniform_int(0, a, rng.get());
        REQUIRE(m.sign() <= a.sign());
        REQUIRE(0 <= m);
        REQUIRE(m <= a);
        REQUIRE(0 <= m.sign());
        REQUIRE(m.sign() <= 2);
        while (m.sign())
            div(m, static_cast<long>(10), /*out*/m);
    }
}

TEST_CASE("mul") {
    REQUIRE(integer(0) * integer(0) == integer(0));
    REQUIRE(integer(5) * integer(0) == integer(0));
    REQUIRE(integer(0) * integer(2) == integer(0));

    REQUIRE(integer(5) * integer(2) == integer(10));
    REQUIRE(integer(-5) * integer(2) == integer(-10));
    REQUIRE(integer(-5) * integer(-2) == integer(10));
    REQUIRE(integer(5) * integer(-2) == integer(-10));
}

TEST_CASE("*=") {
    integer a(5);
    a *= integer(3);
    REQUIRE(a == integer(15));
    a *= 3;
    REQUIRE(a == 45);
    int b = a;
    REQUIRE(b == 45);

    integer c(3);
    c *= c;
    REQUIRE(c == 9);

    integer d(3), e("1000000000000000000000000000000000000");
    d *= e;
    REQUIRE(d.abs.words.size() == 2);
    REQUIRE(d == integer("3000000000000000000000000000000000000"));
}

TEST_CASE("/") {
    REQUIRE(integer(0) / integer(7) == integer(0));
    REQUIRE(integer(7) / integer(7) == integer(1));
    REQUIRE(integer(7) / integer(3) == integer(2));
    REQUIRE(integer(7) / integer(8) == integer(0));

    REQUIRE(integer(-7) / integer(3) == integer(-2));
    REQUIRE(integer(-7) / integer(-3) == integer(2));
    REQUIRE(integer(7) / integer(-3) == integer(-2));
}

TEST_CASE("%") {
    REQUIRE(integer(0) % integer(7) == integer(0));
    REQUIRE(integer(7) % integer(7) == integer(0));
    REQUIRE(integer(7) % integer(3) == integer(1));
    REQUIRE(integer(7) % integer(8) == integer(7));

    REQUIRE(integer(-7) % integer(3) == integer(-1));
    REQUIRE(integer(-7) % integer(-3) == integer(1));
    REQUIRE(integer(7) % integer(-3) == integer(-1));
}

TEST_CASE("mod") {
    REQUIRE(mod(0_i, 7u) == 0);
    REQUIRE(mod(1_i, 7u) == 1);
    REQUIRE(mod(7_i, 7u) == 0);
    REQUIRE(mod(8_i, 7u) == 1);
#if 0
    REQUIRE(integer(14).mod(7u) == 0);
    REQUIRE(integer(-1).mod(7u) == 6);
    REQUIRE(integer(-6).mod(7u) == 1);
    REQUIRE(integer(-7).mod(7u) == 0);

    REQUIRE(integer(0).mod((ulong)7) == 0);
    REQUIRE(integer(1).mod((ulong)7) == 1);
    REQUIRE(integer(7).mod((ulong)7) == 0);
    REQUIRE(integer(8).mod((ulong)7) == 1);
    REQUIRE(integer(14).mod((ulong)7) == 0);
    REQUIRE(integer(-1).mod((ulong)7) == 6);
    REQUIRE(integer(-6).mod((ulong)7) == 1);
    REQUIRE(integer(-7).mod((ulong)7) == 0);

    REQUIRE(integer(5).mod(3u) == 2);
    REQUIRE(integer(5).mod((ulong)3) == 2);
#endif
}

// TODO randomized long division test against cpp_int for big integers!

TEST_CASE("is_x") {
    integer z = 0;
    REQUIRE(z.is_int());
    REQUIRE(z.is_uint());
    REQUIRE(z.is_long());
    REQUIRE(z.is_ulong());

    integer o = 1;
    REQUIRE(o.is_int());
    REQUIRE(o.is_uint());
    REQUIRE(o.is_long());
    REQUIRE(o.is_ulong());

    integer n = -1;
    REQUIRE(n.is_int());
    REQUIRE(!n.is_uint());
    REQUIRE(n.is_long());
    REQUIRE(!n.is_ulong());

    integer p = std::numeric_limits<uint>::max();
    REQUIRE(!p.is_int());
    REQUIRE(p.is_uint());
    REQUIRE(p.is_long());
    REQUIRE(p.is_ulong());

    integer q = (ulong)std::numeric_limits<uint>::max() + 1;
    REQUIRE(!q.is_int());
    REQUIRE(!q.is_uint());
    REQUIRE(q.is_long());
    REQUIRE(q.is_ulong());

    integer a = std::numeric_limits<long>::max();
    REQUIRE(!a.is_int());
    REQUIRE(!a.is_uint());
    REQUIRE(a.is_long());
    REQUIRE(a.is_ulong());

    integer b = (ulong)std::numeric_limits<long>::max() + 1;
    REQUIRE(!b.is_int());
    REQUIRE(!b.is_uint());
    REQUIRE(!b.is_long());
    REQUIRE(b.is_ulong());

    integer c = std::numeric_limits<ulong>::max();
    REQUIRE(!c.is_int());
    REQUIRE(!c.is_uint());
    REQUIRE(!c.is_long());
    REQUIRE(c.is_ulong());

    integer d = c + 1;
    REQUIRE(!d.is_int());
    REQUIRE(!d.is_uint());
    REQUIRE(!d.is_long());
    REQUIRE(!d.is_ulong());
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

TEST_CASE("abs") {
    REQUIRE(abs(integer(0)) == integer(0));
    REQUIRE(abs(integer(10)) == integer(10));
    REQUIRE(abs(integer(-3)) == integer(3));
}

TEST_CASE("pow") {
    REQUIRE(pow(integer(2), 3) == 8);
    REQUIRE(pow(integer(10), 30) == integer("1000000000000000000000000000000"));
}

TEST_CASE("factorial") {
    integer a(1);
    for (int i = 2; i <= 50; i++) {
        integer b = a;
        for (int j = 1; j < i; j++)
            b += a;
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
    REQUIRE(integer(-1).num_bits() == 1);
    REQUIRE(integer(2).num_bits() == 2);
    REQUIRE(integer(3).num_bits() == 2);
    REQUIRE(integer(4).num_bits() == 3);
    REQUIRE(integer(15).num_bits() == 4);
    REQUIRE(integer(16).num_bits() == 5);
    REQUIRE(integer(-16).num_bits() == 5);
}

#include <bit>

TEST_CASE("popcount") {
    for (uint i: {0, -1, 4, 31231, -3123121})
        REQUIRE(integer(i).popcount() == std::popcount(i));
}

TEST_CASE("<<") {
    integer(1) << 64; // regression test
}

TEST_CASE("literal") {
    integer a = 1'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890'234'567'890_i;
    REQUIRE(a.str() == "1234567890234567890234567890234567890234567890234567890234567890234567890234567890");
    integer b = -5_i;
    REQUIRE(b == -5);
}

constexpr std::string os() {
    std::ostringstream s;
    s << -15_i;
    return s.str();
}

TEST_CASE("constexpr ostream") {
    REQUIRE(os() == "-15");
}
