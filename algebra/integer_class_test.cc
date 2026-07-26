#include "algebra/integer_class.h"
#include "algebra/__test.h"

TEST_CASE("add_product") {
    integer a;

    a = 1;
    add_product(a, 2_i, 3_i);
    REQUIRE(a == 7);
    add_product(a, 2_i, 0_i);
    REQUIRE(a == 7);
    add_product(a, 0_i, 3_i);
    REQUIRE(a == 7);
    add_product(a, 1_i, 3_i);
    REQUIRE(a == 10);
    add_product(a, 2_i, 1_i);
    REQUIRE(a == 12);

    a = -1;
    add_product(a, 2_i, 3_i);
    REQUIRE(a == 5);

    a = 1;
    add_product(a, -2_i, 3_i);
    REQUIRE(a == -5);

    a = -1;
    add_product(a, -2_i, 3_i);
    REQUIRE(a == -7);
}

TEST_CASE("ctor") {
    integer a = -1;
    integer b(a.abs);
    REQUIRE(b == 1);
    REQUIRE(b.sign() == 1);
    integer c(a);
    REQUIRE(c == -1);
    REQUIRE(c.sign() == -1);
}

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
    // INT_MIN is the case that used to negate a signed int that already held INT_MIN. The runtime
    // REQUIREs below cannot catch that reliably (whether the optimizer exploits the overflow
    // depends on inlining), but undefined behaviour in a constant expression is ill-formed, so
    // these static_asserts fail to compile if the negation is done in a signed type.
    static_assert(static_cast<int>(integer(std::numeric_limits<int>::min())) == std::numeric_limits<int>::min());
    REQUIRE(static_cast<int>(integer(0)) == 0);
    REQUIRE(static_cast<int>(integer(1)) == 1);
    REQUIRE(static_cast<int>(integer(-1)) == -1);
    int a = std::numeric_limits<int>::max();
    int b = std::numeric_limits<int>::min();
    REQUIRE(static_cast<int>(integer(a)) == a);
    integer e(b);
    REQUIRE(e.sign() == -1);
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
    static_assert(static_cast<long>(integer(std::numeric_limits<long>::min())) == std::numeric_limits<long>::min());
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
    REQUIRE(integer(-7) % integer(-3) == integer(-1));
    REQUIRE(integer(7) % integer(-3) == integer(1));
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
    REQUIRE(z.is_int32());
    REQUIRE(z.is_uint32());
    REQUIRE(z.is_int64());
    REQUIRE(z.is_uint64());

    integer o = 1;
    REQUIRE(o.is_int32());
    REQUIRE(o.is_uint32());
    REQUIRE(o.is_int64());
    REQUIRE(o.is_uint64());

    integer n = -1;
    REQUIRE(n.is_int32());
    REQUIRE(!n.is_uint32());
    REQUIRE(n.is_int64());
    REQUIRE(!n.is_uint64());

    integer p = std::numeric_limits<uint>::max();
    REQUIRE(!p.is_int32());
    REQUIRE(p.is_uint32());
    REQUIRE(p.is_int64());
    REQUIRE(p.is_uint64());

    integer q = (ulong)std::numeric_limits<uint>::max() + 1;
    REQUIRE(!q.is_int32());
    REQUIRE(!q.is_uint32());
    REQUIRE(q.is_int64());
    REQUIRE(q.is_uint64());

    integer a = std::numeric_limits<long>::max();
    REQUIRE(!a.is_int32());
    REQUIRE(!a.is_uint32());
    REQUIRE(a.is_int64());
    REQUIRE(a.is_uint64());

    integer b = (ulong)std::numeric_limits<long>::max() + 1;
    REQUIRE(!b.is_int32());
    REQUIRE(!b.is_uint32());
    REQUIRE(!b.is_int64());
    REQUIRE(b.is_uint64());

    integer c = std::numeric_limits<ulong>::max();
    REQUIRE(!c.is_int32());
    REQUIRE(!c.is_uint32());
    REQUIRE(!c.is_int64());
    REQUIRE(c.is_uint64());

    integer d = c + 1;
    REQUIRE(!d.is_int32());
    REQUIRE(!d.is_uint32());
    REQUIRE(!d.is_int64());
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

TEST_CASE("increment / decrement") {
    integer a = 5;
    REQUIRE(++a == 6);
    REQUIRE(a == 6);
    REQUIRE(a++ == 6);
    REQUIRE(a == 7);

    integer b = 5;
    REQUIRE(--b == 4);
    REQUIRE(b == 4);
    REQUIRE(b-- == 4);
    REQUIRE(b == 3);

    integer c = -5;
    REQUIRE(++c == -4);
    REQUIRE(--c == -5);

    integer z = 0;
    REQUIRE(++z == 1);
    REQUIRE(--z == 0);
    REQUIRE(--z == -1);
    REQUIRE(++z == 0);
    REQUIRE(z.sign() == 0);

    // carry across a word boundary
    integer w2_64 = 1;
    w2_64 <<= 64;
    integer w = w2_64;
    --w;
    ++w;
    REQUIRE(w == w2_64);
    ++w;
    --w;
    REQUIRE(w == w2_64);

    integer n = w2_64;
    n.negate();
    integer nw = n;
    ++nw;
    REQUIRE(nw == n + 1);
    --nw;
    REQUIRE(nw == n);
}

TEST_CASE("mod integer") {
    REQUIRE(mod(integer(7), integer(5)) == 2);
    REQUIRE(mod(integer(-7), integer(5)) == 3);
    REQUIRE(mod(integer(12), integer(5)) == 2);

    integer big = 1;
    big <<= 200;
    big += 12345;
    integer d = 1;
    d <<= 100;
    d += 7;
    integer q = big / d;
    REQUIRE(mod(big, d) == big - q * d);
}

TEST_CASE("static_cast<uint8_t> range") {
    auto u8 = [](const integer& a) { return static_cast<uint8_t>(a); };
    auto u16 = [](const integer& a) { return static_cast<uint16_t>(a); };

    REQUIRE(u8(integer(0)) == 0);
    REQUIRE(u8(integer(200)) == 200);
    REQUIRE(u8(integer(255)) == 255);
    REQUIRE_THROWS(u8(integer(256)));
    REQUIRE_THROWS(u8(integer(300)));
    REQUIRE_THROWS(u8(integer(70000)));
    REQUIRE_THROWS(u8(integer(-1)));

    REQUIRE(u16(integer(65535)) == 65535);
    REQUIRE_THROWS(u16(integer(65536)));
    REQUIRE_THROWS(u16(integer(-1)));
}

TEST_CASE("parse with base") {
    REQUIRE(integer("ff", 16) == 255);
    REQUIRE(integer("FF", 16) == 255);
    REQUIRE(integer("-ff", 16) == -255);
    REQUIRE(integer("101", 2) == 5);
    REQUIRE(integer("-101", 2) == -5);
    REQUIRE(integer("777", 8) == 511);
    REQUIRE(integer("123") == 123);
    REQUIRE(integer("-123") == -123);
    REQUIRE(integer(std::string_view("ff"), 16) == 255);
}

TEST_CASE("mod is euclidean") {
    // mod() returns a value in [0, abs(b))
    REQUIRE(mod(-1_i, 7u) == 6);
    REQUIRE(mod(-6_i, 7u) == 1);
    REQUIRE(mod(-7_i, 7u) == 0);
    REQUIRE(mod(-14_i, 7u) == 0);
    REQUIRE(mod(7_i, 7u) == 0);
    REQUIRE(mod(0_i, 7u) == 0);

    REQUIRE(mod(integer(-10), integer(5)) == 0);
    REQUIRE(mod(integer(-10), integer(3)) == 2);
    REQUIRE(mod(integer(10), integer(3)) == 1);
    REQUIRE(mod(integer(-10), integer(-3)) == 2);
    REQUIRE(mod(integer(10), integer(-3)) == 1);
    REQUIRE(mod(integer(-9), integer(3)) == 0);
    REQUIRE(mod(integer(0), integer(3)) == 0);

    // large operands
    integer big = 1;
    big <<= 200;
    REQUIRE(mod(big, big) == 0);
    integer nbig = -big;
    REQUIRE(mod(nbig, big) == 0);
    REQUIRE(mod(nbig - 1, big) == big - 1);
}

TEST_CASE("is_int8 / is_int16 / is_int32 / is_int64 boundaries") {
    // the negative side reaches one further than the positive side
    REQUIRE(integer(127).is_int8());
    REQUIRE(!integer(128).is_int8());
    REQUIRE(integer(-128).is_int8());
    REQUIRE(!integer(-129).is_int8());

    REQUIRE(integer(INT16_MAX).is_int16());
    REQUIRE(!integer(int32_t(INT16_MAX) + 1).is_int16());
    REQUIRE(integer(INT16_MIN).is_int16());
    REQUIRE(!integer(int32_t(INT16_MIN) - 1).is_int16());

    REQUIRE(integer(INT32_MAX).is_int32());
    REQUIRE(!integer(int64_t(INT32_MAX) + 1).is_int32());
    REQUIRE(integer(INT32_MIN).is_int32());
    REQUIRE(!integer(int64_t(INT32_MIN) - 1).is_int32());

    REQUIRE(integer(INT64_MAX).is_int64());
    REQUIRE(!(integer(INT64_MAX) + 1).is_int64());
    REQUIRE(integer(INT64_MIN).is_int64());
    REQUIRE(!(integer(INT64_MIN) - 1).is_int64());

    // zero fits everywhere, including a zero produced by subtraction
    integer z = 5;
    z -= 5;
    REQUIRE(z.is_zero());
    REQUIRE(z.is_int8());
    REQUIRE(z.is_uint8());
    REQUIRE(z.is_int128());
}

TEST_CASE("is_int128 boundaries") {
    const integer two127 = (integer(1) << 127);

    // INT128_MIN == -2**127 is representable, +2**127 is not
    REQUIRE((-two127).is_int128());
    REQUIRE(!two127.is_int128());
    REQUIRE(static_cast<int128_t>(-two127) == std::numeric_limits<int128_t>::min());

    // one below INT128_MIN is not representable, even though its high word matches
    REQUIRE(!(-(two127 + 1)).is_int128());
    REQUIRE(!(-(two127 + 12345)).is_int128());
    REQUIRE_THROWS(static_cast<int128_t>(-(two127 + 1)));

    // just inside the range on both sides
    REQUIRE((two127 - 1).is_int128());
    REQUIRE(static_cast<int128_t>(two127 - 1) == std::numeric_limits<int128_t>::max());
    REQUIRE((-(two127 - 1)).is_int128());

    // is_uint128 accepts anything up to two words, but nothing negative
    REQUIRE(two127.is_uint128());
    REQUIRE((two127 * 2 - 1).is_uint128());
    REQUIRE(!(-two127).is_uint128());
    REQUIRE(!(two127 * 2).is_uint128());
}

TEST_CASE("is_one") {
    REQUIRE(integer(1).is_one());
    REQUIRE(!integer(-1).is_one());
    REQUIRE(!integer(0).is_one());
    REQUIRE(!integer(2).is_one());
    REQUIRE(!((integer(1) << 64) + 1).is_one());
}

TEST_CASE("plus and minus a signed builtin operand") {
    // a - b with b < 0 must add abs(b), and the sign of a mixed-sign magnitude difference
    // must only be flipped for subtraction
    for (long x : {-12L, -5L, -1L, 0L, 1L, 5L, 12L})
        for (long y : {-12L, -5L, -1L, 0L, 1L, 5L, 12L}) {
            integer a = x;
            a += y;
            REQUIRE(a == x + y);

            integer b = x;
            b -= y;
            REQUIRE(b == x - y);

            REQUIRE(integer(x) + y == x + y);
            REQUIRE(integer(x) - y == x - y);
            REQUIRE(y + integer(x) == y + x);
            REQUIRE(y - integer(x) == y - x);
        }

    // the cases that used to be wrong, spelled out
    integer a = 10;
    a -= -3;
    REQUIRE(a == 13);
    integer b = -10;
    b -= -3;
    REQUIRE(b == -7);
    integer c = 0;
    c -= -5;
    REQUIRE(c == 5);
    integer d = -3;
    d += 5;
    REQUIRE(d == 2);

    // multi word, so the magnitude subtraction is not a single word operation
    const integer big = integer(1) << 80;
    integer e = big;
    e -= -5;
    REQUIRE(e == big + 5);
    integer f = -big;
    f -= -5;
    REQUIRE(f == -(big - 5));
    integer g = -big;
    g += 5;
    REQUIRE(g == -(big - 5));
}
