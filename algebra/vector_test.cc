#include "algebra/vector.h"
#include "algebra/rational.h"
#include "algebra/__test.h"

using V2 = Vec2<rational>;
using V3 = Vec3<rational>;

TEST_CASE("abs of a vector") {
    // an exact element type
    REQUIRE(abs(V2(rational(-3, 2), rational(4))) == V2(rational(3, 2), rational(4)));
    REQUIRE(abs(V3(rational(-1), rational(0), rational(-2))) == V3(rational(1), rational(0), rational(2)));

    // a floating point element type: unqualified abs() must not pick up ::abs(int) and truncate
    const Vec2<double> d = abs(Vec2<double>(-2.5, 1.5));
    REQUIRE(d.x == 2.5);
    REQUIRE(d.y == 1.5);
    const Vec2<float> f = abs(Vec2<float>(-0.25f, -4.5f));
    REQUIRE(f.x == 0.25f);
    REQUIRE(f.y == 4.5f);
}

TEST_CASE("cross 2d") {
    REQUIRE(cross(V2{1, 0}, V2{0, 1}) == 1);
    REQUIRE(cross(V2{0, 1}, V2{1, 0}) == -1);
    REQUIRE(cross(V2{2, 3}, V2{4, 5}) == 2 * 5 - 3 * 4);
    REQUIRE(cross(V2{2, 3}, V2{2, 3}) == 0);
}

TEST_CASE("cross 3d") {
    const V3 x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1};
    REQUIRE(cross(x, y) == z);
    REQUIRE(cross(y, z) == x);
    REQUIRE(cross(z, x) == y);
    REQUIRE(cross(y, x) == -z);
    REQUIRE(cross(z, y) == -x);
    REQUIRE(cross(x, z) == -y);

    const V3 a{1, 2, 3}, b{4, 5, 6};
    REQUIRE(cross(a, b) == V3{-3, 6, -3});
    REQUIRE(cross(b, a) == V3{3, -6, 3});
    REQUIRE(is_zero(cross(a, a)));

    // the result is orthogonal to both operands
    Random rng(13);
    for (int i = 0; i < 30; i++) {
        auto r = [&] { return rational(rng.Uniform<int>(-9, 9), rng.Uniform<int>(1, 5)); };
        const V3 p{r(), r(), r()}, q{r(), r(), r()};
        const V3 c = cross(p, q);
        REQUIRE(dot(c, p) == 0);
        REQUIRE(dot(c, q) == 0);
    }
}

TEST_CASE("dot and lerp") {
    REQUIRE(dot(V3{1, 2, 3}, V3{4, 5, 6}) == 32);
    REQUIRE(dot2(V2{3, 4}) == 25);
    REQUIRE(lerp(V2{0, 0}, V2{2, 4}, rational(1, 2)) == V2{1, 2});
    REQUIRE(lerp(V2{1, 1}, V2{3, 5}, rational(0)) == V2{1, 1});
    REQUIRE(lerp(V2{1, 1}, V2{3, 5}, rational(1)) == V2{3, 5});
}

TEST_CASE("indexing") {
    Vec2<rational> a{1, 2};
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    a[0] = 7;
    REQUIRE(a.x == 7);

    Vec3<rational> b{1, 2, 3};
    REQUIRE(b[0] == 1);
    REQUIRE(b[1] == 2);
    REQUIRE(b[2] == 3);
    b[2] = 9;
    REQUIRE(b.z == 9);

    Vec4<rational> c{1, 2, 3, 4};
    REQUIRE(c[0] == 1);
    REQUIRE(c[1] == 2);
    REQUIRE(c[2] == 3);
    REQUIRE(c[3] == 4);
    c[3] = 8;
    REQUIRE(c.w == 8);

    REQUIRE(Vec3<rational>::dim == 3);

    // usable in a constant expression
    static_assert([]{ Vec2<int> v{3, 4}; return v[0] * 10 + v[1]; }() == 34);
    static_assert([]{ Vec4<int> v{1, 2, 3, 4}; v[2] = 9; return v.z; }() == 9);
}

TEST_CASE("is_zero") {
    REQUIRE(is_zero(V3(0, 0, 0)));
    REQUIRE(!is_zero(V3(0, 0, 1)));
    REQUIRE(!is_zero(V3(-1, 0, 0)));
    REQUIRE(is_zero(V2(0, 0)));
    REQUIRE(!is_zero(V2(1, 0)));
}

TEST_CASE("argmax_abs") {
    REQUIRE(argmax_abs(V3(1, 2, 3)) == 2);
    REQUIRE(argmax_abs(V3(-4, 2, 3)) == 0);
    REQUIRE(argmax_abs(V3(1, -5, 3)) == 1);
    REQUIRE(argmax_abs(V3(0, 0, 0)) == 0);
    REQUIRE(argmax_abs(V3(2, 2, 2)) == 0); // ties keep the lowest index
    REQUIRE(argmax_abs(V3(2, -2, 1)) == 0);
    REQUIRE(argmax_abs(V2(1, -7)) == 1);
}

TEST_CASE("div_colinear") {
    // k such that b * k == a
    REQUIRE(div_colinear(V3(2, 4, 6), V3(1, 2, 3)) == 2);
    REQUIRE(div_colinear(V3(1, 2, 3), V3(2, 4, 6)) == 1/2_q);
    REQUIRE(div_colinear(V3(-1, -2, -3), V3(1, 2, 3)) == -1);
    REQUIRE(div_colinear(V3(0, 0, 0), V3(1, 2, 3)) == 0);
    // the largest component of b is used, so a zero component of b is not divided by
    REQUIRE(div_colinear(V3(0, 6, 0), V3(0, 3, 0)) == 2);
}

TEST_CASE("order / strict_order / loose_order - scalars") {
    // order is strictly monotonic, or all three equal
    REQUIRE(order<rational>(1, 2, 3));
    REQUIRE(order<rational>(3, 2, 1));
    REQUIRE(order<rational>(2, 2, 2));
    REQUIRE(!order<rational>(1, 1, 2));
    REQUIRE(!order<rational>(1, 2, 2));
    REQUIRE(!order<rational>(1, 3, 2));

    REQUIRE(strict_order<rational>(1, 2, 3));
    REQUIRE(strict_order<rational>(3, 2, 1));
    REQUIRE(!strict_order<rational>(2, 2, 2));
    REQUIRE(!strict_order<rational>(1, 1, 2));

    REQUIRE(loose_order<rational>(1, 2, 3));
    REQUIRE(loose_order<rational>(3, 2, 1));
    REQUIRE(loose_order<rational>(2, 2, 2));
    REQUIRE(loose_order<rational>(1, 1, 2));
    REQUIRE(loose_order<rational>(1, 2, 2));
    REQUIRE(!loose_order<rational>(1, 3, 2));
}

TEST_CASE("order / strict_order / loose_order - Vec2") {
    // both components must satisfy it
    REQUIRE(order(V2(0, 0), V2(1, 0), V2(2, 0)));
    REQUIRE(!order(V2(0, 0), V2(1, 1), V2(2, 1)));
    REQUIRE(loose_order(V2(0, 0), V2(1, 1), V2(2, 1)));
    REQUIRE(!strict_order(V2(0, 0), V2(1, 0), V2(2, 0))); // y is not strictly increasing
    REQUIRE(strict_order(V2(0, 0), V2(1, 1), V2(2, 2)));
    REQUIRE(!loose_order(V2(0, 0), V2(3, 1), V2(2, 2))); // x turns back
}

TEST_CASE("same_sign") {
    REQUIRE(same_sign<rational>(2, 5));
    REQUIRE(same_sign<rational>(-2, -5));
    REQUIRE(same_sign<rational>(0, 0));
    REQUIRE(!same_sign<rational>(0, 1));
    REQUIRE(!same_sign<rational>(1, 0));
    REQUIRE(!same_sign<rational>(-1, 1));

    REQUIRE(same_sign(V3(1, -2, 0), V3(7, -1, 0)));
    REQUIRE(!same_sign(V3(1, -2, 0), V3(7, -1, 1)));
}

TEST_CASE("min and minimize") {
    REQUIRE(min<rational>(2, 5) == 2);
    REQUIRE(min<rational>(5, 2) == 2);
    rational a = 5;
    minimize(a, rational(2));
    REQUIRE(a == 2);
    minimize(a, rational(7)); // no change
    REQUIRE(a == 2);
}

TEST_CASE("format Vec") {
    // components are separated by a single space, each formatted with its own formatter
    REQUIRE(std::format("{}", V2(1, 2)) == "1 2");
    REQUIRE(std::format("{}", V3(1, 2, 3)) == "1 2 3");
    REQUIRE(std::format("{}", Vec4<rational>(1, 2, 3, 4)) == "1 2 3 4");
    REQUIRE(std::format("{}", V2(rational(1, 2), rational(-3, 4))) == "1/2 -3/4");
    REQUIRE(std::format("{}", Vec2<int>(3, 4)) == "3 4");
    REQUIRE(std::format("{}", V2(2, 1)) == "2 1");

    // and the rest of the format string is not disturbed
    REQUIRE(std::format("[{}]", V2(1, 2)) == "[1 2]");
    REQUIRE(std::format("{}|{}", V2(1, 2), V2(3, 4)) == "1 2|3 4");
    REQUIRE(std::format("{}, {}", V3(1, 2, 3), 5) == "1 2 3, 5");
}

TEST_CASE("format Vec into other output iterators") {
    std::string s;
    std::format_to(std::back_inserter(s), "{}", V2(1, 2));
    REQUIRE(s == "1 2");

    char buf[32] = {};
    auto r = std::format_to_n(buf, sizeof(buf), "{}", V2(1, 2));
    REQUIRE(r.size == 3);
    REQUIRE(std::string(buf, r.out) == "1 2");

    // writing into a plain pointer: the formatter has to write through the iterator it is
    // given and return the position past the last character it wrote
    char raw[16] = {};
    auto end = std::format_to(raw, "{}", V3(1, 2, 3));
    REQUIRE(std::string(raw, end) == "1 2 3");

    char raw2[16] = {};
    auto end2 = std::format_to(raw2, "{};", V2(rational(1, 2), 3));
    REQUIRE(std::string(raw2, end2) == "1/2 3;");
}

TEST_CASE("format Vec has no format spec") {
    // std::format checks the specifier at compile time, so this needs the runtime interface
    Vec2<int> v{1, 2};
    REQUIRE(std::vformat("{}", std::make_format_args(v)) == "1 2");
    REQUIRE(std::vformat("{:}", std::make_format_args(v)) == "1 2");
    REQUIRE_THROWS_AS(std::vformat("{:d}", std::make_format_args(v)), std::format_error);
}

TEST_CASE("formatting a vector into a counted buffer") {
    // the Vec formatter writes each component through ctx.out() and drops what format_to() returns
    const V2 v(rational(1, 2), rational(3));
    char buf[64] = {};
    const auto r = std::format_to_n(buf, sizeof(buf) - 1, "{}", v);
    REQUIRE(std::string(buf, r.out) == std::format("{}", v));
    REQUIRE(std::string(buf, r.out) == "1/2 3");

    const V3 w(rational(-1), rational(0), rational(7, 2));
    char buf3[64] = {};
    const auto r3 = std::format_to_n(buf3, sizeof(buf3) - 1, "{}", w);
    REQUIRE(std::string(buf3, r3.out) == std::format("{}", w));
}
