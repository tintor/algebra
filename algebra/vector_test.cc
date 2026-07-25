#include "algebra/vector.h"
#include "algebra/rational.h"
#include "algebra/__test.h"

using V2 = Vec2<rational>;
using V3 = Vec3<rational>;

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
