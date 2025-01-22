#include "algebra/segment_segment_squared_distance.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"
using namespace algebra;
using std::array;

TEST_CASE("segment_segment_squared_distance - line vs line") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {3,3,8}, {3,-3,8}) == 64);
}

TEST_CASE("segment_segment_squared_distance - parallel") {
    Vec3<rational> a(0,0,0), b(5,0,0), d(1,2,0);
    REQUIRE(segment_segment_squared_distance(a, b, a + d, b + d) == 4);
}

TEST_CASE("segment_segment_squared_distance - colinear") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {7,0,0}, {10,0,0}) == 4);
}

TEST_CASE("segment_segment_squared_distance - identical") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {0,0,0}, {5,0,0}) == 0);
}

TEST_CASE("segment_segment_squared_distance - overlap") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {5,0,0}, {4,0,0}) == 0);
}

TEST_CASE("segment_segment_squared_distance - vertex vs vertex") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {6,1,0}, {6,9,0}) == 2);
}

TEST_CASE("segment_segment_squared_distance - vertex vs line") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {1,2,0}, {1,9,0}) == 4);
}

TEST_CASE("failing") {
    Vec3<rational> a={-97, -68, 74};
    Vec3<rational> b={54, 73, -33};
    Vec3<rational> c={-6, -59, -74};
    Vec3<rational> d={-77, 61, -35};
    const auto e = segment_segment_squared_distance(a, b, c, d);
    REQUIRE(segment_segment_squared_distance(c, d, b, a) == e);
}

using Vec3Q = Vec3<rational>;

void REQUIRE_segment_segment_squared_distance(Vec3Q a, Vec3Q b, Vec3Q c, Vec3Q d, rational e, const array<Vec3Q, 4>& orig) {
    auto o = segment_segment_squared_distance(a, b, c, d);
    if (o != e) {
        std::print("segment_segment_squared_distance({} | {} | {} | {}) = {}\n", a, b, c, d, o);
        std::print("segment_segment_squared_distance({} | {} | {} | {}) = {}\n", orig[0], orig[1], orig[2], orig[3], e);
    }
    REQUIRE(o == e);
}

void test_segment_segment_squared_distance(Vec3Q a, Vec3Q b, Vec3Q c, Vec3Q d, rational e, const array<Vec3Q, 4>& orig) {
    REQUIRE_segment_segment_squared_distance(a, b, c, d, e, orig);
    REQUIRE_segment_segment_squared_distance(b, a, c, d, e, orig);
    REQUIRE_segment_segment_squared_distance(a, b, d, c, e, orig);
    REQUIRE_segment_segment_squared_distance(b, a, d, c, e, orig);

    REQUIRE_segment_segment_squared_distance(c, d, a, b, e, orig);
    REQUIRE_segment_segment_squared_distance(c, d, b, a, e, orig);
    REQUIRE_segment_segment_squared_distance(d, c, a, b, e, orig);
    REQUIRE_segment_segment_squared_distance(d, c, b, a, e, orig);
}

void test_segment_segment_squared_distance2(Vec3<rational> a, Vec3<rational> b, Vec3<rational> c, Vec3<rational> d, Random& rng) {
    const auto e = segment_segment_squared_distance(a, b, c, d);
    array<Vec3Q, 4> orig = {a, b, c, d};

    for (int i = 0; i < 6; i++) {
        auto conv = [i](Vec3<rational> v) {
            if (i == 0) return v;
            if (i == 1) return xzy(v);
            if (i == 2) return yxz(v);
            if (i == 3) return yzx(v);
            if (i == 4) return zxy(v);
            if (i == 5) return zyx(v);
            return v;
        };
        test_segment_segment_squared_distance(conv(a), conv(b), conv(c), conv(d), e, orig);
    }

    const auto s = rng.SampleQ();
    if (s != 0)
        test_segment_segment_squared_distance(a * s, b * s, c * s, d * s, e * s * s, orig);

    const auto t = rng.SampleVec3Q();
    test_segment_segment_squared_distance(a + t, b + t, c + t, d + t, e, orig);
}

TEST_CASE("segment_segment_squared_distance - stress") {
    Random rng(0);
    for (int i = 0; i < 1000; i++) {
        const auto a = rng.SampleVec3Q();
        const auto b = rng.SampleVec3Q();
        const auto c = rng.SampleVec3Q();
        const auto d = rng.SampleVec3Q();

        test_segment_segment_squared_distance2(a, b, c, d, rng);

        // touching
        test_segment_segment_squared_distance2(a, b, a, b, rng);
        test_segment_segment_squared_distance2(a, b, a, c, rng);
        const auto ab4 = lerp(a, b, 0.4_q);
        const auto ab6 = lerp(a, b, 0.6_q);
        test_segment_segment_squared_distance2(a, b, ab4, c, rng);

        // degenerate
        test_segment_segment_squared_distance2(a, a, b, c, rng);
        test_segment_segment_squared_distance2(a, a, b, b, rng);
        test_segment_segment_squared_distance2(a, a, a, c, rng);

        // colinear
        test_segment_segment_squared_distance2(a, b, b, ab4, rng);
        test_segment_segment_squared_distance2(a, ab4, ab6, b, rng);
        test_segment_segment_squared_distance2(a, ab6, ab4, b, rng);
        test_segment_segment_squared_distance2(a, b, ab4, ab6, rng);
    }
}
