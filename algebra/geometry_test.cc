#include "algebra/geometry.h"
#include "algebra/rational_func.h"
#include "algebra/__test.h"
using namespace algebra;

VEC_OP_vs_sv(+, rational, rational, std::integral auto)
VEC_OP_vs_sv(-, rational, rational, std::integral auto)
VEC_OP_vs_sv(*, rational, rational, std::integral auto)
VEC_OP_vs_sv(/, rational, rational, std::integral auto)

TEST_CASE("ccw - colinear") {
    Vec2<rational> a(4, 1), b(5, 2), c(7, 4);
    REQUIRE(ccw(a, b, c) == 0);
}

TEST_CASE("ccw - non-colinear") {
    Vec2<rational> a(4, 1), b(5, 3), c(7, 4);
    auto e = ccw(a, b, c);
    REQUIRE(e != 0);
    REQUIRE(ccw(b, c, a) == e);
    REQUIRE(ccw(c, a, b) == e);
    REQUIRE(ccw(a, c, b) == -e);
    REQUIRE(ccw(b, a, c) == -e);
    REQUIRE(ccw(c, b, a) == -e);
}

TEST_CASE("segment_vs_segment_intersection_single_point - cross (vertical | horisontal)") {
    rational m, n;
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,1}, {1,-1}));
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,1}, {1,-1}, &m, &n));
    REQUIRE(m == 1/3_q);
    REQUIRE(n == 1/2_q);
}

TEST_CASE("segment_vs_segment_intersection_single_point - cross (diagonal)") {
    rational m, n;
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0,0}, {9,9}, {9,0}, {0,9}, &m, &n));
    REQUIRE(m == 1/2_q);
    REQUIRE(n == 1/2_q);
}

TEST_CASE("segment_vs_segment_intersection_single_point - disjoint (cross)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,3}, {1,1}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - disjoint (collinear)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {3 + pow(10_q, 20), 0}, {5,0}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - parallel") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {0,2}, {3,2}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - touch (collinear)") {
    rational m, n;
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {3,0}, {5,0}, &m, &n));
    REQUIRE(m == 1);
    REQUIRE(n == 0);
}

TEST_CASE("segment_vs_segment_intersection_single_point - touch (T)") {
    rational m, n;
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,0}, {7,7}, &m, &n));
    REQUIRE(m == 1/3_q);
    REQUIRE(n == 0);
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (identical)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {0,0}, {3,0}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (touch)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0,0}, {3,0}, {3,0}, {2,0}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (contained)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0,0}, {10,0}, {3,0}, {5,0}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (side by side)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0,0}, {5,0}, {4,0}, {10,0}));
}

TEST_CASE("point_segment_squared_distance") {
    // degenerate
    REQUIRE(point_segment_squared_distance<rational>({5,0,0}, {5,0,0}, {5,0,0}) == 0);
    REQUIRE(point_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {5,0,0}) == 25);

    // colinear
    REQUIRE(point_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {8,0,0}) == 25);
    REQUIRE(point_segment_squared_distance<rational>({11,0,0}, {5,0,0}, {8,0,0}) == 9);

    // zero
    REQUIRE(point_segment_squared_distance<rational>({5,0,0}, {5,0,0}, {8,0,0}) == 0);
    REQUIRE(point_segment_squared_distance<rational>({6,0,0}, {5,0,0}, {8,0,0}) == 0);
    REQUIRE(point_segment_squared_distance<rational>({8,0,0}, {5,0,0}, {8,0,0}) == 0);

    // normal
    REQUIRE(point_segment_squared_distance<rational>({6,1,0}, {5,0,0}, {8,0,0}) == 1);
    REQUIRE(point_segment_squared_distance<rational>({7,-2,0}, {5,0,0}, {8,0,0}) == 4);

}

TEST_CASE("segment_segment_squared_distance - line vs line") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {3,3,8}, {3,-3,8}) == 64);
}

TEST_CASE("segment_vs_segment_squared_distance - parallel") {
    Vec3<rational> a(0,0,0), b(5,0,0), d(1,2,0);
    REQUIRE(segment_segment_squared_distance(a, b, a + d, b + d) == 4);
}

TEST_CASE("segment_vs_segment_squared_distance - colinear") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {7,0,0}, {10,0,0}) == 4);
}

TEST_CASE("segment_vs_segment_squared_distance - identical") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {0,0,0}, {5,0,0}) == 0);
}

TEST_CASE("segment_vs_segment_squared_distance - overlap") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {5,0,0}, {4,0,0}) == 0);
}

TEST_CASE("segment_vs_segment_squared_distance - vertex vs vertex") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {6,1,0}, {6,9,0}) == 2);
}

TEST_CASE("segment_vs_segment_squared_distance - vertex vs line") {
    REQUIRE(segment_segment_squared_distance<rational>({0,0,0}, {5,0,0}, {1,2,0}, {1,9,0}) == 4);
}

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
        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(_rng);
        }
        if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min, max);
            return dist(_rng);
        }
    }

    rational SampleQ() { return rational{Uniform<int>(-99, 99)}; } //, Uniform<int>(1, 1000)}; }
    Vec3<rational> SampleVec3Q() { return {SampleQ(), SampleQ(), SampleQ()}; }

private:
    std::mt19937_64 _rng;
};

TEST_CASE("failing1") {
    Vec3<rational> a={2,0,0};
    Vec3<rational> b={5,0,0};
    Vec3<rational> c={6,1,0};
    Vec3<rational> d={6,-1,0};

    const auto e = segment_segment_squared_distance(a, b, c, d);
    REQUIRE(segment_segment_squared_distance(b, a, c, d) == e);
    REQUIRE(segment_segment_squared_distance(a, b, d, c) == e);
    REQUIRE(segment_segment_squared_distance(b, a, d, c) == e);

    REQUIRE(segment_segment_squared_distance(c, d, a, b) == e);
    REQUIRE(segment_segment_squared_distance(c, d, b, a) == e);
    REQUIRE(segment_segment_squared_distance(d, c, a, b) == e);
    REQUIRE(segment_segment_squared_distance(d, c, b, a) == e);
}

TEST_CASE("failing") {
    Vec3<rational> a={-97, -68, 74};
    Vec3<rational> b={54, 73, -33};
    Vec3<rational> c={-6, -59, -74};
    Vec3<rational> d={-77, 61, -35};
    const auto e = segment_segment_squared_distance(a, b, c, d);
    REQUIRE(segment_segment_squared_distance(c, d, b, a) == e);
}

void test_segment_segment_squared_distance(Vec3<rational> a, Vec3<rational> b, Vec3<rational> c, Vec3<rational> d, rational e) {
    std::print("\na={}\n", a);
    std::print("b={}\n", b);
    std::print("c={}\n", c);
    std::print("d={}\n", d);
    std::print("e={}\n", e);
    REQUIRE(segment_segment_squared_distance(a, b, c, d) == e);
    REQUIRE(segment_segment_squared_distance(b, a, c, d) == e);
    REQUIRE(segment_segment_squared_distance(a, b, d, c) == e);
    REQUIRE(segment_segment_squared_distance(b, a, d, c) == e);

    REQUIRE(segment_segment_squared_distance(c, d, a, b) == e);
    REQUIRE(segment_segment_squared_distance(c, d, b, a) == e);
    REQUIRE(segment_segment_squared_distance(d, c, a, b) == e);
    REQUIRE(segment_segment_squared_distance(d, c, b, a) == e);

}

void test_segment_segment_squared_distance2(Vec3<rational> a, Vec3<rational> b, Vec3<rational> c, Vec3<rational> d, Random& rng) {
    const auto e = segment_segment_squared_distance(a, b, c, d);

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
        test_segment_segment_squared_distance(conv(a), conv(b), conv(c), conv(d), e);
    }

    const auto s = rng.SampleQ();
    if (s != 0)
        test_segment_segment_squared_distance(a * s, b * s, c * s, d * s, e * s * s);

    const auto t = rng.SampleVec3Q();
    test_segment_segment_squared_distance(a + t, b + t, c + t, d + t, e);
}

#if 0
TEST_CASE("segment_vs_segment_squared_distance - stress") {
    Random rng(0);
    for (int i = 0; i < 10000'000; i++) {
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
#endif
