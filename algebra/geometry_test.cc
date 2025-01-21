#include "algebra/geometry.h"
#include "algebra/rational_func.h"
#include "algebra/__test.h"
using namespace algebra;

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
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {1_q, 1_q}, {1_q, -1_q}));
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {1_q, 1_q}, {1_q, -1_q}, &m, &n));
    REQUIRE(m == 1/3_q);
    REQUIRE(n == 1/2_q);
}

TEST_CASE("segment_vs_segment_intersection_single_point - cross (diagonal)") {
    rational m, n;
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {9_q, 9_q}, {9_q, 0_q}, {0_q, 9_q}, &m, &n));
    REQUIRE(m == 1/2_q);
    REQUIRE(n == 1/2_q);
}

TEST_CASE("segment_vs_segment_intersection_single_point - disjoint (cross)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {1_q, 3_q}, {1_q, 1_q}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - disjoint (collinear)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {3_q + pow(10_q, 20), 0_q}, {5_q, 0_q}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - parallel") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {0_q, 2_q}, {3_q, 2_q}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - touch (collinear)") {
    rational m, n;
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {3_q, 0_q}, {5_q, 0_q}, &m, &n));
    REQUIRE(m == 1_q);
    REQUIRE(n == 0_q);
}

TEST_CASE("segment_vs_segment_intersection_single_point - touch (T)") {
    rational m, n;
    REQUIRE(segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {1_q, 0_q}, {7_q, 7_q}, &m, &n));
    REQUIRE(m == 1/3_q);
    REQUIRE(n == 0_q);
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (identical)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {0_q, 0_q}, {3_q, 0_q}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (touch)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {3_q, 0_q}, {3_q, 0_q}, {2_q, 0_q}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (contained)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {10_q, 0_q}, {3_q, 0_q}, {5_q, 0_q}));
}

TEST_CASE("segment_vs_segment_intersection_single_point - overlap (side by side)") {
    REQUIRE(!segment_vs_segment_intersection_single_point<rational>({0_q, 0_q}, {5_q, 0_q}, {4_q, 0_q}, {10_q, 0_q}));
}

TEST_CASE("segment_vs_segment_squared_distance - line vs line") {
    REQUIRE(segment_segment_squared_distance<rational>({0_q, 0_q, 0_q}, {5_q, 0_q, 0_q}, {3_q, 3_q, 8_q}, {3_q, -3_q, 8_q}) == 64_q);
}

TEST_CASE("segment_vs_segment_squared_distance - parallel") {
    REQUIRE(segment_segment_squared_distance<rational>({0_q, 0_q, 0_q}, {5_q, 0_q, 0_q}, {1_q, 2_q, 0_q}, {4_q, 2_q, 0_q}) == 4_q);
}

TEST_CASE("segment_vs_segment_squared_distance - colinear") {
    REQUIRE(segment_segment_squared_distance<rational>({0_q, 0_q, 0_q}, {5_q, 0_q, 0_q}, {7_q, 0_q, 0_q}, {10_q, 0_q, 0_q}) == 4_q);
}
