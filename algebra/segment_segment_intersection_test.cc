#include "algebra/segment_segment_intersection.h"
#include "algebra/rational_vector.h"
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
