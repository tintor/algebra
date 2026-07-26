#include "algebra/segment_segment_intersection.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"

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

TEST_CASE("segment_segment_intersection_single_point - cross (vertical | horisontal)") {
    rational m, n;
    REQUIRE(segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,1}, {1,-1}));
    REQUIRE(segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,1}, {1,-1}, &m, &n));
    REQUIRE(m == 1/3_q);
    REQUIRE(n == 1/2_q);
}

TEST_CASE("segment_segment_intersection_single_point - cross (diagonal)") {
    rational m, n;
    REQUIRE(segment_segment_intersection_single_point<rational>({0,0}, {9,9}, {9,0}, {0,9}, &m, &n));
    REQUIRE(m == 1/2_q);
    REQUIRE(n == 1/2_q);
}

TEST_CASE("segment_segment_intersection_single_point - disjoint (cross)") {
    REQUIRE(!segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,3}, {1,1}));
}

TEST_CASE("segment_segment_intersection_single_point - disjoint (collinear)") {
    REQUIRE(!segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {3 + pow(10_q, 20), 0}, {5,0}));
}

TEST_CASE("segment_segment_intersection_single_point - parallel") {
    REQUIRE(!segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {0,2}, {3,2}));
}

TEST_CASE("segment_segment_intersection_single_point - touch (collinear)") {
    rational m, n;
    REQUIRE(segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {3,0}, {5,0}, &m, &n));
    REQUIRE(m == 1);
    REQUIRE(n == 0);
}

TEST_CASE("segment_segment_intersection_single_point - touch (T)") {
    rational m, n;
    REQUIRE(segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {1,0}, {7,7}, &m, &n));
    REQUIRE(m == 1/3_q);
    REQUIRE(n == 0);
}

TEST_CASE("segment_segment_intersection_single_point - overlap (identical)") {
    REQUIRE(!segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {0,0}, {3,0}));
}

TEST_CASE("segment_segment_intersection_single_point - overlap (touch)") {
    REQUIRE(!segment_segment_intersection_single_point<rational>({0,0}, {3,0}, {3,0}, {2,0}));
}

TEST_CASE("segment_segment_intersection_single_point - overlap (contained)") {
    REQUIRE(!segment_segment_intersection_single_point<rational>({0,0}, {10,0}, {3,0}, {5,0}));
}

TEST_CASE("segment_segment_intersection_single_point - overlap (side by side)") {
    REQUIRE(!segment_segment_intersection_single_point<rational>({0,0}, {5,0}, {4,0}, {10,0}));
}

// SegmentParams documents that M = A + (B - A) * s and N = C + (D - C) * t are the two endpoints
// of the overlap, so s parametrizes AB and t parametrizes CD. Checking that contract directly,
// against the point returning overload, which computes the same two endpoints.
static void require_overlap(Vec2<rational> a, Vec2<rational> b, Vec2<rational> c, Vec2<rational> d,
                            Vec2<rational> expected_m, Vec2<rational> expected_n) {
    const auto r = segment_segment_intersection_param<rational>(a, b, c, d);
    REQUIRE(std::holds_alternative<SegmentParams<rational>>(r));
    const auto p = std::get<SegmentParams<rational>>(r);
    REQUIRE(a + (b - a) * p.s == expected_m);
    REQUIRE(c + (d - c) * p.t == expected_n);

    // the two overloads must agree on the endpoints
    const auto q = segment_segment_intersection<rational>(a, b, c, d);
    REQUIRE(std::holds_alternative<std::pair<Vec2<rational>, Vec2<rational>>>(q));
    const auto pair = std::get<std::pair<Vec2<rational>, Vec2<rational>>>(q);
    REQUIRE(pair.first == expected_m);
    REQUIRE(pair.second == expected_n);
}

TEST_CASE("segment_segment_intersection_param - overlap (CD inside AB)") {
    require_overlap({0,0}, {10,0}, {2,0}, {8,0}, {2,0}, {8,0});
}

TEST_CASE("segment_segment_intersection_param - overlap (AB inside CD)") {
    require_overlap({2,0}, {8,0}, {0,0}, {10,0}, {2,0}, {8,0});
}

TEST_CASE("segment_segment_intersection_param - overlap (side by side)") {
    require_overlap({0,0}, {5,0}, {4,0}, {10,0}, {4,0}, {5,0});
}

TEST_CASE("segment_segment_intersection_param - overlap (identical)") {
    require_overlap({0,0}, {3,0}, {0,0}, {3,0}, {0,0}, {3,0});
}

TEST_CASE("segment_segment_intersection_param - overlap (identical, reversed)") {
    require_overlap({0,0}, {3,0}, {3,0}, {0,0}, {0,0}, {3,0});
}

TEST_CASE("segment_segment_intersection_param - overlap (diagonal)") {
    require_overlap({0,0}, {6,6}, {2,2}, {9,9}, {2,2}, {6,6});
}

TEST_CASE("segment_segment_intersects") {
    REQUIRE(segment_segment_intersects<rational>({0,0}, {3,0}, {0,2}, {3,2}) == 0); // parallel
    REQUIRE(segment_segment_intersects<rational>({0,0}, {3,0}, {1,1}, {1,-1}) == 1); // cross
    REQUIRE(segment_segment_intersects<rational>({0,0}, {3,0}, {3,0}, {5,0}) == 1); // touch
    REQUIRE(segment_segment_intersects<rational>({0,0}, {5,0}, {4,0}, {10,0}) == 2); // overlap
}

TEST_CASE("segment_segment_intersection - degenerate segments") {
    // both degenerate and equal
    const auto r = segment_segment_intersection<rational>({1,1}, {1,1}, {1,1}, {1,1});
    REQUIRE(std::holds_alternative<Vec2<rational>>(r));
    REQUIRE(std::get<Vec2<rational>>(r) == Vec2<rational>(1, 1));
    // both degenerate and different
    REQUIRE(segment_segment_intersects<rational>({1,1}, {1,1}, {2,2}, {2,2}) == 0);
    // AB degenerate, on CD
    REQUIRE(segment_segment_intersects<rational>({2,2}, {2,2}, {0,0}, {5,5}) == 1);
    // AB degenerate, off CD
    REQUIRE(segment_segment_intersects<rational>({2,3}, {2,3}, {0,0}, {5,5}) == 0);
    // CD degenerate, on AB
    REQUIRE(segment_segment_intersects<rational>({0,0}, {5,5}, {2,2}, {2,2}) == 1);
    // CD degenerate, off AB
    REQUIRE(segment_segment_intersects<rational>({0,0}, {5,5}, {2,3}, {2,3}) == 0);
}
