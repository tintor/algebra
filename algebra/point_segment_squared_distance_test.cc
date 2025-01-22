#include "algebra/point_segment_squared_distance.h"
#include "algebra/rational_func.h"
#include "algebra/__test.h"

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
