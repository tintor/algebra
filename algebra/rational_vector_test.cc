#include "algebra/rational_vector.h"
#include "algebra/__test.h"

TEST_CASE("*") {
    qvec3 a = {1, 2, 3};
    REQUIRE(a * 2 == qvec3(2, 4, 6));
    REQUIRE(a * 2_n == qvec3(2, 4, 6));
    REQUIRE(a * 2_i == qvec3(2, 4, 6));
    REQUIRE(a * 1/2_q == qvec3(1/2_q, 1, 3/2_q));
}
