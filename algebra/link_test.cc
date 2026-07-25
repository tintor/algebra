// Verifies that every header can be included from more than one translation unit,
// i.e. that no header defines a non-inline function or variable at namespace scope.
// The second translation unit is link_test2.cc.
#include "algebra/geometry.h"
#include "algebra/rational_vector.h"
#include "algebra/real.h"
#include "algebra/xrational.h"
#include "algebra/expr.h"
#include "algebra/dual.h"
#include "algebra/__test.h"

int second_translation_unit();

TEST_CASE("linked with a second translation unit") {
    REQUIRE(second_translation_unit() == 42);
}
