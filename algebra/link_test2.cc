#include "algebra/geometry.h"
#include "algebra/rational_vector.h"
#include "algebra/real.h"
#include "algebra/xrational.h"
#include "algebra/expr.h"
#include "algebra/dual.h"

int second_translation_unit() {
    algebra::natural a = 1;
    a <<= 100;
    return (a > 0u) ? 42 : 0;
}
