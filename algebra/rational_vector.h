#pragma once
#include "algebra/rational_func.h"
#include "algebra/xrational.h"
#include "algebra/vector.h"

namespace algebra {

// Vec<rational> x rational_like
VEC_OP_vs_sv(+, rational, rational, rational_like auto)
VEC_OP_vs_sv(-, rational, rational, rational_like auto)
VEC_OP_vs_sv(*, rational, rational, rational_like auto)
VEC_OP_vs_sv(/, rational, rational, rational_like auto)

using qvec2 = Vec2<rational>;
using qvec3 = Vec3<rational>;
using qvec4 = Vec4<rational>;

// Vec<xrational> x xrational_like
VEC_OP_vs_sv(+, xrational, xrational, xrational_like auto)
VEC_OP_vs_sv(-, xrational, xrational, xrational_like auto)
VEC_OP_vs_sv(*, xrational, xrational, xrational_like auto)
VEC_OP_vs_sv(/, xrational, xrational, xrational_like auto)

using xvec2 = Vec2<xrational>;
using xvec3 = Vec3<xrational>;
using xvec4 = Vec4<xrational>;

}
