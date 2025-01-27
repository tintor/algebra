#pragma once
#include "algebra/rational_func.h"
#include "algebra/xrational.h"
#include "algebra/vector.h"

namespace algebra {

VEC_OP_vs_sv(+, rational, rational, integral auto)
VEC_OP_vs_sv(-, rational, rational, integral auto)
VEC_OP_vs_sv(*, rational, rational, integral auto)
VEC_OP_vs_sv(/, rational, rational, integral auto)

VEC_AOP_vs(+=, rational, integral auto)
VEC_AOP_vs(-=, rational, integral auto)
VEC_AOP_vs(*=, rational, integral auto)
VEC_AOP_vs(/=, rational, integral auto)

using qvec2 = Vec2<rational>;
using qvec3 = Vec3<rational>;
using qvec4 = Vec4<rational>;

VEC_OP_vs_sv(+, xrational, xrational, rational_like auto)
VEC_OP_vs_sv(-, xrational, xrational, rational_like auto)
VEC_OP_vs_sv(*, xrational, xrational, rational_like auto)
VEC_OP_vs_sv(/, xrational, xrational, rational_like auto)

VEC_AOP_vs(+=, xrational, rational_like auto)
VEC_AOP_vs(-=, xrational, rational_like auto)
VEC_AOP_vs(*=, xrational, rational_like auto)
VEC_AOP_vs(/=, xrational, rational_like auto)

using xvec2 = Vec2<xrational>;
using xvec3 = Vec3<xrational>;
using xvec4 = Vec4<xrational>;

}
