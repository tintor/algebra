#pragma once
#include "algebra/rational_func.h"
#include "algebra/vector.h"

namespace algebra {

// Vec<rational> x integral
VEC_OP_vs_sv(+, rational, rational, std::integral auto)
VEC_OP_vs_sv(-, rational, rational, std::integral auto)
VEC_OP_vs_sv(*, rational, rational, std::integral auto)
VEC_OP_vs_sv(/, rational, rational, std::integral auto)

// Vec<rational> x rational
VEC_OP_vs_sv(+, rational, rational, rational)
VEC_OP_vs_sv(-, rational, rational, rational)
VEC_OP_vs_sv(*, rational, rational, rational)
VEC_OP_vs_sv(/, rational, rational, rational)

// Vec<rational> x integer
VEC_OP_vs_sv(+, rational, rational, integer)
VEC_OP_vs_sv(-, rational, rational, integer)
VEC_OP_vs_sv(*, rational, rational, integer)
VEC_OP_vs_sv(/, rational, rational, integer)

// Vec<rational> x natural
VEC_OP_vs_sv(+, rational, rational, natural)
VEC_OP_vs_sv(-, rational, rational, natural)
VEC_OP_vs_sv(*, rational, rational, natural)
VEC_OP_vs_sv(/, rational, rational, natural)

using qvec2 = Vec2<rational>;
using qvec3 = Vec3<rational>;
using qvec4 = Vec4<rational>;

}
