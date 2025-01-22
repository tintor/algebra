#pragma once
#include "algebra/vector.h"
#include "algebra/solve_linear.h"

namespace algebra {

template<typename T>
constexpr T segment_segment_squared_distance(const Vec3<T>& pa, const Vec3<T>& pb, const Vec3<T>& qa, const Vec3<T>& qb) {
    const auto A = pb - pa, B = qb - qa, C = pa - qa;
    const auto aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);

    // line/line
    T d = aa * bb - ab * ab;
    T s = ab * bc - bb * ac;
    T t = aa * bc - ab * ac;
    if (d < 0) {
        negate(d);
        negate(s);
        negate(t);
    }
    if (d > 0 && 0 <= s && 0 <= t && s <= d && t <= d) {
        auto e = C;
        e += A * (s / d);
        e -= B * (t / d);
        return dot2(e);
    }

    // vertex/vertex
    auto CpA = C + A;
    auto CmB = C - B;
    T m = dot2(C);
    minimize(m, dot2(CpA));
    minimize(m, dot2(CmB));
    minimize(m, dot2(CpA - B));

    // line/vertex
    if (aa > 0 && 0 <= -ac && -ac <= aa)
        minimize(m, dot2(C - A * (ac / aa)));
    if (aa > 0 && 0 <= ab - ac && ab - ac <= aa)
        minimize(m, dot2(CmB + A * ((ab - ac) / aa)));
    if (bb > 0 && 0 <= bc && bc <= bb)
        minimize(m, dot2(B * (bc / bb) - C));
    if (bb > 0 && 0 <= ab + bc && ab + bc <= bb)
        minimize(m, dot2(B * ((ab + bc) / bb) - CpA));

    return m;
}

}
