#pragma once
#include "algebra/vector.h"
namespace algebra {

template<typename T>
constexpr T point_segment_squared_distance(const Vec3<T>& p, const Vec3<T>& a, const Vec3<T>& b) {
    const auto B = b - a;
    const auto P = p - a;
    T d = dot(B, B);
    T s = dot(P, B);
    if (d < 0) {
        d = -d;
        s = -s;
    }
    return (d > 0 && 0 <= s && s <= d) ? dot2(B * (s / d) - P) : min(dot2(P), dot2(b - p));
}

}
