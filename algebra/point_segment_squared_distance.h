#pragma once
#include "algebra/vector.h"
namespace algebra {

// Note: T needs exact division (rational), the projection is not integral in general.
template<typename T>
constexpr T point_segment_squared_distance(const Vec3<T>& p, const Vec3<T>& a, const Vec3<T>& b) {
    const auto B = b - a;
    const auto P = p - a;
    const T d = dot(B, B); // a sum of squares, never negative
    const T s = dot(P, B);
    return (d > 0 && 0 <= s && s <= d) ? dot2(B * (s / d) - P) : min(dot2(P), dot2(b - p));
}

}
