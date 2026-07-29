#pragma once
#include "algebra/solve_linear.h"
#include <variant>

namespace algebra {

constexpr void negate(std::floating_point auto& a) { a = -a; }

// Twice the signed area of the triangle, with the opposite sign to the usual convention and to
// signed_area2() in polygon2.h: this is negative for a counter clockwise triple. Only its zero is
// used here, which is the same either way.
template<typename T>
constexpr T ccw(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c) {
    return (a.x + b.x) * (a.y - b.y) + (b.x + c.x) * (b.y - c.y) + (c.x + a.x) * (c.y - a.y);
}

template<typename T>
constexpr bool abs_greater(const T& a, const T& b) {
    return abs(a) > abs(b);
}

template<typename T>
struct PointParams {
    T s, t;
};

template<typename T>
struct SegmentParams {
    T s, t;
};

// monospace: disjoint
// PointParams: intersection is single point, returns M such that M = A + (B - A) * s = C + (D - C) * t
// SegmentParams: intersection is line segment, returns (M, N) such that M = A + (B - A) * s and N = C + (D - C) * t
template<typename T>
std::variant<None, PointParams<T>, SegmentParams<T>> segment_segment_intersection_param(
        const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, const Vec2<T>& d) {
    T s, t, det;
    // if not parallel AND not degenerate
    if (__solve_linear(a - c, b - a, c - d, s, t, det)) {
        if (det < 0) {
            negate(s);
            negate(t);
            negate(det);
        }
        if (s >= 0 && s <= det && t >= 0 && t <= det)
            return PointParams<T>{s / det, t / det};
        return None();
    }

    // both degenerate
    if (a == b && c == d) {
        if (a == c)
            return PointParams<T>{T(0), T(0)};
        return None();
    }

    // AB degenerate
    if (a == b) {
        if (ccw(a, c, d) == 0 && loose_order(c, a, d))
            return PointParams<T>{T(0), div_colinear(a - c, d - c)};
        return None();
    }

    // CD degenerate
    if (c == d) {
        if (ccw(c, a, b) == 0 && loose_order(a, c, b))
            return PointParams<T>{div_colinear(c - a, b - a), T(0)};
        return None();
    }

    // if collinear
    if (ccw(a, b, c) == 0) {
        const int i = argmax_abs(a - b);
        T A = a[i];
        T B = b[i];
        T C = c[i];
        T D = d[i];

        const bool swap_ab = A > B;
        const bool swap_cd = C > D;
        if (swap_ab)
            std::swap(A, B);
        if (swap_cd)
            std::swap(C, D);

        if (B < C || D < A)
            return None();
        if (B == C)
            return PointParams{T(swap_ab ? 0 : 1), T(swap_cd ? 1 : 0)};
        if (D == A)
            return PointParams{T(swap_ab ? 1 : 0), T(swap_cd ? 0 : 1)};
        // overlap: [m, n] are its endpoints, m the lower and n the upper one along axis i. Each is
        // an endpoint of AB or of CD, whichever lies inside the other segment, so it can only be
        // expressed as 0 or 1 in that segment's own parameter. s must parametrize AB and t must
        // parametrize CD (see the comment above this function), so both need converting.
        const Vec2<T> m = (A > C) ? (swap_ab ? b : a) : (swap_cd ? d : c);
        const Vec2<T> n = (B < D) ? (swap_ab ? a : b) : (swap_cd ? c : d);
        return SegmentParams{div_colinear(m - a, b - a), div_colinear(n - c, d - c)};
    }

    // parallel (non-colinear)
    return None();
}

// The same intersection as points rather than parameters. Every branch of the parameter version
// already identifies the answer exactly, so this evaluates its result instead of repeating it: for an
// exact T the two agree by construction.
template<typename T>
std::variant<None, Vec2<T>, std::pair<Vec2<T>, Vec2<T>>> segment_segment_intersection(
        const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, const Vec2<T>& d) {
    const auto r = segment_segment_intersection_param(a, b, c, d);
    if (const PointParams<T>* p = std::get_if<PointParams<T>>(&r))
        return lerp(a, b, p->s);
    if (const SegmentParams<T>* s = std::get_if<SegmentParams<T>>(&r))
        return std::pair{lerp(a, b, s->s), lerp(c, d, s->t)};
    return None();
}

// returns (s, t), such that intersection point M = A + (B - A) * s, and M = C + (D - C) * t, or false if no unique solution
// Note: returns false in case of full or partial overlap
template<typename T>
constexpr bool segment_segment_intersection_single_point(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, const Vec2<T>& d, T* s = nullptr, T* t = nullptr) {
    const auto res = segment_segment_intersection_param(a, b, c, d);
    if (!std::holds_alternative<PointParams<T>>(res))
        return false;
    if (s)
        *s = std::get<PointParams<T>>(res).s;
    if (t)
        *t = std::get<PointParams<T>>(res).t;
    return true;
}

// 0 - disjoint
// 1 - intersection is single point
// 2 - intersection is non-degenerate line segment
template<typename T>
int segment_segment_intersects(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, const Vec2<T>& d) {
    return segment_segment_intersection(a, b, c, d).index();
}

}
