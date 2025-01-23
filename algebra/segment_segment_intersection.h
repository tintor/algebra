#pragma once
#include "algebra/solve_linear.h"
#include <variant>

namespace algebra {

constexpr void negate(std::floating_point auto& a) { a = -a; }

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
    if (__solve_linear(a - c, a - b, c - d, s, t, det)) {
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
        // overlap
        return SegmentParams{
            (A > C) ? T(swap_ab ? 1 : 0) : T(swap_cd ? 1 : 0),
            (B < D) ? T(swap_ab ? 0 : 1) : T(swap_cd ? 0 : 1)};
    }

    // parallel (non-colinear)
    return None();
}

template<typename T>
std::variant<None, Vec2<T>, std::pair<Vec2<T>, Vec2<T>>> segment_segment_intersection(
        const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, const Vec2<T>& d) {
    T s, t, det;
    // if not parallel AND not degenerate
    if (__solve_linear(a - c, a - b, c - d, s, t, det)) {
        if (det < 0) {
            negate(s);
            negate(t);
            negate(det);
        }
        if (s >= 0 && s <= det && t >= 0 && t <= det)
            return a - (a - b) * (s / det);
        return None();
    }

    // both degenerate
    if (a == b && c == d) {
        if (a == c)
            return a;
        return None();
    }

    // AB degenerate
    if (a == b) {
        if (ccw(a, c, d) == 0 && loose_order(c, a, d))
            return a;
        return None();
    }

    // CD degenerate
    if (c == d) {
        if (ccw(c, a, b) == 0 && loose_order(a, c, b))
            return c;
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
            return swap_ab ? a : b;
        if (D == A)
            return swap_ab ? b : a;
        // overlap
        return std::pair{
            (A > C) ? (swap_ab ? b : a) : (swap_cd ? d : c),
            (B < D) ? (swap_ab ? a : b) : (swap_cd ? c : d)};
    }

    // parallel (non-colinear)
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
int segment_vs_segment_intersects(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, const Vec2<T>& d) {
    return segment_segment_intersection(a, b, c, d).index();
}

}
