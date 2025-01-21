#pragma once
#include "algebra/vector.h"
#include "algebra/solve_linear.h"
#include <concepts>
#include <variant>
#include <stdexcept>

namespace algebra {

template<typename T>
constexpr void minimize(T& a, const T& b) { if (b < a) a = b; }

template<typename T>
constexpr T segment_segment_squared_distance(const Vec3<T>& pa, const Vec3<T>& pb, const Vec3<T>& qa, const Vec3<T>& qb) {
    auto A = pb - pa, B = qb - qa, C = pa - qa;
    auto aa = dot(A, A), bb = dot(B, B), ab = dot(A, B), ac = dot(A, C), bc = dot(B, C);

    // line/line
    T d = aa * bb - ab * ab;
    T s = ab * bc - bb * ac;
    T t = aa * bc - ab * ac;
    if ((d > 0 && 0 <= s && s <= d && 0 <= t && t <= d) || (d < 0 && d <= s && s <= 0 && d <= t && t <= 0))
        return dot2(C + A * (s / d) - B * (t / d));

    // vertex/vertex
    auto CpA = C + A;
    auto CmB = C - B;
    T m = dot2(C);
    minimize(m, dot2(CpA));
    minimize(m, dot2(CmB));
    minimize(m, dot2(CpA - B));

    // line/vertex
    if (aa > 0 && 0 >= ac && -ac <= aa)
        minimize(m, dot2(C - A * (ac / aa)));
    if (aa > 0 && 0 <= ab - ac && ab - ac <= aa)
        minimize(m, dot2(CmB + A * ((ab - ac) / aa)));
    if (bb > 0 && 0 <= bc && bc <= bb)
        minimize(m, dot2(B * (bc / bb) - C));
    if (bb > 0 && 0 <= ab + bc && (ab + bc) <= bb)
        minimize(m, dot2(B * ((ab + bc) / bb) - CpA));

    return m;
}

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

struct None {};
struct Any {};

// monospace: disjoint
// PointParams: intersection is single point, returns M such that M = A + (B - A) * s = C + (D - C) * t
// SegmentParams: intersection is line segment, returns (M, N) such that M = A + (B - A) * s and N = C + (D - C) * t
template<typename T>
std::variant<None, PointParams<T>, SegmentParams<T>> segment_vs_segment_intersection_param(
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
std::variant<None, Vec2<T>, std::pair<Vec2<T>, Vec2<T>>> segment_vs_segment_intersection(
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
constexpr bool segment_vs_segment_intersection_single_point(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, const Vec2<T>& d, T* s = nullptr, T* t = nullptr) {
    const auto res = segment_vs_segment_intersection_param(a, b, c, d);
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
    return segment_vs_segment_intersection(a, b, c, d).index();
}

template<typename T>
struct Line3 {
    Vec3<T> orig, dir;
};

// plane equation is: f(x) = (n * x + d) / sqrt(den)
// it is sufficient that T is integer (instead of rational)
template<typename T>
struct Plane3 {
    Vec3<T> n; // n != 0
    T d;
    T den; // n is not divided to avoid imprecise irrational numbers
};

template<typename T>
constexpr Vec3<T> operator*(const Plane3<T>& a, const Plane3<T>& b) { return {a.x * b.x, a.y * b.y, a.z * b.z}; }

template<typename T>
constexpr bool operator==(const Plane3<T>& a, const Plane3<T>& b) {
    if (a.den == 1 && b.den == 1)
        return a.n == b.n && a.d == b.d;
    Vec3<T> ans = a.n * a.n;
    Vec3<T> bns = b.n * b.n;
    return same_sign(a.n, b.n) && same_sign(a.d, b.d) && ans * b.den == bns * a.den && a.d * a.d * b.den == b.d * b.d * a.den;
}

// are there T A and T B such that aA + bB = 0
template<typename T>
constexpr bool are_parallel(const Vec3<T>& a, const Vec3<T>& b) {
    if (is_zero(a) || is_zero(b))
        return true;
    return a.x * b.y == b.x * a.y && a.x * b.z == b.x * a.z && a.y * b.z == b.y * a.z;
}

// A + B*x = 0
template<typename T>
constexpr std::variant<None, T, Any> solve_linear(const Vec3<T> a, const Vec3<T>& b) {
    if (is_zero(a) && is_zero(b))
        return Any();
    for (int i = 0; i < 3; i++)
        if (b[i] != 0)
            return -a[i] / b[i];
    return None();
}

// A + B*x = 0
template<typename T>
constexpr std::variant<None, T, Any> solve_linear(const Vec2<T> a, const Vec2<T>& b) {
    if (is_zero(a) && is_zero(b))
        return Any();
    for (int i = 0; i < 2; i++)
        if (b[i] != 0)
            return -a[i] / b[i];
    return None();
}

template<typename T>
constexpr void simplify(Vec3<T>& a) {
}

// result can be:
// - empty
// - line
// - plane
template<typename T>
std::variant<None, Line3<T>, Plane3<T>> plane_intersection(const Plane3<T>& a, const Plane3<T>& b) {
    if (are_parallel(a.n, b.n))
        return (a == b || a == -b) ? a : None();

    // result is a line
    Vec3<T> dir = cross(a.n, b.n);
    simplify(dir.x, dir.y, dir.z);

    const T ad = dot(a.n, dir);
    const T bd = dot(b.n, dir);
    const auto origin = solve_linear(a.d * bd - b.d * ad, a.n * bd - b.n * ad);
    return Line3<T>{std::get<Vec3<T>>(origin), dir};
}

// result can be:
// - empty
// - point
// - line
template<typename T>
std::variant<None, Vec3<T>, Line3<T>> line_plane_intersection(const Line3<T>& p, const Plane3<T>& q) {
    auto res = solve_linear(dot(q.n, p.origin) + q.d, dot(p.dir, q.n));
    if (std::holds_alternative<None>(res))
        return None();
    if (std::holds_alternative<Any>(res))
        return p;
    return p.orig + p.dir * std::get<T>(res);
}

// result can be:
// - empty
// - point
// - line
// - plane
template<typename T>
std::variant<None, Vec3<T>, Line3<T>, Plane3<T>> plane_intersection(const Plane3<T>& a, const Plane3<T>& b, const Plane3<T>& c) {
    Vec3<T> D{a.d, b.d, c.d};
    Vec3<T> X{a.n.x, b.n.x, c.n.x};
    Vec3<T> Y{a.n.y, b.n.y, c.n.y};
    Vec3<T> Z{a.n.z, b.n.z, c.n.z};
    Vec3<T> m;
    T det;
    if (__solve_linear(D, X, Y, Z, m.x, m.y, m.z, det))
        return m / D;

    const bool ab = are_parallel(a.n, b.n);
    const bool bc = are_parallel(b.n, c.n);

    // all three planes parallel
    if (ab && bc)
        return ((a == b || a == -b) && (c == b || c == -b)) ? a : None();

    if (ab)
        return (a == b || a == -b) ? plane_intersection(a, c) : None();
    if (bc)
        return (b == c || b == -c) ? plane_intersection(a, b) : None();
    if (are_parallel(a, c))
        return (a == c || a == -c) ? plane_intersection(a, b) : None();

    auto res = plane_intersection(a, b); // res must be a line
    return line_plane_intersection(std::get<Line3<T>>(res), c);
}

}
