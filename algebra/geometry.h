#pragma once
#include <concepts>
#include <stdexcept>
#include <variant>
#include <print>

namespace algebra {

template<int D, typename T>
struct Vec {
};

template<typename T>
struct Vec<2, T> {
    static constexpr int length = 2;
    T x, y;
    constexpr const T& operator[](int i) const { return (&x)[i]; }
    constexpr T& operator[](int i) { return (&x)[i]; }
};

template<typename T>
struct Vec<3, T> {
    static constexpr int length = 3;
    T x, y, z;
    constexpr const T& operator[](int i) const { return (&x)[i]; }
    constexpr T& operator[](int i) { return (&x)[i]; }
};

template<typename T>
struct Vec<4, T> {
    static constexpr int length = 4;
    T x, y, z, w;
    constexpr const T& operator[](int i) const { return (&x)[i]; }
    constexpr T& operator[](int i) { return (&x)[i]; }
};

template<typename T> using Vec2 = Vec<2, T>;
template<typename T> using Vec3 = Vec<3, T>;
template<typename T> using Vec4 = Vec<4, T>;

// Consider: Variant of Vec3<rational> where x, y, z share the same denominator

#define VEC_OP(OP) \
template<int D, typename T> \
constexpr Vec<D, T> operator OP(const Vec<D, T>& a, const Vec<D, T>& b) { \
    Vec<D, T> c; \
    for (int i = 0; i < D; i++) \
        c[i] = a[i] OP b[i]; \
    return c; \
} \
template<int D, typename T> \
constexpr Vec<D, T> operator OP(const Vec<D, T>& a, const T& b) { \
    Vec<D, T> c; \
    for (int i = 0; i < D; i++) \
        c[i] = a[i] OP b; \
    return c; \
} \
template<int D, typename T> \
constexpr Vec<D, T> operator OP(const T& a, const Vec<D, T>& b) { \
    Vec<D, T> c; \
    for (int i = 0; i < D; i++) \
        c[i] = a OP b[i]; \
    return c; \
}

VEC_OP(+)
VEC_OP(-)
VEC_OP(*)
VEC_OP(/)

template<int D, typename T>
constexpr bool operator==(const Vec<D, T>& a, const Vec<D, T>& b) {
    for (int i = 0; i < D; i++)
        if (a[i] != b[i])
            return false;
    return true;
}

template<int D, typename T>
constexpr bool is_zero(const Vec<D, T>& a) {
    for (int i = 0; i < D; i++)
        if (a[i] != 0)
            return false;
    return true;
}

template<int D, typename T>
constexpr Vec<D, T> abs(const Vec<D, T>& a) {
    Vec<D, T> c;
    for (int i = 0; i < D; i++)
        c[i] = abs(a[i]);
    return c;
}

template<int D, typename T>
constexpr int argmax(const Vec<D, T>& a) {
    int m = 0;
    for (int i = 1; i < D; i++)
        if (a[i] > a[m])
            m = i;
    return m;
}

template<int D, typename T>
constexpr T dot(const Vec<D, T>& a, const Vec<D, T>& b) {
    T c = 0;
    for (int i = 0; i < D; i++)
        c += a[i] * b[i];
    return c;
}

template<typename T>
constexpr auto dot2(const T& a) { return dot(a, a); }

#define SWIZZLE2(A, B) template<int D, typename T> constexpr Vec<2, T> A ## B(const Vec<D, T>& v) { return {v.A, v.B}; }
#define SWIZZLE3(A, B, C) template<int D, typename T> constexpr Vec<3, T> A ## B ## C(const Vec<D, T>& v) { return {v.A, v.B, v.C}; }
#define SWIZZLE4(A, B, C, D) template<int D, typename T> constexpr Vec<4, T> A ## B ## C ## D(const Vec<D, T>& v) { return {v.A, v.B, v.C, v.D}; }

SWIZZLE2(x, y)
SWIZZLE2(x, z)
SWIZZLE2(y, x)
SWIZZLE2(y, z)
SWIZZLE2(z, x)
SWIZZLE2(z, y)

SWIZZLE3(x, y, z)
SWIZZLE3(x, y, w)
SWIZZLE3(x, z, w)
SWIZZLE3(y, z, x)
SWIZZLE3(y, z, w)
SWIZZLE3(z, x, y)

template<typename T> constexpr T cross(const Vec2<T>& a, const Vec2<T>& b) { return {a.x * b.y - a.y * b.x}; }
template<typename T> constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) { return {cross(yz(a), yz(b)), cross(xz(a), xz(b)), cross(xy(a), xy(b))}; }

template<typename T>
constexpr void minimize(T& a, const T& b) { if (b < a) a = b; }

template<typename T>
constexpr bool order(const T& a, const T& b, const T& c) {
    return (a < b && b < c) || (a > b && b > c) || (a == b && b == c);
}

template<typename T>
constexpr bool strict_order(const T& a, const T& b, const T& c) {
    return (a < b && b < c) || (a > b && b > c);
}

template<typename T>
constexpr bool loose_order(const T& a, const T& b, const T& c) {
    return (a <= b && b <= c) || (a >= b && b >= c);
}

template<typename T>
constexpr bool order(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c) {
    return order(a.x, b.x, c.x) && order(a.y, b.y, c.y);
}

template<typename T>
constexpr bool strict_order(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c) {
    return strict_order(a.x, b.x, c.x) && strict_order(a.y, b.y, c.y);
}

template<typename T>
constexpr bool loose_order(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c) {
    return loose_order(a.x, b.x, c.x) && loose_order(a.y, b.y, c.y);
}

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

template<typename T>
constexpr T determinant(const Vec2<T>& a, const Vec2<T>& b) {
    return a.x*b.y - a.y*b.x;
}

template<typename T>
constexpr T determinant(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c) {
    return a.x*b.y*c.z + b.x*c.y*a.z + c.x*a.y*b.z - c.x*b.y*a.z - b.x*a.y*c.z - a.x*c.y*b.z;
}

// A + sB + tC = 0
// return (s, t) or false if no unique solution
template<typename T>
constexpr bool solve_linear(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, T* s, T* t) {
    auto det = determinant(b, c);
    if (det == 0)
        return false;
    if (s)
        *s = determinant(a, c) / det;
    if (t)
        *t = determinant(a, b) / det;
    return true;
}

// similar to above, but with delayed division
template<typename T>
constexpr bool __solve_linear(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, T& s, T& t, T& det) {
    det = determinant(b, c);
    if (det == 0)
        return false;
    s = determinant(a, c);
    t = determinant(a, b);
    return true;
}

// A + sB + tC + rD = 0
// return (s, t, r) or false if no unique solution
template<typename T>
constexpr bool solve_linear(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c, const Vec3<T>& d, T* s, T* t, T* r) {
    auto det = determinant(b, c, d);
    if (det == 0)
        return false;
    if (s)
        *s = determinant(a, c, d) / det;
    if (t)
        *t = determinant(a, b, d) / det;
    if (r)
        *r = determinant(a, b, c) / det;
    return true;
}

// similar to above, but with delayed division
template<typename T>
constexpr bool __solve_linear(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c, const Vec3<T>& d, T& s, T& t, T& r, T& det) {
    det = determinant(b, c, d);
    if (det == 0)
        return false;
    s = determinant(a, c, d);
    t = determinant(a, b, d);
    r = determinant(a, b, c);
    return true;
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

// return k such that B*k = A (assuming B != 0)
template<typename T>
constexpr T div_colinear(const Vec2<T>& a, const Vec2<T>& b) {
    return abs_greater(b.x, b.y) ? (a.x / b.x) : (a.y / b.y);
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
        const int i = argmax(abs(a - b));
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
        const int i = argmax(abs(a - b));
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
std::variant<None, PointParams<T>, SegmentParams<T>> segment_vs_segment_intersection_params(
        const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c, const Vec3<T>& d) {
    throw std::runtime_error("not implemented");
}

template<typename T>
std::variant<None, Vec3<T>, std::pair<Vec3<T>, Vec3<T>>> segment_vs_segment_intersection(
        const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c, const Vec3<T>& d) {
    throw std::runtime_error("not implemented");
}

template<typename T>
int segment_vs_segment_intersects(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c, const Vec3<T>& d) {
    throw std::runtime_error("not implemented");
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
constexpr bool same_sign(const T& a, const T& b) {
    if (sign(a) > 0)
        return sign(b) > 0;
    if (sign(a) < 0)
        return sign(b) < 0;
    return sign(b) == 0;
}

template<typename T>
constexpr bool same_sign(const Vec3<T>& a, const Vec3<T>& b) {
    return same_sign(a.x, b.x) && same_sign(a.y, b.y) && same_sign(a.z, b.z);
}

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
