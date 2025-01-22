#pragma once
#include "algebra/vector.h"
#include "algebra/solve_linear.h"
#include "algebra/point_segment_squared_distance.h"
#include "algebra/segment_segment_squared_distance.h"
#include "algebra/segment_segment_intersection.h"
#include <concepts>
#include <stdexcept>

namespace algebra {

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
    simplify(dir.x, dir.y, dir.z); // TODO generalize this (for floats it would be normalization, for bigints division by gcd)

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
