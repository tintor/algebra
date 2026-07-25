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

// same set of points, opposite orientation
template<typename T>
constexpr Plane3<T> operator-(const Plane3<T>& a) { return {-a.n, -a.d, a.den}; }

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
// - point
// - line (when the line lies inside the plane)
template<typename T>
std::variant<None, Vec3<T>, Line3<T>> line_plane_intersection(const Line3<T>& p, const Plane3<T>& q) {
    // q.n * (p.orig + p.dir * t) + q.d == 0
    const T den = dot(p.dir, q.n);
    const T num = dot(q.n, p.orig) + q.d;
    if (den == 0) {
        if (num == 0)
            return p; // whole line is inside the plane
        return None();
    }
    return p.orig + p.dir * (-num / den);
}

// result can be:
// - empty
// - line
// - plane
// Note: the point on the resulting line generally has non-integer coordinates,
// so T needs exact division (rational) even though Plane3 itself works with integers.
template<typename T>
std::variant<None, Line3<T>, Plane3<T>> plane_intersection(const Plane3<T>& a, const Plane3<T>& b) {
    if (are_parallel(a.n, b.n)) {
        if (a == b || a == -b)
            return a;
        return None();
    }

    // result is a line
    Vec3<T> dir = cross(a.n, b.n);
    // point on both planes (note: this uses dir before it is scaled down below)
    const Vec3<T> orig = (cross(b.n, dir) * -a.d + cross(dir, a.n) * -b.d) / dot(dir, dir);
    simplify(dir.x, dir.y, dir.z); // TODO generalize this (for floats it would be normalization, for bigints division by gcd)
    return Line3<T>{orig, dir};
}

template<typename T>
constexpr std::variant<None, Vec3<T>, Line3<T>, Plane3<T>> __widen(const std::variant<None, Line3<T>, Plane3<T>>& a) {
    if (const Line3<T>* e = std::get_if<Line3<T>>(&a))
        return *e;
    if (const Plane3<T>* e = std::get_if<Plane3<T>>(&a))
        return *e;
    return None();
}

template<typename T>
constexpr std::variant<None, Vec3<T>, Line3<T>, Plane3<T>> __widen(const std::variant<None, Vec3<T>, Line3<T>>& a) {
    if (const Vec3<T>* e = std::get_if<Vec3<T>>(&a))
        return *e;
    if (const Line3<T>* e = std::get_if<Line3<T>>(&a))
        return *e;
    return None();
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
    // D + X*m.x + Y*m.y + Z*m.z == 0 is the three plane equations
    if (__solve_linear(D, X, Y, Z, m.x, m.y, m.z, det))
        return m / det;

    const bool ab = are_parallel(a.n, b.n);
    const bool bc = are_parallel(b.n, c.n);

    // all three planes parallel
    if (ab && bc) {
        if ((a == b || a == -b) && (c == b || c == -b))
            return a;
        return None();
    }

    if (ab) {
        if (a == b || a == -b)
            return __widen(plane_intersection(a, c));
        return None();
    }
    if (bc) {
        if (b == c || b == -c)
            return __widen(plane_intersection(a, b));
        return None();
    }
    if (are_parallel(a.n, c.n)) {
        if (a == c || a == -c)
            return __widen(plane_intersection(a, b));
        return None();
    }

    // normals are linearly dependent, but no two planes are parallel
    auto res = plane_intersection(a, b); // res must be a line
    return __widen(line_plane_intersection(std::get<Line3<T>>(res), c));
}

}
