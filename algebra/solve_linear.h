#pragma once
#include "algebra/vector.h"
#include <variant>

namespace algebra {

template<typename T>
constexpr T determinant(const Vec2<T>& a, const Vec2<T>& b) {
    return a.x*b.y - a.y*b.x;
}

template<typename T>
constexpr T determinant(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c) {
    return a.x*b.y*c.z + b.x*c.y*a.z + c.x*a.y*b.z - c.x*b.y*a.z - b.x*a.y*c.z - a.x*c.y*b.z;
}

struct None {};
struct Any {};

// TODO: Is distinction of None and Any important? solve_linear below `returns false` in both cases
// A + B*x = 0
template<int D, typename T>
constexpr std::variant<None, T, Any> solve_linear(const Vec<D, T> a, const Vec<D, T>& b) {
    if (is_zero(a) && is_zero(b))
        return Any();
    for (int i = 0; i < D; i++)
        if (b[i] != 0)
            return -a[i] / b[i];
    return None();
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

// A + sB + tC = 0
// return (s, t) or false if no unique solution
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

// A + sB + tC + rD = 0
// return (s, t, r) or false if no unique solution
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

}
