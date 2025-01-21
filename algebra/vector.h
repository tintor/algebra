#pragma once

namespace algebra {

template<int D, typename T>
struct _Vec {
    static constexpr int dim = D;
    constexpr const T& operator[](int i) const { return reinterpret_cast<const T*>(this)[i]; }
    constexpr T& operator[](int i) { return reinterpret_cast<T*>(this)[i]; }
};

template<int D, typename T>
struct Vec {
};

template<typename T>
struct Vec<2, T> : public _Vec<2, T> {
    T x, y;
    Vec() {}
    Vec(const T& x, const T& y) : x(x), y(y) {}
};

template<typename T>
struct Vec<3, T> : public _Vec<3, T> {
    T x, y, z;
    Vec() {}
    Vec(const T& x, const T& y, const T& z) : x(x), y(y), z(z) {}
};

template<typename T>
struct Vec<4, T> : public _Vec<4, T> {
    T x, y, z, w;
    Vec() {}
    Vec(const T& x, const T& y, const T& z, const T& w) : x(x), y(y), z(z), w(w) {}
};

template<typename T> using Vec2 = Vec<2, T>;
template<typename T> using Vec3 = Vec<3, T>;
template<typename T> using Vec4 = Vec<4, T>;

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
constexpr int argmax_abs(const Vec<D, T>& a) {
    int m = 0;
    for (int i = 1; i < D; i++)
        if (abs_greater(a[i], a[m]))
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

}
