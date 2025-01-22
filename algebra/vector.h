#pragma once
#include <format>

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

#define __VEC_OP_vv(OP, TC, TA, TB) \
constexpr auto operator OP(const algebra::Vec<D, TA>& a, const algebra::Vec<D, TB>& b) { \
    algebra::Vec<D, TC> c; \
    for (int i = 0; i < D; i++) \
        c[i] = a[i] OP b[i]; \
    return c; \
}

#define __VEC_OP_vs(OP, TC, TA, TB) \
constexpr auto operator OP(const algebra::Vec<D, TA>& a, const TB& b) { \
    algebra::Vec<D, TC> c; \
    for (int i = 0; i < D; i++) \
        c[i] = a[i] OP b; \
    return c; \
}

#define __VEC_OP_sv(OP, TC, TA, TB) \
constexpr auto operator OP(const TA& a, const algebra::Vec<D, TB>& b) { \
    algebra::Vec<D, TC> c; \
    for (int i = 0; i < D; i++) \
        c[i] = a OP b[i]; \
    return c; \
}

#define VEC_OP(OP) \
template<int D, typename T> __VEC_OP_vv(OP, T, T, T) \
template<int D, typename T> __VEC_OP_vs(OP, T, T, T) \
template<int D, typename T> __VEC_OP_sv(OP, T, T, T)

#define VEC_OP_vs(OP, TC, TA, TB) template<int D> __VEC_OP_vs(OP, TC, TA, TB)
#define VEC_OP_sv(OP, TC, TA, TB) template<int D> __VEC_OP_sv(OP, TC, TA, TB)
#define VEC_OP_vs_sv(OP, TC, TA, TB) VEC_OP_vs(OP, TC, TA, TB) VEC_OP_sv(OP, TC, TB, TA)

VEC_OP(+)
VEC_OP(-)
VEC_OP(*)
VEC_OP(/)

#define VEC_AOP(OP) \
template<int D, typename T> \
constexpr algebra::Vec<D, T>& operator OP(algebra::Vec<D, T>& a, const algebra::Vec<D, T>& b) { \
    for (int i = 0; i < D; i++) \
        a[i] OP b[i]; \
    return a; \
} \
template<int D, typename T> \
constexpr algebra::Vec<D, T>& operator OP(algebra::Vec<D, T>& a, const T& b) { \
    for (int i = 0; i < D; i++) \
        a[i] OP b; \
    return a; \
}

VEC_AOP(+=)
VEC_AOP(-=)
VEC_AOP(*=)
VEC_AOP(/=)

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

template<int D, typename T>
constexpr auto operator-(const Vec<D, T>& a) {
    Vec<D, T> c;
    for (int i = 0; i < D; i++)
        c[i] = -a[i];
    return c;
}

template<typename T>
constexpr auto dot2(const T& a) { return dot(a, a); }

template<int D, typename T>
constexpr Vec<D, T> lerp(const Vec<D, T>& a, const Vec<D, T>& b, const T& t) {
    Vec<D, T> c = b;
    c -= a;
    c *= t;
    c += a;
    return c;
}

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
SWIZZLE3(x, z, y)
SWIZZLE3(x, z, w)
SWIZZLE3(y, x, z)
SWIZZLE3(y, z, x)
SWIZZLE3(y, z, w)
SWIZZLE3(z, x, y)
SWIZZLE3(z, y, x)

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

template<typename T>
constexpr bool same_sign(const T& a, const T& b) {
    if (sign(a) > 0)
        return sign(b) > 0;
    if (sign(a) < 0)
        return sign(b) < 0;
    return sign(b) == 0;
}

template<int D, typename T>
constexpr bool same_sign(const Vec<D, T>& a, const Vec<D, T>& b) {
    for (int i = 0; i < D; i++)
        if (!same_sign(a[i], b[i]))
            return false;
    return true;
}

// return k such that B*k = A (assuming k exists and is unique)
template<int D, typename T>
constexpr T div_colinear(const Vec<D, T>& a, const Vec<D, T>& b) {
    const int i = argmax_abs(b);
    return a[i] / b[i];
}

template<typename T>
constexpr T min(const T& a, const T& b) { return (a < b) ? a : b; }

template<typename T>
constexpr void minimize(T& a, const T& b) { if (b < a) a = b; }

constexpr uint64_t hash_fn_64bit(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdllu;
    return k;
}

}

template <int D, typename T>
struct std::formatter<algebra::Vec<D, T>, char> {
    constexpr auto parse(auto& ctx) {
        return ctx.begin();
    }

    constexpr auto format(const algebra::Vec<D, T>& a, auto& ctx) const {
        for (int i = 0; i < D; i++) {
            if (i > 0)
                std::format_to(ctx.out(), " ");
            std::format_to(ctx.out(), "{}", a[i]);
        }
        return ctx.out();
    }
};

template<int D, typename T>
struct std::hash<algebra::Vec<D, T>> {
    constexpr size_t operator()(const algebra::Vec<D, T>& a) const {
        uint64_t seed = 0;
        for (int i = 0; i < D; i++)
            seed = algebra::hash_fn_64bit(seed ^ std::hash<T>()(a[i]));
        return seed;
    }
};
