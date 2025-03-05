#pragma once
#include <format>

namespace algebra {

template<typename T>
struct dual {
    T real;
    T dual;
};

template<typename T>
dual<T> operator+(const dual<T>& a, const dual<T>& b) { return {a.real + b.real, a.dual + b.dual}; }
template<typename T>
dual<T> operator-(const dual<T>& a, const dual<T>& b) { return {a.real - b.real, a.dual - b.dual}; }
template<typename T>
dual<T> operator*(const dual<T>& a, const dual<T>& b) { return {a.real * b.real, a.real * b.dual + b.real * a.dual}; }
template<typename T>
dual<T> operator/(const dual<T>& a, const dual<T>& b) { return {a.real / b.real, (a.dual * b.real - b.dual * a.real) / (b.real * b.real)}; }
template<typename T>
dual<T> sqrt(const dual<T>& a) { auto s = sqrt(a.real); return {s, a.dual / s / 2}; }
template<typename T>
dual<T> pow(const dual<T>& a, const T& k) { return {pow(a.real, k), k * a.dual * pow(a.real, k - T(1))}; }
template<typename T>
dual<T> exp(const dual<T>& a) { auto s = exp(a.real); return {s, s * a.dual}; }
template<typename T>
dual<T> log(const dual<T>& a) { return {log(a.real), a.dual / a.real}; }
template<typename T>
dual<T> sin(const dual<T>& a) { return {sin(a.real), a.dual * cos(a.real)}; }
template<typename T>
dual<T> cos(const dual<T>& a) { return {cos(a.real), -a.dual * sin(a.real)}; }
template<typename T>
dual<T> tan(const dual<T>& a) { return {tan(a.real), a.dual / (cos(a.real) * cos(a.real))}; }
template<typename T>
dual<T> atan(const dual<T>& a) { return {atan(a.real), a.dual / (T(1) + a.real * a.real)}; }
template<typename T>
dual<T> abs(const dual<T>& a) { return {abs(a.real), signum(a.real) * a.dual}; }

}

template<typename T>
struct std::formatter<algebra::dual<T>, char> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }

    constexpr auto format(const algebra::dual<T>& a, auto& ctx) const {
        using namespace algebra;
        auto it = ctx.out();
        if (a.real == T(0)) {
            if (a.dual == T(1))
                return std::format_to(it, "eps");
            if (a.dual == T(0))
                return std::format_to(it, "0");
            return std::format_to(it, "{}*eps", a.dual);
        }
        it = std::format_to(it, "{}", a.real);
        if (a.dual > T(0))
            it = std::format_to(it, "+{}*eps", a.dual);
        if (T(0) > a.dual)
            it = std::format_to(it, "{}*eps", a.dual);
        return it;
    }
};
