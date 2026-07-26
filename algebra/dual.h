#pragma once
#include <cmath>
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
dual<T> sqrt(const dual<T>& a) { using std::sqrt; auto s = sqrt(a.real); return {s, a.dual / s / 2}; }
template<typename T>
dual<T> pow(const dual<T>& a, const T& k) { using std::pow; return {pow(a.real, k), k * a.dual * pow(a.real, k - T(1))}; }
template<typename T>
dual<T> exp(const dual<T>& a) { using std::exp; auto s = exp(a.real); return {s, s * a.dual}; }
template<typename T>
dual<T> log(const dual<T>& a) { using std::log; return {log(a.real), a.dual / a.real}; }
template<typename T>
dual<T> sin(const dual<T>& a) { using std::sin; using std::cos; return {sin(a.real), a.dual * cos(a.real)}; }
template<typename T>
dual<T> cos(const dual<T>& a) { using std::sin; using std::cos; return {cos(a.real), -a.dual * sin(a.real)}; }
template<typename T>
dual<T> tan(const dual<T>& a) { using std::tan; using std::cos; return {tan(a.real), a.dual / (cos(a.real) * cos(a.real))}; }
template<typename T>
dual<T> atan(const dual<T>& a) { using std::atan; return {atan(a.real), a.dual / (T(1) + a.real * a.real)}; }
template<typename T>
dual<T> abs(const dual<T>& a) {
    using std::abs;
    // spelled out with comparisons instead of signum(), which is declared in integer_backend.h
    // and so is not visible for a floating point T. abs is not differentiable at 0.
    if (a.real < T(0))
        return {abs(a.real), -a.dual};
    if (T(0) < a.real)
        return {abs(a.real), a.dual};
    return {abs(a.real), T(0)};
}

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
