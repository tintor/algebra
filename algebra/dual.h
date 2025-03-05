#pragma once

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
dual<T> pow(const dual<T>& a, const T& y) { return {pow(a.real, y), y * a.dual * pow(a.real, y - T(1))); }
template<typename T>
dual<T> sin(const dual<T>& a) { return {sin(a.real), a.dual * cos(a.real)}; }
template<typename T>
dual<T> cos(const dual<T>& a) { return {cos(a.real), -a.dual * sin(a.real)); }
template<typename T>
dual<T> tan(const dual<T>& a) { return {tan(a.real), a.dual / (cos(a.real) * cos(a.real))); }
template<typename T>
dual<T> atan(const dual<T>& a) { return {atan(a.real), a.dual / (T(1) + a.real * a.real)); }

}
