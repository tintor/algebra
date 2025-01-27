#pragma once
#include "algebra/rational.h"
#include "algebra/natural_func.h"
#include "algebra/integer_func.h"

namespace algebra {

struct xrational;
template<> struct IsNumberClass<xrational> : std::true_type {};
template<typename T> concept xrational_like = rational_like<T> || std::same_as<T, xrational>;

// represents number of form: rational * sqrt(root)
// it is closed under: multiplication and division
// squaring always produces rational
//
struct xrational {
    rational base;
    natural root; // must be positive!

    xrational() : base(0), root(1) { }
    xrational(const xrational& a) : base(a.base), root(a.root) { }
    xrational(rational_like auto base) : base(std::move(base)), root(1) { }
    xrational(rational base, natural root) : base(std::move(base)), root(std::move(root)) {
        if (this->root == 0)
            throw std::runtime_error("root must be positive");
        simplify();
    }

    xrational& operator=(const xrational& a) {
        base = a.base;
        root = a.root;
        return *this;
    }

    xrational& operator=(const rational_like auto& a) {
        base = a;
        root = 1;
        return *this;
    }

    void simplify() {
        if (base.num.sign() == 0) {
            root = 1;
            return;
        }
        natural w = 1, r = 1;
        exact_sqrt(root, w, r);
        if (w != 1)
            base *= integer(w);
        root = std::move(r);
    }

    void negate() { base.negate(); }

    void invert() {
        base.invert();
        if (root != 1)
            base /= root;
    }
};

constexpr void negate(xrational& a) { a.negate(); }

constexpr xrational operator-(const xrational& a) { return {-a.base, a.root}; }

constexpr void __addition_compatible(const xrational& a, const xrational& b) {
    if (a.root != b.root)
        throw std::runtime_error("adding xrationals with different roots");
}

constexpr xrational operator+(const xrational& a, const xrational& b) {
    __addition_compatible(a, b);
    return {a.base + b.base, a.root};
}

constexpr xrational operator+(const xrational& a, const rational_like auto& b) { return {a.base + b, a.root}; }

constexpr xrational& operator+=(xrational& a, const xrational& b) {
    __addition_compatible(a, b);
    a.base += b.base;
    return a;
}

constexpr xrational& operator+=(xrational& a, const rational_like auto& b) { a.base += b; return a; }

constexpr xrational operator-(const xrational& a, const xrational& b) {
    __addition_compatible(a, b);
    return {a.base - b.base, a.root};
}

constexpr xrational operator-(const xrational& a, const rational_like auto& b) { return {a.base - b, a.root}; }

constexpr xrational& operator-=(xrational& a, const xrational& b) {
    __addition_compatible(a, b);
    a.base -= b.base;
    return a;
}

constexpr xrational& operator-=(xrational& a, const rational_like auto& b) { a.base -= b; return a; }

constexpr xrational sqr(const xrational& a) {
    natural num, den;
    mul(a.base.num.abs, a.base.num.abs, num);
    mul(a.base.den.abs, a.base.den.abs, den);
    mul(num, a.root);
    return {rational{num, den}};
}

constexpr xrational operator*(const xrational& a, const xrational& b) {
    if (&a == &b || (a.base == b.base && a.root == b.root))
        return sqr(a);
    if (a.root == b.root)
        return {a.base * b.base * integer(a.root)};

    // TODO maybe use gcd(a.root, b.root) to avoid difficult factorization
    natural whole = a.base.num.abs * b.base.num.abs;
    natural root = 1;
    exact_sqrt(a.root, whole, root);
    exact_sqrt(b.root, whole, root);
    integer w = whole;
    if ((a.base.sign() < 0) != (b.base.sign() < 0))
        w.negate();
    return {rational{std::move(w), a.base.den.abs * b.base.den.abs}, std::move(root)};
}

constexpr xrational operator*(const xrational& a, const rational_like auto& b) { return {a.base * b, a.root}; }

constexpr xrational& operator*=(xrational& a, const xrational& b) {
    a = a * b;
    return a;
}

constexpr xrational& operator*=(xrational& a, const rational_like auto& b) {
    a.base *= b;
    return a;
}

constexpr xrational operator/(const xrational& a, const xrational& b) {
    if (&a == &b || (a.base == b.base && a.root == b.root))
        return xrational{rational{1}};
    if (a.root == b.root)
        return xrational{a.base / b.base};

    rational e = a.base;
    e /= b.base;
    e /= integer(b.root);
    return {std::move(e), a.root * b.root};
}

constexpr xrational operator/(const xrational& a, const rational_like auto& b) { return {a.base / b, a.root}; }

constexpr xrational& operator/=(xrational& a, const xrational& b) {
    a = a / b;
    return a;
}

constexpr xrational& operator/=(xrational& a, const rational_like auto& b) {
    a.base /= b;
    return b;
}

struct bit_range {
    uint64_t min, max;
    constexpr bit_range(uint64_t min, uint64_t max) : min(min), max(max) { }
    constexpr bit_range(const natural& a) { min = max = a.num_bits(); }
};

constexpr bit_range operator*(bit_range a, bit_range b) {
    return {a.min + b.min - 1, a.max + b.max};
}

constexpr bit_range operator+(bit_range a, bit_range b) {
    return {std::min(a.min, b.min), std::max(a.max, b.max) + 1};
}

constexpr bool operator==(const xrational& a, const xrational& b) {
    if (signum(a.base.num) != signum(b.base.num))
        return false;
    if (a.root == b.root)
        return a.base == b.base;

    natural p = a.base.num.abs * a.base.num.abs;
    p *= a.root;
    p *= b.base.den.abs;
    p *= b.base.den.abs;
    natural q = b.base.num.abs * b.base.num.abs;
    q *= b.root;
    q *= a.base.den.abs;
    q *= a.base.den.abs;
    return p == q;
}

constexpr bool operator==(const xrational& a, const rational_like auto& b) { return a.root == 1 && a.base == b; }
constexpr bool operator==(const rational_like auto& a, const xrational& b) { return b.root == 1 && a == b.base; }

constexpr bool operator<(const xrational& a, const xrational& b) {
    if (signum(a.base.num) != signum(b.base.num))
        return signum(a.base.num) < signum(b.base.num);
    if (a.root == b.root)
        return a.base < b.base;
    if (a.base == b.base)
        return a.root < b.root;

    // TODO could use interval arithmetic with doubles first to find bounds for p and q
    natural p = a.base.num.abs * a.base.num.abs;
    p *= a.root;
    p *= b.base.den.abs;
    p *= b.base.den.abs;
    natural q = b.base.num.abs * b.base.num.abs;
    q *= b.root;
    q *= a.base.den.abs;
    q *= a.base.den.abs;
    return p < q;
}

constexpr bool operator<(const xrational& a, const rational_like auto& b) { return a < xrational(b); }
constexpr bool operator<(const rational_like auto& a, const xrational& b) { return xrational(a) < b; }

constexpr xrational sqrt(const xrational& a) {
    if (a.base.sign() < 0)
        throw std::runtime_error("sqrt() of negative");
    if (a.root != 1)
        throw std::runtime_error("sqrt of xrational with root");
    natural whole = 1, root = 1;
    exact_sqrt(a.base.num.abs, whole, root);
    exact_sqrt(a.base.den.abs, whole, root);
    return {rational{whole, a.base.den}, root};
}

constexpr xrational& operator<<=(xrational& a, std_int auto b) { a.base <<= b; return a; }
ALGEBRA_SHIFT_OP(xrational)

namespace literals {
constexpr auto operator""_x(const char* s) { return xrational(rational(s)); }
}

}

template <>
struct std::formatter<algebra::xrational, char> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }

    constexpr auto format(const algebra::xrational& a, auto& ctx) const {
        if (a.root != 1) {
            if (a.base == 1) {
                std::format_to(ctx.out(), "sqrt({})", a.root);
                return ctx.out();
            }
            if (a.base.num == 1 && a.base.den != 1) {
                std::format_to(ctx.out(), "sqrt({})/{}", a.root, a.base.den);
                return ctx.out();
            }
            if (a.base.num == -1 && a.base.den != 1) {
                std::format_to(ctx.out(), "-sqrt({})/{}", a.root, a.base.den);
                return ctx.out();
            }
            std::format_to(ctx.out(), "{}*sqrt({})", a.base, a.root);
            return ctx.out();
        }
        std::format_to(ctx.out(), "{}", a.base);
        return ctx.out();
    }
};

constexpr std::ostream& operator<<(std::ostream& os, const algebra::xrational& a) { return os << std::format("{}", a); }

template<>
struct std::hash<algebra::xrational> {
    constexpr size_t operator()(const algebra::xrational& a) const {
        uint64_t seed = 0;
        seed = algebra::integer_backend::hash_fn_64bit(seed ^ std::hash<algebra::integer>()(a.base.num));
        seed = algebra::integer_backend::hash_fn_64bit(seed ^ std::hash<algebra::integer>()(a.base.den));
        seed = algebra::integer_backend::hash_fn_64bit(seed ^ std::hash<algebra::natural>()(a.root));
        return seed;
    }
};
