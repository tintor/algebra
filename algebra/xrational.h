#pragma once
#include "algebra/rational.h"
#include "algebra/natural_func.h"
#include "algebra/integer_func.h"

namespace algebra {

// represents number of form: rational * sqrt(root)
// it is closed under: multiplication and division
// squaring always produces rational
//
struct xrational {
    rational base;
    natural root; // must be positive!

    xrational(rational base, natural root = 1) : base(std::move(base)), root(std::move(root)) {
        if (this->root == 0)
            throw std::runtime_error("root must be positive");
        simplify();
    }
    xrational(std::integral auto a) : base(a), root(1) { }
    xrational(integer a) : base(std::move(a)), root(1) { }

    xrational& operator=(std::integral auto a) {
        base = a;
        root = 1;
        return *this;
    }
    xrational& operator=(integer a) {
        base.num = std::move(a);
        base.den = 1;
        root = 1;
        return *this;
    }
    xrational& operator=(rational a) {
        base = std::move(a);
        root = 1;
        return *this;
    }

    void simplify() {
        natural w = 1, r = 1;
        exact_sqrt(root, w, r);
        if (w != 1)
            base *= integer(w);
        root = std::move(r);
    }
};

constexpr xrational operator+(const xrational& a, const xrational& b) {
    if (a.root == b.root)
        return {a.base + b.base, a.root};
    throw std::runtime_error("adding xrationals with different roots");
}

constexpr xrational operator-(const xrational& a, const xrational& b) {
    if (a.root == b.root)
        return {a.base - b.base, a.root};
    throw std::runtime_error("subtracting xrationals with different roots");
}

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

template<typename T> concept xrational_like = rational_like<T> || std::same_as<T, xrational>;
constexpr bool operator>(const xrational_like auto& a, const xrational_like auto& b) { return a < b; }
constexpr bool operator>=(const xrational_like auto& a, const xrational_like auto& b) { return !(a < b); }
constexpr bool operator<=(const xrational_like auto& a, const xrational_like auto& b) { return b >= a; }
constexpr bool operator!=(const xrational_like auto& a, const xrational_like auto& b) { return !(a == b); }

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

constexpr xrational& operator<<=(xrational& a, std::integral auto b) {
    a.base *= b;
    return a;
}

ALGEBRA_SHIFT_OP(xrational)

namespace literals {
constexpr auto operator""_x(const char* s) { return xrational(rational(s)); }
}

}

template <>
struct std::formatter<algebra::xrational, char> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }

    constexpr auto format(const algebra::xrational& a, auto& ctx) const {
        if (a.base == 1 && a.root != 1) {
            std::format_to(ctx.out(), "sqrt({})", a.root);
            return ctx.out();
        }
        if (a.base.num == 1 && a.base.den != 1 && a.root != 1) {
            std::format_to(ctx.out(), "sqrt({})/{}", a.root, a.base.den);
            return ctx.out();
        }
        if (a.base.num == -1 && a.base.den != 1 && a.root != 1) {
            std::format_to(ctx.out(), "-sqrt({})/{}", a.root, a.base.den);
            return ctx.out();
        }
        std::format_to(ctx.out(), "{}", a.base);
        if (a.root != 1)
            std::format_to(ctx.out(), "*sqrt({})", a.root);
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
