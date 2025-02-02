#pragma once
#include "algebra/natural.h"
#include "algebra/rational.h"
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
    natural root; // must be positive! Note: root is not fully simplified! It might have some square factors! Full factorization is expensive.

    constexpr xrational() : base(0), root(1) { }
    constexpr xrational(const xrational& a) : base(a.base), root(a.root) { }
    constexpr xrational(rational_like auto base) : base(std::move(base)), root(1) { }
    constexpr xrational(rational base, natural root) : base(std::move(base)), root(std::move(root)) {
        if (this->root == 0)
            throw std::runtime_error("root must be positive");
        simplify();
    }

    constexpr xrational& operator=(const xrational& a) {
        base = a.base;
        root = a.root;
        return *this;
    }

    constexpr xrational& operator=(const rational_like auto& a) {
        base = a;
        root = 1;
        return *this;
    }

    constexpr void simplify() {
        if (base.num.sign() == 0) {
            root = 1;
            return;
        }

        // Remove all squares of 2 from root
        bool simplify = false;
        auto z = root.num_trailing_zeros();
        if (z > 1) {
            root >>= z;
            base.num <<= z / 2;
            simplify = true;
        }

        // Remove all squares of 3 from root
        // first test is with mod9() as div(root, 9, q) will allocate memory for q!
        natural q;
        if (root > 1 && root.mod9() == 0) {
            do {
                if (div(root, 9, root))
                    break;
                std::swap(root, q);
                if (base.den.mod3() == 0)
                    base.den /= 3;
                else
                    base.num *= 3;
            } while (root > 1);
        }

        if (root > 1 && exact_sqrt(root, q)) {
            base.num *= q;
            root = 1;
            simplify = true;
        }

        if (simplify)
            base.simplify();
    }

    constexpr bool is_negative() const { return base.is_negative(); }
    constexpr void negate() { base.negate(); }

    constexpr void invert() {
        base.invert();
        if (root != 1)
            base /= root;
    }
};

constexpr void negate(xrational& a) { a.negate(); }
constexpr auto signum(const xrational& a) { return signum(a.base); }

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
    integer num, den;
    mul(a.base.num, a.base.num, num);
    mul(a.base.den, a.base.den, den);
    mul(num, a.root);
    return {rational{num, den}};
}

constexpr xrational operator*(const xrational& a, const xrational& b) {
    if (&a == &b || (a.base == b.base && a.root == b.root))
        return sqr(a);
    if (a.root == b.root)
        return {a.base * b.base * integer(a.root)};

    // TODO maybe use gcd(a.root, b.root) to avoid difficult factorization
    natural whole = abs(a.base.num * b.base.num).to_natural();
    natural root = 1;
    exact_sqrt(a.root, whole, root);
    exact_sqrt(b.root, whole, root);
    integer w = whole;
    if ((a.base.sign() < 0) != (b.base.sign() < 0))
        w.negate();
    return {rational{std::move(w), a.base.den * b.base.den}, std::move(root)};
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

constexpr bool operator==(const xrational& a, const xrational& b) {
    if (signum(a) != signum(b))
        return false;
    if (a.root == b.root)
        return a.base == b.base;

    natural e = gcd(a.root, b.root);
    if (e == 1)
        return false;

    natural as, bs, ae, be;

    div(a.root, e, ae, /*dummy*/as);
    if (!__exact_sqrt1(ae, as))
        return false;
    div(b.root, e, be, /*dummy*/bs);
    if (!exact_sqrt(be, bs) || !__exact_sqrt2(ae, as))
        return false;

    bit_range ra = bit_range(as.num_bits()) * bit_range(b.base.den.num_bits()) * bit_range(a.base.num.num_bits());
    bit_range rb = bit_range(bs.num_bits()) * bit_range(a.base.den.num_bits()) * bit_range(b.base.num.num_bits());
    if (ra < rb || rb < ra)
        return false;

    as.words.reserve_bits(ra.max);
    bs.words.reserve_bits(rb.max);

    as *= b.base.den.abs;
    as *= a.base.num.abs;
    bs *= a.base.den.abs;
    bs *= a.base.num.abs;
    return as == bs;
}

constexpr bool operator==(const xrational& a, const rational_like auto& b) { return a.root == 1 && a.base == b; }
constexpr bool operator==(const rational_like auto& a, const xrational& b) { return b.root == 1 && a == b.base; }

constexpr bool __less_abs(const xrational& a, const xrational& b) {
    if (a.root == b.root)
        return a.base < b.base;
    if (a.base == b.base)
        return a.root < b.root;

    natural aa, bb;
    aa.words.reserve_bits((a.base.num.num_bits() + b.base.den.num_bits()) * 2 + a.root.num_bits());
    bb.words.reserve_bits((b.base.num.num_bits() + a.base.den.num_bits()) * 2 + b.root.num_bits());

    mul(a.base.num.abs, b.base.den.abs, aa);
    mul(b.base.num.abs, a.base.den.abs, bb);

    aa *= aa;
    aa *= a.root;

    bb *= bb;
    bb *= b.root;
    return aa < bb;
}

constexpr bool operator<(const xrational& a, const xrational& b) {
    if (signum(a) != signum(b))
        return signum(a) < signum(b);
    return (a > 0) ? __less_abs(a, b) : __less_abs(b, a);
}

constexpr bool operator<(const xrational& a, const rational_like auto& b) {
    // TODO optimize this allocation
    return a.root.is_one() ? a.base < b : a < xrational(b);
}
constexpr bool operator<(const rational_like auto& a, const xrational& b) {
    // TODO optimize this allocation
    return b.root.is_one() ? a < b.base : xrational(a) < b;
}

constexpr xrational sqrt(const xrational& a) {
    if (a.base.sign() < 0)
        throw std::runtime_error("sqrt() of negative");
    if (a.root != 1)
        throw std::runtime_error("sqrt of xrational with root");
    natural whole = 1, root = 1;
    exact_sqrt(a.base.num.to_natural(), whole, root);
    exact_sqrt(a.base.den.to_natural(), whole, root);
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

namespace algebra {

constexpr xrational abs(xrational a) {
    if (a.base.is_negative())
        a.negate();
    return a;
}

constexpr xrational pow(const xrational& base, integer exp) {
    if (exp == 0)
        return 1;

    const bool invert = exp < 0;
    xrational result = 1;
    xrational _base = base;
    while (exp) {
        if (exp % 2)
            result *= _base;
        _base *= _base;
        exp >>= 1;
    }
    if (invert)
        result.invert();
    return result;
}

}
