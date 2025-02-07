#pragma once
#include "algebra/rational.h"
#include <memory>

namespace algebra {

template<typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> c;
    c.reserve(a.size() + b.size());
    c = a;
    std::copy(b.begin(), b.end(), std::back_inserter(c));
    return c;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& a, const T& b) {
    std::vector<T> c;
    c.reserve(a.size() + 1);
    c = a;
    c.push_back(b);
    return c;
}

template<typename T>
std::vector<T> operator+(const T& a, const std::vector<T>& b) {
    std::vector<T> c;
    c.reserve(1 + b.size());
    c.push_back(a);
    std::copy(b.begin(), b.end(), std::back_inserter(c));
    return c;
}

template<typename T>
std::vector<T> subvec(const std::vector<T>& a, size_t start) {
    std::vector<T> c;
    c.reserve(a.size() - start);
    std::copy(a.begin() + start, a.end(), std::back_inserter(c));
    return c;
}

struct expr;
using expr_ptr = std::shared_ptr<expr>;

struct expr {
    virtual ~expr() {}
    virtual int sign() const = 0;
};

template<typename T>
constexpr const T* dcast(const expr_ptr& a) { return dynamic_cast<const T*>(a.get()); }
template<typename T>
constexpr const T* dcast(const expr* a) { return dynamic_cast<const T*>(a); }

struct expr_matrix : public expr {
    int rows, cols;
    std::vector<expr_ptr> data;

    virtual constexpr int sign() const {
        if (rows == 0 && cols == 0)
            return data[0]->sign();
        throw std::runtime_error("indeterminate sign");
    }
};

struct expr_var : public expr {
    std::string name;
    virtual constexpr int sign() const { throw std::runtime_error("indeterminate sign"); }
};

struct expr_rel : public expr {
    expr_ptr left;
    bool less;
    bool eq;
    bool greater;
    expr_ptr right;
    virtual constexpr int sign() const { throw std::runtime_error("unimplemented"); }
};

struct expr_pi : public expr {
    virtual constexpr int sign() const { return 1; }
};

struct expr_e : public expr {
    virtual constexpr int sign() const { return 1; }
};

struct expr_sin : public expr {
    expr_ptr value;
    constexpr expr_sin(expr_ptr value) : value(value) {}
    virtual constexpr int sign() const;
};

struct expr_cos : public expr {
    expr_ptr value;
    constexpr expr_cos(expr_ptr value) : value(value) {}
    virtual constexpr int sign() const;
};

struct expr_integer : public expr {
    integer value;
    constexpr expr_integer(integer value) : value(std::move(value)) {}
    virtual constexpr int sign() const { return value.sign(); }
};

constexpr bool is_integer(expr_ptr a) { return dcast<expr_integer>(a); }
constexpr const integer& integer_value(expr_ptr a) { return dcast<expr_integer>(a)->value; }

struct expr_rational : public expr {
    rational value;
    constexpr expr_rational(rational value) : value(std::move(value)) {}
    virtual constexpr int sign() const { return value.sign(); }
};

constexpr bool is_rational(expr_ptr a) { return dcast<expr_rational>(a) || dcast<expr_integer>(a); }
constexpr rational rational_value(expr_ptr a) {
    if (const expr_integer* i = dcast<expr_integer>(a))
        return i->value;
    return dcast<expr_rational>(a)->value;
}

struct expr_power : public expr {
    expr_ptr base;
    rational exp;
    constexpr expr_power(expr_ptr base, rational exp) : base(base), exp(std::move(exp)) {}
    virtual constexpr int sign() const {
        int b = base->sign();
        if (b > 0)
            return 1;
        if (b == 0 && exp.sign() > 0)
            return 0;
        if (b < 0 && exp.is_integer())
            return exp.num.is_even();
        throw std::runtime_error("indeterminate sign of pow(negative, non-integer)");
    }
};

constexpr bool is_power(expr_ptr a) { return dcast<expr_power>(a); }
constexpr const auto& power_base(expr_ptr a) { return dcast<expr_power>(a)->base; }
constexpr const auto& power_exp(expr_ptr a) { return dcast<expr_power>(a)->exp; }
constexpr bool is_sqrt(expr_ptr a) { using namespace algebra::literals; return is_power(a) && power_exp(a) == 1/2_q; }
constexpr bool is_cbrt(expr_ptr a) { using namespace algebra::literals; return is_power(a) && power_exp(a) == 1/3_q; }

struct expr_sum : public expr {
    std::vector<expr_ptr> values;
    mutable std::optional<int> _sign;
    constexpr expr_sum(std::vector<expr_ptr> v) : values(std::move(v)) {}
    virtual constexpr int sign() const;
};

constexpr bool is_sum(expr_ptr a) { return dcast<expr_sum>(a); }
constexpr const auto& sum_values(expr_ptr a) { return dcast<expr_sum>(a)->values; }

struct expr_negation : public expr {
    expr_ptr value;
    constexpr expr_negation(expr_ptr value) : value(value) {}
    virtual constexpr int sign() const { return -value->sign(); }
};

constexpr bool is_negation(expr_ptr a) { return dcast<expr_negation>(a); }
constexpr expr_ptr negation_value(expr_ptr a) { return dcast<expr_negation>(a)->value; }

struct expr_product : public expr {
    std::vector<expr_ptr> values;
    mutable std::optional<int> _sign;
    constexpr expr_product(std::vector<expr_ptr> v) : values(std::move(v)) {}
    virtual constexpr int sign() const {
        if (_sign != std::nullopt) return *_sign;
        int sign = 1;
        for (const auto& e: values) {
            int s = e->sign();
            if (s < 0) sign = -sign;
            if (s == 0) return 0;
        }
        _sign = sign;
        return sign;
    }
};

constexpr bool is_product(expr_ptr a) { return dcast<expr_product>(a); }
constexpr const auto& product_values(expr_ptr a) { return dcast<expr_product>(a)->values; }
constexpr bool is_product_rx(expr_ptr a) { return is_product(a) && product_values(a).size() == 2 && is_rational(product_values(a)[0]); }

// ==========

constexpr expr_ptr make_integer(const integer& a) { return std::make_shared<expr_integer>(a); }
constexpr expr_ptr make_rational(const rational& a) {
    if (a.is_integer())
        return make_integer(a.num);
    return std::make_shared<expr_rational>(a);
}

namespace literals {
constexpr auto operator""_e(const char* s) { return make_rational(rational(s)); }
}

const expr_ptr ZERO_EXPR = make_integer(0);
const expr_ptr ONE_EXPR = make_integer(1);
const expr_ptr E_EXPR = std::make_shared<expr_e>();
const expr_ptr PI_EXPR = std::make_shared<expr_pi>();

constexpr expr_ptr operator+(expr_ptr a, expr_ptr b);
constexpr expr_ptr operator-(expr_ptr a, expr_ptr b);
constexpr expr_ptr operator*(expr_ptr a, expr_ptr b);
constexpr expr_ptr operator/(expr_ptr a, expr_ptr b);

constexpr expr_ptr make_product(std::vector<expr_ptr> v);
constexpr expr_ptr pow(expr_ptr a, const rational& b);

constexpr expr_ptr operator+(std_int auto a, expr_ptr b) { return make_integer(a) + b; }
constexpr expr_ptr operator+(expr_ptr a, std_int auto b) { return a + make_integer(b); }

constexpr expr_ptr operator-(std_int auto a, expr_ptr b) { return make_integer(a) - b; }
constexpr expr_ptr operator-(expr_ptr a, std_int auto b) { return a - make_integer(b); }

constexpr expr_ptr operator*(std_int auto a, expr_ptr b) { return make_integer(a) * b; }
constexpr expr_ptr operator*(expr_ptr a, std_int auto b) { return a * make_integer(b); }

constexpr expr_ptr operator/(std_int auto a, expr_ptr b) { return make_integer(a) / b; }
constexpr expr_ptr operator/(expr_ptr a, std_int auto b) { return a / make_integer(b); }

constexpr bool identical(expr_ptr a, expr_ptr b) {
    auto ai = dcast<expr_integer>(a);
    if (ai) {
        auto bi = dcast<expr_integer>(b);
        return bi && ai->value == bi->value;
    }

    auto ar = dcast<expr_rational>(a);
    if (ar) {
        auto br = dcast<expr_rational>(b);
        return br && ar->value == br->value;
    }

    auto aw = dcast<expr_power>(a);
    if (aw) {
        auto bw = dcast<expr_power>(b);
        return bw && identical(aw->base, bw->base) && aw->exp == bw->exp;
    }

    auto as = dcast<expr_sum>(a);
    if (as) {
        auto bs = dcast<expr_sum>(b);
        if (!bs || as->values.size() != bs->values.size())
            return false;
        for (size_t i = 0; i < as->values.size(); i++)
            if (!identical(as->values[i], bs->values[i]))
                return false;
        return true;
    }

    auto ap = dcast<expr_product>(a);
    if (ap) {
        auto bp = dcast<expr_product>(b);
        if (!bp || ap->values.size() != bp->values.size())
            return false;
        for (size_t i = 0; i < ap->values.size(); i++)
            if (!identical(ap->values[i], bp->values[i]))
                return false;
        return true;
    }

    auto an = dcast<expr_negation>(a);
    if (an) {
        auto bn = dcast<expr_negation>(b);
        return bn && an->value == bn->value;
    }

    if (dcast<expr_e>(a))
        return dcast<expr_e>(b);
    if (dcast<expr_pi>(a))
        return dcast<expr_pi>(b);

    throw std::runtime_error("unreachable");
}

constexpr std::optional<int> safe_sign(expr_ptr a);

constexpr expr_ptr make_sum(std::vector<expr_ptr> v) {
    // move all rationals to the front of sum and combine them
    rational a;
    int r = 0;
    int w = 0;
    while (r < v.size()) {
        if (is_rational(v[r]))
            a += rational_value(v[r]);
        else
            v[w++] = v[r];
        r++;
    }
    v.resize(w);
    if (a.num.sign())
        v.insert(v.begin(), make_rational(a));

    // a + b + a -> 2*a + b
    for (size_t i = 0; i < v.size(); i++) {
        int count = 1;
        size_t j = v.size() - 1;
        while (j > i) {
            if (safe_sign(v[i] - v[j]) == 0) {
            //if (identical(v[i], v[j])) {
                count += 1;
                v[j] = v.back();
                v.pop_back();
            }
            j--;
        }
        if (count > 1)
            v[i] = count * v[i];
    }

    // e + b - e -> b
    for (size_t i = 0; i < v.size(); i++) {
        bool matched = false;
        for (size_t j = i + 1; j < v.size(); j++) {
            if (safe_sign(v[i] + v[j]) == 0) {
                v[j] = v.back();
                v.pop_back();
                matched = true;
                break;
            }
        }
        if (matched) {
            v[i] = v.back();
            v.pop_back();
            i--;
        }
    }

    if (v.empty())
        return ZERO_EXPR;
    if (v.size() == 1)
        return v[0];
    if (v.size() == 2)
        return v[0] + v[1];
    return make_shared<expr_sum>(std::move(v));
}

// TODO a + 2 * a -> 3 * a
constexpr expr_ptr operator+(expr_ptr a, expr_ptr b) {
    if (is_rational(a) && is_rational(b))
        return make_rational(rational_value(a) + rational_value(b));

    if (safe_sign(a) == 0)
        return b;
    if (safe_sign(b) == 0)
        return a;

    if (is_sum(a) && is_sum(b))
        return make_sum(sum_values(a) + sum_values(b));
    if (is_sum(a))
        return make_sum(sum_values(a) + b);
    if (is_sum(b))
        return make_sum(a + sum_values(b));

    if (a == b)
        return 2 * a;

    if (is_negation(a) && negation_value(a) == b)
        return ZERO_EXPR;
    if (is_negation(b) && negation_value(b) == a)
        return ZERO_EXPR;

    if (is_product_rx(a) && is_product_rx(b) && product_values(a)[1] == product_values(b)[1])
        return make_rational(rational_value(product_values(a)[0]) + rational_value(product_values(b)[0])) * product_values(a)[1];
    if (is_product_rx(a) && product_values(a)[1] == b)
        return make_rational(rational_value(product_values(a)[0]) + 1) * b;
    if (is_product_rx(b) && product_values(b)[1] == a)
        return make_rational(rational_value(product_values(b)[0]) + 1) * a;

    return std::make_shared<expr_sum>(std::vector<expr_ptr>{a, b});
}

constexpr expr_ptr operator-(expr_ptr a) {
    if (is_rational(a))
        return make_rational(-rational_value(a));
    if (is_negation(a))
        return negation_value(a);
    if (is_sum(a)) {
        std::vector<expr_ptr> v = sum_values(a);
        for (auto& e: v)
            e = -e;
        return make_sum(v);
    }
    if (is_product(a) && is_rational(product_values(a)[0]))
        return make_rational(-rational_value(product_values(a)[0])) * make_product(subvec(product_values(a), 1));
    return std::make_shared<expr_negation>(a);
}

constexpr expr_ptr make_product(std::vector<expr_ptr> v) {
    // a * 2 * b * 3 -> 6 * a * b
    rational a = 1;
    int r = 0;
    int w = 0;
    while (r < v.size()) {
        if (is_rational(v[r]))
            a *= rational_value(v[r]);
        else
            v[w++] = v[r];
        r++;
    }
    v.resize(w);
    if (a != 1)
        v.insert(v.begin(), make_rational(a));

    // a * b * a -> a^2 * b
    for (size_t i = 0; i < v.size(); i++) {
        int count = 1;
        size_t j = v.size() - 1;
        while (j > i) {
            if (identical(v[i], v[j])) {
                count += 1;
                v[j] = v.back();
                v.pop_back();
            }
            j--;
        }
        if (count > 1)
            v[i] = pow(v[i], count);
    }

    if (v.empty())
        return ONE_EXPR;
    if (v.size() == 1)
        return v[0];
    if (v.size() == 2)
        return v[0] * v[1];
    return std::make_shared<expr_product>(std::move(v));
}

constexpr expr_ptr operator*(expr_ptr a, expr_ptr b) {
    if (is_rational(a) && is_rational(b))
        return make_rational(rational_value(a) * rational_value(b));

    if (is_rational(b))
        return b * a;

    if (is_negation(a) && is_negation(b))
        return negation_value(a) * negation_value(b);
    if (is_negation(a))
        return -(negation_value(a) * b);
    if (is_negation(b))
        return -(a * negation_value(b));

    if (is_rational(a) && rational_value(a) == -1)
        return -b;
    if (safe_sign(a) == 0 || safe_sign(b) == 0)
        return ZERO_EXPR;
    if (a == b)
        return pow(a, 2);

    if (is_power(a) && is_power(b) && power_base(a) == power_base(b))
        return pow(power_base(a), power_exp(a) + power_exp(b));
    if (is_power(a) && power_base(a) == b)
        return pow(b, power_exp(a) + 1);
    if (is_power(b) && power_base(b) == a)
        return pow(a, power_exp(b) + 1);

    if (is_product(a) && is_product(b))
        return make_product(product_values(a) + product_values(b));
    if (is_product(a))
        return make_product(product_values(a) + b);
    if (is_product(b))
        return make_product(a + product_values(b));

    return std::make_shared<expr_product>(std::vector<expr_ptr>{a, b});
}

constexpr expr_ptr sin(expr_ptr a) {
    return std::make_shared<expr_sin>(a);
}

constexpr expr_ptr cos(expr_ptr a) {
    return std::make_shared<expr_cos>(a);
}

constexpr expr_ptr pow(expr_ptr a, const rational& b) {
    using namespace algebra::literals;
    if (b == 1)
        return a;
    if (b == 0)
        return make_rational(1_q);
    if (b > 0 && b.is_even() && is_negation(a))
        return pow(-a, b);
    if (b == 2 && is_sum(a) && sum_values(a).size() == 2) {
        auto p = sum_values(a)[0];
        auto q = sum_values(a)[1];
        return make_sum(std::vector<expr_ptr>{pow(p, 2), make_product(std::vector<expr_ptr>{2_e, p, q}), pow(q, 2)});
    }
    if (is_integer(a) && b == 1/2_q && !integer_value(a).is_negative()) {
        const natural v = integer_value(a).abs;
        if (is_possible_square(v)) {
            const natural m = isqrt(v);
            if (m * m == v)
                return make_integer(m);
        }
    }
    if (is_rational(a) && b.is_integer())
        return make_rational(pow(rational_value(a), b.num));
    if (is_power(a))
        return pow(power_base(a), power_exp(a) * b);
    if (is_product(a)) {
        std::vector<expr_ptr> v = product_values(a);
        for (auto& e: v)
            e = pow(e, b);
        return make_product(v);
    }
    return std::make_shared<expr_power>(a, b);
}

constexpr expr_ptr sqrt(expr_ptr a) { using namespace algebra::literals; return pow(a, 1/2_q); }
constexpr expr_ptr cbrt(expr_ptr a) { using namespace algebra::literals; return pow(a, 1/3_q); }

constexpr expr_ptr operator-(expr_ptr a, expr_ptr b) {
    if (is_rational(a) && is_rational(b))
        return make_rational(rational_value(a) - rational_value(b));
    return a + -b;
}

constexpr expr_ptr operator/(expr_ptr a, expr_ptr b) {
    using namespace algebra::literals;
    if (b->sign() == 0)
        throw std::runtime_error("division by zero");
    if (is_rational(a) && is_rational(b))
        return make_rational(rational_value(a) / rational_value(b));
    return a * pow(b, -1_q);
}

#define EXPR_CMP(OP) \
constexpr bool operator OP(expr_ptr a, expr_ptr b) { return (a - b)->sign() OP 0; } \
constexpr bool operator OP(std_int auto a, expr_ptr b) { return make_rational(a) OP b; } \
constexpr bool operator OP(expr_ptr a, std_int auto b) { return a OP make_rational(b); }

EXPR_CMP(<)
EXPR_CMP(>)
EXPR_CMP(<=)
EXPR_CMP(>=)
EXPR_CMP(==)
EXPR_CMP(!=)

constexpr bool needs_parenthesis(expr_ptr a) {
    return !is_sqrt(a) && !is_cbrt(a) && !is_integer(a) && !dcast<expr_pi>(a) && !dcast<expr_e>(a) && !dcast<expr_sin>(a) && !dcast<expr_cos>(a);
}

}

template <>
struct std::formatter<const algebra::expr*, char> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }

    void format_expression(const algebra::expr* a, auto& ctx) const {
        using namespace algebra;
        using namespace algebra::literals;
        if (auto integer = dcast<expr_integer>(a)) {
            std::format_to(ctx.out(), "{}", integer->value);
        } else if (auto rational = dcast<expr_rational>(a)) {
            std::format_to(ctx.out(), "{}", rational->value);
        } else if (auto power = dcast<expr_power>(a)) {
            if (power->exp == 1/2_q) {
                std::format_to(ctx.out(), "sqrt(");
                format_expression(power->base.get(), ctx);
                std::format_to(ctx.out(), ")");
            } else if (power->exp == 1/3_q) {
                std::format_to(ctx.out(), "cbrt(");
                format_expression(power->base.get(), ctx);
                std::format_to(ctx.out(), ")");
            } else {
                if (!needs_parenthesis(power->base)) {
                    format_expression(power->base.get(), ctx);
                } else {
                    std::format_to(ctx.out(), "(");
                    format_expression(power->base.get(), ctx);
                    std::format_to(ctx.out(), ")");
                }
                if (power->exp.is_integer())
                    std::format_to(ctx.out(), "^{}", power->exp);
                else
                    std::format_to(ctx.out(), "^({})", power->exp);
            }
        } else if (auto neg = dcast<expr_negation>(a)) {
            auto v = neg->value;
            if (!needs_parenthesis(v)) {
                std::format_to(ctx.out(), "-");
                format_expression(v.get(), ctx);
            } else {
                std::format_to(ctx.out(), "-(");
                format_expression(v.get(), ctx);
                std::format_to(ctx.out(), ")");
            }
        } else if (auto sum = dcast<expr_sum>(a)) {
            for (int i = 0; i < sum->values.size(); i++) {
                auto b = sum->values[i];
                if (is_negation(b)) {
                    if (i > 0)
                        std::format_to(ctx.out(), " - ");
                    else
                        std::format_to(ctx.out(), "-");
                    format_expression(negation_value(b).get(), ctx);
                } else if (is_rational(b) && b->sign() < 0) {
                    if (i > 0)
                        std::format_to(ctx.out(), " - {}", abs(rational_value(b)));
                    else
                        std::format_to(ctx.out(), "{}", rational_value(b));
                } else {
                    if (i > 0)
                        std::format_to(ctx.out(), " + ");
                    format_expression(b.get(), ctx);
                }
            }
        } else if (auto product = dcast<expr_product>(a)) {
            for (int i = 0; i < product->values.size(); i++) {
                if (i > 0)
                    std::format_to(ctx.out(), "*");
                auto b = product->values[i];
                if (!needs_parenthesis(b))
                    format_expression(b.get(), ctx);
                else {
                    std::format_to(ctx.out(), "(");
                    format_expression(b.get(), ctx);
                    std::format_to(ctx.out(), ")");
                }
            }
        } else if (dcast<expr_pi>(a)) {
            std::format_to(ctx.out(), "π");
        } else if (dcast<expr_e>(a)) {
            std::format_to(ctx.out(), "e");
        } else if (auto c = dcast<expr_cos>(a)) {
            std::format_to(ctx.out(), "cos(");
            format_expression(c->value.get(), ctx);
            std::format_to(ctx.out(), ")");
        } else if (auto s = dcast<expr_sin>(a)) {
            std::format_to(ctx.out(), "sin(");
            format_expression(s->value.get(), ctx);
            std::format_to(ctx.out(), ")");
        } else {
            std::format_to(ctx.out(), "###");
        }
    }

    auto format(const algebra::expr* a, auto& ctx) const {
        if (!a)
            return std::format_to(ctx.out(), "null");
        format_expression(a, ctx);
        return ctx.out();
    }
};

constexpr std::ostream& operator<<(std::ostream& os, const algebra::expr* a) { return os << std::format("{}", a); }

template <>
struct std::formatter<algebra::expr_ptr, char> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(algebra::expr_ptr a, auto& ctx) const {
        return std::format_to(ctx.out(), "{}", const_cast<const algebra::expr*>(a.get()));
    }
};

constexpr std::ostream& operator<<(std::ostream& os, algebra::expr_ptr a) { return os << std::format("{}", a); }

namespace algebra {

class unknown_sign_error : public std::runtime_error {
public:
    unknown_sign_error(std::string m) : std::runtime_error(m) {}
};

template<typename T>
struct interval {
    T min, max;
};

}

template <typename T>
struct std::formatter<algebra::interval<T>, char> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }

    auto format(const algebra::interval<T>& a, auto& ctx) const {
        return std::format_to(ctx.out(), "[{}, {}]", a.min, a.max);
    }
};

namespace algebra {

template<typename T>
constexpr interval<T> operator-(const interval<T>& a) {
    return interval{-a.max, -a.min};
}

template<typename T>
constexpr interval<T> operator+(const interval<T>& a, const interval<T>& b) {
    return interval{a.min + b.min, a.max + b.max};
}

template<typename T>
constexpr interval<T> operator*(const interval<T>& a, const interval<T>& b) {
    auto p = a.min * b.min;
    auto q = a.min * b.max;
    if (p > q) std::swap(p, q);
    auto r = a.max * b.min;
    auto s = a.max * b.max;
    if (r > s) std::swap(r, s);

    return interval{(p < r) ? p : r, (q > s) ? q : s};
}

template<typename T>
constexpr interval<T>& operator+=(interval<T>& a, const interval<T>& b) {
    a = a + b;
    return a;
}

template<typename T>
constexpr interval<T>& operator*=(interval<T>& a, const interval<T>& b) {
    a = a * b;
    return a;
}

template<typename T>
constexpr interval<T> pow(const interval<T>& a, const integer& b) {
    interval<T> e{1, 1};
    for (int i = 0; i < abs(b); i++)
        e *= a;
    if (b < 0) {
        if (e.max >= 0 && e.min <= 0)
            throw std::runtime_error("division by zero");
        e = interval<T>{T(1) / e.max, T(1) / e.min};
    }
    return e;
}

// TODO make this a virtual function
constexpr std::optional<interval<rational>> bounds(expr_ptr a) {
    using namespace algebra::literals;
    if (is_rational(a))
        return interval<rational>{rational_value(a), rational_value(a)};
    if (dcast<expr_e>(a))
        return interval<rational>{2.7_q, 2.8_q};
    if (dcast<expr_pi>(a))
        return interval<rational>{3.1_q, 3.2_q};
    if (dcast<expr_sin>(a) || dcast<expr_cos>(a))
        return interval<rational>{-1, 1};
    if (is_sum(a)) {
        interval<rational> b{0, 0};
        for (expr_ptr e: sum_values(a)) {
            auto bo = bounds(e);
            if (bo == std::nullopt)
                return std::nullopt;
            b += *bo;
        }
        return b;
    }
    if (is_product(a)) {
        interval<rational> b{1, 1};
        for (expr_ptr e: product_values(a)) {
            auto bo = bounds(e);
            if (bo == std::nullopt)
                return std::nullopt;
            b *= *bo;
        }
        return b;
    }
    if (is_negation(a)) {
        auto bo = bounds(negation_value(a));
        if (bo == std::nullopt)
            return std::nullopt;
        return -*bo;
    }
    if (is_power(a) && power_exp(a).is_integer()) {
        auto bo = bounds(power_base(a));
        if (bo == std::nullopt)
            return std::nullopt;
        auto b = *bo;
        const integer& e = power_exp(a).num;
        if (e < 0 && b.max > 0 && b.min < 0)
            return std::nullopt;
        return pow(b, e);
    }
    return std::nullopt;
}

constexpr std::optional<int> safe_sign(expr_ptr a) {
    try {
        return a->sign();
    } catch(unknown_sign_error) {
        return std::nullopt;
    }
}

constexpr int expr_sin::sign() const {
    throw unknown_sign_error(std::format("unknown sign of {}", static_cast<const expr*>(this)));
}

constexpr int expr_cos::sign() const {
    throw unknown_sign_error(std::format("unknown sign of {}", static_cast<const expr*>(this)));
}

constexpr int expr_sum::sign() const {
    std::vector<expr_ptr> positive, negative;
    bool sign_unknown = false;
    for (const auto& e: values) {
        auto s = safe_sign(e);
        if (s == std::nullopt) {
            sign_unknown = true;
            break;
        }
        if (*s > 0)
            positive.push_back(e);
        if (*s < 0)
            negative.push_back(e);
    }
    if (sign_unknown) {
        auto bo = bounds(values[0]);
        if (bo == std::nullopt)
            throw unknown_sign_error(std::format("unknown sign and bounds of {}", values[0]));

        interval<rational> b = *bo;
        for (int i = 1; i < values.size(); i++) {
            bo = bounds(values[i]);
            if (bo == std::nullopt)
                throw unknown_sign_error(std::format("unknown sign and bounds of {}", values[i]));
            b += *bo;
        }

        if (b.min > 0)
            return 1;
        if (b.max < 0)
            return -1;
        throw unknown_sign_error(std::format("unknown sign of {}", static_cast<const expr*>(this)));
    }

    if (positive.size() > 0 && negative.size() == 0) return 1;
    if (positive.size() == 0 && negative.size() == 0) return 0;
    if (positive.size() == 0 && negative.size() > 0) return -1;

    auto p = make_sum(positive);
    auto n = make_sum(negative);

    auto p_bounds = bounds(p);
    auto n_bounds = bounds(n);
    if (p_bounds != std::nullopt && n_bounds != std::nullopt) {
        interval<rational> e = *p_bounds + *n_bounds;
        if (e.min > 0)
            return 1;
        if (e.max < 0)
            return -1;
    }

    auto a = pow(p, 2);
    auto b = pow(n, 2);
    if (is_rational(a) && is_rational(b)) {
        if (rational_value(a) == rational_value(b)) return 0;
        if (rational_value(a) > rational_value(b)) return 1;
        return -1;
    }
    if ((is_power(a) && power_exp(a) > 1) || (is_power(b) && power_exp(b) > 1))
        throw unknown_sign_error(std::format("unknown sign of {}", static_cast<const expr*>(this)));
    return (a - b)->sign();
}

}
