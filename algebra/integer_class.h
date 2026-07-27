#pragma once
#include "algebra/natural.h"

namespace algebra {

struct integer;
template<> struct IsNumberClass<integer> : std::true_type {};

struct integer {
    // The whole value: the magnitude in the words, the sign in the sign of the size, which is how
    // integer_backend already stores it. Magnitude arithmetic goes through the cnatural/vnatural/
    // inatural conversions below, or through __magnitude() when it has to happen in place.
    integer_backend words;

    constexpr integer() {}
    // integer_backend's constructor already puts the sign in the sign of the size
    constexpr integer(std_int auto a) : words(a) { }
    constexpr integer(integer&& o) : words(std::move(o.words)) { }
    constexpr integer(natural&& o) : words(std::move(o.words)) { words.set_negative(false); }
    constexpr integer(const integer& o) : words(o.words) { }
    constexpr integer(const natural& o) : words(o.words) { words.set_negative(false); }

    constexpr void operator=(std_int auto a) { words = a; }
    constexpr void operator=(integer&& o) { words = std::move(o.words); }
    constexpr void operator=(natural&& o) { words = std::move(o.words); words.set_negative(false); }
    constexpr void operator=(const integer& o) { words = o.words; }
    constexpr void operator=(const natural& o) { words = o.words; words.set_negative(false); }

    // The low word, which is only meaningful when the value is not zero.
    constexpr uint64_t low_word() const { return words.size() ? words[0] : 0; }

    constexpr auto sign() const { return words.sign(); }
    constexpr bool is_negative() const { return sign() < 0; }
    constexpr bool is_even() const { return (low_word() % 2) == 0; }
    constexpr bool is_odd() const { return low_word() % 2; }
    constexpr bool is_one() const { return words.size() == 1 && words[0] == 1 && sign() >= 0; }
    constexpr bool is_zero() const { return words.size() == 0; }

    // TODO remove this and deprecate natural
    constexpr natural to_natural() const { natural m; m.words = words; m.words.set_negative(false); return m; }

    constexpr bool is_int8() const {
        if (words.size() > 1)
            return false;
        if (words[0] <= 127)
            return true;
        return sign() < 0 && words[0] == 128;
    }

    constexpr bool is_int16() const {
        if (words.size() > 1)
            return false;
        if (words[0] <= INT16_MAX)
            return true;
        return sign() < 0 && words[0] == static_cast<uint64_t>(INT16_MAX) + 1;
    }

    constexpr bool is_int32() const {
        if (words.size() > 1)
            return false;
        if (words[0] <= INT32_MAX)
            return true;
        return sign() < 0 && words[0] == static_cast<uint64_t>(INT32_MAX) + 1;
    }

    constexpr bool is_int64() const {
        if (words.size() > 1)
            return false;
        if (words[0] <= INT64_MAX)
            return true;
        return sign() < 0 && words[0] == static_cast<uint64_t>(INT64_MAX) + 1;
    }

    constexpr bool is_int128() const {
        if (words.size() > 2)
            return false;
        if (words.size() < 2)
            return true;

        uint64_t w = words[1];
        if ((w & (uint64_t(1) << 63)) == 0)
            return true;
        // only INT128_MIN itself reaches past the positive range, and its magnitude is exactly
        // 2**127, so the low word has to be zero as well
        return sign() < 0 && w == uint64_t(1) << 63 && words[0] == 0;
    }

    constexpr bool is_uint8() const { return sign() >= 0 && words.size() <= 1 && low_word() <= UINT8_MAX; }
    constexpr bool is_uint16() const { return sign() >= 0 && words.size() <= 1 && low_word() <= UINT16_MAX; }
    constexpr bool is_uint32() const { return sign() >= 0 && words.size() <= 1 && low_word() <= UINT32_MAX; }
    constexpr bool is_uint64() const { return sign() >= 0 && words.size() <= 1; }
    constexpr bool is_uint128() const { return sign() >= 0 && words.size() <= 2; }

    constexpr integer(std::string_view s, unsigned base = 10)
            : words(natural((s.size() && s[0] == '-') ? s.substr(1) : s, base).words) {
        if (s.size() && s[0] == '-')
            words.negate();
    }
    constexpr integer(const char* s, unsigned base = 10) : integer(std::string_view(s), base) {}

    constexpr operator char() const {
        Check(is_int8(), "integer -> int8 overflow");
        return is_negative() ? -static_cast<int>(words[0]) : words[0];
    }
    constexpr operator uint8_t() const {
        Check(is_uint8(), "integer -> uint8 overflow");
        return words[0];
    }
    constexpr operator short() const {
        Check(is_int16(), "integer -> int16 overflow");
        return is_negative() ? -static_cast<int>(words[0]) : words[0];
    }
    constexpr operator uint16_t() const {
        Check(is_uint16(), "integer -> uint16 overflow");
        return words[0];
    }
    constexpr operator int() const {
        Check(is_int32(), "integer -> int32 overflow");
        // negate in unsigned: for INT_MIN the magnitude is 2**31, which does not fit in an int,
        // so negating after the cast is signed overflow. Converting an out of range unsigned
        // value to a signed type is well defined (modular) since C++20, same as operator int128_t.
        if (is_negative())
            return static_cast<int>(-static_cast<uint32_t>(words[0]));
        return static_cast<int>(words[0]);
    }
    constexpr operator unsigned() const {
        Check(is_uint32(), "integer -> uint32 overflow");
        return words[0];
    }
    constexpr operator long() const {
        Check(is_int64(), "integer -> int64 overflow");
        // see operator int() for why the negation is done in unsigned
        if (is_negative())
            return static_cast<long>(-words[0]);
        return static_cast<long>(words[0]);
    }
    constexpr operator unsigned long() const {
        Check(is_uint64(), "integer -> uint64 overflow");
        return words[0]; // is_uint64() already rejected a negative value
    }
    constexpr operator unsigned long long() const {
        Check(is_uint64(), "integer -> uint64 overflow");
        return words[0]; // is_uint64() already rejected a negative value
    }
    static_assert(sizeof(long) == 8);
    static_assert(sizeof(long long) == 8);
    constexpr operator int128_t() const {
        Check(is_int128(), "integer -> int128 overflow");
        if (sign() == 2) return (uint128_t(words[1]) << 64) | words[0];
        if (sign() == 1) return uint128_t(words[0]);
        if (sign() == -1) return -uint128_t(words[0]);
        if (sign() == -2) return -((uint128_t(words[1]) << 64) | words[0]);
        return 0;
    }
    constexpr operator uint128_t() const {
        Check(is_uint128(), "integer -> unt128 overflow");
        if (sign() == 2) return (uint128_t(words[1]) << 64) | words[0];
        if (sign() == 1) return words[0];
        return 0;
    }

    constexpr int str_size_upper_bound(unsigned base = 10) const { return is_negative() + to_natural().str_size_upper_bound(base); }
    constexpr int str(char* buffer, int buffer_size, unsigned base = 10, bool upper = true) const {
        int result = 0;
        if (is_negative()) {
            if (buffer_size <= 0)
                throw std::runtime_error("buffer too small");
            buffer_size -= 1;
            *buffer++ = '-';
            result = 1;
        }
        return result + to_natural().str(buffer, buffer_size, base, upper);
    }

    constexpr std::string str(unsigned base = 10, bool upper = true) const {
        std::string s;
        s.resize(str_size_upper_bound(base));
        s.resize(str(s.data(), s.size(), base, upper));
        return s;
    }
    constexpr std::string hex() const { return str(16); }

    constexpr void negate() { words.negate(); }

    constexpr size_t popcount() const {
        if (!is_negative())
            return to_natural().popcount();

        size_t c = 0;
        uint64_t carry = 1;
        for (int i = 0; i < words.size(); i++) {
            uint128_t w = (uint128_t)words[i] + carry;
            carry = w >> 64;
            c += std::popcount(~static_cast<uint64_t>(w));
        }
        return c;
    }

    constexpr int size_of() const { return words.size() * 8; }

    constexpr operator bool() const { return sign(); }

    constexpr operator cnatural() const { return {words.data(), words.size()}; }
    constexpr operator vnatural() { return {{words.data(), words.size()}, words.capacity()}; }
    constexpr operator inatural() { return {words.data(), words.size()}; }

    // Borrows the magnitude as a natural for work that has to happen in place. natural is exactly
    // one integer_backend, so the swap is a pointer exchange and not a copy; the sign is restored
    // when the scope ends, and a result of zero comes back positive.
    class magnitude_ref {
        integer_backend& _words;
        bool _negative;
    public:
        natural value;
        constexpr magnitude_ref(integer_backend& w) : _words(w), _negative(w.sign() < 0) {
            value.words.swap(w);
            value.words.set_negative(false);
        }
        // Not copyable or movable: two of these over one backend would each restore in turn, and
        // the second would put the pre-operation value back, losing the result silently. magnitude()
        // returns a prvalue, so copy elision means nothing here needs to be copied or moved.
        magnitude_ref(const magnitude_ref&) = delete;
        magnitude_ref(magnitude_ref&&) = delete;
        magnitude_ref& operator=(const magnitude_ref&) = delete;
        magnitude_ref& operator=(magnitude_ref&&) = delete;
        constexpr ~magnitude_ref() {
            _words.swap(value.words);
            _words.set_negative(_negative && _words.size() != 0);
        }
        constexpr natural* operator->() { return &value; }
        constexpr natural& operator*() { return value; }
    };
    constexpr magnitude_ref magnitude() { return magnitude_ref(words); }

    constexpr integer& operator++() {
        // natural::operator++/-- work on the magnitude and preserve the sign of words
        if (is_negative())
            --*magnitude();
        else
            ++*magnitude();
        return *this;
    }

    constexpr integer& operator--() {
        if (is_negative())
            ++*magnitude();
        else if (words.size() > 0)
            --*magnitude();
        else
            *this = -1;
        return *this;
    }

    constexpr integer operator++(int) { integer a = *this; operator++(); return a; }
    constexpr integer operator--(int) { integer a = *this; operator--(); return a; }

    constexpr bool bit(int64_t i) const {
        const size_t w = i / 64;
        return w < static_cast<size_t>(words.size()) && (words[w] & (uint64_t(1) << (i % 64)));
    }
    constexpr auto num_bits() const { return algebra::num_bits(static_cast<cnatural>(*this)); }
    constexpr auto num_trailing_zeros() const { return algebra::num_trailing_zeros(static_cast<cnatural>(*this)); }

    template<std::floating_point T>
    constexpr operator T() const {
        auto a = static_cast<T>(to_natural());
        return (sign() < 0) ? -a : a;
    }

    constexpr void swap(integer& o) { words.swap(o.words); }

    constexpr uint64_t mod2() const {
        return low_word() % 2;
    }

    constexpr uint64_t mod3() const {
        uint64_t m = algebra::mod3(static_cast<cnatural>(*this));
        if (is_negative())
            m = (m * 2) % 3;
        return m;
    }

    constexpr uint64_t mod4() const {
        uint64_t m = low_word() % 4;
        if (is_negative())
            m = (m * 3) % 4;
        return m;
    }

    constexpr uint64_t mod5() const {
        uint64_t m = algebra::mod5(static_cast<cnatural>(*this));
        if (is_negative())
            m = (m * 4) % 5;
        return m;
    }
};

constexpr bool operator==(const integer& a, const integer& b) {
    return a.is_negative() == b.is_negative() && __equal(static_cast<cnatural>(a), static_cast<cnatural>(b));
}

constexpr bool operator==(const integer& a, std_int auto b) {
    if (b < 0)
        return a.is_negative() && __equal_u(static_cast<cnatural>(a), abs_unsigned(b));
    return !a.is_negative() && __equal_u(static_cast<cnatural>(a), make_unsigned(b));
}

constexpr void negate(integer& a) { a.negate(); }
constexpr integer operator-(integer a) { a.negate(); return a; }

template<bool plus>
constexpr integer& __add(integer& a, const integer& b) {
    integer a_copy = a;

    if (a.words.size() < b.words.size())
        a.words.resize(b.words.size());
    bool a_neg = a.is_negative();
    const uint64_t carry = __add_and_return_carry(a, a_neg, b, plus == b.is_negative());
    if (carry)
        a.words.push_back(carry);
    else
        a.words.normalize();
    a.words.set_negative(a_neg);
    return a;
}

constexpr integer& operator+=(integer& a, const integer& b) { return __add<true>(a, b); }
constexpr integer& operator-=(integer& a, const integer& b) { return __add<false>(a, b); }

template<bool plus>
constexpr integer& __add(integer& a, std_int auto b) {
    auto ub = make_unsigned(b);
    if (b < 0) // a + b == a - abs(b) and a - b == a + abs(b)
        return __add<!plus>(a, static_cast<decltype(ub)>(~ub + 1));

    if (plus == !a.is_negative()) {
        *a.magnitude() += ub;
        return a;
    }

    if (!__less_u(static_cast<cnatural>(a), ub)) {
        *a.magnitude() -= ub;
        return a;
    }
    // here plus == a.is_negative(), so the result of the magnitude subtraction is
    // negative only for (a >= 0) - b
    a = ub - static_cast<decltype(ub)>(a.to_natural());
    if constexpr (!plus)
        a.negate();
    return a;
}

constexpr integer& operator+=(integer& a, std_int auto b) { return __add<true>(a, b); }
constexpr integer& operator-=(integer& a, std_int auto b) { return __add<false>(a, b); }

constexpr integer operator+(integer a, std_int auto b) { return a += b; }
constexpr integer operator+(std_int auto a, integer b) { return b += a; }
constexpr integer operator+(integer a, const integer& b) { return a += b; }

constexpr integer operator-(integer a, const integer& b) { return a -= b; }
constexpr integer operator-(integer a, std_int auto b) { return a -= b; }
constexpr integer operator-(std_int auto a, integer b) { b -= a; return -b; }

constexpr void mul(const integer& a, const integer& b, integer& c) {
    c = a.to_natural() * b.to_natural();
    c.words.set_negative(a.is_negative() != b.is_negative());
}

constexpr void mul(integer& a, const integer& b) {
    const bool negative = a.is_negative() != b.is_negative();
    // b may be a itself; take its magnitude before borrowing, since a borrow empties a.words
    const natural bn = b.to_natural();
    *a.magnitude() *= bn;
    a.words.set_negative(negative);
}

constexpr integer operator*(const integer& a, const integer& b) {
    integer c;
    mul(a, b, c);
    return c;
}

constexpr integer operator*(const integer& a, const natural& b) {
    integer c;
    c = a.to_natural() * b;
    c.words.set_negative(a.is_negative());
    return c;
}
constexpr integer operator*(const natural& a, const integer& b) { return b * a; }

// TODO avoid memory allocation here for int128!
constexpr integer operator*(const integer& a, std_int auto b) {
    integer c;
    mul(a, integer(b), c);
    return c;
}

constexpr integer& operator*=(integer& a, const integer& b) {
    mul(a, b);
    return a;
}

constexpr integer& operator*=(integer& a, const natural& b) {
    const bool negative = a.is_negative();
    *a.magnitude() *= b;
    a.words.set_negative(negative);
    return a;
}

constexpr integer& operator*=(integer& a, std_int auto b) {
    const bool negative = a.is_negative() != (b < 0);
    *a.magnitude() *= abs_unsigned(b);
    a.words.set_negative(negative);
    return a;
}

constexpr std::string str(const integer& a) {
    return a.is_negative() ? "-" + str(a.to_natural()) : str(a.to_natural());
}

constexpr std::string stre(const integer& a) {
    return a.is_negative() ? "-" + stre(a.to_natural()) : stre(a.to_natural());
}

// The magnitude products, done on integer's own backend. These mirror natural's add_product and
// sub_product, which are themselves thin wrappers over the same kernels, so nothing is copied and
// no natural temporary is built -- which matters because sub_product sits on the division path.
constexpr void __magnitude_add_product(integer& a, cnatural b, cnatural c) {
    if (b.size == 0 || c.size == 0)
        return;
    const int A = a.words.size();
    a.words.resize(std::max<int>(A, b.size + c.size) + 1);
    vnatural va{{a.words.data(), A}, a.words.capacity()};
    if (b.size < c.size)
        __add_product(va, b, c);
    else
        __add_product(va, c, b);
    a.words.downsize(va.size);
}

constexpr void __magnitude_sub_product(integer& a, cnatural b, cnatural c) {
    if (b.size == 0 || c.size == 0)
        return;
    inatural ia{a.words.data(), a.words.size()};
    const bool ok = (b.size < c.size) ? __sub_product(ia, b, c) : __sub_product(ia, c, b);
    Check(ok, "sub_product() assumes A >= B * C");
    a.words.downsize(ia.size);
}

constexpr void __magnitude_add_product(integer& a, cnatural b, const uint64_t c) {
    if (b.size == 0 || c == 0)
        return;
    const int A = a.words.size();
    a.words.resize(std::max<int>(A, b.size + 1) + 1);
    vnatural va{{a.words.data(), A}, a.words.capacity()};
    __add_product(va, b, c);
    a.words.downsize(va.size);
}

constexpr void __magnitude_sub_product(integer& a, cnatural b, const uint64_t c) {
    if (b.size == 0 || c == 0)
        return;
    inatural ia{a.words.data(), a.words.size()};
    Check(__sub_product(ia, b, c), "sub_product() assumes A >= B * c");
    a.words.downsize(ia.size);
}

constexpr void __magnitude_complement(integer& a) {
    __complement(inatural{a.words.data(), a.words.size()});
    a.words.normalize();
}

template<bool plus>
constexpr void __add_product(integer& a, const integer& b, const integer& c) {
    if (b.words.empty() || c.words.empty())
        return;
    const bool a_negative = a.is_negative();
    const bool bc_negative = b.is_negative() != c.is_negative();

    const cnatural cb = b, cc = c;
    if ((plus && a_negative == bc_negative) || (!plus && a_negative != bc_negative)) {
        __magnitude_add_product(a, cb, cc);
        a.words.set_negative(a_negative);
    } else if (a.num_bits() > b.num_bits() + c.num_bits()) {
        __magnitude_sub_product(a, cb, cc);
        a.words.set_negative(a_negative);
    } else {
        const int m = mul_max_size(cb, cc);
        a.words.resize(m + 1);
        a.words[m] = 1;
        __magnitude_sub_product(a, cb, cc);
        if (a.words.size() > m) {
            a.words[m] -= 1;
            a.words.normalize();
        } else {
            a.words.resize(m);
            __magnitude_complement(a);
            a.words.set_negative(!a_negative);
        }
    }
}

constexpr void add_product(integer& a, const integer& b, const integer& c) { __add_product<true>(a, b, c); }
constexpr void sub_product(integer& a, const integer& b, const integer& c) { __add_product<false>(a, b, c); }

constexpr int __signum(const integer& a) {
    if (a.is_negative())
        return -1;
    return a.is_zero() ? 0 : 1;
}

template<bool plus>
constexpr void __add_product(integer& a, const integer& b, const uint64_t cu, const bool c_negative) {
    if (b.words.empty() || cu == 0)
        return;
    const bool a_negative = a.is_negative();
    const bool bc_negative = b.is_negative() != c_negative;

    if ((plus && a_negative == bc_negative) || (!plus && a_negative != bc_negative)) {
        __magnitude_add_product(a, static_cast<cnatural>(b), cu);
        a.words.set_negative(a_negative);
    } else if (a.num_bits() > b.num_bits() + num_bits(cu)) {
        __magnitude_sub_product(a, static_cast<cnatural>(b), cu);
        a.words.set_negative(a_negative);
    } else {
        const int m = mul_max_size(static_cast<cnatural>(b), {&cu, 1});
        a.words.resize(m + 1);
        a.words[m] = 1;
        __magnitude_sub_product(a, static_cast<cnatural>(b), cu);
        if (a.words.size() > m) {
            a.words[m] -= 1;
            a.words.normalize();
        } else {
            a.words.resize(m);
            __magnitude_complement(a);
            a.words.set_negative(!a_negative);
        }
    }
}

constexpr void add_product(integer& a, const integer& b, const std_int auto c) { __add_product<true>(a, b, abs_unsigned(c), c < 0); static_assert(sizeof(c) <= 8); }
constexpr void sub_product(integer& a, const integer& b, const std_int auto c) { __add_product<false>(a, b, abs_unsigned(c), c < 0); static_assert(sizeof(c) <= 8); }

constexpr void div(const integer& a, const integer& b, integer& quot, integer& rem) {
    if (b == 1) {
        quot = a;
        rem = 0;
        return;
    }
    if (b == -1) {
        quot = -a;
        rem = 0;
        return;
    }
    const bool a_negative = a.is_negative();
    const bool negative = a.is_negative() != b.is_negative();
    // Take the operand magnitudes before borrowing the outputs: quot or rem may be the same object
    // as a or b, and a borrow leaves the borrowed value's words empty while it is alive.
    const natural an = a.to_natural();
    const natural bn = b.to_natural();
    {
        auto q = quot.magnitude();
        auto r = rem.magnitude();
        div(an, bn, *q, *r);
    }
    quot.words.set_negative(negative);
    rem.words.set_negative(a_negative);
}

constexpr integer operator/(const integer& a, const integer& b) {
    integer quot, rem;
    div(a, b, quot, rem);
    return quot;
}

// TODO generalize for any std_int
constexpr int64_t div(const integer& a, int64_t b, integer& quot) {
    if (b == 1) {
        quot = a;
        return 0;
    }
    if (b == -1) {
        quot = -a;
        return 0;
    }
    int64_t rem;
    {
        const natural an = a.to_natural(); // quot may alias a, see div() above
        auto q = quot.magnitude();
        rem = div(an, abs_unsigned(b), *q);
    }
    if (!quot.is_zero())
        quot.words.set_negative(a.is_negative() != (b < 0));
    return a.is_negative() ? -rem : rem;
}

constexpr integer operator/(const integer& a, const std_int auto b) {
    integer c = a.to_natural() / abs_unsigned(b);
    c.words.set_negative(a.is_negative() != (b < 0));
    return c;
}

constexpr integer& operator/=(integer& a, const integer& b) {
    integer rem;
    div(a, b, a, rem);
    return a;
}
constexpr integer& operator/=(integer& a, const std_int auto b) {
    const bool negative = a.is_negative();
    *a.magnitude() /= abs_unsigned(b);
    a.words.set_negative(negative != (b < 0));
    return a;
}

constexpr integer operator%(const integer& a, const integer& divisor) {
    integer quotient, remainder;
    div(a, divisor, quotient, remainder);
    return remainder;
}

// TODO generalize for any std_int
constexpr int64_t operator%(const integer& a, int64_t b) {
    uint64_t m = a.to_natural() % abs_unsigned(b);
    return (a.sign() >= 0) ? m : -static_cast<int64_t>(m);
}

constexpr int operator%(const integer& a, int b) { return a % (int64_t)b; }
constexpr int64_t operator%(const integer& a, unsigned b) { return a % (int64_t)b; }

// Note: return type is integer instead of uint64_t, as it can be negative (can't fit into int64_t either)
constexpr integer operator%(const integer& a, uint64_t b) {
    integer c = a.to_natural() % b;
    if (a.is_negative())
        c.negate();
    return c;
}

constexpr integer mod(const integer& a, const integer& b) {
    // Note: mod() on the magnitudes would resolve back to this function, since natural
    // converts implicitly to integer. Call the natural kernel explicitly.
    natural r;
    mod(a.to_natural(), b.to_natural(), /*out*/r);
    if (a.is_negative() && !r.words.empty()) {
        // result is in [0, abs(b)) range
        natural e = b.to_natural();
        e -= r;
        return e;
    }
    return r;
}

constexpr uint64_t mod(const integer& a, uint64_t b) {
    uint64_t m = a.to_natural() % b;
    return (a.is_negative() && m) ? (b - m) : m;
}

// TODO this version doesn't use uint128_t %, is it faster?
// TODO if it is, move it to natural.h instead, it doesn't belong here
constexpr unsigned mod(const integer& a, uint32_t b) {
    if (b == 0)
        throw std::runtime_error("division by zero");
    uint64_t m = (uint64_t(1) << 32) % b;
    m = (m * m) % b;

    uint64_t acc = 0;
    for (auto i = a.words.size(); i-- > 0;) {
        acc *= m;
        acc += a.words[i] % b;
        acc %= b;
    }
    if (a.is_negative()) {
        acc *= b - 1;
        acc %= b;
    }
    return acc;
}

constexpr integer& operator%=(integer& a, const integer& b) { a = a % b; return a; }
// TODO issue temporary memory allocation for cent / ucent
constexpr integer& operator%=(integer& a, std_int auto b) { a = a % integer(b); return a; }

constexpr bool operator<(const integer& a, const integer& b) {
    if (a.is_negative())
        return !b.is_negative() || __less(static_cast<cnatural>(b), static_cast<cnatural>(a));
    return !b.is_negative() && __less(static_cast<cnatural>(a), static_cast<cnatural>(b));
}
// TODO issue temporary memory allocation for cent / ucent
constexpr bool operator<(const integer& a, std_int auto b) { return a < integer(b); }
constexpr bool operator<(std_int auto a, const integer& b) { return integer(a) < b; }

constexpr integer operator~(integer a) {
    if (a.sign() >= 0) {
        a += 1;
        a.negate();
        return a;
    }
    a.negate();
    a -= 1;
    return a;
}

namespace literals {
constexpr auto operator""_i(const char* s) { return integer(s); }
}

constexpr void operator<<=(integer& a, int64_t i) {
    const bool negative = a.is_negative();
    *a.magnitude() <<= i;
    a.words.set_negative(negative);
}

ALGEBRA_SHIFT_OP(integer)

static_assert(sizeof(integer) == 16);

}

template <>
struct std::formatter<algebra::integer, char> : public std::formatter<algebra::natural, char> {
    constexpr auto format(const algebra::integer& a, auto& ctx) const { return format_padded(a, ctx); }
};

#if 0
template <>
struct std::formatter<algebra::neg_integer, char> {
    constexpr auto parse(auto& ctx) {
        auto it = ctx.begin();
        if (it == ctx.end() || *it != '}')
            throw std::format_error("Invalid format specifier for integer.");
        return it;
    }

    constexpr auto format(const algebra::neg_integer& a, auto& ctx) const {
        if (a.a->sign() >= 0)
            return std::format_to(ctx.out(), "{}", a.a->abs);
        return std::format_to(ctx.out(), "-{}", a.a->abs);
    }
};
#endif

constexpr std::ostream& operator<<(std::ostream& os, const algebra::integer& a) { return os << a.str(); }

template<>
struct std::hash<algebra::integer> {
    constexpr size_t operator()(const algebra::integer& a) const {
        // the backend hash includes the sign in the size, which is what makes -x hash apart from x
        return std::hash<algebra::integer_backend>()(a.words);
    }
};
