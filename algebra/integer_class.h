#pragma once
#include "algebra/kernels.h"
#include "algebra/integer_backend.h"

namespace algebra {

struct integer;
template<> struct IsNumberClass<integer> : std::true_type {};

// The magnitude of an integer, as a value. The sign aware code reaches the in place magnitude
// operations through the __abs_ functions instead, which is what keeps the two layers apart.
constexpr integer abs(const integer& a);



struct integer {

    // The whole value: the magnitude in the words, the sign in the sign of the size, which is how
    // integer_backend already stores it. Magnitude arithmetic goes through the cnatural/vnatural/
    // inatural conversions below, or through __magnitude() when it has to happen in place.
    integer_backend words;

    constexpr integer(std::initializer_list<uint64_t> a) : words(a) { } // words, least significant first
    constexpr integer() {}
    // integer_backend's constructor already puts the sign in the sign of the size
    constexpr integer(std_int auto a) : words(a) { }
    constexpr integer(integer&& o) : words(std::move(o.words)) { }
    constexpr integer(const integer& o) : words(o.words) { }

    constexpr void operator=(std_int auto a) { words = a; }
    constexpr void operator=(integer&& o) { words = std::move(o.words); }
    constexpr void operator=(const integer& o) { words = o.words; }

    // The low word, which is only meaningful when the value is not zero.
    constexpr uint64_t low_word() const { return words.size() ? words[0] : 0; }

    constexpr auto sign() const { return words.sign(); }
    constexpr bool is_negative() const { return sign() < 0; }
    constexpr bool is_even() const { return (low_word() % 2) == 0; }
    constexpr bool is_odd() const { return low_word() % 2; }
    constexpr bool is_one() const { return words.size() == 1 && words[0] == 1 && sign() >= 0; }
    constexpr bool is_zero() const { return words.size() == 0; }

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

    constexpr void __abs_parse(std::string_view s, unsigned base);
    constexpr integer(std::string_view s, unsigned base = 10) {
        const bool negative = s.size() && s[0] == '-';
        __abs_parse(negative ? s.substr(1) : s, base);
        if (negative)
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

    constexpr int str_size_upper_bound(unsigned base = 10) const { return is_negative() + __abs_str_size_upper_bound_member(base); }
    constexpr int str(char* buffer, int buffer_size, unsigned base = 10, bool upper = true) const {
        int result = 0;
        if (is_negative()) {
            if (buffer_size <= 0)
                throw std::runtime_error("buffer too small");
            buffer_size -= 1;
            *buffer++ = '-';
            result = 1;
        }
        return result + __abs_str_buffer(buffer, buffer_size, base, upper);
    }

    constexpr std::string str(unsigned base = 10, bool upper = true) const {
        std::string s;
        s.resize(str_size_upper_bound(base));
        s.resize(str(s.data(), s.size(), base, upper));
        return s;
    }
    constexpr std::string hex() const { return str(16); }

    // ---- the magnitude layer, ported from integer ----
    // These work on the words and leave the sign alone. The sign aware code reaches them through
    // the __abs_ free functions, never by an operator name.
    constexpr void __abs_mul_add(uint64_t a, uint64_t b) {
        if (a == 0) {
            words = b;
            return;
        }
        uint64_t carry = __mul_add_return_carry(*this, a, b);
        if (carry)
            words.push_back(carry);
    }

    constexpr integer& __abs_add_word(const uint64_t b) {
        // note: for an empty (zero) value the carry is b itself, not 1
        const uint64_t carry = __add_and_return_carry(*this, b);
        if (carry)
            words.push_back(carry);
        return *this;
    }

    constexpr integer& __abs_add_word(const uint128_t b) {
        uint128_t carry = __add_and_return_carry(*this, b);
        if (carry) {
            if (carry >> 64)
                words.push_back(carry, carry >> 64);
            else
                words.push_back(carry);
        }
        return *this;
    }

    constexpr integer& __abs_sub_word(uint64_t b) {
        // note: spelled out instead of *this >= b, because operator<(integer, unsigned) is
        // declared further down this header and is not visible from inside the class
        Check(words.size() > 1 || words[0] >= b, "integer can't be negative");
        inatural a = *this;
        __sub(a, b);
        words.downsize(a.size);
        return *this;
    }

    constexpr integer& __abs_sub_word(uint128_t b) {
        if (b <= UINT64_MAX)
            return __abs_sub_word(static_cast<uint64_t>(b));
        Check(words.size() > 2 || (words.size() == 2 && unsafe_u128() >= b), "integer can't be negative");
        inatural a = *this;
        __sub(a, b);
        words.downsize(a.size);
        return *this;
    }

    constexpr integer& __abs_increment() {
        if (__increment_and_return_carry(*this))
            words.push_back(1);
        return *this;
    }

    constexpr integer& __abs_decrement() {
        Check(!words.empty(), "decrementing zero integer");
        inatural a = *this;
        __decrement(a);
        words.downsize(a.size);
        return *this;
    }

    constexpr int64_t __abs_popcount_member() const {
        int64_t c = 0;
        for (int i = 0; i < words.size(); i++)
            c += std::popcount(words[i]);
        return c;
    }

    constexpr int __abs_str_buffer(char* buffer, int buffer_size, unsigned base, bool upper) const;
    constexpr int __abs_str_size_upper_bound_member(unsigned base) const;

    constexpr uint128_t unsafe_u128() const {
        uint128_t a = words[0];
        if (words.size() > 1)
            a |= static_cast<uint128_t>(words[1]) << 64;
        return a;
    }

    // the magnitude modulo a builtin, which fits in the builtin
    constexpr uint128_t __abs_mod_word(const uint128_t b) const {
        Check(b > 0, "division by zero");
        if (b <= UINT64_MAX)
            return __mod(static_cast<cnatural>(*this), static_cast<uint64_t>(b));
        return __mod(static_cast<cnatural>(*this), b);
    }

    constexpr void set_zero() { words.set_zero(); }
    constexpr void negate() { words.negate(); }

    constexpr size_t popcount() const {
        if (!is_negative())
            return __abs_popcount_member();

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


    constexpr integer& operator++();
    constexpr integer& operator--();
    constexpr integer operator++(int) { integer a = *this; operator++(); return a; }
    constexpr integer operator--(int) { integer a = *this; operator--(); return a; }

    constexpr bool bit(int64_t i) const {
        const size_t w = i / 64;
        return w < static_cast<size_t>(words.size()) && (words[w] & (uint64_t(1) << (i % 64)));
    }
    constexpr auto num_bits() const { return algebra::num_bits(static_cast<cnatural>(*this)); }
    constexpr auto num_trailing_zeros() const { return algebra::num_trailing_zeros(static_cast<cnatural>(*this)); }

    template<std::floating_point T> constexpr T __abs_to_float() const;

    template<std::floating_point T>
    constexpr operator T() const {
        const T a = __abs_to_float<T>();
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

    constexpr uint64_t mod6() const {
        uint64_t m = algebra::mod6(static_cast<cnatural>(*this));
        if (is_negative())
            m = (m * 5) % 6;
        return m;
    }

    constexpr uint64_t mod7() const {
        uint64_t m = algebra::mod7(static_cast<cnatural>(*this));
        if (is_negative())
            m = (m * 6) % 7;
        return m;
    }

    constexpr uint64_t mod9() const {
        uint64_t m = algebra::mod9(static_cast<cnatural>(*this));
        if (is_negative())
            m = (m * 8) % 9;
        return m;
    }

    constexpr uint64_t mod10() const {
        uint64_t m = algebra::mod10(static_cast<cnatural>(*this));
        if (is_negative())
            m = (m * 9) % 10;
        return m;
    }

    constexpr uint64_t mod8() const {
        uint64_t m = low_word() % 8;
        if (is_negative())
            m = (m * 7) % 8;
        return m;
    }
};

class magnitude_ref {
    integer_backend& _words;
    bool _negative;
public:
    integer value;
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
    constexpr integer* operator->() { return &value; }
    constexpr integer& operator*() { return value; }
};
constexpr magnitude_ref magnitude(integer& a) { return magnitude_ref(a.words); }

constexpr void __mul(const integer& a, const integer& b, integer& q);
constexpr void __mul(const integer& a, uint128_t b, integer& out);
constexpr void __mul(const uint64_t* a, const int A, uint64_t b, uint64_t carry, integer& out);
constexpr void __mul(const integer& a, uint64_t b, uint64_t carry, integer& out);

constexpr void __square(integer& a);
constexpr void __abs_mul(cnatural a, cnatural b, integer& out);
constexpr void __abs_mul(integer& a, const integer& b);
constexpr integer& __abs_mul(integer& a, std_unsigned_int auto b);
constexpr integer __low_words(const integer& a, int n);
constexpr void __divide_2n1n(const integer& a, const integer& b, int n, integer& q, integer& r);
constexpr void __divide_3n2n(const integer& a, const integer& b, int n, integer& q, integer& r);
constexpr void divide_bz(const integer& a, const integer& d, integer& q, integer& r);
constexpr uint64_t __abs_div(const integer& a, uint64_t b, integer& q);
constexpr void __abs_div(cnatural a, cnatural b, integer& q, integer& r);
constexpr void __abs_mod(cnatural a, cnatural b, integer& r);
constexpr void __abs_mod(integer& a, const integer& b);
constexpr integer& __abs_shl(integer& a, int64_t b);
constexpr integer& __abs_div(integer& a, std_int auto b);
constexpr std::string __abs_stre(const integer& a);

constexpr integer& __abs_add(integer& a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a.__abs_add_word(static_cast<uint64_t>(b));
    static_assert(sizeof(b) <= 16);
    return a.__abs_add_word(static_cast<uint128_t>(b));
}

constexpr integer& __abs_sub(integer& a, const std_unsigned_int auto b) {
    if (b <= UINT64_MAX)
        return a.__abs_sub_word(static_cast<uint64_t>(b));
    static_assert(sizeof(b) <= 16);
    return a.__abs_sub_word(static_cast<uint128_t>(b));
}

constexpr integer& __abs_inc(integer& a) { return a.__abs_increment(); }
constexpr integer& __abs_dec(integer& a) { return a.__abs_decrement(); }
constexpr size_t __abs_popcount(const integer& a) { return a.__abs_popcount_member(); }
constexpr int __abs_str_size_upper_bound(const integer& a, unsigned base) { return a.__abs_str_size_upper_bound_member(base); }
constexpr int __abs_str(const integer& a, char* buffer, int buffer_size, unsigned base, bool upper) {
    return a.__abs_str_buffer(buffer, buffer_size, base, upper);
}
constexpr std::string __abs_str(const integer& a) { return str(static_cast<cnatural>(a)); }


constexpr integer& integer::operator++() {
    // the magnitude increment and decrement preserve the sign of words, and the borrow brings a
    // result of zero back positive
    if (is_negative())
        __abs_dec(*magnitude(*this));
    else
        __abs_inc(*magnitude(*this));
    return *this;
}

constexpr integer& integer::operator--() {
    if (is_negative())
        __abs_inc(*magnitude(*this));
    else if (words.size() > 0)
        __abs_dec(*magnitude(*this));
    else
        *this = -1;
    return *this;
}

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
        __abs_add(*magnitude(a), ub);
        return a;
    }

    if (!__less_u(static_cast<cnatural>(a), ub)) {
        __abs_sub(*magnitude(a), ub);
        return a;
    }
    // here plus == a.is_negative(), so the result of the magnitude subtraction is
    // negative only for (a >= 0) - b
    a = ub - static_cast<decltype(ub)>(abs(a));
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
    const bool negative = a.is_negative() != b.is_negative();
    // views for the operands and a borrow for the result, so nothing is copied. c may alias a or b,
    // so the views are taken before the borrow empties c.
    const cnatural ca = a, cb = b;
    {
        auto m = magnitude(c);
        __abs_mul(ca, cb, *m);
    }
    c.words.set_negative(negative);
}

constexpr void mul(integer& a, const integer& b) {
    const bool negative = a.is_negative() != b.is_negative();
    if (&a == &b) {
        auto m = magnitude(a);
        __abs_mul(*m, *m); // squaring: the magnitude mul() takes the &a == &b path
    } else {
        // b is a distinct object, so a view of it stays valid while a's magnitude is borrowed
        const cnatural cb = b;
        integer out;
        __abs_mul(static_cast<cnatural>(a), cb, out);
        a.words = std::move(out.words);
    }
    a.words.set_negative(negative);
}

constexpr integer operator*(const integer& a, const integer& b) {
    integer c;
    mul(a, b, c);
    return c;
}

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


constexpr integer& operator*=(integer& a, std_int auto b) {
    const bool negative = a.is_negative() != (b < 0);
    __abs_mul(*magnitude(a), abs_unsigned(b));
    a.words.set_negative(negative);
    return a;
}

constexpr std::string str(const integer& a) {
    const std::string m = __abs_str(abs(a));
    return a.is_negative() ? "-" + m : m;
}

constexpr std::string stre(const integer& a) {
    const std::string m = __abs_stre(abs(a));
    return a.is_negative() ? "-" + m : m;
}

// The magnitude products, done on integer's own backend. These mirror integer's add_product and
// sub_product, which are themselves thin wrappers over the same kernels, so nothing is copied and
// no integer temporary is built -- which matters because sub_product sits on the division path.
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
    // Each output borrows the magnitude of its target below, and magnitude_ref restores on scope
    // exit; two borrows over one backend would restore in turn, so the second would put the
    // pre-operation value back and the result would be lost without a word of warning.
    Check(&quot != &rem, "div() needs separate quotient and remainder");
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
    const integer an = abs(a);
    const integer bn = abs(b);
    {
        auto q = magnitude(quot);
        auto r = magnitude(rem);
        __abs_div(an, bn, *q, *r);
    }
    quot.words.set_negative(negative);
    rem.words.set_negative(a_negative);
}

constexpr integer operator/(const integer& a, const integer& b) {
    integer quot, rem;
    div(a, b, quot, rem);
    return quot;
}

// The magnitude divided by an unsigned builtin of any width, returning the remainder.
// A divisor above UINT64_MAX cannot go through the single word kernel, which would keep only its
// low word and silently divide by the wrong number, so it takes the general kernel instead. The
// remainder is below the divisor and so fits into T.
template<std_unsigned_int T>
constexpr T __abs_div_word(const integer& a, T b, integer& q) {
    Check(b != 0, "division by zero");
    if constexpr (sizeof(T) > 8) {
        if (b > UINT64_MAX) {
            integer quot, rem;
            __abs_div(static_cast<cnatural>(a), static_cast<cnatural>(integer(b)), quot, rem);
            q = std::move(quot);
            return static_cast<T>(rem);
        }
    }
    return __abs_div(a, static_cast<uint64_t>(b), q);
}

// The divisor's signedness decides its own overload. A single int64_t parameter would take a
// uint64_t above INT64_MAX by conversion and read it as negative, which negated the quotient.
template<std_signed_int T>
constexpr T div(const integer& a, T b, integer& quot) {
    if (b == 1) {
        quot = a;
        return 0;
    }
    if (b == -1) {
        quot = -a;
        return 0;
    }
    make_unsigned_t<T> rem; // an int128_t divisor leaves a remainder that does not fit a word
    {
        const integer an = abs(a); // quot may alias a, see div() above
        auto q = magnitude(quot);
        rem = __abs_div_word(an, abs_unsigned(b), *q);
    }
    if (!quot.is_zero())
        quot.words.set_negative(a.is_negative() != (b < 0));
    const T r = static_cast<T>(rem);
    return a.is_negative() ? -r : r;
}

template<std_unsigned_int T>
constexpr T div(const integer& a, T b, integer& quot) {
    if (b == 1) {
        quot = a;
        return 0;
    }
    T rem;
    {
        const integer an = abs(a); // quot may alias a, see div() above
        auto q = magnitude(quot);
        rem = __abs_div_word(an, b, *q);
    }
    if (!quot.is_zero())
        quot.words.set_negative(a.is_negative());
    return rem;
}

constexpr integer operator/(const integer& a, const std_int auto b) {
    // the magnitude division by name: abs(a) is an integer, so operator/ here would recurse
    integer c = a;
    c.words.set_negative(false);
    __abs_div(c, abs_unsigned(b));
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
    __abs_div(*magnitude(a), abs_unsigned(b));
    a.words.set_negative(negative != (b < 0));
    return a;
}

constexpr integer operator%(const integer& a, const integer& divisor) {
    // the remainder alone: going through div() would build a quotient only to discard it
    const bool negative = a.is_negative();
    integer r;
    __abs_mod(static_cast<cnatural>(a), static_cast<cnatural>(divisor), r);
    integer c = std::move(r);
    c.words.set_negative(negative && !c.is_zero());
    return c;
}

// TODO generalize for any std_int
constexpr int64_t operator%(const integer& a, int64_t b) {
    Check(b != 0, "division of integer by zero");
    uint64_t m = __mod(static_cast<cnatural>(a), abs_unsigned(b));
    return (a.sign() >= 0) ? m : -static_cast<int64_t>(m);
}

constexpr int operator%(const integer& a, int b) { return a % (int64_t)b; }
constexpr int64_t operator%(const integer& a, unsigned b) { return a % (int64_t)b; }

// Note: return type is integer instead of uint64_t, as it can be negative (can't fit into int64_t either)
constexpr integer operator%(const integer& a, uint64_t b) {
    Check(b > 0, "division of integer by zero");
    integer c = __mod(static_cast<cnatural>(a), b);
    if (a.is_negative())
        c.negate();
    return c;
}

constexpr integer mod(const integer& a, const integer& b) {
    // Note: mod() on the magnitudes would resolve back to this function, since integer
    // converts implicitly to integer. Call the integer kernel explicitly.
    integer r;
    __abs_mod(abs(a), abs(b), /*out*/r);
    if (a.is_negative() && !r.words.empty()) {
        // result is in [0, abs(b)) range
        integer e = abs(b);
        e -= r;
        return e;
    }
    return r;
}

constexpr uint64_t mod(const integer& a, uint64_t b) {
    const uint64_t m = __mod(static_cast<cnatural>(a), b); // by name: abs(a) % b would recurse
    return (a.is_negative() && m) ? (b - m) : m;
}

// TODO this version doesn't use uint128_t %, is it faster?
// TODO if it is, move it to integer.h instead, it doesn't belong here
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
    __abs_shl(*magnitude(a), i);
    a.words.set_negative(negative);
}

ALGEBRA_SHIFT_OP(integer)

constexpr void __mul(const integer& a, const integer& b, integer& q) {
    if (&a != &q)
        q.set_zero();
    auto A = a.words.size();
    auto B = b.words.size();
    q.words.resize(A + B);
    vnatural vq = q;
    __mul({a.words.data(), A}, {b.words.data(), B}, vq, /*init*/false);
    q.words.downsize(vq.size);
}

constexpr void __mul(const integer& a, uint128_t b, integer& out) {
    if (&a != &out)
        out.set_zero();
    auto as = a.words.size();
    out.words.resize(a.words.size() + 2);

    for (auto i = as; i-- > 0;) {
        const auto aw = a.words[i];
        if (aw == 0)
            continue;

        uint64_t* ow = &out.words[i];
        *ow = 0;

        auto acc = __mulq(aw, b);
        acc += *ow;
        *ow = acc;
        acc >>= 64;

        acc += __mulq(aw, b >> 64);
        acc += *++ow;
        *ow = acc;
        acc >>= 64;

        if (acc) {
            acc += *++ow;
            *ow = acc;
            acc >>= 64;

            if (acc) {
                acc += *++ow;
                *ow = acc;
            }
        }
    }
    out.words.normalize();
}

constexpr void __mul(const uint64_t* a, const int A, uint64_t b, uint64_t carry, integer& out) {
    if (b == 0) {
        out.words.set_zero();
        if (carry)
            out.words.push_back(carry);
        return;
    }
    if (a != out.words.data())
        out.words.reset(A, /*initialize*/false);
    for (int i = 0; i < A; ++i) {
        auto acc = __mulq(a[i], b) + carry;
        out.words[i] = acc;
        carry = acc >> 64;
    }
    if (carry)
        out.words.push_back(carry);
}

constexpr void __mul(const integer& a, uint64_t b, uint64_t carry, integer& out) {
    __mul(a.words.data(), a.words.size(), b, carry, out);
}

// ---- the magnitude layer definitions ----

// bitwise operations on magnitudes; integer has no two's complement bitwise layer of its own
template<bool GROW, typename Op>
constexpr integer& __abs_bitwise_op(integer& a, const integer& b, Op op) {
    const int bs = b.words.size();
    if constexpr (GROW) {
        if (bs > a.words.size())
            a.words.resize(bs);
        for (int i = 0; i < bs; i++)
            a.words[i] = op(a.words[i], b.words[i]);
    } else {
        const int n = std::min<int>(a.words.size(), bs);
        if (n == 0) {
            a.words.set_zero(); // downsize(0) alone would leave a stale word behind
            return a;
        }
        a.words.downsize(n);
        for (int i = 0; i < n; i++)
            a.words[i] = op(a.words[i], b.words[i]);
    }
    a.words.normalize();
    return a;
}

constexpr integer& __abs_or(integer& a, const integer& b) {
    return __abs_bitwise_op<true>(a, b, [](uint64_t x, uint64_t y) { return x | y; });
}
constexpr integer& __abs_and(integer& a, const integer& b) {
    return __abs_bitwise_op<false>(a, b, [](uint64_t x, uint64_t y) { return x & y; });
}
constexpr integer& __abs_xor(integer& a, const integer& b) {
    return __abs_bitwise_op<true>(a, b, [](uint64_t x, uint64_t y) { return x ^ y; });
}

constexpr integer __abs_or_copy(integer a, const integer& b) { return __abs_or(a, b); }
constexpr integer __abs_and_copy(integer a, const integer& b) { return __abs_and(a, b); }
constexpr integer __abs_xor_copy(integer a, const integer& b) { return __abs_xor(a, b); }

constexpr integer& __abs_or(integer& a, uint64_t b) {
    if (a.words.size() == 0)
        a = b;
    else
        a.words[0] |= b;
    return a;
}

constexpr integer& __abs_and(integer& a, uint64_t b) {
    if (a.words.size() == 0)
        return a;
    a.words.downsize(1);
    a.words[0] &= b;
    a.words.normalize();
    return a;
}

constexpr integer& __abs_xor(integer& a, uint64_t b) {
    if (a.words.size() == 0) {
        a = b;
    } else {
        a.words[0] ^= b;
        a.words.normalize();
    }
    return a;
}

constexpr integer __abs_or_copy(integer a, uint64_t b) { return __abs_or(a, b); }
constexpr integer __abs_and_copy(integer a, uint64_t b) { return __abs_and(a, b); }
constexpr integer __abs_xor_copy(integer a, uint64_t b) { return __abs_xor(a, b); }



constexpr void square(integer& a) {
    if (a.words.size() == 0)
        return;
    if (a.words.size() == 1) {
        uint128_t p = __mulq(a.words[0], a.words[0]);

        a.words.reset_one_without_init();
        a.words[0] = p;

        uint64_t high = p >> 64;
        if (high)
            a.words.push_back(high);
        return;
    }
    __square(a);
}


// q must not alias a or b
constexpr void mul_karatsuba(const integer& a, const integer& b, integer& q) {
    auto A = a.words.size();
    auto B = b.words.size();
    // TODO this case will disappear when a.words.empty() is removed!
    if (A == 0 || B == 0) {
        q.set_zero();
        return;
    }
    if (A == 1) {
        if (B == 1) {
            q = __mulq(a.words[0], b.words[0]);
            return;
        }
        if (a.words[0] == 1) {
            q = b;
            return;
        }
    }
    if (b.words[0] == 1 && B == 1) {
        q = a;
        return;
    }

    if (is_power_of_two(static_cast<cnatural>(a))) { // the kernel by name: integer converts to both cnatural and uint64_t
        const size_t z = (A - 1) * 64 + std::countr_zero(a.words[A - 1]); // = a.num_trailing_zeros() but O(1)
        const size_t bits = b.num_bits() + z;
        const size_t words = (bits + 63) / 64;
        q.words.reset(words, /*init*/false); // preallocates memory!
        // TODO this can be done directly without moving
        q = b;
        q <<= z;
        return;
    }
    if (is_power_of_two(static_cast<cnatural>(b))) {
        const size_t z = (B - 1) * 64 + std::countr_zero(b.words[B - 1]); // = b.num_trailing_zeros() but O(1)
        const size_t bits = a.num_bits() + z;
        const size_t words = (bits + 63) / 64;
        q.words.reset(words, /*init*/false); // preallocates memory!
        // TODO this can be done directly without moving
        q = a;
        q <<= z;
        return;
    }

    int Q = mul_max_size(a, b);
    q.words.reset(Q);
    vnatural vq = q;
    if (std::min(A, B) <= 2 || std::max(A, B) < KARATSUBA_LIMIT) {
        __add_product(vq, static_cast<cnatural>(a), static_cast<cnatural>(b));
    } else {
        const int W = 4 * std::max(A, B);
        if (W <= 1024) {
            uint64_t w[1024];
            __mul_karatsuba_rec(a, b, vq, w, w + W);
        } else {
            auto w = new uint64_t[W];
            __mul_karatsuba_rec(a, b, vq, w, w + W);
            delete[] w;
        }
    }
    q.words.resize(vq.size);
    // mul_max_size() is an upper bound and __add_product() does not normalize, so the top word can
    // still be zero here. An unnormalized value compares unequal to the same value from operator*.
    q.words.normalize();
}

constexpr integer mul_karatsuba(const integer& a, const integer& b) {
    integer q;
    mul_karatsuba(a, b, q);
    return q;
}


template<std::floating_point T>
constexpr T integer::__abs_to_float() const {
    // the magnitude only; operator T() applies the sign
    const int exponent = static_cast<int>(num_bits()) - std::numeric_limits<T>::digits;
    if (exponent >= std::numeric_limits<T>::max_exponent)
        return std::numeric_limits<T>::infinity();
    if (exponent <= 0)
        return words[0];
    const auto m = extract_u64(*this, exponent);
    return std::ldexp(static_cast<T>(m), exponent);
}


constexpr void __abs_mod(integer& a, const integer& b) {
    const int A = a.words.size();
    const int B = b.words.size();

    Check(B != 0, "division by zero");
    if (B <= 1) {
        a = __mod(a, b.words[0]);
        return;
    }

    uint64_t* r = a.words.data() + A;
    int R = 0;
    for (auto i = A; i-- > 0;) {
        r -= 1;
        R += 1;
        inatural ir {r, R};
        const uint64_t w = __saturated_div(ir, b);
        // r -= b * w. Unlike the four other call sites this one dropped the result, so an estimate
        // that was too large would leave a silently corrupted remainder behind instead of failing.
        Check(__sub_product(ir, b, w), "__saturated_div() overestimated the quotient");
        R = ir.size;
    }
    a.words.resize(R);
}

// After the shift operators, since the division helpers below shift their operands.

constexpr void __square(integer& a) {
    if (a.words.size() == 2) {
        auto carry = __mulq(a.words[0], a.words[0]);
        uint64_t b0 = carry;
        uint64_t b1 = carry >> 64;

        uint128_t pq = __mulq(a.words[0], a.words[1]);
        carry = b1 + pq; // can't use pq*2 here due to dword overflow
        b1 = carry;
        uint64_t b2 = carry >> 64;

        carry = b1 + pq;
        b1 = carry;
        carry >>= 64;
        carry += b2;
        b2 = carry;

        uint64_t b3 = carry >> 64;
        carry = b2 + __mulq(a.words[1], a.words[1]);
        b2 = carry;
        b3 += carry >> 64;

        if (b3) {
            a.words.reset(4, /*initialize*/false);
            a.words[0] = b0;
            a.words[1] = b1;
            a.words[2] = b2;
            a.words[3] = b3;
            return;
        }
        if (b2) {
            a.words.reset(3, /*initialize*/false);
            a.words[0] = b0;
            a.words[1] = b1;
            a.words[2] = b2;
            return;
        }
        a.words.reset(2, /*initialize*/false);
        a.words[0] = b0;
        a.words[1] = b1;
        return;
    }

    // a*a == sum(a[i]*a[i] << 2i) + 2 * sum(a[i]*a[j] << (i+j)) for i<j,
    // which needs half the multiplications of a full a*a
    const int n = a.words.size();
    integer r;
    r.words.reset(n * 2); // zero initialized

    // cross products, each counted once
    for (int i = 0; i < n; i++) {
        const uint64_t ai = a.words[i];
        if (ai == 0)
            continue;
        uint128_t carry = 0;
        for (int j = i + 1; j < n; j++) {
            carry += r.words[i + j];
            carry += __mulq(ai, a.words[j]);
            r.words[i + j] = carry;
            carry >>= 64;
        }
        for (int k = i + n; carry; k++) {
            carry += r.words[k];
            r.words[k] = carry;
            carry >>= 64;
        }
    }

    // double them (the sum of cross products is less than a*a/2, so this fits)
    uint64_t shifted = 0;
    for (int k = 0; k < n * 2; k++) {
        const uint64_t w = r.words[k];
        r.words[k] = (w << 1) | shifted;
        shifted = w >> 63;
    }

    // add the squares of the individual words
    uint128_t carry = 0;
    for (int i = 0; i < n; i++) {
        carry += r.words[i * 2];
        carry += __mulq(a.words[i], a.words[i]);
        r.words[i * 2] = carry;
        carry >>= 64;
        carry += r.words[i * 2 + 1];
        r.words[i * 2 + 1] = carry;
        carry >>= 64;
    }

    r.words.normalize();
    a = std::move(r);
}

constexpr void __abs_mul(cnatural a, cnatural b, integer& out) {
    if (a.size == 0 || b.size == 0) {
        out.set_zero();
        return;
    }

    if (a.size == 1) {
        if (b.size == 1) {
            uint128_t p = __mulq(a[0], b[0]);
            out.words.reset_one_without_init();
            out.words[0] = p;
            uint64_t high = p >> 64;
            if (high)
                out.words.push_back(high);
            return;
        }
        __mul(b.words, b.size, a[0], /*carry*/0, out);
        return;
    }

    if (b.size == 1) {
        __mul(a.words, a.size, b[0], /*carry*/0, out);
        return;
    }

    out.set_zero();
    out.words.resize(a.size + b.size);
    vnatural vq = out;
    __mul(a, b, vq, /*init*/false);
    out.words.downsize(vq.size);
}

constexpr void __abs_mul(integer& a, const integer& b) {
    if (a.words.size() == 0 || b.words.size() == 0) {
        a.words.set_zero();
        return;
    }

    if (b.words.size() == 1) {
        a.__abs_mul_add(b.words[0], /*carry*/0);
        return;
    }

    if (a.words.size() == 1) {
        __mul(b, a.words[0], /*carry*/0, a);
        return;
    }

    if (&a == &b)
        __square(a);
    else
        __mul(a, b, a);
}

constexpr integer& __abs_mul(integer& a, std_unsigned_int auto b) {
    if (b <= UINT64_MAX) {
        a.__abs_mul_add(static_cast<uint64_t>(b), 0);
    } else {
        static_assert(sizeof(b) == 16 || sizeof(b) <= 8);
        __mul(a, static_cast<uint128_t>(b), a);
    }
    return a;
}

constexpr integer __low_words(const integer& a, int n) {
    const int s = std::min<int>(n, a.words.size());
    integer c;
    c.words.reset(s, /*initialize*/false);
    for (int i = 0; i < s; i++)
        c.words[i] = a.words[i];
    c.words.normalize();
    return c;
}

constexpr void __divide_2n1n(const integer& a, const integer& b, int n, integer& q, integer& r) {
    if ((n & 1) || n <= BZ_BASE_CASE_SIZE) {
        __abs_div(a, b, /*out*/q, /*out*/r);
        return;
    }
    const int h = n / 2;

    // A = [a1 a2 a3 a4] and B = [b1 b2], each part of h words
    integer q1, r1;
    __divide_3n2n(a >> (64 * h), b, h, /*out*/q1, /*out*/r1); // [a1 a2 a3] / B

    integer t = r1 << (64 * h);
    t += __low_words(a, h); // [r1 a4]

    integer q2;
    __divide_3n2n(t, b, h, /*out*/q2, /*out*/r);

    q = q1 << (64 * h);
    q += q2;
}

constexpr void __divide_3n2n(const integer& a, const integer& b, int n, integer& q, integer& r) {
    // A = [a1 a2 a3] and B = [b1 b2], each part of n words
    const integer b1 = b >> (64 * n);
    const integer a12 = a >> (64 * n);

    integer r1;
    if ((a12 >> (64 * n)) < b1) { // a1 < b1
        __divide_2n1n(a12, b1, n, /*out*/q, /*out*/r1);
    } else {
        // quotient would not fit into n words, so use the largest value that does
        q = 1;
        q <<= 64 * n;
        q -= 1u;
        // r1 = [a1 a2] - b1 * q = [a1 a2] - (b1 << 64*n) + b1
        r1 = a12;
        r1 -= b1 << (64 * n);
        r1 += b1;
    }

    // r = [r1 a3] - q * b2, with at most two corrections
    r = r1 << (64 * n);
    r += __low_words(a, n);
    const integer qb2 = q * __low_words(b, n);
    while (r < qb2) {
        q -= 1u;
        r += b;
    }
    r -= qb2;
}

constexpr void divide_bz(const integer& a, const integer& d, integer& q, integer& r) {
    Check(!d.words.empty(), "division by zero");
    const int n = d.words.size();
    if (n <= BZ_BASE_CASE_SIZE || a.words.size() <= n) {
        __abs_div(a, d, /*out*/q, /*out*/r);
        return;
    }

    // normalize, so that the top bit of the divisor is set
    const int shift = std::countl_zero(d.words.back());
    const integer an = a << shift;
    const integer dn = d << shift;

    // divide chunk by chunk, most significant chunk first
    const int chunks = (an.words.size() + n - 1) / n;
    q.set_zero();
    r.set_zero();
    integer qi, t;
    for (int i = chunks; i-- > 0;) {
        t = r << (64 * n);
        t += __low_words(an >> (64 * n * i), n);
        __divide_2n1n(t, dn, n, /*out*/qi, /*out*/r);
        q <<= 64 * n;
        q += qi;
    }
    r >>= shift;
}

constexpr uint64_t __abs_div(const integer& a, uint64_t b, integer& q) {
    Check(b != 0, "division by zero");
    if (&a != &q)
        q.words.reset(a.words.size());
    vnatural vq = q;
    uint64_t r = __div(a, b, vq);
    q.words.downsize(vq.size);
    return r;
}

constexpr void __abs_div(cnatural a, cnatural b, integer& q, integer& r) {
    Check(b.size != 0, "division by zero");
    if (b.size <= 1) {
        if (a.words != q.words.data())
            q.words.reset(a.size);
        vnatural vq = q;
        r = __div(a, b[0], vq);
        q.words.resize(vq.size);
        return;
    }
    if (__less(a, b)) {
        r.words.resize(a.size);
        std::copy(a.words, a.words + a.size, r.words.data());
        q.set_zero(); // update Q after R in case &A == &Q
        return;
    }

    // NOTE max word size of R is b.word.size + 1
    const int Q = a.size - b.size + 1; // b.size >= 2 here, so this is less than a.size
    if (a.words != q.words.data())
        q.words.reset(Q, /*initialize*/false);

#if 0
    if (Q == 1) {
        q.words.reset(1);
        const uint64_t w = __saturated_div(a, b);
        r.words.reset(a.size - 1, /*initialize*/false);
        std::copy(a.words + 1, a.words + a.size, r.words.data());
        vnatural vr = r;
        sub_product(vr, b, w);
        r.words.downsize(vr.size);
        q.words[0] = w;
        return;
    }
#endif

    r.set_zero();
    r.words.reset(a.size - Q, /*initialize*/false);
    std::copy(a.words + Q, a.words + a.size, r.words.data());

    for (int i = Q; i-- > 0;) {
        if (r.words.size() || a[i])
            r.words.insert_first_word(a[i]);

        const uint64_t w = __saturated_div(r, b);
        vnatural vr = r;
        Check(__sub_product(vr, b, w), "__saturated_div() overestimated the quotient");
        r.words.downsize(vr.size);
        q.words[i] = w;
    }
    q.words.downsize(Q);
    q.words.normalize();
}

constexpr void __abs_mod(cnatural a, cnatural b, integer& r) {
    Check(b.size != 0, "division by zero");
    if (b.size <= 1) {
        r = __mod(a, b[0]);
        return;
    }
    r.words.set_zero();
    r.words.reserve(b.size + 1); // the remainder never exceeds this, so the loop below never grows it
    for (auto i = a.size; i-- > 0;) {
        if (r.words.size() || a[i])
            r.words.insert_first_word(a[i]);
        const uint64_t q = __saturated_div(r, b);
        inatural ir = r;
        Check(__sub_product(ir, b, q), "mod() remainder went negative");
        r.words.downsize(ir.size);
    }
}

constexpr integer& __abs_shl(integer& a, int64_t b) {
    if (b > 0) {
        if (a.words.size() == 0)
            return a;
        auto word_shift = b / 64;
        auto bit_shift = b % 64;

        if (bit_shift) {
            uint64_t carry = 0;
            for (int i = 0; i < a.words.size(); ++i) {
                auto current = a.words[i];
                a.words[i] = (current << bit_shift) | carry;
                carry = current >> (64 - bit_shift);
            }
            if (carry)
                a.words.push_back(carry);
        }
        a.words.insert_first_n_words(word_shift);
        return a;
    }
    if (b < 0) {
        b = -b; // TODO undefined behavior for INT64_MIN, but that is huge amount of shift
        auto word_shift = b / 64;
        auto bit_shift = b % 64;

        if (word_shift >= a.words.size()) {
            a.words.set_zero();
            return a;
        }

        a.words.erase_first_n_words(word_shift);
        if (bit_shift != 0) {
            uint64_t carry = 0;
            for (auto idx = a.words.size(); idx-- > 0;) {
                auto current = a.words[idx];
                a.words[idx] = (current >> bit_shift) | carry;
                carry = (current << (64 - bit_shift));
            }
        }
        a.words.normalize();
    }
    return a;
}

constexpr integer& __abs_div(integer& a, std_int auto b) {
    Check(b > 0, "division of integer by zero or negative number");
    if constexpr (sizeof(b) > 8) {
        // the single word kernel below would keep only the low word of the divisor
        if (abs_unsigned(b) > UINT64_MAX) {
            integer quot, rem;
            __abs_div(static_cast<cnatural>(a), static_cast<cnatural>(integer(b)), quot, rem);
            a = std::move(quot);
            return a;
        }
    }
    __abs_div(a, static_cast<uint64_t>(b), a);
    return a;
}

constexpr std::string __abs_stre(const integer& a) {
    std::string s = "[";
    for (int i = 0; i < a.words.size(); i++)
        s += std::format(" {}", a.words[i]);
    s += " ]";
    return s;
}

constexpr int integer::__abs_str_buffer(char* buffer, int buffer_size, unsigned base, const bool upper) const {
    // emits the digits of the magnitude only; the caller writes any sign
    Check(base >= 2, "str() with base less than 2");
    Check(base <= 36, "str() with base greater than 36");
    char* p = buffer;
    const char* end = buffer + buffer_size;

    const char A = (upper ? 'A' : 'a') - 10;
    if (words.size() == 0) {
        Check(p < end, "buffer too small");
        *p++ = '0';
    } else if (words.size() == 1) {
        uint64_t a = low_word();
        if (base == 10 && buffer_size >= 20 * words.size()) {
            while (a) {
                *p++ = '0' + int(a % 10);
                a /= 10;
            }
        } else if (base == 16 && buffer_size >= 16 * words.size()) {
            uint64_t w = low_word();
            while (w) {
                int c = w % 16;
                *p++ = (c < 10) ? ('0' + c) : (A + c);
                w /= 16;
            }
        } else {
            while (a) {
                const int c = a % base;
                Check(p < end, "buffer too small");
                *p++ = (c < 10) ? ('0' + c) : (A + c);
                a /= base;
            }
        }
    } else {
        if (base == 10 && buffer_size >= 20 * words.size()) {
            // divide by 10**19 (the largest power of 10 that fits into uint64_t) and split
            // the remainder with 64-bit arithmetic: 19x fewer full precision divisions
            constexpr uint64_t CHUNK = 10'000'000'000'000'000'000ull;
            integer a = abs(*this);
            while (a) {
                uint64_t rem = __abs_div(a, CHUNK, /*out*/a);
                if (a) {
                    for (int i = 0; i < 19; i++) {
                        *p++ = '0' + int(rem % 10);
                        rem /= 10;
                    }
                } else {
                    while (rem) {
                        *p++ = '0' + int(rem % 10);
                        rem /= 10;
                    }
                }
            }
        } else if (base == 16 && buffer_size >= 16 * words.size()) {
            for (int i = 0; i < words.size(); i++) {
                auto w = words[i];
                if (i == words.size() - 1) {
                    while (w) {
                        int c = w % 16;
                        *p++ = (c < 10) ? ('0' + c) : (A + c);
                        w /= 16;
                    }
                } else {
                    for (int j = 0; j < 16; j++) {
                        int c = w % 16;
                        *p++ = (c < 10) ? ('0' + c) : (A + c);
                        w /= 16;
                    }
                }
            }
        } else {
            // same chunking as above, with the largest power of base that fits into uint64_t
            uint64_t chunk = base;
            int chunk_digits = 1;
            while (chunk <= UINT64_MAX / base) {
                chunk *= base;
                chunk_digits += 1;
            }
            integer n = abs(*this);
            while (n) {
                uint64_t rem = __abs_div(n, chunk, /*out*/n);
                for (int i = 0; (n && i < chunk_digits) || (!n && rem); i++) {
                    const int c = rem % base;
                    Check(p < end, "buffer too small");
                    *p++ = (c < 10) ? ('0' + c) : (A + c);
                    rem /= base;
                }
            }
        }
    }
    std::reverse(buffer, p);
    return p - buffer;
}

constexpr int integer::__abs_str_size_upper_bound_member(unsigned base) const {
    if (words.size() == 0)
        return 1;
    int m;
    switch (base) {
    case 2: m = 64; break;
    case 3: m = 41; break;
    case 4: m = 32; break;
    case 5: m = 28; break;
    case 6: m = 25; break;
    case 7: m = 23; break;
    case 8: m = 22; break;
    case 9: m = 21; break;
    case 10: m = 20; break;
    case 11: m = 19; break;
    case 12: m = 18; break;
    case 13: m = 18; break;
    case 14: m = 17; break;
    case 15: m = 17; break;
    case 16: m = 16; break;
    default: m = 16; // enough for any base above 16, avoids a dependency on log_upper()
    }
    return words.size() * m;
}

constexpr void integer::__abs_parse(std::string_view s, unsigned base) {
    const char* p = s.data();
    const char* end = s.data() + s.size();
    Check(p < end, "expecting digit instead of end of string");

    uint64_t acc = 0;
    unsigned count = 0;
    if (base == 10) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            Check('0' <= c && c <= '9', "expecting 0-9 for base 10");
            acc = acc * 10 + c - '0';
            count += 1;
            if (count == 19) {
                const uint64_t m = 10'000'000'000'000'000'000ull;
                __abs_mul_add(m, acc);
                acc = 0;
                count = 0;
            }
        }
        if (count)
            __abs_mul_add(pow(10ull, count), acc);
        words.normalize();
        return;
    }

    // Appends the low `bits` bits of acc. The chunk widths below are not all whole words, so the
    // shift has to be by bits: insert_first_word() would shift a 63 bit octal chunk by 64.
    const auto append = [this](uint64_t acc, unsigned bits) {
        *this <<= bits;
        if (words.size() == 0)
            words.push_back(acc);
        else
            words[0] |= acc;
    };

    if (base == 2) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            Check(c == '0' || c == '1', "expecting 0-1 for base 2");
            acc = acc * 2 + c - '0';
            count += 1;
            if (count == 64) {
                append(acc, count);
                acc = 0;
                count = 0;
            }
        }
    } else if (base == 8) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            Check('0' <= c && c <= '7', "expecting 0-7 for base 8");
            acc = acc * 8 + c - '0';
            count += 3;
            if (count == 63) {
                append(acc, count);
                acc = 0;
                count = 0;
            }
        }
    } else if (base == 16) {
        while (p < end) {
            if (*p == '\'') {
                p++;
                continue;
            }
            char c = *p++;
            int d;
            if ('0' <= c && c <= '9')
                d = c - '0';
            else if ('a' <= c && c <= 'f')
                d = c - 'a' + 10;
            else if ('A' <= c && c <= 'F')
                d = c - 'A' + 10;
            else
                Fail("expecting 0-9 or A-F for base 16");
            acc = acc * 16 + d;
            count += 4;
            if (count == 64) {
                append(acc, count);
                acc = 0;
                count = 0;
            }
        }
    } else
        Fail("unsupported base");
    if (count)
        append(acc, count);
    words.normalize();
}


constexpr integer abs(const integer& a) { integer m; m.words = a.words; m.words.set_negative(false); return m; }

// Borrows the magnitude for work that has to happen in place. The value is exactly
// one integer_backend, so the swap is a pointer exchange and not a copy; the sign is restored
// when the scope ends, and a result of zero comes back positive.

static_assert(sizeof(integer) == 16);

}

template <>
struct std::formatter<algebra::integer, char> {
    int width = 0;
    char fill = ' ';
    char align = '>';
    unsigned base = 10;
    bool upper = false;

    constexpr auto parse(auto& ctx) {
        auto it = ctx.begin();
        if (it != ctx.end() && it + 1 != ctx.end() && (it[1] == '>' || it[1] == '<' || it[1] == '^')) {
            fill = *it++;
            align = *it++;
        }
        while (it != ctx.end() && '0' <= *it && *it <= '9')
            width = width * 10 + *it++ - '0';
        if (it != ctx.end()) {
            if (*it == 'b' || *it == 'B') {
                base = 2;
                it++;
            } else if (*it == 'o') {
                base = 8;
                it++;
            } else if (*it == 'd') {
                base = 10;
                it++;
            } else if (*it == 'x' || *it == 'X') {
                base = 16;
                upper = (*it++ == 'X');
            }
        }
        algebra::Check(it != ctx.end() && *it == '}', "Invalid format specifier for integer.");
        return it;
    }

    // shared by std::formatter<integer>
    constexpr auto format_padded(const auto& a, auto& ctx) const {
        auto it = ctx.out();
        const int bound = a.str_size_upper_bound(base);

        char c_array[100];
        std::string str;
        char* buffer;
        if (bound <= 100) {
            buffer = c_array;
        } else {
            str.resize(bound);
            buffer = str.data();
        }

        const int n = a.str(buffer, bound, base, upper);
        int pre = 0;
        int post = 0;
        if (width > n) {
            if (align == '>')
                pre = width - n;
            if (align == '<')
                post = width - n;
            if (align == '^') {
                pre = (width - n) / 2;
                post = width - pre - n;
            }
        }
        for (int i = 0; i < pre; i++)
            *it++ = fill;
        for (int i = 0; i < n; i++)
            *it++ = buffer[i];
        for (int i = 0; i < post; i++)
            *it++ = fill;
        return it;
    }

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
