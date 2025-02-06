#pragma once
#include "algebra/util.h"
#include <string>
#include <vector>

namespace algebra {

struct cnatural {
    const uint64_t* const words;
    const int size;

    uint64_t operator[](int i) const { return words[i]; }
    uint64_t back() const { return words[size - 1]; }
    int countl_zero() const { return std::countl_zero(back()); }
};

struct inatural {
    uint64_t* const words; // pointer address can't change, but words can
    int size;

    uint64_t operator[](int i) const { return words[i]; }
    uint64_t& operator[](int i) { return words[i]; }
    uint64_t back() { return words[size - 1]; }
    uint64_t& back() const { return words[size - 1]; }
    void normalize() {
        while (size > 0 && back() == 0)
            size -= 1;
    }
    operator cnatural() const { return {words, size}; }
};

struct vnatural : public inatural {
    const int capacity;

    void push_back(uint64_t a) {
        Check(size < capacity);
        words[size++] = a;
    }
};

template<typename T, int S>
class maybe_stack {
public:
    constexpr maybe_stack(int size) {
        _t = (size > S) ? new T[size] : _stack;
    }
    constexpr ~maybe_stack() {
        if (_t != _stack)
            delete[] _t;
    }
    constexpr operator T*() { return _t; }
private:
    T _stack[S];
    T* _t;
};

constexpr int64_t add_max_size(cnatural a, cnatural b) {
    if (a.size > b.size)
        return (a.back() == UINT64_MAX) ? a.size + 1 : a.size;
    if (a.size < b.size)
        return (b.back() == UINT64_MAX) ? b.size + 1 : b.size;
    return can_overflow_with_carry(a.back(), b.back()) ? a.size + 1 : a.size;
}

constexpr int64_t mul_max_size(cnatural a, cnatural b) {
    return (a.size && b.size) ? a.size + b.size - 1 + (127 - a.countl_zero() - b.countl_zero()) / 64 : 0;
}

constexpr int64_t div_max_size(cnatural a, cnatural b) {
    return (a.size >= b.size) ? a.size - b.size + (64 + b.countl_zero() - a.countl_zero()) / 64 : 0;
}

constexpr int64_t num_bits(cnatural a) { return a.size ? 64 * a.size - a.countl_zero() : 0; }

constexpr int64_t num_trailing_zeros(cnatural a) {
    for (int i = 0; i < a.size; i++)
        if (a[i])
            return std::countr_zero(a[i]) + i * 64;
    return 0;
}

constexpr bool is_power_of_two(uint64_t a) { return a && (a & (a - 1)) == 0; }

constexpr bool is_power_of_two(cnatural a) {
    if (a.size == 0)
        return false;
    auto e = a.words + a.size - 1;
    if (*e & (*e - 1))
        return false;
    while (--e >= a.words)
        if (*e)
            return false;
    return true;
}

constexpr bool __less(cnatural a, cnatural b) {
    if (a.size > b.size)
        return false;
    if (a.size < b.size)
        return true;
    for (auto i = a.size; i-- > 0;) {
        if (a[i] > b[i])
            return false;
        if (a[i] < b[i])
            return true;
    }
    return false;
}

// return A < B * c
constexpr bool __less_a_bc_scalar(cnatural a, cnatural b, uint64_t c) {
    if (c == 0 || a.size > b.size + 1)
        return false;
    if (a.size < b.size)
        return true;

    auto a_bits = num_bits(a);
    auto bc_bits = num_bits(b) + num_bits(c);
    if (bc_bits < a_bits)
        return false;
    if (a_bits < bc_bits - 1)
        return true;

    // A == B || A == B + 1
    int cmp = 0;
    uint128_t carry = 0;
    for (auto i = 0; i < b.size; i++) {
        carry += __mulq(b[i], c);
        uint64_t bc = carry;
        carry >>= 64;

        if (a[i] > bc) cmp = 1;
        if (a[i] < bc) cmp = -1;
    }

    const uint64_t aa = (a.size > b.size) ? a[b.size] : 0;
    const uint64_t bc = carry;
    return (aa == bc) ? (cmp < 0) : (aa < bc);
}

constexpr void __mul(cnatural a, cnatural b, vnatural& q, bool init = true);

// unsafe, because if ignores words beyond the first two
constexpr uint128_t __unsafe_u128(cnatural a) {
    uint128_t c = a[0];
    if (a.size > 1)
        c |= static_cast<uint128_t>(a[1]) << 64;
    return c;
}

// returns A < B * C
constexpr bool __less_a_bc(cnatural a, cnatural b, cnatural c) {
    if (b.size == 0 || c.size == 0)
        return false;
    if (a.size == 0)
        return true;
    if (b.size == 1 && b[0] == 1)
        return __less(a, c);
    if (c.size == 1 && c[0] == 1)
        return __less(a, b);
    if (a.size <= 2 && b.size == 1 && c.size == 1)
        return __unsafe_u128(a) < __mulq(b[0], c[0]);

    const auto aa = num_bits(a);
    const auto bc = num_bits(b) + num_bits(c);
    if (aa < bc - 1)
        return true;
    if (bc < aa)
        return false;

    maybe_stack<uint64_t, 1024 / 8> w(b.size + c.size);
    vnatural v {{w, 0}, b.size + c.size};
    __mul(b, c, v);
    return __less(a, v);
}

// returns A * B < C
constexpr bool __less_ab_c(cnatural a, cnatural b, cnatural c) {
    if (c.size == 0)
        return false;
    if (a.size == 0 || b.size == 0)
        return true;
    if (a.size == 1 && a[0] == 1)
        return __less(b, c);
    if (b.size == 1 && b[0] == 1)
        return __less(a, c);
    if (a.size == 1 && b.size == 1 && c.size <= 2)
        return __mulq(a[0], b[0]) < __unsafe_u128(c);

    const auto ab = num_bits(a) + num_bits(b);
    const auto cc = num_bits(c);
    if (ab < cc)
        return true;
    if (cc < ab - 1)
        return false;

    maybe_stack<uint64_t, 1024 / 8> w(a.size + b.size);
    vnatural v {{w, 0}, a.size + b.size};
    __mul(a, b, v);
    return __less(v, c);
}

constexpr bool __sub_product(inatural& a, cnatural b, cnatural c);

// returns A * B < C * D
constexpr bool __less_ab_cd(cnatural a, cnatural b, cnatural c, cnatural d) {
    if (c.size == 0 || d.size == 0)
        return false;
    if (a.size == 0 || b.size == 0)
        return true;
    if (a.size == 1 && a[0] == 1)
        return __less_a_bc(b, c, d);
    if (b.size == 1 && b[0] == 1)
        return __less_a_bc(a, c, d);
    if (c.size == 1 && c[0] == 1)
        return __less_ab_c(a, b, d);
    if (d.size == 1 && d[0] == 1)
        return __less_ab_c(a, b, c);
    if (a.size == 1 && b.size == 1 && c.size == 1 && d.size == 1)
        return __mulq(a[0], b[0]) < __mulq(c[0], d[0]);

    const auto ab = num_bits(a) + num_bits(b);
    const auto cd = num_bits(c) + num_bits(d);
    if (ab < cd - 1)
        return true;
    if (cd < ab - 1)
        return false;

    if (a.size + b.size <= c.size + d.size) {
        maybe_stack<uint64_t, 1024 / 8> w(a.size + b.size);
        vnatural v {{w, 0}, a.size + b.size};
        __mul(a, b, v);
        return !__sub_product(v, c, d);
    }

    maybe_stack<uint64_t, 1024 / 8> w(c.size + d.size);
    vnatural v {{w, 0}, c.size + d.size};
    __mul(c, d, v);
    return __sub_product(v, a, b) && v.size > 0;
}

// returns signum(A * B - C * D)
constexpr int __det_ab_cd(cnatural a, cnatural b, cnatural c, cnatural d) {
    if (a.size == 0 || b.size == 0)
        return (c.size == 0 || d.size == 0) ? 0 : -1;
    if (c.size == 0 || d.size == 0)
        return 1;
    // TODO
    /*if (a.size == 1 && a[0] == 1)
        return __less_a_bc(b, c, d);
    if (b.size == 1 && b[0] == 1)
        return __less_a_bc(a, c, d);
    if (c.size == 1 && c[0] == 1)
        return __less_ab_c(a, b, d);
    if (d.size == 1 && d[0] == 1)
        return __less_ab_c(a, b, c);*/
    if (a.size == 1 && b.size == 1 && c.size == 1 && d.size == 1) {
        const auto ab = __mulq(a[0], b[0]);
        const auto cd = __mulq(c[0], d[0]);
        if (ab < cd) return -1;
        if (ab > cd) return 1;
        return 0;
    }

    const auto ab = num_bits(a) + num_bits(b);
    const auto cd = num_bits(c) + num_bits(d);
    if (ab < cd - 1)
        return -1;
    if (cd < ab - 1)
        return 1;

    if (a.size + b.size <= c.size + d.size) {
        maybe_stack<uint64_t, 1024 / 8> w(a.size + b.size);
        vnatural v {{w, 0}, a.size + b.size};
        __mul(a, b, v);
        if (!__sub_product(v, c, d))
            return -1;
        return (v.size == 0) ? 0 : 1;
    }

    maybe_stack<uint64_t, 1024 / 8> w(c.size + d.size);
    vnatural v {{w, 0}, c.size + d.size};
    __mul(c, d, v);
    if (!__sub_product(v, a, b))
        return 1;
    return (v.size == 0) ? 0 : -1;
}

constexpr bool __equal(cnatural a, cnatural b) {
    if (a.size != b.size)
        return false;
    for (auto i = a.size; i-- > 0;)
        if (a[i] != b[i])
            return false;
    return true;
}

constexpr bool __increment_and_return_carry(inatural a) {
    for (int i = 0; i < a.size; ++i)
        if (++a[i])
            return false;
    return true;
}

constexpr void __increment(vnatural& a) {
    for (int i = 0; i < a.size; ++i)
        if (++a[i])
            return;
    a.push_back(1);
}

// assumes A != 0
constexpr void __decrement(inatural& a) {
    for (int i = 0; i < a.size; ++i)
        if (a[i]--)
            return;
    Check(a.size > 0);
    a.size--;
}

// returns ~a + 1 (except for 0, which returns 0)
constexpr void __complement(inatural a) {
    for (int i = 0; i < a.size; ++i) {
        a[i] = ~a[i];
        if (++a[i]) {
            i += 1;
            while (i < a.size) {
                a[i] = ~a[i];
                i += 1;
            }
            return;
        }
    }
}

constexpr int __normalized_size(cnatural a) {
    int s = a.size;
    while (s > 0 && a[s - 1] == 0)
        s -= 1;
    return s;
}

// supports q == a
// b != q
constexpr void __mul(cnatural a, cnatural b, vnatural& q, bool init) {
    Check(a.size + b.size <= q.capacity);
    q.size = a.size + b.size;
    if (init)
        for (int i = a.size; i < q.size; i++)
            q[i] = 0;
    for (int i = a.size; i-- > 0;) {
        const uint64_t w = a[i];
        if (w == 0)
            continue;

        uint64_t* qi = q.words + i;
        auto acc = __mulq(w, b[0]);
        *qi = acc;
        acc >>= 64;
        int j = 1;

        while (j < b.size) {
            acc += __mulq(w, b[j]);
            acc += qi[j];
            qi[j] = acc;
            acc >>= 64;
            j += 1;
        }

        if (acc) {
            acc += qi[j];
            qi[j] = acc;
            acc >>= 64;
            if (acc) {
                j += 1;
                acc += qi[j];
                qi[j] = acc;
            }
        }
    }
    q.normalize();
}

constexpr void __add(vnatural& a, const uint64_t b, int shift = 0) {
    if (b == 0)
        return;
    Check(a.capacity >= shift);
    while (a.size < shift)
        a[a.size++] = 0;
    uint128_t carry = b;
    for (int i = shift; i < a.size; ++i) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
        if (carry == 0)
            return;
    }
    a.push_back(carry);
}

constexpr uint64_t __add_and_return_carry(inatural a, const uint64_t b) {
    if (b == 0)
        return 0;
    uint128_t carry = b;
    for (int i = 0; i < a.size; ++i) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
        if (carry == 0)
            return 0;
    }
    return carry;
}

constexpr uint128_t __add_and_return_carry(inatural a, const uint128_t b) {
    if (b == 0)
        return 0;
    uint128_t carry = b;
    for (int i = 0; i < a.size; ++i) {
        uint128_t e = carry + a[i];
        a[i] = e;
        if (e >= carry)
            carry = e >> 64;
        else
            carry = (e >> 64) + UINT64_MAX + 1;
        if (carry == 0)
            return 0;
    }
    return carry;
}

constexpr uint64_t __add_and_return_carry(inatural a, cnatural b, int shift = 0) {
    if (b.size == 0)
        return 0;
    Check(a.size >= b.size + shift);
    uint128_t carry = 0;
    for (int i = 0; i < b.size; ++i) {
        carry += a[i + shift];
        carry += b[i];
        a[i + shift] = carry;
        carry >>= 64;
    }
    if (carry == 0)
        return 0;
    for (int i = b.size + shift; i < a.size; ++i) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
        if (carry == 0)
            return 0;
    }
    return carry;
}

constexpr void __add(vnatural& a, cnatural b, int shift = 0) {
    if (b.size == 0)
        return;
    Check(a.capacity >= b.size + shift);
    while (a.size < b.size + shift)
        a[a.size++] = 0;
    a.push_back(__add_and_return_carry(a, b, shift));
}

// assuming a >= b
constexpr void __sub(inatural& a, cnatural b) {
    int i = 0;
    for (; i < b.size; ++i) {
        if (a[i] < b[i])
            goto borrow;
        no_borrow:
        a[i] -= b[i];
    }
    goto end;

    for (; i < b.size; ++i) {
        a[i]--;
        if (~a[i] && a[i] >= b[i])
            goto no_borrow;
        borrow:
        a[i] -= b[i];
    }
    for (; i < a.size; ++i)
        if (a[i]--)
            break;

    end:
    a.normalize();
}

// assuming a >= b
constexpr void __sub(inatural& a, const uint64_t b) {
    const bool decrement = a[0] < b;
    a[0] -= b;
    if (!decrement) {
        if (a[0] == 0)
            a.size = 0;
        return;
    }
    int i = 0;
    while (true)
        if (a[++i]--)
            return;
    a.size--;
}

constexpr void __sub(inatural& a, const uint128_t b) {
    const bool subtract = a[0] < b;
    a[0] -= static_cast<uint64_t>(b);
    if (!subtract) {
        if (a[0] == 0)
            a.size = 0;
        return;
    }
    const bool decrement = a[1] < (b >> 64);
    a[1] -= static_cast<uint64_t>(b >> 64);
    if (!decrement) {
        if (a[1] == 0)
            a.size = 1;
        return;
    }
    int i = 1;
    while (true)
        if (a[++i]--)
            return;
    a.size--;
}

// A = A * b + c
constexpr uint64_t __mul_add_return_carry(inatural a, const uint64_t b, const uint64_t c) {
    uint128_t carry = c;
    for (int i = 0; i < a.size; ++i) {
        carry += __mulq(a[i], b);
        a[i] = carry;
        carry >>= 64;
    }
    return carry;
}

constexpr uint64_t __mod(cnatural a, const uint64_t b) {
    uint128_t acc = 0;
    for (int i = a.size; i-- > 0;) {
        acc <<= 64;
        acc |= a[i];
        acc = __divq_mod(acc, b);
    }
    return acc;
}

constexpr uint128_t __mod(cnatural a, const uint128_t b) {
    // compute m = (2**128) % b
    uint128_t m = UINT128_MAX % b;
    m += 1;
    if (m == b)
        m = 0;

    uint128_t res = 0;
    for (int i = 0; i < a.size; i += 2) {
        res = mul_mod(res, m, b);
        uint128_t acc = a[i];
        if (i + 1 < a.size)
            acc |= static_cast<uint128_t>(a[i + 1]) << 64;
        res = add_mod(res, acc % b, b);
    }
    return res;
}

constexpr int mod10(cnatural a) {
    uint64_t m = 0;
    int i = a.size;
    while (i > 0) {
        m = m * 6 + a[--i] % 10;
        m %= 10;
    }
    return m;
}

constexpr int KARATSUBA_LIMIT = 32;

// long mul  16384               takes 235/230ms*
// karatsuba 16384 with 8 limit  takes ****ms (fails, likely due to small W buffer size)
// karatsuba 16384 with 16 limit takes  33ms* (regresses size 16 compared to long mul)
// karatsuba 16384 with 32 limit takes  26ms*
// karatsuba 16384 with 64 limit takes  26ms*

constexpr void swap(cnatural& a, cnatural& b) {
    char buffer[sizeof(cnatural)];
    std::memcpy(buffer, &a, sizeof(cnatural));
    std::memcpy(&a, &b, sizeof(cnatural));
    std::memcpy(&b, buffer, sizeof(cnatural));
}

// assuming a.size >= b.size
constexpr void __mul_karatsuba_rec(cnatural a, cnatural b, vnatural& q, uint64_t* w, const uint64_t* we) {
    if (a.size < b.size)
        swap(a, b);
    if (a.size < KARATSUBA_LIMIT) {
        __mul(a, b, q);
        return;
    }

    const int m = (a.size + 1) / 2; // m >= 2

    // AA and BB are stored in Q, while R and P are stored in W
    cnatural a0 {a.words, m};
    cnatural a1 {a.words + m, a.size - m};

    if (b.size <= m) {
        vnatural r {{w, 0}, a1.size + b.size};
        w += r.capacity;
        Check(w <= we);
        __mul_karatsuba_rec(a1, b, r, w, we); // r = a1 * b
        __mul_karatsuba_rec(a0, b, q, w, we); // q = a0 * b
        __add(q, r, m); // q = a * b
    } else {
        vnatural aa {{q.words, m}, m + 1};
        std::copy(a.words, a.words + m, aa.words);
        __add(aa, a1); // aa += a1

        // Top of b0 might be 0 (un-normalized), __add() needs to be robust to that!
        cnatural b0 {b.words, m};
        cnatural b1 {b.words + m, b.size - m};

        vnatural bb {{q.words + aa.capacity, m}, q.capacity - aa.capacity};
        std::copy(b.words, b.words + m, bb.words);
        __add(bb, b1); // bb = b0 + b1

        vnatural r {{w, 0}, aa.size + bb.size + 1};
        w += r.capacity;
        Check(w <= we);
        __mul_karatsuba_rec(aa, bb, r, w, we); // r = aa * bb
        // TODO ^ How is this working? AA and BB are stored in Q and nested __mul_karatsuba_rec call will overwrite them? Tests are pasing with 256 words.

        vnatural p {{w, 0}, a1.size + b1.size};
        w += p.capacity;
        Check(w <= we);
        __mul_karatsuba_rec(a1, b1, p, w, we); // p = a1 * b1
        __mul_karatsuba_rec(a0, b0, q, w, we); // q = a0 * b0

        __sub(r, p);
        __sub(r, q);
        __add(q, p, m + m);
        __add(q, r, m);
    }
}

// supports a == q
constexpr uint64_t __div(cnatural a, uint64_t b, vnatural& q) {
    Check(a.size <= q.capacity);
    q.size = a.size;
    uint128_t acc = 0;
    for (auto i = a.size; i-- > 0;) {
        acc <<= 64;
        acc |= a[i];
        uint64_t r;
        __divq(acc, b, q[i], r);
        acc = r;
    }
    q.normalize();
    return acc;
}

// A += B * c
constexpr void __add_product(vnatural& a, cnatural b, const uint64_t c, int shift = 0) {
    while (a.size < b.size + shift)
        a.push_back(0);
    uint128_t carry = 0;
    for (int i = 0; i < b.size; i++) {
        carry += a[i + shift];
        carry += __mulq(b[i], c);
        a[i + shift] = carry;
        carry >>= 64;
    }
    if (carry == 0)
        return;
    for (int i = b.size + shift; i < a.size; i++) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
        if (carry == 0)
            return;
    }
    a.push_back(carry);
    carry >>= 64;
    if (carry)
        a.push_back(carry);
}

// A += B * C
constexpr void __add_product(vnatural& a, cnatural b, cnatural c) {
    if (c.size == 0)
        return;
    for (int i = b.size; i-- > 0;)
        if (b[i])
            __add_product(a, c, b[i], i);
}

// Assumes A >= B * c (returns false otherwise)
// A -= (B * c) << (shift * 64)
constexpr bool __sub_product(inatural& a, cnatural b, const uint64_t c, int shift = 0) {
    if (a.size < b.size + shift)
        return false;
    uint128_t carry = 0; // for bc
    for (int i = 0; i < b.size; i++) {
        carry += __mulq(b[i], c);
        uint64_t bc = carry;
        carry >>= 64;

        carry += a[i + shift] < bc;
        a[i + shift] -= bc;
    }
    for (int i = b.size + shift; carry; i++) {
        uint64_t bc = carry;
        carry >>= 64;

        if (i >= a.size)
            return false;
        carry += a[i] < bc;
        a[i] -= bc;
    }
    a.normalize();
    return true;
}

// Assumes A >= B * C (returns false otherwise)
// A -= B * C
constexpr bool __sub_product(inatural& a, cnatural b, cnatural c) {
    for (int i = b.size; i-- > 0; )
        if (b[i] != 0)
            if (!__sub_product(a, c, b[i], i))
                return false;
    return true;
}

constexpr std::string str(cnatural a) {
    if (a.size == 0)
        return "0";
    std::vector<uint64_t> aa;
    aa.resize(a.size);
    std::copy(a.words, a.words + a.size, aa.data());
    vnatural b {{aa.data(), a.size}, a.size};

    std::string s;
    while (b.size)
        s += '0' + __div(b, 10, b);
    std::reverse(s.begin(), s.end());
    return s;
}

// returns static_cast<ucent>((a >> e) & UINT128_MAX) - without memory allocation
constexpr uint128_t extract_u128(cnatural a, int64_t e) {
    Check(e >= 0);
    const auto i = e / 64;
    const auto b = e % 64;

    if (i >= a.size)
        return 0;
    if (b == 0) {
        uint128_t res = a[i];
        if (i + 1 >= a.size)
            return res;
        res |= static_cast<uint128_t>(a[i + 1]) << 64;
        return res;
    }

    uint128_t res = a[i] >> b;
    if (i + 1 >= a.size)
        return res;
    res |= static_cast<uint128_t>(a[i + 1]) << (64 - b);
    if (i + 2 >= a.size)
        return res;
    res |= static_cast<uint128_t>(a[i + 2]) << (128 - b);
    return res;
}

// returns (a >> e).words[0] - without memory allocation
constexpr uint64_t extract_u64(cnatural a, int64_t e) {
    Check(e >= 0);
    const auto word_shift = e / 64;
    const auto bit_shift = e % 64;

    if (word_shift >= a.size)
        return 0;
    if (bit_shift == 0)
        return a[word_shift];

    uint64_t res = a[word_shift] >> bit_shift;
    if (word_shift + 1 >= a.size)
        return res;
    res |= a[word_shift + 1] << (64 - bit_shift);
    return res;
}

// returns largest uint64_t Q such that A >= B * Q (assuming B != 0)
constexpr uint64_t __saturated_div(cnatural a, cnatural b) {
    if (__less(a, b))
        return 0;
    if (__equal(a, b))
        return 1;
    if (a.size == 1) // since a > b  ==>  B == 1
        return a[0] / b[0];
    if (a.size == 2)
        return std::min<uint128_t>(__unsafe_u128(a) / __unsafe_u128(b), UINT64_MAX);

    const auto e = num_bits(a) - 128;
    uint128_t q = extract_u128(a, e);

    const uint128_t bq = extract_u128(b, e);
    if (bq)
        q /= bq;

    if (q > UINT64_MAX)
        q = UINT64_MAX;

    const uint64_t c = q;
    // if (a >= b * c)
    if (!__less_a_bc_scalar(a, b, c))
        return c;
    if (!__less_a_bc_scalar(a, b, c - 1))
        return c - 1;
    return c - 2;
}

}
