#pragma once
#include "algebra/util.h"
#include <string>
#include <vector>

namespace algebra {

constexpr int64_t add_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    if (A > B)
        return (a[A - 1] == UINT64_MAX) ? A + 1 : A;
    if (A < B)
        return (b[B - 1] == UINT64_MAX) ? B + 1 : B;
    return can_overflow_with_carry(a[A - 1], b[A - 1]) ? A + 1 : A;
}

constexpr int64_t mul_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    return (A && B) ? A + B - 1 + (127 - std::countl_zero(a[A - 1]) - std::countl_zero(b[B - 1])) / 64 : 0;
}

constexpr int64_t div_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    return (A >= B) ? A - B + (64 + std::countl_zero(b[B - 1]) - std::countl_zero(a[A - 1])) / 64 : 0;
}

constexpr int64_t num_bits(const uint64_t* a, const size_t A) { return A ? 64 * A - std::countl_zero(a[A - 1]) : 0; }

constexpr int64_t num_trailing_zeros(const uint64_t* a, const int A) {
    for (int i = 0; i < A; i++)
        if (a[i])
            return std::countr_zero(a[i]) + i * 64;
    return 0;
}

constexpr bool is_power_of_two(const uint64_t* a, const int A) {
    if (A == 0)
        return false;
    auto e = a + A - 1;
    if (*e & (*e - 1))
        return false;
    while (--e >= a)
        if (*e)
            return false;
    return true;
}

constexpr bool __less(const uint64_t* a, const int A, const uint64_t* b, const int B) {
    if (A > B)
        return false;
    if (A < B)
        return true;
    for (auto i = A; i-- > 0;) {
        if (a[i] > b[i])
            return false;
        if (a[i] < b[i])
            return true;
    }
    return false;
}

// return A < B * c
constexpr bool __less_a_bc(const uint64_t* a, const int A, const uint64_t* b, const int B, uint64_t c) {
    if (c == 0 || A > B + 1)
        return false;
    if (A < B)
        return true;

    // A == B || A == B + 1
    int cmp = 0;
    uint128_t carry = 0;
    for (auto i = 0; i < B; i++) {
        carry += __mulq(b[i], c);
        uint64_t bc = carry;
        carry >>= 64;

        if (a[i] > bc) cmp = 1;
        if (a[i] < bc) cmp = -1;
    }

    const uint64_t aa = (A > B) ? a[B] : 0;
    const uint64_t bc = carry;
    return (aa == bc) ? (cmp < 0) : (aa < bc);
}

constexpr bool __equal(const uint64_t* a, const int A, const uint64_t* b, const int B) {
    if (A != B)
        return false;
    for (auto i = A; i-- > 0;)
        if (a[i] != b[i])
            return false;
    return true;
}

constexpr bool __increment_and_return_carry(uint64_t* a, const int A) {
    for (int i = 0; i < A; ++i)
        if (++a[i])
            return false;
    return true;
}

constexpr void __increment(uint64_t* a, int& A) {
    for (int i = 0; i < A; ++i)
        if (++a[i])
            return;
    a[A++] = 1;
}

// assumes A != 0
constexpr void __decrement(uint64_t* a, int& A) {
    for (int i = 0; i < A; ++i)
        if (a[i]--)
            return;
    A--;
}

// returns ~a + 1 (except for 0, which returns 0)
constexpr void __complement(uint64_t* a, const int A) {
    for (int i = 0; i < A; ++i) {
        a[i] = ~a[i];
        if (++a[i]) {
            i += 1;
            while (i < A) {
                a[i] = ~a[i];
                i += 1;
            }
            return;
        }
    }
}

// `q` needs to have capacity of at least A + B!
// supports q == a
// b != q
constexpr void __mul(const uint64_t* a, const int A, const uint64_t* b, const int B, uint64_t* q, int& Q, bool init = true) {
    // TODO if a != q it might be possible to swap a and b depending on their sizes! maybe have A be smaller than B

    Q = A + B;
    if (init)
        for (int i = A; i < Q; i++)
            q[i] = 0;
    for (int i = A; i-- > 0;) {
        const uint64_t w = a[i];
        if (w == 0)
            continue;

        uint64_t* qi = q + i;
        auto acc = __mulq(w, *b);
        *qi = acc;
        acc >>= 64;
        int j = 1;

        while (j < B) {
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
    while (Q && !q[Q - 1])
        Q -= 1;
}

constexpr void __add(uint64_t* a, int& A, const uint64_t b, int shift = 0) {
    if (b == 0)
        return;
    while (A < shift)
        a[A++] = 0;
    uint128_t carry = b;
    for (int i = shift; i < A; ++i) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
        if (carry == 0)
            return;
    }
    a[A++] = carry;
}

constexpr uint128_t __add_and_return_carry(uint64_t* a, const int A, const uint128_t b) {
    uint128_t carry = b;
    for (int i = 0; i < A; ++i) {
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

constexpr void __add(uint64_t* a, int& A, const uint64_t* b, const int B, int shift = 0) {
    if (B == 0)
        return;
    while (A < B + shift)
        a[A++] = 0;
    uint128_t carry = 0;
    for (int i = 0; i < B; ++i) {
        carry += a[i + shift];
        carry += b[i];
        a[i + shift] = carry;
        carry >>= 64;
    }
    if (carry == 0)
        return;
    for (int i = B + shift; i < A; ++i) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
        if (carry == 0)
            return;
    }
    a[A++] = carry;
}

// assuming a >= b
constexpr void __sub(uint64_t* a, int& A, const uint64_t* b, const int B) {
    int i = 0;
    for (; i < B; ++i) {
        if (a[i] < b[i])
            goto borrow;
        no_borrow:
        a[i] -= b[i];
    }
    goto end;

    for (; i < B; ++i) {
        a[i]--;
        if (~a[i] && a[i] >= b[i])
            goto no_borrow;
        borrow:
        a[i] -= b[i];
    }
    for (; i < A; ++i)
        if (a[i]--)
            break;

    end:
    while (A && !a[A - 1])
        A -= 1;
}

// assuming a >= b
constexpr void __sub(uint64_t* a, int& A, const uint64_t b) {
    const bool decrement = *a < b;
    *a -= b;
    if (!decrement) {
        if (*a == 0)
            A = 0;
        return;
    }
    int i = 0;
    while (true)
        if (a[++i]--)
            return;
    A--;
}

constexpr void __sub(uint64_t* a, int& A, const uint128_t b) {
    const bool subtract = *a < b;
    *a -= static_cast<uint64_t>(b);
    if (!subtract) {
        if (*a == 0)
            A = 0;
        return;
    }
    const bool decrement = a[1] < (b >> 64);
    a[1] -= static_cast<uint64_t>(b >> 64);
    if (!decrement) {
        if (a[1] == 0)
            A = 1;
        return;
    }
    int i = 1;
    while (true)
        if (a[++i]--)
            return;
    A--;
}

// A = A * b + c
constexpr uint64_t __mul_add_return_carry(uint64_t* a, const int A, uint64_t b, uint64_t c) {
    uint128_t carry = c;
    for (int i = 0; i < A; ++i) {
        carry += __mulq(a[i], b);
        a[i] = carry;
        carry >>= 64;
    }
    return carry;
}

constexpr uint64_t __mod(const uint64_t* a, const int A, const uint64_t b) {
    Check(b != 0, "division by zero");
    uint128_t acc = 0;
    for (int i = A; i-- > 0;) {
        acc <<= 64;
        acc |= a[i];
        acc = __divq_mod(acc, b);
    }
    return acc;
}

constexpr uint128_t __mod(const uint64_t* a, const int A, const uint128_t b) {
    // compute m = (2**128) % b
    uint128_t m = UINT128_MAX % b;
    m += 1;
    if (m == b)
        m = 0;

    uint128_t res = 0;
    for (int i = 0; i < A; i += 2) {
        res = mul_mod(res, m, b);
        uint128_t acc = a[i];
        if (i + 1 < A)
            acc |= static_cast<uint128_t>(a[i + 1]) << 64;
        res = add_mod(res, acc % b, b);
    }
    return res;
}

constexpr int mod10(const uint64_t* a, const int A) {
    uint64_t m = 0;
    int i = A;
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

// assuming a.size >= b.size
constexpr void __mul_karatsuba_rec(const uint64_t* a, int A, const uint64_t* b, int B, uint64_t* q, int& Q, uint64_t* w, const uint64_t* we) {
    if (A < B) {
        std::swap(a, b);
        std::swap(A, B);
    }
    if (A < KARATSUBA_LIMIT) {
        __mul(a, A, b, B, q, Q);
        return;
    }

    const int m = (A + 1) / 2; // m >= 2

    // AA and BB are stored in Q, while R and P are stored in W
    const uint64_t* a1 = a + m;
    const int A1 = A - m;

    int R;
    if (B <= m) {
        auto r = w;
        w += A1 + B;
        __mul_karatsuba_rec(a1, A1, b, B, r, R, w, we); // r = a1 * b
        __mul_karatsuba_rec(a, m, b, B, q, Q, w, we); // q = a0 * b
        __add(q, Q, r, R, m); // q = a * b
    } else {
        uint64_t* aa = q;
        std::copy(a, a + m, aa);
        int AA = m;
        __add(aa, AA, a1, A1); // aa = a0 + a1

        const uint64_t* b1 = b + m;
        const int B1 = B - m;

        uint64_t* bb = aa + AA;
        std::copy(b, b + m, bb);
        int BB = m;
        __add(bb, BB, b1, B1); // bb = b0 + b1

        auto r = w;
        w += AA + BB + 1;
        __mul_karatsuba_rec(aa, AA, bb, BB, r, R, w, we); // r = aa * bb
        // TODO ^ How is this working? AA and BB are stored in Q and nested __mul_karatsuba_rec call will overwrite them? Tests are pasing with 256 words.

        uint64_t* p = w;
        w += A1 + B1;
        int P;
        __mul_karatsuba_rec(a1, A1, b1, B1, p, P, w, we);
        __mul_karatsuba_rec(a, m, b, m, q, Q, w, we);

        __sub(r, R, p, P);
        __sub(r, R, q, Q);
        __add(q, Q, p, P, m + m);
        __add(q, Q, r, R, m);
    }
}

// supports a == q
constexpr uint64_t __div(const uint64_t* a, const int A, uint64_t b, uint64_t* q, int& Q) {
    Check(b != 0, "division by zero");
    Q = A;
    uint128_t acc = 0;
    for (auto i = A; i-- > 0;) {
        acc <<= 64;
        acc |= a[i];
        uint64_t r;
        __divq(acc, b, q[i], r);
        acc = r;
    }
    while (Q > 0 && q[Q - 1] == 0)
        Q -= 1;
    return acc;
}

// A += B * c
constexpr void __add_product(uint64_t* a, int& A, const uint64_t* b, const int B, const uint64_t c, int shift = 0) {
    while (A < B + shift)
        a[A++] = 0;
    uint128_t carry = 0;
    for (int i = 0; i < B; i++) {
        carry += a[i + shift];
        carry += __mulq(b[i], c);
        a[i + shift] = carry;
        carry >>= 64;
    }
    if (carry == 0)
        return;
    for (int i = B + shift; i < A; i++) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
        if (carry == 0)
            return;
    }
    a[A++] = carry;
    carry >>= 64;
    if (carry)
        a[A++] = carry;
}

// A += B * C
constexpr void __add_product(uint64_t* a, int& A, const uint64_t* b, const int B, const uint64_t* c, const int C) {
    if (C == 0)
        return;
    for (int i = B; i-- > 0;)
        if (b[i])
            __add_product(a, A, c, C, b[i], i);
}

// Assumes A >= B * c
// A -= B * c
constexpr int __sub_product(uint64_t* a, int A, const uint64_t* b, const int B, const uint64_t c) {
    uint128_t carry = 0; // for bc
    uint64_t borrow = 0; // for a
    for (int i = 0; i < B; i++) {
        carry += __mulq(b[i], c);
        uint64_t bc = carry;
        carry >>= 64;

        int128_t diff = (int128_t)a[i] - bc - borrow;
        if (diff < 0) {
            a[i] = diff + ((uint128_t)1 << 64);
            borrow = 1;
        } else {
            a[i] = diff;
            borrow = 0;
        }
    }
    for (int i = B; carry || borrow; i++) {
        uint64_t bc = carry;
        carry >>= 64;

        int128_t diff = (int128_t)a[i] - bc - borrow;
        if (diff < 0) {
            a[i] = diff + ((uint128_t)1 << 64);
            borrow = 1;
        } else {
            a[i] = diff;
            borrow = 0;
        }
    }
    while (A && !a[A - 1])
        A -= 1;
    return A;
}

// Assumes A >= B * c
// A -= B * c
constexpr int __sub_product(uint64_t* a, int A, const uint64_t* b, const int B, const uint64_t* c, const int C) {
    for (int i = 0; i < B; i++)
        if (b[i] != 0)
            A = i + __sub_product(a + i, A - i, c, C, b[i]);
    return A;
}

constexpr std::string str(uint64_t* a, int A) {
    if (A == 0)
        return "0";
    std::string s;
    while (A)
        s += '0' + __div(a, A, 10, a, A);
    std::reverse(s.begin(), s.end());
    return s;
}

constexpr std::string str(const uint64_t* a, int A) {
    if (A == 0)
        return "0";
    std::vector<uint64_t> aa;
    aa.resize(A);
    std::copy(a, a + A, aa.data());
    return str(aa.data(), A);
}


}
