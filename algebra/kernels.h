#pragma once
#include "algebra/util.h"
#include <string>
#include <vector>

namespace algebra {

constexpr size_t add_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    if (A > B)
        return (a[A - 1] == UINT64_MAX) ? A + 1 : A;
    if (A < B)
        return (b[B - 1] == UINT64_MAX) ? B + 1 : B;
    return can_overflow_with_carry(a[A - 1], b[A - 1]) ? A + 1 : A;
}

constexpr size_t mul_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    return (A && B) ? A + B - 1 + (127 - std::countl_zero(a[A - 1]) - std::countl_zero(b[B - 1])) / 64 : 0;
}

constexpr size_t div_max_size(const uint64_t* a, const size_t A, const uint64_t* b, const int B) {
    return (A >= B) ? A - B + (64 + std::countl_zero(b[B - 1]) - std::countl_zero(a[A - 1])) / 64 : 0;
}

constexpr size_t num_bits(const uint64_t* a, const size_t A) { return A ? 64 * A - std::countl_zero(a[A - 1]) : 0; }

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

constexpr bool __equal(const uint64_t* a, const int A, const uint64_t* b, const int B) {
    if (A != B)
        return false;
    for (auto i = A; i-- > 0;)
        if (a[i] != b[i])
            return false;
    return true;
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

// shift must be >= 0
constexpr void __add(uint64_t* a, int& A, const uint64_t* b, const int B, int shift = 0) {
    while (A < B + shift)
        a[A++] = 0;

    uint128_t acc = 0;
    for (int i = 0; i < B; ++i) {
        acc += a[shift + i];
        acc += b[i];
        a[shift + i] = acc;
        acc >>= 64;
    }
    for (int i = shift + B; i < A; ++i) {
        acc += a[i];
        a[i] = acc;
        acc >>= 64;
    }
    if (acc)
        a[A++] = acc;
}

// assuming a >= b
constexpr void __sub(uint64_t* a, int& A, const uint64_t* b, const int B) {
    uint64_t borrow = 0;
    for (int i = 0; i < B; ++i) {
        int128_t diff = (int128_t)a[i] - b[i] - borrow;
        if (diff < 0) {
            a[i] = diff + ((uint128_t)1 << 64);
            borrow = 1;
        } else {
            a[i] = diff;
            borrow = 0;
        }
    }
    for (int i = B; borrow && i < A; ++i) {
        int128_t diff = (int128_t)a[i] - borrow;
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
        __add(aa, AA, a1, A1, 0); // aa = a0 + a1

        const uint64_t* b1 = b + m;
        const int B1 = B - m;

        uint64_t* bb = aa + AA;
        int BB = m;
        std::copy(b, b + m, bb);
        __add(bb, BB, b1, B1, 0); // bb = b0 + b1

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
        __add(q, Q, p, P, m * 2);
        __add(q, Q, r, R, m);
    }
}

// `Q` must be initialized to buffer capacity
// `q` must be all zero
constexpr void __mini_mul(const uint64_t* a, const int A, const uint64_t* b, const int B, uint64_t* q, int& Q) {
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

// A += B * C
constexpr void __add_product(uint64_t* a, int& A, const uint64_t* b, const int B, const uint64_t* c, const int C) {
    for (int i = B; i-- > 0;) {
        const uint64_t w = b[i];
        if (w == 0)
            continue;

        uint128_t acc = 0;
        uint64_t* qi = a + i;
        int j = 0;
        while (j < C) {
            acc += __mulq(w, c[j]);
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
                acc >>= 64;
                Check(acc == 0);
            }
        }
    }
    while (A && !a[A - 1])
        A -= 1;
}

// A += B * c
constexpr void __add_product(uint64_t* a, int& A, const uint64_t* b, const int B, const uint64_t c) {
    uint128_t carry = 0;
    for (int i = 0; i < B; i++) {
        carry += a[i];
        carry += __mulq(b[i], c);
        a[i] = carry;
        carry >>= 64;
    }
    if (carry == 0)
        return;
    for (int i = B; i < A; i++) {
        carry += a[i];
        a[i] = carry;
        carry >>= 64;
    }
    if (carry)
        a[A++] = carry;
}

// Assumes A >= B * c
// A -= B * c
constexpr void __sub_product(uint64_t* a, int& A, const uint64_t* b, const int B, const uint64_t c) {
    Check(A >= B);
    uint128_t carry = 0; // for bc
    uint64_t borrow = 0; // for a
    for (int i = 0; i < B; ++i) {
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
    for (int i = B; borrow && i < A; ++i) {
        int128_t diff = (int128_t)a[i] - borrow;
        if (diff < 0) {
            a[i] = diff + ((uint128_t)1 << 64);
            borrow = 1;
        } else {
            a[i] = diff;
            borrow = 0;
        }
    }

    Check(!borrow);
    while (A && !a[A - 1])
        A -= 1;
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
