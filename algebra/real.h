#pragma once
#include "algebra/real_class.h"

namespace algebra {

// returns result * (base ** exp)
template<int B>
real<B> pow(real<B> base, int64_t exp, real<B> result = 1) {
    if (exp < 0)
        return result / pow(base, -exp, real<B>(1));

    if (base == B) {
        static_assert(sizeof(result.exp) == 4);
        const int64_t e = static_cast<int64_t>(result.exp) + exp;
        if (e < std::numeric_limits<int>::min() || e > std::numeric_limits<int>::max())
            throw std::runtime_error("exp overflow");
        return {result.num, static_cast<int>(e)};
    }
    if (exp == 0)
        return result;
    if (exp == 1)
        return result * base;

    if (exp & 1)
        result *= base;
    exp >>= 1;
    while (exp) {
        base *= base;
        if (exp & 1)
            result *= base;
        exp >>= 1;
    }
    return result;
}

template<int B>
real<B> abs(real<B> a) {
    if (a.num.sign() < 0)
        a.num.negate();
    return a;
}

}
