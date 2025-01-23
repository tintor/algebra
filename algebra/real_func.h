#pragma once
#include "algebra/real.h"

namespace algebra {

template<int B>
real<B> pow(real<B> base, int64_t exp, real<B> result = 1) {
    if (exp < 0)
        return 1 / pow(base, -exp, result);

    if (base == B) {
        static_assert(sizeof(result.exp) == 4);
        if (static_cast<int>(static_cast<int64_t>(result.exp) + exp) != result.exp + static_cast<int>(exp))
            throw std::runtime_error("exp overflow");
        return {result.num, result.exp + static_cast<int>(exp)};
    }
    if (exp == 0)
        return 1;
    if (exp == 1)
        return base;

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
