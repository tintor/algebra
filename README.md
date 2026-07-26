# Algebra for Modern C++

```
#include "algebra/algebra.h"
#include <print>
using namespace algebra;
using namespace algebra::literals;

int main(int argc, char* argv[]) {
    int e = 2;
    integer i = 7_i + e;
    rational r = 5/6_q;
    rational a = r * i;
    std::print("{} | {:.2f}\n", a, a); // prints 15/2 7.50
    std::print("{:.20}\n", sqrt(2_q,  8)); // prints 1.41421356237309504880

    decimal d = 1.1_d; // stored exactly (unlike float and double which can't represent this value)
    std::print("{}\n", d); // prints 1.1
    return 0;
}
```

### Features
- header-only and no dependencies
- full `constexpr` and `std::format()` support
- arbitrary precision and compact algebraic data types
- `natural` / `integer` / `rational` / `real<>` / `decimal` classes behave similarly to built-in `int` and `float` types (except for overflow)
- no heap allocation for integer values in `[-UINT64, UINT64]` range
- all types support casting to and from all built-in integer and floating point types
- no silent overflow / failures (std::runtime_error is thrown)
- output using `std::format()` / `std::print()` / `std::ostream` / `.str()`
- `real` allows more compact and efficient representation than `rational`, but requires rounding
- `real<2>` is similar to built-in `float` and `double`, but with arbitrary long mantissa, and 32-bit exponent
- `decimal` alias for `real<10>`
- `sizeof(natural)` is 16 bytes, `sizeof(integer)` is 32 bytes, `sizeof(rational)` is 64 bytes, while `std::vector<>` is 24 bytes

### Limitations
- multiplication and division currently use `O(N^2)` algorithms where N is number of 64-bit words used
  (`mul_karatsuba()` and `divide_bz()` are available, but are not used by the operators yet)
- `__int128` support requires GNU extensions (`-std=gnu++23`, or clang, which accepts it in either mode)

### Headers

| header | contents |
| --- | --- |
| `algebra/algebra.h` | includes everything below |
| `algebra/natural.h` | `natural` and functions on it (also pulls in `natural_class.h`) |
| `algebra/integer.h` | `integer` and functions on it (also pulls in `integer_class.h`) |
| `algebra/rational.h` | `rational` and functions on it (also pulls in `rational_class.h`) |
| `algebra/real.h` | `real<Base>`, `decimal` (also pulls in `real_class.h`) |
| `algebra/xrational.h` | `xrational`: `rational * sqrt(natural)` |
| `algebra/expr.h` | symbolic expressions (`expr`, `expr_ptr`) |
| `algebra/vector.h` | `Vec<D, T>` with `Vec2` / `Vec3` / `Vec4` aliases |
| `algebra/rational_vector.h` | `qvec2/3/4` and `xvec2/3/4` aliases plus mixed-type vector operators |
| `algebra/solve_linear.h` | small linear systems and determinants |
| `algebra/geometry.h` | `Line3`, `Plane3` and their intersections |
| `algebra/dual.h` | `dual<T>` dual numbers for forward mode automatic differentiation |
| `algebra/kernels.h`, `algebra/util.h` | low level word array kernels (names starting with `__` are internal) |

---

## `class natural`

Arbitrary precision non-negative integer.

Overloaded operators:
- arithmetic 	`+` `-` `*` `/` `%` `+=` `-=` `*=` `/=` `%=` `++` `--`
- relational `<` `>` `<=` `>=` `==` `!=`
- shift `<<` `>>` `<<=` `>>=`
- bitwise `~` `|` `&` `^` `|=` `&=` `^=`

Subtraction below zero, division by zero and assignment of a negative value all throw `std::runtime_error`.

#### `integer_backend natural::words`
- Allows low level access to vector of individual words of this number.
#### `natural::natural()`
- Initializes to `0` value.
#### `natural::natural(std::integral auto a)`
#### `natural::natural(natural&& o)`
#### `natural::natural(const natural& o)`
#### `natural::natural(std::string_view s, unsigned base = 10)`
- Supported bases are 2, 8, 10 and 16. `'` is allowed as a digit separator.
#### `natural::natural(const char* s, uint32_t base = 10)`
#### `void natural::swap(natural& o)`
#### `void natural::set_zero()`
#### `uint64_t natural::low_word() const`
- Least significant word, or 0 when the value is zero.
#### `size_t natural::num_trailing_zeros() const`
- Returns number of trailing zeros in binary representation.
#### `bool natural::is_even() const`
- Same as `(a & 1) == 0`, but avoids temporary allocation for result of `&`
#### `bool natural::is_odd() const`
- Same as `(a & 1) == 1`, but avoids temporary allocation for result of `&`
#### `bool natural::is_one() const`
#### `bool natural::is_uint8() const`
#### `bool natural::is_uint16() const`
#### `bool natural::is_uint32() const`
#### `bool natural::is_uint64() const`
#### `bool natural::is_uint128() const`
#### `void natural::mul_add(uint64_t a, uint64_t carry)`
- `*this = *this * a + carry`
#### `std::string natural::str(uint32_t base = 10, bool upper = true) const`
#### `std::string natural::hex() const`
- Same as `natural::str(16)`
#### `int natural::str_size_upper_bound(uint32_t base = 10) const`
#### `int natural::str(char* buffer, int buffer_size, uint32_t base = 10, bool upper = true) const`
- Writes into a caller provided buffer and returns the number of characters written.
#### `int64_t natural::num_bits() const`
#### `bool natural::bit(int64_t i) const`
#### `int64_t natural::popcount() const`
#### `int64_t natural::size_of() const`
- Number of bytes used by the words of this number.
#### `int mod2() const`
#### `int mod3() const`
#### `int mod4() const`
#### `int mod5() const`
#### `int mod6() const`
#### `int mod7() const`
#### `int mod8() const`
#### `int mod9() const`
#### `int mod10() const`
- Remainder modulo a small constant, without a division.

## `class integer`

Arbitrary precision signed integer. The sign is stored in `abs.words`.

#### `natural integer::abs`

#### `integer::integer()`
#### `integer::integer(std::integral auto a)`
#### `integer::integer(integer&& o)`
#### `integer::integer(natural&& o)`
#### `integer::integer(const integer& o)`
#### `integer::integer(const natural& o)`
#### `integer::integer(std::string_view s, unsigned base = 10)`
#### `integer::integer(const char* s, unsigned base = 10)`

#### `void integer::operator=(std::integral auto a)`
#### `void integer::operator=(integer&& o)`
#### `void integer::operator=(natural&& o)`
#### `void integer::operator=(const integer& o)`
#### `void integer::operator=(const natural& o)`

#### `int integer::sign() const`
- Negative, zero or positive; the magnitude of the returned value is the word count.
#### `bool integer::is_negative() const`
#### `bool integer::is_even() const`
#### `bool integer::is_odd() const`
#### `bool integer::is_one() const`
#### `bool integer::is_zero() const`
#### `natural integer::to_natural() const`
#### `bool integer::is_int8() const`
#### `bool integer::is_int16() const`
#### `bool integer::is_int32() const`
#### `bool integer::is_int64() const`
#### `bool integer::is_int128() const`
#### `bool integer::is_uint8() const`
#### `bool integer::is_uint16() const`
#### `bool integer::is_uint32() const`
#### `bool integer::is_uint64() const`
#### `bool integer::is_uint128() const`

#### `std::string integer::str(unsigned base = 10, bool upper = true) const`
#### `std::string integer::hex() const`
#### `int integer::str_size_upper_bound(unsigned base = 10) const`
#### `int integer::str(char* buffer, int buffer_size, unsigned base = 10, bool upper = true) const`

#### `void integer::negate()`
#### `size_t integer::popcount() const`
- Number of set bits in two's complement representation.
#### `int integer::size_of() const`
#### `auto integer::num_bits() const`
#### `auto integer::num_trailing_zeros() const`
#### `bool integer::bit(size_t i) const`
#### `void integer::swap(integer& o)`
#### `uint64_t integer::mod2() const`
#### `uint64_t integer::mod3() const`
#### `uint64_t integer::mod4() const`
#### `uint64_t integer::mod5() const`
- Non-negative remainder, also for negative values.

## `class rational`

Always kept in lowest terms with a positive denominator.

#### `integer rational::num`
#### `integer rational::den`
#### `rational::rational()`
#### `rational::rational(integer a)`
#### `rational::rational(natural a)`
#### `rational::rational(integer a, integer b)`
- Initializes rational as a/b, and simplifies by removing common divisor.
#### `static rational rational::normalized(integer num, integer den)`
- Same as `rational(num, den)`, but assuming they are already simplified.
#### `rational::rational(std::integral auto a)`
#### `rational::rational(std::integral auto a, std::integral auto b)`
#### `rational::rational(float x)`
#### `rational::rational(double x)`
- Exact conversion.
#### `rational::rational(std::string_view s)`
- Accepts `123`, `-1/2`, `1.25` and `1e-3` forms.
#### `rational::rational(const std::string& s)`
#### `rational::rational(const char* s)`

#### `void rational::simplify()`
- You can use `.simplify()` after directly modifying `.num` and `.den` fields, to remove common factors from them.
- It throws exception if `den` is zero.
- Note that `rational` is automatically simplified after all arithmetic operations.
#### `void rational::invert()`
- Swap `num` and `den` in-place. Throws exception if `num` is zero.
#### `void rational::negate()`
- Same as `a = -a`, but performed in-place without memory allocation.
#### `std::string rational::str() const`
#### `int rational::sign() const`
#### `bool rational::is_integer() const`
#### `bool rational::is_even() const`
#### `bool rational::is_odd() const`
#### `bool rational::is_negative() const`
#### `bool rational::is_zero() const`
#### `rational::operator float() const`
#### `rational::operator double() const`

`std::format()` supports `{:.N}` and `{:.Nf}`, which round to N digits after the decimal point.

## `class real<int Base>`

`num * Base**exp`, with `decimal` as an alias for `real<10>`.

#### `integer real::num`
#### `int real::exp`
#### `real::real(I a, int exp = 0)`
#### `real::real(integer a, int exp = 0)`
#### `real::real(natural a, int exp = 0)`
#### `real::real(float a)`
#### `real::real(double a)`
#### `real::real(const rational& a)`
- Exact conversion; throws if the denominator is not a power of `Base`.
#### `static real real::round(const rational& a, int digits)`
- Nearest value with `digits` digits after the point (truncated towards zero).
#### `void real::normalize()`
- Moves trailing factors of `Base` from `num` into `exp`.
#### `std::string real::str() const`

## `class xrational`

`rational * sqrt(natural)`. Closed under multiplication and division; addition requires
compatible roots.

#### `rational xrational::base`
#### `natural xrational::root`
- Must be positive. Not fully simplified: it can still contain square factors.
#### `xrational::xrational()`
#### `xrational::xrational(rational_like auto base)`
#### `xrational::xrational(rational base, natural root)`
#### `void xrational::simplify()`
#### `bool xrational::is_rational() const`
#### `bool xrational::is_zero() const`
#### `bool xrational::is_negative() const`
#### `void xrational::negate()`
#### `void xrational::invert()`

## `class expr` / `class expr_ptr`
- `expr_ptr` is an alias for `std::shared_ptr<expr>`.
- Node types: `expr_integer`, `expr_rational`, `expr_power`, `expr_sum`, `expr_product`,
  `expr_negation`, `expr_sin`, `expr_cos`, `expr_pi`, `expr_e`, `expr_var`.
- Constants: `ZERO_EXPR`, `ONE_EXPR`, `PI_EXPR`, `E_EXPR`.

Overloaded operators:
- arithmetic 	`+` `-` `*` `/`
- relational `<` `>` `<=` `>=` `==` `!=` (these compare *values*, by determining the sign of the difference)

---

# Functions

### algebra/util.h
#### `void Check(bool value, const char* message = ...)`
- Throws `std::runtime_error` with the source location when `value` is false.
#### `[[noreturn]] void Fail(const char* message)`
#### `int num_bits(std::unsigned_integral auto)`
#### `uint64_t pow(uint64_t base, unsigned exp)`
#### `uint128_t add_mod(uint128_t a, uint128_t b, uint128_t m)`
#### `uint128_t mul_mod(uint128_t a, uint128_t b, uint128_t m)`
#### `uint64_t pow_mod(uint64_t a, uint64_t n, uint64_t m)`
- All three assume the operands are already in `[0, m)`.

### algebra/natural.h
#### `natural power_of_two(size_t e)`
#### `void mul(const natural&, const natural&, natural& out)`
#### `void mul(natural& a, const natural& b)`
#### `void square(natural& a)`
- `a = a * a`, using half the multiplications of `mul(a, a, out)`.
#### `void mul_karatsuba(const natural& a, const natural& b, natural& q)`
#### `natural mul_karatsuba(const natural& a, const natural& b)`
- Sub-quadratic multiplication. `q` must not alias `a` or `b`.
#### `void add_product(natural& acc, const natural& a, const natural& b)`
#### `void add_product(natural& acc, const natural& a, uint64_t b)`
#### `void sub_product(natural& acc, const natural& a, const natural& b)`
#### `void sub_product(natural& acc, const natural& a, uint64_t b)`
- `acc += a * b` / `acc -= a * b` without memory allocation. Subtraction assumes `acc >= a * b`.
#### `uint64_t div(const natural& dividend, uint64_t divisor, natural& quotient)`
- Returns the remainder.
#### `void div(const natural& dividend, const natural& divisor, natural& quotient, natural& remainder)`
#### `void divide_bz(const natural& a, const natural& d, natural& q, natural& r)`
- Recursive (Burnikel-Ziegler) division; same result as `div()`.
#### `void mod(const natural& dividend, const natural& divisor, natural& remainder)`
#### `void mod(natural& dividend, const natural& divisor)`
#### `uint128_t extract_u128(cnatural a, int64_t e)`
- returns `static_cast<unsigned __int128>(a >> e)` without memory allocation
#### `uint64_t extract_u64(cnatural a, int64_t e)`
- returns `static_cast<uint64_t>(a >> e)` without memory allocation
#### `void uniform_sample_bits(const size_t n, auto& rng, natural& out)`
#### `natural uniform_sample_bits(const size_t n, auto& rng)`
- uniformly sample from `[0, (2**n)-1]`
#### `void uniform_sample(const natural& count, auto& rng, natural& out)`
- uniformly sample from [0, count-1]
#### `natural uniform_sample(const natural& count, auto& rng)`
#### `natural uniform_sample(const natural& min, const natural& max, auto& rng)`
#### `natural pow(natural base, std::integral auto exp)`
#### `natural pow(natural base, const natural& exp)`
#### `auto num_trailing_zeros(std::unsigned_integral auto a)`
#### `natural gcd(natural a, natural b)`
#### `auto gcd(std::integral auto a, std::integral auto b)`
#### `natural lcm(const natural& a, const natural& b)`
#### `uint64_t isqrt(uint64_t x)`
#### `natural isqrt(const natural& x)`
- Largest `q` with `q * q <= x`.
#### `natural isqrt2(const natural& x)`
#### `natural isqrt3(const natural& x)`
- Alternative `isqrt()` implementations, kept for benchmarking.
#### `natural isqrt_hardware(const natural& a)`
- Very fast, but only approximate for large values.
#### `natural iroot(const natural& a, uint32_t n)`
- Largest `q` with `q**n <= a`.
#### `bool exact_sqrt(const natural& a, natural& b)`
- Sets `b` to `sqrt(a)` and returns true when `a` is a perfect square.
#### `void exact_sqrt(natural a, natural& whole, natural& root)`
- Factors `sqrt(a)` into `whole * sqrt(root)`, accumulating into already initialized arguments.
#### `bool is_possible_square(const natural& a)`
- Cheap filter that rejects ~98% of non-squares.
#### `bool is_power_of_two(const natural& a)`
#### `bool is_power_of_three(natural a)`
#### `bool is_prime(uint64_t a)`
- Deterministic Miller-Rabin.
#### `bool is_likely_prime(const natural& n, int rounds)`
- Miller-Rabin with the first `rounds` primes as bases (at most 40).
- It returns false if n is composite and returns true if n is probably prime.
- Higher value of `rounds` indicates more accuracy.
#### `std::vector<std::pair<uint64_t, int>> factorize(uint64_t a)`
#### `std::vector<std::pair<natural, int>> factorize(natural a)`
- Prime factorization as (factor, exponent) pairs.
#### `uint64_t try_fermat_factorize(uint64_t n)`
- Returns a divisor of `n`, or 0 when Fermat's method does not find one quickly.
#### `void add_mod(natural& a, const natural& b, const natural& m)`
#### `void sub_mod(natural& a, const natural& b, const natural& m)`
#### `void mul_mod(const natural& a, const natural& b, const natural& m, natural& out)`
- assume that the operands are in `[0, m-1]` range
#### `void pow_mod(natural a, const natural& b, const natural& m, natural& out)`
#### `natural pow_mod(natural a, const natural& b, const natural& m)`
- returns `(a**b) % m`
#### `void binominal(const natural& n, uint64_t k, natural& out)`
- Binomial coefficient (n over k).
#### `uint64_t log_lower(natural a, uint64_t base)`
#### `uint64_t log_upper(natural a, uint64_t base)`
#### `void invert_bits(natural& a)`
#### `void complement(natural& a)`

### algebra/integer.h
#### `integer abs(integer a)`
#### `bool abs_greater(const integer& a, const integer& b)`
- `abs(a) > abs(b)`, minimizing memory allocation.
#### `int signum(const integer& a)`
#### `void add_product(integer& acc, const integer& a, const integer& b)`
#### `void sub_product(integer& acc, const integer& a, const integer& b)`
#### `void add_product(integer& acc, const integer& a, std::integral auto c)`
#### `void sub_product(integer& acc, const integer& a, std::integral auto c)`
#### `void div(const integer& a, const integer& b, integer& quot, integer& rem)`
#### `int64_t div(const integer& a, int64_t b, integer& quot)`
#### `integer mod(const integer& a, const integer& b)`
#### `void mod(integer& a, const integer& b)`
#### `uint64_t mod(const integer& a, uint64_t b)`
#### `unsigned mod(const integer& a, uint32_t b)`
- All `mod()` overloads return a value in `[0, abs(b))`, unlike `operator%` which truncates towards zero.
#### `integer uniform_sample(const integer& min, const integer& max, auto& rng)`
#### `integer exp2(std::integral auto exp)`
#### `integer pow(integer base, std::integral auto exp)`
#### `integer pow(integer base, const natural& exp)`
#### `integer pow(integer base, std::integral auto exp, integer result)`
- returns `result * (base ** exp)`
#### `bool is_power_of_two(const integer& a)`
#### `integer gcd(const integer& a, const integer& b)`
#### `bool inverse_mod(const natural& a, const natural& m, natural& out)`
- returns x such that `(a * x) mod m == 1`, or false if such number doesn't exist
#### `void binominal_mod(const natural& n, uint64_t k, const natural& m, natural& out)`
#### `void simplify(integer& x, integer& y)`
#### `void simplify(integer& x, integer& y, integer& z)`
- Divides all arguments by their common divisor.
#### `bool less_ab_c(const integer& a, const integer& b, const integer& c)`
- returns `a * b < c` (cheaper than naive multiplication)
#### `bool less_a_bc(const integer& a, const integer& b, const integer& c)`
- returns `a < b * c` (cheaper than naive multiplication)
#### `bool less_ab_cd(const integer& a, const integer& b, const integer& c, const integer& d)`
- returns `a * b < c * d` (cheaper than naive multiplication)

### algebra/rational.h
#### `int approx_log2(const rational& a)`
#### `rational sqrt(const integer& x, unsigned iterations)`
#### `rational sqrt(const rational& x, unsigned iterations)`
- Newton iteration; the number of correct digits roughly doubles per iteration.
#### `rational nth_root(const rational& base, const integer& exp, unsigned iterations)`
#### `rational pow(const rational& base, long exp)`
#### `void pow(const rational& base, const integer& exp, rational& out)`
#### `rational pow(const rational& base, const integer& exp)`
#### `rational pow(const rational& base, const rational& exp, unsigned iterations)`
#### `rational fract(const rational& a)`
- Fractional part, always non-negative.
#### `rational abs(rational a)`
#### `bool abs_greater(const rational& a, const rational& b)`
#### `rational round(const rational& a, unsigned digits, unsigned base = 10)`
- Truncates towards zero to `digits` digits.
#### `integer trunc(const rational& a)`
- round towards 0 to integer
#### `rational PI(unsigned n)`
- https://en.wikipedia.org/wiki/Chudnovsky_algorithm for computing PI
- `n` is both the number of series terms and the number of square root iterations.
#### `rational sin(rational x, unsigned n)`
#### `rational cos(rational x, unsigned n)`
#### `rational exp(rational x, unsigned n)`
- Taylor series with `n` terms.
#### `void simplify(rational& x, rational& y)`
#### `void simplify(rational& x, rational& y, rational& z)`
- Scales all arguments by the same factor, preserving their ratios.

### algebra/real.h
#### `rational to_rational(const real<B>& a)`
#### `integer shift<B>(const integer& a, std::integral auto exp)`
- `a * B**exp`, for `exp >= 0`.
#### `real<B> pow(real<B> base, int64_t exp, real<B> result = 1)`
- returns `result * (base ** exp)`
#### `real<B> abs(real<B> a)`

### algebra/xrational.h
#### `xrational sqr(const xrational& a)`
#### `xrational sqrt(const xrational& a)`
#### `xrational abs(xrational a)`
#### `xrational pow(const xrational& base, integer exp)`
#### `auto signum(const xrational& a)`

### algebra/expr.h
#### `expr_ptr make_integer(const integer& a)`
#### `expr_ptr make_rational(const rational& a)`
#### `expr_ptr make_sum(std::vector<expr_ptr> v)`
#### `expr_ptr make_product(std::vector<expr_ptr> v)`
#### `expr_ptr pow(expr_ptr a, const rational& b)`
#### `expr_ptr sin(expr_ptr a)`
#### `expr_ptr cos(expr_ptr a)`
#### `expr_ptr sqrt(expr_ptr a)`
#### `expr_ptr cbrt(expr_ptr a)`
#### `bool identical(expr_ptr a, expr_ptr b)`
- Structural equality, unlike `operator==` which compares values.
#### `std::optional<int> safe_sign(expr_ptr a)`
- `a->sign()`, or `nullopt` when the sign cannot be determined.
#### `std::optional<interval<rational>> bounds(expr_ptr a)`
- Lower and upper bound of the value, when they can be determined.
#### `bool is_integer/is_rational/is_power/is_sqrt/is_cbrt/is_sum/is_product/is_negation(expr_ptr)`
#### `integer_value / rational_value / power_base / power_exp / sum_values / product_values / negation_value`
- Accessors for the corresponding node type.

### algebra/vector.h
#### `Vec<D, T>` with `Vec2<T>`, `Vec3<T>`, `Vec4<T>` aliases
- Arithmetic operators, `==`, `dot()`, `dot2()`, `cross()`, `lerp()`, `abs()`, `is_zero()`,
  `argmax_abs()`, `div_colinear()`, `same_sign()`, `order()` / `strict_order()` / `loose_order()`,
  and swizzles such as `xy()`, `yz()`, `xzy()`.

### algebra/solve_linear.h
#### `T determinant(const Vec2<T>& a, const Vec2<T>& b)`
#### `T determinant(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c)`
#### `std::variant<None, T, Any> solve_linear(const Vec<D, T>& a, const Vec<D, T>& b)`
- solves `A + B*x = 0`
#### `bool solve_linear(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& c, T* s, T* t)`
- solves `A + sB + tC = 0`, false when there is no unique solution
#### `bool solve_linear(const Vec3<T>& a, ..., T* s, T* t, T* r)`
- solves `A + sB + tC + rD = 0`

### algebra/geometry.h
#### `struct Line3<T> { Vec3<T> orig, dir; }`
#### `struct Plane3<T> { Vec3<T> n; T d; T den; }`
- Plane equation is `f(x) = (n * x + d) / sqrt(den)`.
#### `bool are_parallel(const Vec3<T>& a, const Vec3<T>& b)`
#### `std::variant<None, Vec3<T>, Line3<T>> line_plane_intersection(const Line3<T>&, const Plane3<T>&)`
#### `std::variant<None, Line3<T>, Plane3<T>> plane_intersection(const Plane3<T>&, const Plane3<T>&)`
#### `std::variant<None, Vec3<T>, Line3<T>, Plane3<T>> plane_intersection(const Plane3<T>&, const Plane3<T>&, const Plane3<T>&)`
#### `T point_segment_squared_distance(const Vec3<T>& p, const Vec3<T>& a, const Vec3<T>& b)`
#### `T segment_segment_squared_distance(const Vec3<T>& pa, const Vec3<T>& pb, const Vec3<T>& qa, const Vec3<T>& qb)`
#### `segment_segment_intersection(a, b, c, d)`
- Returns `None`, a point, or a segment. `segment_segment_intersection_param()` returns the
  parameters instead of the points, and `segment_segment_intersects()` returns 0, 1 or 2.

### algebra/dual.h
#### `struct dual<T> { T real, dual; }`
- Dual numbers for forward mode automatic differentiation, with `+` `-` `*` `/` and
  `sqrt` `pow` `exp` `log` `sin` `cos` `tan` `atan` `abs`.

### Literals
`using namespace algebra::literals;`

| literal | type |
| --- | --- |
| `123_n` | `natural` |
| `123_i` | `integer` |
| `1/2_q` | `rational` |
| `1.5_f` | `real<2>` |
| `1.1_d` | `decimal` |
| `2_x` | `xrational` |
| `3_e` | `expr_ptr` |
