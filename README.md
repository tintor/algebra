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
    std::print("{} | {:.2f}\n", a, a); // prints 15/2 | 7.50
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
- `integer` / `rational` / `real<>` / `decimal` classes behave similarly to built-in `int` and `float` types (except for overflow)
- no heap allocation for integer values in `[-UINT64, UINT64]` range
- all types cast to any built-in integer and floating point type, and construct from any built-in
  integer; `rational` and `real<>` also construct from `float` and `double` exactly, while an
  `integer` is built from a floating point value with `round_to_zero()`
- no silent overflow / failures (std::runtime_error is thrown)
- output using `std::format()` / `std::print()` / `std::ostream` / `.str()`
- `real` allows more compact and efficient representation than `rational`, but requires rounding
- `real<2>` is similar to built-in `float` and `double`, but with arbitrary long mantissa, and 32-bit exponent
- `decimal` alias for `real<10>`
- `sizeof(integer)` is 16 bytes and `sizeof(rational)` is 32 bytes, while `std::vector<>` is 24 bytes

### Limitations
- multiplication and division currently use `O(N^2)` algorithms where N is number of 64-bit words used
  (`mul_karatsuba()` and `divide_bz()` are available, but are not used by the operators yet)
- the boolean and buffer operations on 2d regions are quadratic in the number of edges
- `real<Base>` division is not exact: it rounds to a fixed number of digits

### Headers

| header | contents |
| --- | --- |
| `algebra/algebra.h` | includes everything below |
| `algebra/integer.h` | `integer` and functions on it (also pulls in `integer_class.h`) |
| `algebra/rational.h` | `rational` and functions on it (also pulls in `rational_class.h`) |
| `algebra/real.h` | `real<Base>`, `decimal` (also pulls in `real_class.h`) |
| `algebra/xrational.h` | `xrational`: `rational * sqrt(integer)` |
| `algebra/expr.h` | symbolic expressions (`expr`, `expr_ptr`) |
| `algebra/vector.h` | `Vec<D, T>` with `Vec2` / `Vec3` / `Vec4` aliases |
| `algebra/rational_vector.h` | `qvec2/3/4` and `xvec2/3/4` aliases plus mixed-type vector operators |
| `algebra/solve_linear.h` | small linear systems and determinants |
| `algebra/geometry.h` | `Line3`, `Plane3` and their intersections; pulls in the distance and intersection headers below |
| `algebra/point_segment_squared_distance.h` | point to segment squared distance in 3d |
| `algebra/segment_segment_squared_distance.h` | segment to segment squared distance in 3d |
| `algebra/segment_segment_intersection.h` | segment intersection in 2d, as points or as parameters |
| `algebra/polygon2.h` | `MultiPolygon2<T>`: a 2d region as rings plus a complement flag |
| `algebra/polygon2_boolean.h` | union, intersection, difference and symmetric difference of regions |
| `algebra/polygon2_buffer.h` | dilate, erode and buffer by a convex structuring element |
| `algebra/polygon2_arc.h` | `ArcPolygon2<T>`: the same, with circular arc edges |
| `algebra/polygon2_arc_boolean.h` | `ArcRegion<T>`: boolean combinations of arc regions, as a tree |
| `algebra/dual.h` | `dual<T>` dual numbers for forward mode automatic differentiation |
| `algebra/kernels.h`, `algebra/util.h`, `algebra/types.h` | low level word array kernels and 128-bit helpers (names starting with `__` are internal) |

---

## `class integer`

Arbitrary precision signed integer. The magnitude lives in `words` and the sign in the sign of its
word count, so a value of one word or less needs no heap allocation at all.

Overloaded operators:
- arithmetic 	`+` `-` `*` `/` `%` `+=` `-=` `*=` `/=` `%=` `++` `--`
- relational `<` `>` `<=` `>=` `==` `!=`
- shift `<<` `>>` `<<=` `>>=`
- bitwise `~`

Division by zero throws `std::runtime_error`, and so does a conversion to a built-in type that would
not fit.

#### `integer_backend integer::words`
- Allows low level access to the vector of individual words of this number.

#### `integer::integer()`
#### `integer::integer(std_int auto a)`
#### `integer::integer(integer&& o)`
#### `integer::integer(const integer& o)`
#### `integer::integer(std::initializer_list<uint64_t> a)`
- The words, least significant first.
#### `integer::integer(std::string_view s, unsigned base = 10)`
- Accepts a leading `-`, and `'` as a digit separator. Bases 2, 8, 10 and 16.
#### `integer::integer(const char* s, unsigned base = 10)`

#### `void integer::operator=(std_int auto a)`
#### `void integer::operator=(integer&& o)`
#### `void integer::operator=(const integer& o)`

#### `int integer::sign() const`
- Negative, zero or positive; the magnitude of the returned value is the word count.
#### `bool integer::is_negative() const`
#### `bool integer::is_even() const`
#### `bool integer::is_odd() const`
#### `bool integer::is_one() const`
#### `bool integer::is_zero() const`
#### `uint64_t integer::low_word() const`
- The least significant word, which is only meaningful when the value is not zero.
#### `void integer::set_zero()`
#### `void integer::negate()`
- Same as `a = -a`, but in place and without memory allocation.
#### `void integer::swap(integer& o)`

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
- Whether the value fits into that built-in type, which is what the corresponding cast requires.

#### `std::string integer::str(unsigned base = 10, bool upper = true) const`
#### `std::string integer::hex() const`
#### `int integer::str_size_upper_bound(unsigned base = 10) const`
#### `int integer::str(char* buffer, int buffer_size, unsigned base = 10, bool upper = true) const`
- Writes into a caller provided buffer and returns the number of characters written.

#### `size_t integer::popcount() const`
- Number of set bits in the two's complement representation.
#### `int integer::size_of() const`
- Number of bytes used by the words of this number.
#### `auto integer::num_bits() const`
#### `auto integer::num_trailing_zeros() const`
#### `bool integer::bit(int64_t i) const`
- Bit `i` of the magnitude. Note that this is not the two's complement bit that `popcount()` counts.

#### `uint64_t integer::mod2() const`
#### `uint64_t integer::mod3() const`
#### `uint64_t integer::mod4() const`
#### `uint64_t integer::mod5() const`
#### `uint64_t integer::mod6() const`
#### `uint64_t integer::mod7() const`
#### `uint64_t integer::mod8() const`
#### `uint64_t integer::mod9() const`
#### `uint64_t integer::mod10() const`
- Remainder modulo a small constant, without a division. Non-negative, also for a negative value.

## `class rational`

Always kept in lowest terms with a positive denominator.

#### `integer rational::num`
#### `integer rational::den`
#### `rational::rational()`
#### `rational::rational(integer a)`
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

For `integer`, `std::format()` supports a fill and alignment (`{:*>10}`, `{:>10}`, `{:^10}`), a width,
and a base: `b`/`B` for binary, `o` for octal, `d` for decimal, `x`/`X` for hexadecimal.

## `class real<int Base>`

`num * Base**exp`, with `decimal` as an alias for `real<10>`.

#### `integer real::num`
#### `int real::exp`
#### `real::real(I a, int exp = 0)`
#### `real::real(integer a, int exp = 0)`
#### `real::real(float a)`
#### `real::real(double a)`
#### `real::real(const rational& a)`
- Exact conversion; throws if the denominator is not a power of `Base`.
#### `static real real::round(const rational& a, int digits)`
- The nearest value with `digits` digits after the point, with halves going away from zero.
#### `void real::normalize()`
- Moves trailing factors of `Base` from `num` into `exp`.
#### `std::string real::str() const`

## `class xrational`

`rational * sqrt(integer)`. Closed under multiplication and division; addition requires
compatible roots.

#### `rational xrational::base`
#### `integer xrational::root`
- Must be positive. Not fully simplified: it can still contain square factors.
#### `xrational::xrational()`
#### `xrational::xrational(rational_like auto base)`
#### `xrational::xrational(rational base, integer root)`
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

## `class MultiPolygon2<T = rational>`

A region of the plane, as a set of rings plus a complement flag. `Ring2<T>` is a
`std::vector<Vec2<T>>`; the closing edge from `back()` to `front()` is implicit, so a ring never
repeats its first vertex.

Membership uses the nonzero winding rule, flipped by `complement`. An outer boundary winds counter
clockwise and a hole winds clockwise, so a hole cancels the shell containing it and nesting to any
depth works without tracking which ring is whose hole. The flag is what makes inversion exact and
free: the complement of a bounded region is unbounded and has no finite ring representation.

`T` needs exact arithmetic *and* division for the predicates to hold, which is why the default is
`rational` rather than `integer`.

#### `std::vector<Ring2<T>> MultiPolygon2::rings`
#### `bool MultiPolygon2::complement`
#### `MultiPolygon2()`
- The empty region. `~MultiPolygon2<T>()` is the whole plane.
#### `MultiPolygon2(Ring2<T> ring)`
#### `MultiPolygon2(std::vector<Ring2<T>> rings, bool complement = false)`
#### `bool MultiPolygon2::is_empty() const`
#### `bool MultiPolygon2::is_whole_plane() const`

Overloaded operators: `~` (complement), `|` `&` `-` `^` (boolean operations, in
`polygon2_boolean.h`), and `==`, which compares the rings structurally rather than as point sets.

## `class ArcPolygon2<T = rational>`

The same, with edges that are line segments or circular arcs. A ring is a `std::vector<ArcVertex<T>>`
(`ArcRing2<T>`), where each vertex carries the *bulge* of the edge leaving it:

#### `Vec2<T> ArcVertex::p`
#### `T ArcVertex::bulge`
- `tan(theta/4)` for the arc's included angle `theta`, and 0 for a straight edge. A positive bulge
  puts the arc on the left of the edge, a negative one on the right.

That choice is what keeps everything rational: a rational bulge with rational endpoints gives a
rational centre and squared radius, so no coordinate is ever irrational. A bulge cannot describe a
full circle, since `theta == 2*pi` needs `tan(pi/2)`; `circle_ring()` uses two half circle edges.

Members and operators mirror `MultiPolygon2`: `rings`, `complement`, `is_empty()`,
`is_whole_plane()`, `~` and `==`.

## `class ArcRegion<T = rational>`

A boolean combination of arc regions, kept as a tree and evaluated on demand. Arc regions have no
explicit boolean result: two arcs meet at `cx +- sqrt(r*r - dy*dy)`, which is not rational, and
cutting a further arc at such a point nests the radicals. Membership *is* exactly computable, so
`contains()` on a combination is the combination of `contains()` on its operands, and every leaf test
is the exact rational predicate from `polygon2_arc.h`.

That gives exact union, intersection, difference, symmetric difference and complement over arc
regions, closed under further combination, with no tolerance anywhere. What it does not give is a
ring list or an exact area for a combination -- writing those down is the step that needs the
irrational points.

#### `ArcRegion(ArcPolygon2<T> p)` / `ArcRegion(ArcRing2<T> ring)`
- A leaf.
#### `size_t ArcRegion::leaf_count() const`
- How many arc regions the combination rests on.

Overloaded operators: `|` `&` `-` `^` and `~`. Note that `~region` is a strict negation of
`contains()`, unlike `~polygon`, which flips the complement flag and so leaves the boundary belonging
to both sides.

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

### algebra/integer.h
#### `integer power_of_two(size_t e)`
#### `integer exp2(std_int auto exp)`
- `2**e`. `exp2()` throws for a negative exponent.
#### `void mul(const integer& a, const integer& b, integer& out)`
#### `void mul(integer& a, const integer& b)`
#### `void square(integer& a)`
- `a = a * a`, using half the multiplications of `mul(a, a, out)`.
#### `void mul_karatsuba(const integer& a, const integer& b, integer& q)`
#### `integer mul_karatsuba(const integer& a, const integer& b)`
- Sub-quadratic multiplication. `q` must not alias `a` or `b`.
#### `void add_product(integer& acc, const integer& a, const integer& b)`
#### `void sub_product(integer& acc, const integer& a, const integer& b)`
#### `void add_product(integer& acc, const integer& a, std_int auto c)`
#### `void sub_product(integer& acc, const integer& a, std_int auto c)`
- `acc += a * b` / `acc -= a * b` without memory allocation.
#### `void div(const integer& a, const integer& b, integer& quot, integer& rem)`
- `quot` and `rem` have to be different objects; either may alias `a` or `b`.
#### `T div(const integer& a, T b, integer& quot)`
- Returns the remainder, for a signed or unsigned built-in `b` of any width.
#### `void divide_bz(const integer& a, const integer& d, integer& q, integer& r)`
- Recursive (Burnikel-Ziegler) division; same result as `div()`.
#### `T operator%(const integer& a, T b)`
- Truncates towards zero, so the remainder carries the sign of the dividend. The result type is the
  widest one that holds it: signed for a signed divisor, and `integer` for an unsigned one, where a
  negative remainder fits neither the divisor's type nor its signed counterpart. Divisors of every
  width up to 128 bits have an overload.
#### `integer mod(const integer& a, const integer& b)`
#### `void mod(integer& a, const integer& b)`
#### `uint64_t mod(const integer& a, uint64_t b)`
#### `unsigned mod(const integer& a, uint32_t b)`
- All `mod()` overloads return a value in `[0, abs(b))`, unlike `operator%` which truncates towards
  zero. The `integer&` overload replaces its argument in place, so it is chosen for a non const
  lvalue: spell the operand `const` (or use the return value) to get the value form.
#### `integer abs(const integer& a)`
#### `bool abs_greater(const integer& a, const integer& b)`
- `abs(a) > abs(b)`, minimizing memory allocation.
#### `int signum(const integer& a)`
#### `void invert_bits(integer& a)`
#### `void complement(integer& a)`
- Bitwise complement of the magnitude, and its two's complement. Both reject a negative value.

#### `void uniform_sample_bits(const size_t n, auto& rng, integer& out)`
#### `integer uniform_sample_bits(const size_t n, auto& rng)`
- uniformly sample from `[0, (2**n)-1]`
#### `void uniform_sample(const integer& count, auto& rng, integer& out)`
#### `integer uniform_sample(const integer& count, auto& rng)`
- uniformly sample from `[0, count-1]`; `count` has to be positive
#### `integer uniform_sample(const integer& min, const integer& max, auto& rng)`
- uniformly sample from `[min, max]`

#### `integer pow(integer base, std_int auto exp)`
#### `integer pow(integer base, const integer& exp)`
#### `integer pow(integer base, std_int auto exp, integer result)`
- returns `result * (base ** exp)`
#### `integer gcd(integer a, integer b)`
#### `auto gcd(std_int auto a, std_int auto b)`
- Of the magnitudes, so the sign of either argument does not matter.
#### `integer lcm(const integer& a, const integer& b)`
- `abs(a * b) / gcd(a, b)`, with the sign of `a * b`.
#### `uint64_t isqrt(std_unsigned_int auto x)`
#### `integer isqrt(const integer& x)`
- Largest `q` with `q * q <= x`.
#### `integer isqrt2(const integer& x)`
#### `integer isqrt3(const integer& x)`
- Alternative `isqrt()` implementations, kept for benchmarking.
#### `integer isqrt_hardware(const integer& a)`
- Very fast, but only approximate for large values.
#### `integer iroot(const integer& a, uint32_t n)`
- Largest `q` with `q**n <= a`. Throws for `n == 0`, which is not a root.
#### `bool exact_sqrt(const integer& a, integer& b)`
- Sets `b` to `sqrt(a)` and returns true when `a` is a perfect square.
#### `void exact_sqrt(integer a, integer& whole, integer& root)`
- Factors `sqrt(a)` into `whole * sqrt(root)`, accumulating into already initialized arguments.
#### `bool is_possible_square(const integer& a)`
- Cheap filter that rejects ~98% of non-squares.
#### `bool is_power_of_two(const integer& a)`
#### `bool is_power_of_three(const integer& n)`
#### `std::pair<int, int> mod63_65(const integer& a)`
- `a % 63` and `a % 65`, in one pass and without a division.
#### `void round_to_zero(std::floating_point auto a, integer& b)`
- The value truncated towards zero. Throws for nan and infinity.

#### `bool is_prime(uint64_t a)`
- Deterministic Miller-Rabin.
#### `bool is_likely_prime(const integer& n, int rounds)`
- Miller-Rabin with the first `rounds` primes as bases (at most 40).
- It returns false if n is composite and returns true if n is probably prime.
- Higher value of `rounds` indicates more accuracy.
#### `std::vector<std::pair<uint64_t, int>> factorize(std_unsigned_int auto a)`
#### `std::vector<std::pair<integer, int>> factorize(integer a)`
- Prime factorization as (factor, exponent) pairs.
#### `uint64_t try_fermat_factorize(uint64_t n)`
- A divisor of `n` strictly between 1 and `n`, or 0 when Fermat's method does not find one quickly.
  A prime `n` reports 0, since the difference of squares it factors into is `n = a*a - b*b` with
  `a - b == 1`.

#### `void add_mod(integer& a, const integer& b, const integer& m)`
#### `void sub_mod(integer& a, const integer& b, const integer& m)`
#### `void mul_mod(const integer& a, const integer& b, const integer& m, integer& out)`
- assume that the operands are in `[0, m-1]` range
#### `void pow_mod(integer a, const integer& b, const integer& m, integer& out)`
#### `integer pow_mod(integer a, const integer& b, const integer& m)`
- returns `(a**b) % m`
#### `bool inverse_mod(const integer& a, const integer& m, integer& out)`
- returns x such that `(a * x) mod m == 1`, or false if such number doesn't exist
#### `void binominal(const integer& n, uint64_t k, integer& out)`
- Binomial coefficient (n over k).
#### `void binominal_mod(const integer& n, uint64_t k, const integer& m, integer& out)`
- The same coefficient reduced modulo `m`.
#### `uint64_t log_lower(const integer& n, uint64_t base)`
#### `uint64_t log_upper(const integer& n, uint64_t base)`
- The base has to be at least two.

#### `void simplify(integer& x, integer& y)`
#### `void simplify(integer& x, integer& y, integer& z)`
- Divides all arguments by their common divisor.
#### `bool less_ab_c(const integer& a, const integer& b, const integer& c)`
- returns `a * b < c` (cheaper than naive multiplication)
#### `bool less_a_bc(const integer& a, const integer& b, const integer& c)`
- returns `a < b * c` (cheaper than naive multiplication)
#### `bool less_ab_cd(const integer& a, const integer& b, const integer& c, const integer& d)`
- returns `a * b < c * d` (cheaper than naive multiplication)

### algebra/kernels.h
#### `uint128_t extract_u128(cwords a, int64_t e)`
- returns `static_cast<unsigned __int128>(a >> e)` without memory allocation
#### `uint64_t extract_u64(cwords a, int64_t e)`
- returns `static_cast<uint64_t>(a >> e)` without memory allocation
- `integer` converts to `cwords`, a read only view of its words.

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
- The part beyond the integer part, so `trunc(a) + fract(a) == a`. The sign follows `a`, the same way
  `std::modf` splits a floating point value.
#### `rational abs(rational a)`
#### `bool abs_greater(const rational& a, const rational& b)`
#### `integer round_to_nearest(const rational& a)`
- The nearest integer, with halves going away from zero.
#### `rational round(const rational& a, unsigned digits, unsigned base = 10)`
- The nearest multiple of `base**-digits`, with halves going away from zero.
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
- solves `A + B*x = 0`, returning `None` when there is no solution and `Any` when every `x` is one
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

### algebra/polygon2.h
#### `T signed_area2(const Ring2<T>& ring)`
- Twice the signed area, which stays integral when `T` is. Positive is counter clockwise.
#### `T signed_area(const Ring2<T>& ring)`
#### `T signed_area(const MultiPolygon2<T>& a)`
- Throws for an unbounded region, which has no finite area.
#### `bool is_ccw(const Ring2<T>& ring)`
#### `void reverse(Ring2<T>& ring)`
- Reverses the orientation of a ring, turning a shell into a hole and back.
#### `bool on_boundary(const Ring2<T>& ring, const Vec2<T>& p)`
#### `bool on_boundary(const MultiPolygon2<T>& a, const Vec2<T>& p)`
#### `int winding_number(const Ring2<T>& ring, const Vec2<T>& p)`
#### `int winding_number(const MultiPolygon2<T>& a, const Vec2<T>& p)`
- How many times the region wraps counter clockwise around `p`. Undefined on the boundary, so
  callers test `on_boundary()` first.
#### `bool contains(const MultiPolygon2<T>& a, const Vec2<T>& p)`
- Closed region membership: the boundary belongs to the region on either side of a complement, so it
  is contained by both `a` and `~a`.
#### `void bounding_box(const MultiPolygon2<T>& a, Vec2<T>& min, Vec2<T>& max)`
- The smallest axis aligned box containing every vertex. Throws for an unbounded or empty region.
#### `void simplify(MultiPolygon2<T>& a)`
- Drops repeated and collinear vertices, and rings that enclose no area (fewer than three vertices,
  or all vertices on one line). The interior is unchanged, but the boundary of a dropped sliver goes
  with it, so a point that was only on such a sliver stops being contained.

### algebra/polygon2_boolean.h
#### `MultiPolygon2<T> boolean_op(BoolOp op, const MultiPolygon2<T>& a, const MultiPolygon2<T>& b)`
- `BoolOp` is `Union`, `Intersection`, `Difference` or `SymmetricDifference`. The operators `|`, `&`,
  `-` and `^` call it.
- Exact for an exact `T`, with no epsilon anywhere: every edge is cut at every crossing, each
  fragment is classified by sampling a point a provably safe step off its midpoint, and the
  surviving fragments are stitched into rings oriented with the interior on their left.
- Cost is quadratic in the number of edges, for the cutting and for the classification.
- The result is not `simplify()`ed: cutting leaves collinear vertices where the inputs met.

### algebra/polygon2_buffer.h
#### `MultiPolygon2<T> buffer(const MultiPolygon2<T>& a, const T& size, Element element = square_element<T>)`
- Positive size grows, negative shrinks, zero is the identity.
- `element` maps a size to a convex structuring element and is any callable, so one that takes more
  than a size goes in as a lambda: `buffer(a, r, [](const T& s) { return polygon_element(s, 16); })`.
- Note what is *not* offered: buffering by a Euclidean distance. Moving an edge out by `r` needs the
  unit normal, i.e. `sqrt(dx*dx + dy*dy)`, and a round join needs a circular arc whose intersections
  with its neighbours are irrational as well. Neither is representable in `rational`, so the shape to
  buffer by is given explicitly instead and the result stays exact.
#### `MultiPolygon2<T> dilate(const MultiPolygon2<T>& a, const Ring2<T>& b)`
#### `MultiPolygon2<T> erode(const MultiPolygon2<T>& a, const Ring2<T>& b)`
- Minkowski sum and its dual with a convex `b` that contains the origin.
#### `Ring2<T> square_element(const T& r)`
- `max(|dx|, |dy|) <= r`, i.e. buffering in the Chebyshev metric.
#### `Ring2<T> diamond_element(const T& r)`
- `|dx| + |dy| <= r`, i.e. the Manhattan metric.
#### `Ring2<T> polygon_element(const T& r, int sides)`
- A convex polygon with `2*sides` vertices inscribed in the circle of radius `r`, with rational
  coordinates from the Pythagorean parametrisation. It is a subset of the disk, so it
  under-approximates a round buffer as closely as wanted.
#### `Ring2<T> reflect(Ring2<T> a)`
- Negates every vertex, giving the `-B` that erosion needs. Note that the ring comes back with the
  opposite orientation; `dilate()` builds convex hulls from the element's vertices and tests the
  origin against a winding number, so neither cares.
#### `Ring2<T> convex_hull(std::vector<Vec2<T>> p)`
- The counter clockwise convex hull, by monotone chain. Exact for an exact `T`.

### algebra/polygon2_arc.h
#### `Vec2<T> arc_midpoint(const Vec2<T>& a, const Vec2<T>& b, const T& bulge)`
#### `Vec2<T> arc_center(const Vec2<T>& a, const Vec2<T>& b, const T& bulge)`
#### `T arc_radius2(const Vec2<T>& a, const Vec2<T>& b, const T& bulge)`
- The midpoint, centre and *squared* radius of an arc. All rational; the radius itself is not.
  `arc_center()` and `arc_radius2()` require a non-zero bulge.
#### `bool on_arc(const Vec2<T>& a, const Vec2<T>& b, const T& bulge, const Vec2<T>& p)`
- Whether `p` lies on the arc (or straight edge) from `a` to `b`, endpoints included.
#### `Ring2<T> chord_ring(const ArcRing2<T>& ring)`
#### `MultiPolygon2<T> chord_polygon(const ArcPolygon2<T>& a)`
- The same vertices with every arc replaced by its chord.
#### `int winding_number(const ArcPolygon2<T>& a, const Vec2<T>& p)`
- The chord polygon's winding number, corrected by one per arc whose circular segment contains `p`.
  Undefined on the boundary and on the chord of an arc.
#### `bool on_boundary(const ArcPolygon2<T>& a, const Vec2<T>& p)`
#### `bool contains(const ArcPolygon2<T>& a, const Vec2<T>& p)`
- Closed membership, as for `MultiPolygon2`. A point on a chord is stepped off first, by a step that
  provably crosses nothing.
#### `T signed_area_chords(const ArcPolygon2<T>& a)`
- The area of the chord polygon. The true area differs by the circular segments, whose area is
  `r*r*(theta - sin theta)/2` and so is not representable, which is why there is no exact
  `signed_area()` for an arc region.
#### `ArcRing2<T> circle_ring(const Vec2<T>& center, const T& radius)`
- A closed circle, as the two half circle edges that a single bulge cannot express.

### algebra/polygon2_arc_boolean.h
#### `bool contains(const ArcRegion<T>& r, const Vec2<T>& p)`
- Evaluates the tree at `p`.
#### `void area_bounds(const ArcRegion<T>& r, const Vec2<T>& min, const Vec2<T>& max, int depth, T& lower, T& undecided)`
- Area by subdividing the given box, since the exact area involves `r*r*(theta - sin theta)/2`.
  Boxes that sample as fully inside or fully outside settle, the rest are split until `depth` runs
  out, and what still straddles the boundary is reported as `undecided`.
- The test per box is five sample points, not an exact containment test, so `lower` and `undecided`
  are a good estimate rather than a proven bound. The box has to contain the region.

### algebra/dual.h
#### `struct dual<T> { T real, dual; }`
- Dual numbers for forward mode automatic differentiation, with `+` `-` `*` `/` and
  `sqrt` `pow` `exp` `log` `sin` `cos` `tan` `atan` `abs`.

### Literals
`using namespace algebra::literals;`

| literal | type |
| --- | --- |
| `123_i` | `integer` |
| `1/2_q` | `rational` |
| `1.5_f` | `real<2>` |
| `1.1_d` | `decimal` |
| `2_x` | `xrational` |
| `3_e` | `expr_ptr` |
