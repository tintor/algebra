# Algebra for Modern C++

```
#include "algebra/rational.h"
using namespace algebra;

int main(int argc, char* argv[]) {
    int e = 2;
    integer i = 7_i + e;
    rational r = 5/6_q;
    rational a = r * i;
    std::print("{} | {:.2f}\n", a, a); // prints 15/2 7.50
    std::print("{:.20}\n", sqrt(2_q,  8)); // prints 1.41421356237309504880
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
- output using `std::format()` / `std::print()` / `std::ostream` / `.str()`
- `real` allows more compact and efficient representation than `rational`, but requires rounding
- `real<2>` is similar to built-in `float` and `double`, but with arbitrary long mantissa, and 32-bit exponent
- `decimal` alias for `real<10>`
- `sizeof(integer)` is 16 bytes, `sizeof(rational)` is 32 bytes, while `std::vector<>` is 24 bytes

### Limitations
- multiplication and division currently use `O(N^2)` algorithms where N is number of 64-bit words used

### Classes
### `class natural`
Overloaded operators:
- arithmetic 	`+` `-` `*` `/` `%` `+=` `-=` `*=` `/=` `%=` `++` `--`
- relational `<` `>` `<=` `>=` `==` `!=`
- shift `<<` `>>` `<<=` `>>=`
- bitwise `~` `|` `&` `^` `|=` `&=` `^=`
#### `integer_backend natural::words`
- Allows low level access to vector of individual words of this number.
#### `natural::natural()`
- Initializes to `0` value.
#### `natural::natural(std::integral auto a)`
#### `natural::natural(natural&& o)`
#### `natural::natural(const natural& o)`
#### `natural::natural(std::string_view s, unsigned base = 10)`
#### `natural::natural(const char* s, uint32_t base = 10)`
#### `void natural::swap(natural& o)`
#### `size_t natural::num_trailing_zeros() const`
- Returns number of trailing zeros in binary representation.
#### `bool natural::is_even() const`
- Same as `(a & 1) == 0`, but avoids temporary allocation for result of `&`
#### `bool natural::is_odd() const`
- Same as `(a & 1) == 1`, but avoids temporary allocation for result of `&`
#### `void natural::mul_add(natural::word a, natural::word carry)`
#### `std::string natural::str(uint32_t base = 10, bool upper = true) const`
#### `std::string natural::hex() const`
- Same as `natual::str(16)`
#### `size_type natural::str_size_upper_bound(uint32_t base = 10) const`
#### `size_type natural::str(char* buffer, int buffer_size, uint32_t base = 10, bool upper = true) const`

#### `size_t natural::num_bits() const`
#### `bool natural::bit(size_t i) const`
#### `size_t natural::popcount() const`
#### `size_t natural::size_of() const`

### `class integer`
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

#### `size_type integer::sign() const`
#### `bool integer::is_negative() const`
#### `bool integer::is_even() const`
#### `bool integer::is_odd() const`
#### `bool integer::is_one() const`
#### `bool integer::is_zero() const`
#### `bool integer::is_int8() const`
#### `bool integer::is_int16() const`
#### `bool integer::is_int32() const`
#### `bool integer::is_int() const`
#### `bool integer::is_int64() const`
#### `bool integer::is_long() const`
#### `bool integer::is_int128() const`
#### `bool integer::is_cent() const`
#### `bool integer::is_uchar() const`
#### `bool integer::is_uint8() const`
#### `bool integer::is_ushort() const`
#### `bool integer::is_uint16() const`
#### `bool integer::is_uint() const`
#### `bool integer::is_uint32() const`
#### `bool integer::is_ulong() const`
#### `bool integer::is_uint64() const`
#### `bool integer::is_ucent() const`
#### `bool integer::is_uint128() const`

#### `std::string integer::str(unsigned base = 10, bool upper = true) const`
#### `std::string integer::hex() const`
#### `int integer::str_size_upper_bound(unsigned base = 10) const`
#### `int integer::str(char* buffer, int buffer_size, unsigned base = 10, bool upper = true) const`

#### `void integer::negate()`
#### `size_t integer::popcount() const`
#### `int integer::size_of() const`
#### `auto integer::num_bits() const`
#### `auto integer::num_trailing_zeros() const`
#### `void integer::swap(integer& o)`

### `class rational`
#### `integer rational::num`
#### `integer rational::den`
`rational()`
`rational(integer a)`
`rational(integer a, integer b)`
`rational(integer a, integer b, int)`
`rational(std::integral auto a)`
`rational(std::integral auto a, std::integral auto b)`
`rational(float x)`
`rational(double x)`
`rational(std::string_view s)`
`rational(const std::string& s)`
`rational(const char* s)`

#### `void rational::simplify()`
- You can use `.simplify()` after directly modifying `.num` and `.den` fields, to remove common factors from them.
- It throws exception if `den` is zero.
- Note that `rational` is automatically simplified after all arithmetic operations.
#### `void rational::invert()`
- Swap `num` and `den` in-place. Throws exception if `num` is zero.
#### `void rational::negate()`
- Same as `a = -a`, but performed in-place without memory allocation.
#### `std::string rational::str() const`
#### `size_type rational::sign() const`
#### `bool rational::is_integer() const`
#### `bool rational::is_even() const`
#### `bool rational::is_odd() const`

### `class real<int Base>`
### `class expr`
### `class expr_ptr`
- Alias for `std::shared_ptr<expr>`
Overloaded operators:
- arithmetic 	`+` `-` `*` `/` `%`
- relational `<` `>` `<=` `>=` `==` `!=`

### Functions

### algebra/natural.h
#### `int num_bits(std::unsigned_integeral auto)`
#### `void mul(const natural&, const natural&, natural&)`
#### `void add_product(natural& acc, const natural& a, const natural& b)`
#### `void sub_product(natural& acc, const natural& a, const natural& b)`
#### `uint64_t div(const natural& dividend, uint64_t divisor, natural& quotient)`
#### `unsigned __int128 extract_128bits(const natural& a, uint e)`
- returns `static_cast<unsigned __int128>(a >> e)` without memory allocation (except for result itself)
#### `uint64_t extract_64bits(const natural& a, uint e)`
- returns `static_cast<uint64_t>(a >> e)` without memory allocation
#### `void div(const natural& dividend, const natural& divisor, natural& quotient, natural& remainder)`
#### `natural pow(natural base, std::integral auto exp)`
#### `natural pow(natural base, const natural& _exp)`
#### `natural uniform_int(const natural& min, const natural& max, auto& rng)`
#### `auto num_trailing_zeros(std::unsigned_integral auto a)`
#### `natural gcd(natural a, natural b)`
#### `natural isqrt(const natural& x)`
#### `constexpr bool is_prime(const uint64_t a)`
#### `bool is_prime(const natural& a)`
#### `bool is_power_of_two(const natural& a)`
#### `uint64_t log_lower(natural a, uint64_t base)`
#### `uint64_t log_upper(natural a, uint64_t base)`

### algebra/integer.h
#### `void add_product(integer& acc, const integer& a, const integer& b)`
#### `void sub_product(integer& acc, const integer& a, const integer& b)`
#### `void div(const integer& a, const integer& b, integer& quot, integer& rem)`
#### `long div(const integer& a, long b, integer& quot)`
#### `uint64_t mod(const integer& a, uint64_t b)`
#### `integer abs(integer a)`
#### `integer uniform_int(const integer& min, const integer& max, auto& rng)`
#### `integer pow(integer base, std::integral auto exp)`
#### `integer pow(integer base, const natural& exp)`

### algebra/rational.h
#### `rational sqrt(const integer& x, unsigned iterations)`
#### `rational sqrt(const rational& x, unsigned iterations)`
#### `rational nth_root(const rational& base, const integer& exp, unsigned iterations)`
#### `rational pow(const rational& base, long exp)`
#### `void pow(const rational& base, const integer& exp, rational& out)`
#### `rational pow(const rational& base, const integer& exp)`
#### `rational pow(const rational& base, const rational& exp, unsigned iterations)`
#### `rational fract(const rational& a)`
#### `rational abs(rational a)`
#### `rational round(const rational& a, unsigned digits, unsigned base = 10)`
#### `integer trunc(const rational& a)`
#### `rational PI(unsigned n)`
- https://en.wikipedia.org/wiki/Chudnovsky_algorithm for computing PI

### algebra/real.h
#### `rational to_rational(const real& a)`

### algebra/expr.h
#### `expr_ptr make_integer(const integer& a)`
#### `expr_ptr make_rational(const rational& a)`
#### `expr_ptr pow(expr_ptr a, const rational& b)`
#### `expr_ptr sin(expr_ptr a)`
#### `expr_ptr cos(expr_ptr a)`
#### `expr_ptr sqrt(expr_ptr a)`
#### `expr_ptr cbrt(expr_ptr a)`
