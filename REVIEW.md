# Code review findings

Review of the whole repository (24 headers, ~11.5k LOC).
Items marked **[confirmed]** were reproduced by compiling and running a probe program
against the headers (g++ 15, `-std=gnu++23 -O1`).

Checkboxes track fix progress. One bug per commit: failing test first, then fix.

---

## 1. Wrong results / crashes — confirmed by running

- [x] **1.1 `natural -= uint64_t` destroys the value when the low word becomes 0.**
  `kernels.h:531-544` — `__sub(inatural&, uint64_t)` sets `a.size = 0` whenever `a[0] == 0`
  after the subtract, regardless of higher words. Same bug in the `uint128_t` overload
  (`:546`), which additionally reads `a[1]` when `size == 1`.
  `(2^64+5) - 5u` → `0`. **[confirmed]**

- [x] **1.2 `natural % uint128_t` is wrong for >= 3 words.**
  `kernels.h:589-605` — the Horner loop runs least-significant-first (`i = 0; i += 2`),
  but Horner requires most-significant-first, so word *i* gets weight `m^(k-i)` instead of `m^i`.
  `(2^200+12345) % (2^100+7)` → `4722366459672795283456`, correct is `12394`. **[confirmed]**

- [x] **1.3 `algebra::pow(uint64_t, unsigned)` returns 0 for base 2, exp >= 32.**
  `util.h:100`: `return 1 << exp;` — `int` shift (UB for exp >= 31).
  Should be `uint64_t(1) << exp`. `pow(2,40)` → `0`. **[confirmed]**

- [x] **1.4 `++integer` / `--integer` change the value by 2.**
  `integer_class.h:194-221` — `++abs` already increments, then
  `__increment_and_return_carry(*this)` increments the *same* storage again
  (`operator inatural()` returns `abs.words`). `++integer(5)` → 7. **[confirmed]**

- [x] **1.5 `mod(integer, integer)` infinite-recurses → segfault.**
  `integer_class.h:540-541`: `mod(a.abs, b.abs)` has no viable 2-arg `mod(natural,natural)`
  overload, so the `natural → integer` conversion selects the enclosing function itself. **[confirmed]**

- [x] **1.6 `inverse_mod` never terminates.**
  `integer.h:127`: `sub_product(r, q, new_r)` should be `sub_product(e, q, new_r)`
  (compare the correct `t` block above it); the following `r.swap(new_r); new_r.swap(e)`
  then restores the old `r`, so the Euclid state never advances. **[confirmed]**

- [x] **1.7 `real<2>(rational)` double-counts the numerator's trailing zeros.**
  `real_class.h:65-67` adds `z` to `exp`, then `normalize()` adds it again.
  `real<2>(rational(4))` → `num=1, exp=4` = **16**; `real<2>(4.0f)` likewise. **[confirmed]**

- [x] **1.8 `cos(rational,n)` and `exp(rational,n)` compute the wrong function.**
  `rational.h:223` starts the term accumulator at `a = x` instead of `1`, so it sums
  `x^(2k+1)/(2k)!`; `rational.h:236-246` omits the leading `1` and divides by `i` instead of
  `i+1`, computing `x*e^x`. Both are accidentally right at `x = 1`.
  `cos(2)` → `-1.832` (want `-0.4161`); `exp(2)` → `14.778` = `2e^2`. **[confirmed]**

- [x] **1.9 `mul_karatsuba` is wrong when the second operand is a power of two.**
  `natural_class.h:479`: `std::countr_zero(a.words[B - 1])` should be `b.words[B - 1]`. **[confirmed]**

- [x] **1.10 `mul_karatsuba` throws on operands with interior zero limbs.**
  Root cause `kernels.h:351-353`: `__mul` skips writing `q[i]` when `a[i] == 0` (`continue`),
  but the `init` loop only zeroes `q[a.size ..]`. Every `__mul` call from the Karatsuba base
  case writes into an uninitialized `W` slice. Existing tests only use dense/all-ones operands.
  `Check failed ... vnatural::push_back`. **[confirmed]**

- [x] **1.11 `divide_bz` produces the wrong quotient/remainder.**
  One concrete cause: `natural_class.h:970` captures `A = a.words.size()` *before* `a <<= shift`
  (`:979`), so if `a` grows a word, `end = min(i, A)` never copies the top word. **[confirmed]**

- [x] **1.12 `std::format("{:.N}", rational)` reads out of bounds when the value has fewer
  integer digits than `N`.** `rational_class.h:451-452`: `s[i + s.size() - r]` with `r > s.size()`;
  needs zero padding. `{:.3}` of `1/1000` → `'0.  1'`. **[confirmed]**

- [x] **1.13 `operator<(uint128_t, natural)` compares `b` with itself.**
  `natural_class.h:390`: `return static_cast<uint64_t>(b) < b;` — should be
  `static_cast<uint64_t>(a) < b`, and it throws for naturals > 64 bits. **[confirmed]**

- [x] **1.14 `__PI` overflows 32-bit unsigned.**
  `rational.h:180`: `545140134*a` with `unsigned a` wraps for `a >= 8`, corrupting every
  Chudnovsky term past ~digit 110. `545140134*8` → `66153776`. **[confirmed]**

- [x] **1.15 Header-only library can't be linked from two TUs.**
  Three non-`inline` definitions at namespace scope: `REAL_FRACT_DIGITS` (`real_class.h:251`),
  `__divide_2n1n` (`natural_class.h:954`), `try_fermat_factorize` (`natural.h:715`).
  `multiple definition of ...` x3. **[confirmed]**

- [x] **1.16 `solve_linear` contradicts its own doc comment.**
  For the documented `A + sB + tC = 0` it returns `-s` (2-D) and `-s`, `-r` (3-D). In-repo
  callers pass sign-flipped arguments and cancel it, so any new caller silently gets a sign
  error. Returns `s = -1` where `s = 1`. **[confirmed]**

- [x] **1.17 `rational(float/double)` returns the absolute value for negative input.**
  `rational_class.h:179-190` — `num` is assigned from the *signed* scaled mantissa and then
  negated again by `if (x < 0) num = -num;`. `rational(-1.5)` → `3/2`.
  Found while fixing 1.7; not covered by any test (all negative literals in the tests go
  through the string parser or unary minus). **[confirmed]**

- [x] **1.18 `natural::operator--` leaves an un-normalized value.**
  `kernels.h:313` — `__decrement()` returns as soon as a word can be decremented without a
  borrow, so `natural(1); --a;` has `size == 1` with word `0`: `a == 0u` is true while
  `bool(a)` is also true, and `2**64 - 1` keeps a leading zero word. Found while fixing 1.4.
  **[confirmed]**

- [x] **1.19 `__add(vnatural&, cnatural, shift)` pushes a zero carry word.**
  `kernels.h:455` called `a.push_back(carry)` unconditionally, so adding into a buffer whose
  size already equals its capacity threw `Check failed ... vnatural::push_back` even when there
  was no carry out. This is what made `mul_karatsuba` throw on sparse operands (the
  `mul_max_size` bound leaves no spare word). Found while fixing 1.10. **[confirmed]**

## 2. Bugs found by reading (not exercised by tests)

- [x] **2.1 Memory leak in `integer_backend`'s move assignment** — `integer_backend.h:95-104`
  overwrites `_words` without `delete[]`ing the current buffer. Every move-assign of a
  heap-sized `natural`/`integer` leaks.

- [x] **2.2 `natural` underflow is silent**, contradicting the README's "no silent overflow /
  failures":
  - `natural(3) - natural(5)` → `18446744073709551614` (`operator-=`/`__sub` have no
    precondition check) **[confirmed]**
  - `natural(3) - 5u` (`natural_class.h:327-329`) skips the check the 128-bit path performs,
    and `__sub`'s borrow loop then walks past `a.size` — an out-of-bounds write.
  - `natural_class.h:340`: `Check(a <= b, "natural can't be negative")` in
    `operator-(unsigned, natural)` has the comparison **backwards** (should be `a >= b`).

- [x] **2.3 Division by zero is a SIGFPE, not the documented exception.**
  In `__div` (`natural_class.h:818`), `mod(const natural&,...)` (`:910`) and
  `mod(natural&,...)` (`:929`), the `b.size <= 1` branch runs `__mod(a, b[0])` *before* the
  `Check(b.size != 0, "division by zero")` on the next line — so that `Check` is unreachable.

- [x] **2.4 `uint8_t(integer)` validates the wrong range** — `integer_class.h:98` calls
  `Check(is_uint16())`. `uint8_t(integer(300))` → `44`. **[confirmed]**

- [x] **2.5 `integer(string_view, base)` ignores `base`** — `integer_class.h:86` forwards only
  the string to `natural`, so `integer("ff", 16)` throws. **[confirmed]**

- [x] **2.6 `rational(int, int)` accepts a zero denominator** — `rational_class.h:29-48` never
  calls `simplify()`, so `rational(1,0)` yields `1/0` while `rational(integer(1), integer(0))`
  throws. **[confirmed]**

- [x] **2.7 `rational("-0.5")` loses the sign** — `rational_class.h:60`: `integer("-0")` is
  `+0`, so the later `if (num >= 0) num += frac` adds instead of subtracts. Returns `1/2`. **[confirmed]**

- [x] **2.8 `a /= a` and `a %= a` on `rational` are wrong** — `rational_class.h:330-335`
  mutates `a.num` before reading `b.num`. `rational(1,2) /= itself` → `1/2` instead of `1`.
  **[confirmed]** (`+=`/`-=` are safe only because the `a.den == b.den` fast path catches
  self-aliasing.)

- [x] **2.9 `pow(rational, negative exp)` hangs** — `rational.h:56-61`: `exp >>= 1` on a
  negative `long` is an arithmetic shift that sticks at `-1`, while `_base *= _base` grows
  without bound. Only `-1` and `-2` are special-cased.

- [x] **2.10 `mod()` returns the modulus instead of 0 for exact multiples** —
  `integer.h:160-163` and `integer_class.h:542-545`, `:551`: `mod(-10, 5)` → `5`.

- [x] **2.11 `pow(integer base, exp, integer result)` returns `1` for `exp == 0`**
  (`integer.h:73`) where every other branch means `result * base^exp` → should return `result`.
  It also tests `base == 2` before the `exp < 0` check.

- [x] **2.12 `simplify(rational&, rational&, rational&)`** — `rational.h:254` passes `z.den`
  where `z.num` is meant.

- [x] **2.13 `xrational::operator/=(xrational&, rational_like)` returns `b`, not `a`**
  (`xrational.h:243`) — ill-formed as soon as it is instantiated.

- [x] **2.14 `operator==(xrational, xrational)`** — `xrational.h:276`: `bs *= a.base.num.abs`
  should be `b.base.num.abs`.

- [x] **2.15 `isqrt3`** — `natural.h:408`: `v -= a` runs on a freshly-constructed (zero)
  `natural`, i.e. an out-of-bounds underflow, before `v` is overwritten. Leftover from `isqrt2`.

- [x] **2.16 `try_fermat_factorize`** (`natural.h:715-736`) — `isqrt(b_sq)` binds the
  `uint64_t` overload, truncating a `uint128_t`; and `b * b == b_sq` computes `b*b` in 64-bit.
  Broken for any `n` where `a^2` exceeds 64 bits.

- [x] **2.17 `cross(Vec3,Vec3)`** — `vector.h:188`: the y component is `a.x*b.z - a.z*b.x`,
  i.e. negated relative to the standard cross product.

- [x] **2.18 `__gcd_inner`** — `natural.h:145` uses `__builtin_ctzl` (64-bit) on `T` that can
  be `unsigned __int128` via `larger_type`.

- [x] **2.19 Non-x86-64 build is broken** — `util.h:60`: `m = a % b;` references an undeclared
  `m` (should be `r`).

- [x] **2.20 `geometry.h` does not compile if instantiated**: `operator*(Plane3,Plane3)` uses
  `a.x/a.y/a.z` on a type with members `n,d,den` (`:27`); `line_plane_intersection` uses
  `p.origin` where the member is `orig` (`:71`); `plane_intersection` calls `solve_linear` with
  scalars and `std::get<Vec3<T>>` on a `variant<None,T,Any>` (`:61-62`), needs a nonexistent
  `operator-(Plane3)` (`:53`), passes `Plane3` to `are_parallel(Vec3,...)` (`:106`), and returns
  `m / D` where Cramer's rule needs `m / det` (`:93`, leaving `det` unused).
  `same_sign` (`vector.h:222`) calls a free `sign()` that the library never defines.

- [x] **2.21 `expr.h`: `==` inside the simplification code silently means pointer identity.**
  Corrected after testing: `EXPR_CMP(==)` defines a *value* comparison as `(a-b)->sign() == 0`,
  which would recurse (`operator+` → `==` → `operator-` → `operator+` → ...), but that
  declaration comes *after* `operator+`/`operator*` in the header, so their `a == b` resolves to
  `std::shared_ptr`'s pointer comparison instead. So there is no recursion today — but the
  meaning of `==` there depends on declaration order, and structurally equal nodes were not
  being simplified (`sqrt(2) + sqrt(2)` built from two nodes stayed a sum). Fixed by using
  `identical()` explicitly, and by making `identical()` total (it threw "unreachable" for
  `sin`/`cos`/variable nodes).

- [x] **2.22 `rational::invert()` does not reject zero.** `rational_class.h:96` swaps `num` and
  `den` unconditionally, producing a rational with a zero denominator, even though the README
  documents "Throws exception if `num` is zero". Found while fixing 2.9 (`pow(0, -3)` returned
  `1/0` instead of throwing). **[confirmed]**

- [x] **2.23 `isqrt(uint64_t)` is wrong near the top of the range and for large perfect
  squares.** `natural.h:206` — for `x` close to `UINT64_MAX` the `std::sqrt` estimate is `2**32`
  and the correction `q * q > x` overflows to 0, so it returned `2**32` (whose square does not
  fit in 64 bits); and because `double` can round a large perfect square down, the truncated
  estimate could be one too small with no upward correction. Found while fixing 2.16.
  **[confirmed]**

## 3. Correct today, but fragile

- **`integer` carries a duplicate word buffer.** `integer_class.h:11` has both `natural abs`
  and `integer_backend words`. `words` is a *deep copy* — so `integer(integer&&)` (`:15`)
  heap-allocates and copies, defeating the move; `sizeof(integer)` is 32, not the README's 16
  (both `static_assert`s are commented out); `integer(std_int)` (`:14`) double-negates `words`
  so the two disagree on sign; and most arithmetic carries a `// TODO update words`, leaving it
  stale. Largest structural issue in the repo.
- **[fixed]** **Debug validation left in a hot path** — `integer_class.h:403-405`:
  `natural orig_a_abs = a.abs;` plus `Check(orig_a_abs - b.abs*c.abs == a.abs)` performs a full
  copy *and* a full multiplication inside `sub_product`, whose stated purpose is "without memory
  allocation". `integer_class.h:433`: `integer aa = a;` is an unused full copy.
- `sub_product(natural&,...)` (`natural_class.h:782-806`) and `__div`'s `__sub_product` call
  (`:863`) discard the `false` return that signals a violated precondition, leaving the operand
  half-modified.
- `natural` relies on the invariant "`words[0] == 0` whenever `size() == 0`" in ~15 comparison
  operators that read `words[0]` before checking `size()`. `operator&=(natural&, uint64_t)`
  (`:1211`) writes `words[0]` with `size() == 0`, which can break it.
- `inatural::back()` (`kernels.h:23-24`): the `const` overload returns a mutable `uint64_t&`,
  the non-`const` one returns by value — exactly inverted.
- `integer_backend::operator[](int i)` (`:141-142`) silently ignores `i` in small-buffer mode,
  so out-of-range indices read/write word 0 instead of failing.
- `integer_backend::resize()` doesn't zero the inline word when growing 0→1, so a stale value
  can resurrect after `pop_back()`.
- `natural_class.h:414`: `uint64_t(b >> 64)` sits in a branch made dead by
  `if constexpr (sizeof(b) <= 8)` *without an `else`* — still compiled, and GCC flags
  `-Wshift-count-overflow`.
- `__saturated_div` (`kernels.h:834-861`) backs off by at most 2 without verifying the third
  candidate; the margin is genuinely ~2 but nothing documents or asserts it.
- `natural::str()` doesn't validate `base` — base 1 loops until the buffer check throws,
  base 0 divides by zero.
- `__diff`/`__sub` (`kernels.h:467-528`) `goto` into the middle of a second `for` loop as a
  hand-rolled borrow state machine; correct, but no comment explains the scheme.
- `pow(base, exp, out)` in `rational.h:96` assigns `out = 1` before copying `base`, so aliasing
  `out` with `base` silently corrupts it.
- `is_power_of_three` (`natural.h:1065`) uses `goto again` that bypasses the `while (a > 1)`
  condition.
- `point_segment_squared_distance.h:11` — `if (d < 0)` on `dot(B,B)` is unreachable. Both
  distance functions truncate silently if instantiated with `T = integer`.
- `SWIZZLE4`'s macro parameter `D` shadows the template parameter `int D` (`vector.h:168`);
  unused today, broken on first use.
- `Vec::operator[]` (`vector.h:9-10`) `reinterpret_cast`s `this` to `T*`, so it is UB in
  principle and can never run in the `constexpr` evaluation it is marked for.

## 4. Duplicated code worth collapsing

- `isqrt`, `isqrt2`, `isqrt3` (`natural.h:311-434`) repeat identical 1-word and 2-word fast
  paths three times; `isqrt2`/`isqrt3` appear to be abandoned experiments (with commented-out
  debug prints).
- `operator|` / `&` / `^` and their `=` forms (`natural_class.h:1151-1251`) are six
  near-identical bodies.
- `mod3/5/6/7/9` live in `natural_class.h:158-221` each marked `// TODO move to kernels`;
  `mod10` already lives in `kernels.h`.
- `std::formatter<integer>` (`integer_class.h:618-652`) copies the entire `format()` body from
  `std::formatter<natural>` instead of delegating.
- `__add(natural&, const uint64_t*, int, int)` (`natural_class.h:419`) duplicates
  `__add_and_return_carry(inatural, cnatural, shift)`; its carry loop also lacks the early exit.
- `hash_fn_64bit` exists twice (`integer_backend.h:183`, `vector.h:250`); `int128_t`/`uint128_t`
  are typedef'd twice (`types.h:16-17`, `util.h:40-41`); `invert_bits` (`natural.h:1081`)
  duplicates `operator~`.
- `abs_greater(rational,rational)` (`rational.h:158-160`) special-cases `a.den == 1` to exactly
  what the general line computes.

## 5. Comments and README that do not match the code

- `kernels.h:671` — the `TODO ... How is this working? AA and BB are stored in Q and nested call
  will overwrite them?` is answerable: the recursive call writes into the `W` slice `r`, and `q`
  is only overwritten after `aa`/`bb` are consumed. Replace with that explanation.
- `util.h:123` — `// TODO this seems to be buggy!` on `mul_mod`: the swap preserves the invariant
  (`a*b == b*a`), and the loop is correct for `a,b in [0,m)`. The actual bug is next door in
  `__mod(cnatural, uint128_t)` (1.2) — the note is likely misattributed.
- `kernels.h:839` — `// since a > b ==> B == 1` should say `b.size == 1`.
- `natural_class.h:835` — `std::min(a.size, a.size - b.size + 1)` where the first argument is
  always larger (`b.size >= 2` here); the `min` and its stale `TODO` are noise.
- `natural_class.h:445` says `mul_karatsuba` has "no support for `&a == &q`" but the
  power-of-two branches also break for `&b == &q` (they `reset` `q` before `q = b`).
- `algebra::round(rational, digits)` (`rational.h:163`) and `real::round` truncate toward zero —
  they do not round.
- README: the opening example includes `"algebra/algebra.h"`, **which does not exist**;
  `sizeof(integer)`/`sizeof(rational)` are 32/64, not 16/32; `is_prime(const natural&)` does not
  exist; `is_likely_prime(n, rounds, rng)` has no `rng` parameter; `extract_128bits`/
  `extract_64bits` are actually `extract_u128`/`extract_u64`; `mul_mod(natural&, const natural&,
  const natural&)` is spelled `__mul_mod`; the documented `size_type` return type does not exist
  (it is `int`); `rational(integer, integer, int)` is private (use `rational::normalized`);
  `real(const rational&, int prec)` is really `real(const rational&)` plus static `real::round`;
  `sin`/`cos`/`exp`/`lcm`/`iroot`/`factorize(natural)`/`divide_bz`/`mul_karatsuba` are
  undocumented.

## 6. Performance

- `PI(n)` is intractable beyond ~n=12: `sqrt(rational, iterations)` doubles the operand size
  every Newton step with no truncation, so `PI(20)` never finishes. `PI(4)` is off by **7.0**
  because the same `n` controls both the series length and the sqrt iterations. Truncate each
  Newton iterate to the working precision and decouple the two parameters.
- `sin`/`cos` reduce with `x %= 2*PI(10)`, dragging a ~140-digit denominator through the series.
- `integer` move constructor/assignment allocate (see section 3).
- `__square`'s general branch (`natural_class.h:634-653`) restarts carry propagation per partial
  product and computes both `i*j` and `j*i`, so it is no faster than long multiplication — the
  usual 2x symmetry saving is not taken.
- `natural::str`/`kernels.h:774` convert one digit per full-precision division; chunking by
  `10^19` would cut the divisions ~19x.
- `lcm` (`natural.h:200`) forms `a*b` before dividing; `(a/gcd)*b` avoids the large intermediate
  (and `lcm(0,0)` currently divides by zero).
- `operator&=`/`^=` resize the left operand to the right's length even though the extra words
  can only be zero.
- `operator<(std_signed_int, const natural b)` (`natural_class.h:374`) takes `natural` **by
  value**; `solve_linear` (`solve_linear.h:23`) takes its first `Vec` by value while the rest are
  by reference.
- `is_prime` rebuilds both base arrays on every call (`natural.h:608-609`); they could be
  `static constexpr`.
- `factorize(uint64_t)` has no `p*p > a` early exit, so a semiprime with a 2^32-ish factor walks
  the full wheel.
- `#include <regex>` in `rational_class.h` for one parse pattern is a heavy compile-time and
  runtime cost for a "no dependencies" header-only library.

## 7. Build / portability

- **Missing includes**, currently masked by clang/libc++ include order: `<cstdint>`
  (`integer_backend.h` — first error under g++), `<cstring>` (`kernels.h` uses `std::memcpy`),
  `<random>` (`natural.h` uses `std::uniform_int_distribution`), `<optional>`, `<array>`.
- **`-std=c++23` fails with GCC.** In strict mode GCC's `__int128` is not `std::integral`, so
  `std_int<uint128_t>` is false and everything taking a `uint128_t` fails to compile. The build
  works only with clang or `-std=gnu++23`. Either support GCC or document clang-only.
- `constexpr` is claimed broadly but several functions can never be constant-evaluated:
  `swap(cnatural&,cnatural&)` uses `std::memcpy` (`kernels.h:625`), so all of Karatsuba is
  runtime-only; `rational(std::string_view)` defines a `static std::regex` inside a `constexpr`
  function (`rational_class.h:54`); `__divide_2n1n` is not marked `constexpr` but `__divide_bz`
  is.
- `.bazelrc` applies `-ffast-math` globally to an *exact* arithmetic library. It does not affect
  the integer paths, but it does affect `std::sqrt`/`frexp`/`ldexp` seeds in `isqrt`,
  `isqrt_hardware`, `round_to_zero` and `rational(float)` — where correctness depends on the FP
  result being within a known error bound.
