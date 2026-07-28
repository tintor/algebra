# Code review findings — round 2

Full read of every header and test in `algebra/` (~19.7k lines) at commit `c1fa9ef`
(branch `rename-natural-test-helpers`). This is a fresh pass; `REVIEW.md` holds the
previous round and is not repeated here, except where an entry there is marked
**[fixed]** but the code says otherwise (section 1, item 1.10).

The whole suite passes as of this review (`bazel test //...`, 29 targets, all green), so
everything below is a gap the tests do not cover.

Items marked **[verified]** were reproduced by compiling and running a program against
the headers; the observed output is quoted. The rest are from reading.

---

## 1. Wrong results, undefined behaviour, crashes

**1.1 `integer(s, 8)` is wrong for octal strings of 42 digits or more** — `integer_class.h:1754`
**[verified]**

The base-2/8/16 parser accumulates a chunk of digits and then calls
`words.insert_first_word(acc)`, which shifts the value left by a full 64 bits. Base 2
chunks at `count == 64` and base 16 at `count == 64`, both correct; base 8 chunks at
`count == 63` (21 octal digits) but still inserts a whole word, so everything parsed
before the second chunk is multiplied by 2^64 instead of 2^63.

```
integer(std::string(42,'1').c_str(), 8)
  got 24305883351495604531780561668432761417
  ref 12152941675747802266549093122563150409   (exactly 2x smaller high part)
```

One chunk is fine (the shift lands on an empty value), which is why
`integer("777", 8) == 511` — the only base-8 assertion in the suite — passes. Fix: use
`count == 64` with a 64/3 split, or shift by `count` bits as the tail path already does.

**1.2 `mul_karatsuba()` can return an unnormalized value, so `==` says it differs from `a*b`** — `integer_class.h:1130-1146` **[verified]**

`Q = mul_max_size(a, b)` is an upper bound that overshoots by one word when the two top
words are both small (`clz(a.back()) + clz(b.back()) >= 64`). `__add_product()` does not
normalize, and `q.words.resize(vq.size)` keeps the full `Q`, so the result carries a zero
top word. `operator==` compares sizes first, so the value is unequal to itself computed
another way, and `std::hash` disagrees too.

```
a = (5*2^30) << 64, b = (3*2^30) << 64
mul_karatsuba: size 4, top word 0        a*b: size 3
same digits, a * b == mul_karatsuba(a, b)  ->  false
```

Fix: `q.words.normalize()` before returning. Note `__abs_div` has the same shape at
`integer_class.h:1456` (`q.words.resize(vq.size)`) but `__div` normalizes internally, so
that one is safe.

**1.3 `div(integer, T)` and `integer / T` silently truncate a 128-bit divisor** — `integer_class.h:773`, `integer_class.h:1567` **[verified]**

`div(const integer&, T, integer&)` passes `abs_unsigned(b)` (a `uint128_t` for
`T = int128_t`) into `__abs_div(const integer&, uint64_t, integer&)`, and
`__abs_div(integer&, std_int auto b)` does an explicit `static_cast<uint64_t>(b)`. Both
paths keep only the low 64 bits of the divisor, with no check.

```
a = 123456789012345678901234567890, d = 2^70 + 1
div(a, (int128_t)d, q)   -> q = a          (divisor truncated to 1)
a / (uint128_t)d         -> q = a
expected                    104571967
```

This is the worst finding: a plain `/` on a public type returns a wrong answer with no
diagnostic. Minimum fix is `static_assert(sizeof(T) <= 8)` on those overloads (there are
already `// TODO generalize for any std_int` comments); better is to route a wide divisor
through the `integer` path.

**1.4 `expr_power::sign()` returns 0 for a negative base with an odd exponent** — `expr.h:143` **[verified]**

```cpp
if (b < 0 && exp.is_integer())
    return exp.num.is_even();      // bool -> 1 or 0; 0 means "zero"
```

`sign((-2)^3)` returns `0`, i.e. "the value is zero", instead of `-1`. Every comparison
operator on `expr_ptr` is `(a - b)->sign() OP 0`, so any comparison whose difference
contains an odd power of a negative base silently answers as if the terms cancelled.
Fix: `return exp.num.is_even() ? 1 : -1;`.

**1.5 `expr_matrix::sign()` indexes an empty vector** — `expr.h:70-72` **[verified]**

```cpp
if (rows == 0 && cols == 0)
    return data[0]->sign();
```

A 0x0 matrix has no elements. Reproduced as an out-of-bounds `std::vector::operator[]`
(caught here by libstdc++ assertions; UB in the `-O3 -DNDEBUG` build). The intent was
almost certainly `rows == 1 && cols == 1`.

**1.6 `round_to_zero()` drops the sign** — `integer.h:499-512` **[verified]**

`if (m < 0) m = -m;` takes the magnitude and nothing puts the sign back:
`round_to_zero(-3.75, b)` gives `b == 3`. The only in-tree caller (`iroot`) passes a
positive value, so this is latent, but it is a public function whose name promises
signed truncation. It also has no `isnan`/`isinf` guard, so
`static_cast<uint64_t>(inf * 2^53)` is UB — reachable from `iroot(a, n)` when
`static_cast<double>(a)` overflows to infinity (a > ~1e308, i.e. ~1024 bits), which is
an ordinary input size for this library.

**1.7 `isqrt(T)` computes garbage for `uint128_t`** — `integer.h:177-188` **[verified]**

The template is constrained to `std_unsigned_int`, which admits `uint128_t`, but the body
clamps to `MAX = UINT32_MAX` "so `q*q` does not overflow" — a 64-bit assumption.
`isqrt(uint128_t(1) << 80)` returns `4294967295` instead of `2^40`. `__isqrt_u128()`
right below it is the correct routine. Constrain `isqrt` to `sizeof(T) <= 8` or dispatch.
`factorize<T>` (`integer.h:324`) has the same exposure and additionally calls
`std::countr_zero(a)`, which does not accept `__int128` at all.

**1.8 `util.h::pow(T base, unsigned exp)` truncates to the width of `T`** — `util.h:115` **[verified]**

The declared return type is `uint64_t`, but `base *= base` happens in `T`:
`pow((uint32_t)100000, 3)` gives `141006540800000` instead of `10^15`
(`pow((uint64_t)100000, 3)` is correct). Compute in `uint64_t`.

**1.9 `util.h::min(signed, unsigned)` returns the wrong element** — `util.h:217-218` **[verified]**

`larger_type<A,B>` picks `A` when the sizes are equal, and the comparison `a < b` is then
done at that width with the usual signed/unsigned conversion:
`min(int32_t(-1), uint32_t(5))` returns `5`. The overload set is also asymmetric — there
is `min(signed, unsigned)` and `max(signed, signed)`, but no `max(signed, unsigned)` or
`min(signed, signed)`. Nothing in the library calls these (dead code today), so the fix
is cheap: delete them or compare correctly.

**1.10 `solve_linear()` still takes its first `Vec` by value** — `solve_linear.h:23`

`REVIEW.md` section 6 marks this **[fixed]**, and `README.md:478` documents the signature
as `const Vec<D,T>&`, but the code is unchanged:
`solve_linear(const Vec<D, T> a, const Vec<D, T>& b)`. For `Vec3<rational>` that is three
heap-allocating copies per call, on the geometry hot path. (The other half of that
REVIEW.md entry, `operator<` taking a `natural` by value, is genuinely gone with the
class.)

**1.11 `uniform_sample(count, rng)` returns garbage for `count == 0`** — `integer.h:40-44` **[verified]**

`Check` only rejects a negative count. For zero, the `is_uint64()` fast path builds
`std::uniform_int_distribution<uint64_t>(0, 0 - 1)`, i.e. the full 64-bit range:

```
uniform_sample(integer(0), rng)  ->  2469588189546311528
```

An empty range should throw. Add `Check(!count.is_zero(), ...)`.

**1.12 `log_lower(n, 1)` / `log_upper(n, 1)` never terminate** — `integer.h:1241-1265`

`a /= base` with `base == 1` leaves `a` unchanged, so the loop hangs. `base == 0` throws
from `__abs_div`, which is fine; base 1 needs the same check.

**1.13 `x / x` returns 1 for a zero `xrational`** — `xrational.h:218` **[verified]**

`if (&a == &b) return xrational{rational{1}};` fires before any zero test, so
`z / z == 1` for `z == 0`. Every other division-by-zero path in the library throws.

---

## 2. Fragile: correct today, one edit away from wrong

**2.1 `__mul(cnatural, cnatural, vnatural&, bool)` reads `b[0]` unconditionally** —
`kernels.h:398`. With `b.size == 0` and `a.size > 0` the first iteration dereferences
`b.words[0]`, which for an empty magnitude is the small-buffer word of an unrelated
object (or garbage). All current callers pre-check for zero; nothing in the function does.

**2.2 The Karatsuba scratch buffer size has no derivation** — `integer_class.h:1136`
(`const int W = 4 * std::max(A, B)`). Working the recursion out gives a peak of about
`4n + 3*log2(n/32)` words, i.e. `4n` is *just* enough only because the base case leaves
slack: at `A = B = 16384` the requirement is 65522 words against a `W` of 65536, 14 words
to spare. The comment at `kernels.h:711-714` records that `KARATSUBA_LIMIT = 8` "fails,
likely due to small W buffer size" — that guess is right, and the same cliff is one
tuning change away at the current limit. `Check(w <= we)` turns it into a thrown
exception rather than corruption, but it should be a derived bound with a comment, not 4x.

**2.3 `bit_range` underflows for a zero operand** — `util.h:205`.
`operator*` computes `{a.min + b.min - 1, ...}`; `bit_range(0)` sets `min = max = 0`, so a
product involving it yields `min == UINT64_MAX` and the `<` disjointness test becomes
nonsense. `xrational::operator==` (`xrational.h:266`) is the only user and its operands
happen to be provably nonzero there.

**2.4 `div(a, b, quot, rem)` with `&quot == &rem` silently corrupts** —
`integer_class.h:741-745`. `magnitude_ref` documents that two live borrows over one
backend lose the result ("the second would put the pre-operation value back"), and
`div()` takes two out-parameters with no check that they differ.

**2.5 `erase_first_n_words(n)` with `n > size()` flips the sign** —
`integer_backend.h:301-307`. `_size += -n` walks past zero and turns a positive value into
a negative one of bogus magnitude. Only `__abs_shl` calls it, and it pre-checks; a `Check`
here is one line.

**2.6 `__sub_product()` leaves its target half-modified when it returns false** —
`kernels.h:832-857`. The contract says "returns false otherwise" but not that `a` is then
garbage. `__abs_mod(integer&, const integer&)` (`integer_class.h:1186`) **ignores the
return value entirely**, unlike the four other call sites, all of which wrap it in
`Check`.

**2.7 `abs()` on a floating-point `Vec` may resolve to `::abs(int)`** — `vector.h:130-135`.
Unqualified `abs` inside `namespace algebra` finds `algebra::abs(integer/rational/Vec)`;
for `T = double` it depends on `<cmath>` having injected `::abs(double)` into the global
namespace. Add `using std::abs;`.

**2.8 `__safe_step()`/`__ray_meets_segment()` skip edges parallel to the ray** —
`polygon2_boolean.h:100-103`, `polygon2_arc.h:196-206`. The justification ("it either misses,
or overlaps and hits an endpoint anyway") relies on the *other* edge at that endpoint not
also being parallel, which `simplify()` normally guarantees but nothing enforces.

**2.9 `contains(ArcPolygon2, p)` can fall through to an undefined winding** —
`polygon2_arc.h:274-284`. It tries two step directions and breaks out of the loop when the
new point is clear; if neither is clear it uses the last one anyway, and applies the
winding formula on a chord where the comment says it does not apply. It should `Check` or
try a third direction.

**2.10 `maybe_stack` is copyable** — `kernels.h:41-55`. A copy would carry a `_t` that
points into the source's inline array and double-`delete[]` a heap one. Delete the copy
and move members.

**2.11 `boolean_op` and `buffer` silently misbehave for an integral `T`** —
`polygon2_boolean.h:195`, `polygon2_arc_boolean.h:131`. `(f.a + f.b) / T(2)` and the
`best / 2` inside `__safe_step` truncate to zero for `T = integer`, putting the sample
point exactly on the boundary. The header comment says "T needs exact arithmetic" but
nothing rejects a ring type that is exact yet not a field.

---

## 3. Duplicated code

**3.1 Four hand-written binary GCDs.** `__gcd_inner` (`integer.h:104`),
`gcd(integer, integer)`'s loop (`integer.h:146-155`), and two more inside
`rational::simplify()` — one on `uint64_t` (`rational_class.h:246-251`) and one on
`integer` (`rational_class.h:263-268`). The last two exist to avoid the sign handling and
the `az - z` re-shift; both could call the first two with a small helper.

**3.2 Three copies of the schoolbook division loop.**
`__abs_div(cnatural, cnatural, integer&, integer&)` (`integer_class.h:1489-1498`),
`__abs_mod(cnatural, cnatural, integer&)` (`integer_class.h:1511-1518`) and
`__abs_mod(integer&, const integer&)` (`integer_class.h:1181-1188`) run the same
`__saturated_div` + `__sub_product` loop with different bookkeeping. The third also
differs in *behaviour* (it ignores the `__sub_product` result, see 2.6).

**3.3 Two `str()` implementations for the magnitude.** `str(cnatural)` (`kernels.h:869`)
divides by 10 once per digit; `integer::__abs_str_buffer` (`integer_class.h:1579`) chunks
by `10^19`. The free function `str(const integer&)` (`integer_class.h:595`) routes to the
slow one *and* copies through `abs(a)`, so `str(a)` is ~19x slower than `a.str()` and
allocates an extra integer. Everything should go through the member.

**3.4 `segment_segment_intersection` duplicates `segment_segment_intersection_param`** —
`segment_segment_intersection.h:117-186` repeats all six branches of the function above it
(degenerate, collinear, overlap, swap bookkeeping) to return points instead of parameters.
One could be a thin wrapper over the other.

**3.5 `factorize(integer)`'s trial-division body is pasted twice** —
`integer.h:958-1016`, once for `p += 2` and once for `p += 4`. ~28 lines each, including
the perfect-square loop and two `is_likely_prime` calls. The `uint64_t` overload
(`integer.h:362-396`) has the same duplication.

**3.6 `abs_greater` is defined twice with different intent** —
`integer.h:1141` / `rational.h:191` ("minimizing memory allocation") and a generic
`abs(a) > abs(b)` template in `segment_segment_intersection.h:23`. The generic one silently
covers any type the specific ones do not, including `xrational`, where it is not defined
at all — so `argmax_abs` / `div_colinear` on an `xvec2` picks up the slow path or fails to
compile depending on the instantiation.

**3.7 `algebra.h` includes `integer.h` twice** (lines 3-4).

**3.8 The `100`-digit scale constant appears twice** in `real_class.h` (`operator/` line
158, `operator/=` line 164) as a bare literal.

---

## 4. Comments and documentation that do not match the code

**4.1 `README.md` documents a class that no longer exists.** 98 occurrences of `natural`,
including a full `## class natural` API section, `algebra/natural.h` in the headers table,
`natural integer::to_natural()`, `rational::rational(natural)`, `real::real(natural, int)`,
`xrational::root` typed as `natural`, and the `123_n` literal. All of that was removed by
commits `a660945`, `91f1500` and `af34c48`.

**4.2 `README.md` sizes are wrong.** It claims `sizeof(integer)` is 32 bytes and
`sizeof(rational)` 64; the code static-asserts 16 and 32.

**4.3 `README.md` has no mention of the polygon feature set.** `polygon2.h`,
`polygon2_arc.h`, `polygon2_boolean.h`, `polygon2_buffer.h`, `polygon2_arc_boolean.h`,
`point_segment_squared_distance.h` and `types.h` are absent from the headers table, and
`MultiPolygon2` / `ArcPolygon2` / `ArcRegion` / `boolean_op` / `buffer` / `contains` /
`winding_number` are undocumented — even though GOAL.md marks both the polygon work and
"fill the missing functions in README.md" as done.

**4.4 `README.md`: "all types support casting to and from all built-in integer and
floating point types."** `integer` has `operator T()` for floating point but no
constructor or assignment from `float`/`double` (`round_to_zero` is the substitute), so
the "from" half is false for the central type.

**4.5 `lcm()`'s comment says "of the magnitudes"; it is not** — `integer.h:164-172`
**[verified]**. `lcm(-4, 6)` returns `-12`: the sign of `a*b` survives. Either take the
magnitudes or fix the comment.

**4.6 `mod(integer&, const integer&)` carries `binominal_mod`'s comment** —
`integer.h:1188` (comment) / `integer.h:1189` reads `// returns (n k) mod m` above an in-place modulo.

**4.7 `__less_a_bc_scalar`'s comment is garbled** — `kernels.h:131` says `// A == B || A == B + 1`
where the invariant being described is about *bit counts* being equal or off by one.

**4.8 `round()` does not round** — `rational.h:198`, and `real<B>::round`
(`real_class.h:104`). Both truncate toward zero. The README says so for the rational one
("Truncates towards zero") and contradicts itself for the real one ("Nearest value ...
(truncated towards zero)"), but the *names* still promise rounding, and
`__reduce_mod_two_pi` relies on the truncating behaviour.
**[verified]** `round(rational("0.19"), 1) == 1/10`.

**4.9 `Random::Uniform` for 128-bit types is biased, and the comment claims otherwise** —
`__test.h:56`. `do { ... } while (false)` runs once and takes `r %= limit`, i.e. plain
modulo bias, under a comment reading "rejection sampling, so the result is uniform". When
`span == U_MAX` the code sets `limit = span`, so `max` can never be produced at all.

**4.10 `__mul_mod`'s comment says "for testing"** — `integer.h:768` — but `pow_mod` and
`is_likely_prime`, the two hottest number-theory routines, both call it.

**4.11 `circle_ring`'s orientation comment is inverted** — `polygon2_arc.h:312` says
"right vertex bulging below, left vertex bulging above". With `perp(v) = {-v.y, v.x}` and
`bulge = -1`, the edge from the right vertex to the left one bulges *above*. The ring is
counter-clockwise as claimed; the per-vertex description is backwards.

**4.12 `__ray_meets_circle` is conservative in a way the comment does not cover** —
`polygon2_arc.h:233`. When both interval ends are strictly *inside* the circle
(`C < 0`, `f(tmax) < 0`) an upward parabola cannot cross zero in between, yet the final
`return B*B >= A*C` is trivially true because `A*C < 0`. The documented conservatism is
"whole circle, not just the arc"; this extra one costs unnecessary halvings in
`__arc_safe_step` and belongs in the comment or, better, in a `C > 0` guard.

**4.13 `area_bounds` documents a bound it does not provide** —
`polygon2_arc_boolean.h:110-117`: "Returns the lower bound and the amount still
undecided, so lower <= true area <= lower + undecided". The decision is five point samples
per box; a box whose five samples are all inside is counted as fully inside. That is a
heuristic, not a bound (the two-level forced subdivision reduces the error without
eliminating it). It also assumes the caller's `min`/`max` box contains the region, with no
check.

**4.14 `ccw()` is clockwise-positive** — `segment_segment_intersection.h:10`. For
`a=(0,0), b=(1,0), c=(0,1)` it returns `-1`, the negation of the usual convention and of
`signed_area2()` in `polygon2.h`. Only `== 0` tests use it today, so nothing is wrong, but
two files in the same library now disagree about the sign of orientation.

**4.15 `fract()` and `trunc()` do not compose** — `rational.h:180`, `rational.h:204`
**[verified]**. `trunc(-3/2) == -1` and `fract(-3/2) == 1/2`, so `trunc + fract == -1/2`,
not `-3/2`. The README does say "always non-negative" for `fract`, and a test pins the
behaviour, but the pair is a trap; at minimum the two should reference each other.

**4.16 `mul_mod(integer, ...)`'s precondition is undocumented** — `integer.h:774`.
`add_mod`/`sub_mod` right above it both say "assumes a and b are in [0, m-1] range";
`mul_mod` needs the same and says nothing, while internally calling `add_mod`.

---

## 5. Interface and architecture

**5.1 Operators on `std::shared_ptr` and `std::vector` are hijacked in namespace `algebra`.**
`expr.h:566-576` defines `==`, `!=`, `<`, `<=`, `>`, `>=` on `expr_ptr` (a
`std::shared_ptr<expr>`), and `expr.h:10-35` defines `+` on `std::vector<T>`. Both are
found by ordinary lookup from any scope with `using namespace algebra`. The comparison
operators are non-templates, so they *beat* `std::shared_ptr`'s own `operator==` — pointer
identity comparison quietly becomes a symbolic value comparison that can throw
`unknown_sign_error`. **[verified]** `a == b` on two distinct `make_integer(2)` nodes
returns `true`, and `std::vector<int> + std::vector<int>` compiles in `main`. Wrapping
`expr_ptr` in a real class type would contain both.

**5.2 `mod(integer&, const integer&)` and `mod(const integer&, const integer&)` overload on
constness with different semantics** — `integer.h:1189` vs `integer_class.h:848`. For a
non-const lvalue the in-place one wins, so `mod(a, b);` mutates `a` while
`mod(as_const(a), b)` returns a value. Different names would remove the hazard.

**5.3 `simplify()` means three different things.** `rational::simplify()` normalizes;
`simplify(integer&, integer&)` divides by the GCD (rescaling a direction vector);
`simplify(MultiPolygon2&)` drops degenerate vertices and rings. All are reachable by
unqualified call from the same scope.

**5.4 Three public `isqrt` variants plus `isqrt_hardware`** — `integer.h:567-665`. The
README says isqrt2/isqrt3 are "kept for benchmarking", but they sit in the public header
with commented-out debug prints (`integer.h:603`, `613`) and duplicated final-correction
logic.

**5.5 `real<B>::operator/` hardcodes 100 digits of precision** — `real_class.h:158`. The
result's accuracy depends on neither operand nor any parameter, and the type has no way to
ask for more. Also `real` is the one arithmetic type here whose `/` is inexact, which no
comment in the header states.

**5.6 `REAL_FRACT_DIGITS` is a mutable global in a header** — `real_class.h:250`. Formatting
of every `real<B>` depends on process-wide state; not thread-safe and invisible at the call
site.

**5.7 `buffer()` takes a plain function pointer for the element** —
`polygon2_buffer.h:135` (`Ring2<T> (*element)(const T&)`). `polygon_element` takes a side
count, so the one element the header advertises as "a finer convex polygon approximates a
disk as closely as wanted" cannot be passed to `buffer()` at all. `std::function` or a
template parameter would fix it.

**5.8 `ZERO_EXPR`/`ONE_EXPR`/`E_EXPR`/`PI_EXPR` are non-inline `const` globals in a header** —
`expr.h:208-211`. Every translation unit gets its own dynamically-initialized copy
(`link_test.cc` only checks that this links, which it does), so pointer identity differs
across TUs and static-initialization order matters. `inline const` fixes it.

**5.9 `expr_sum::_sign` is declared `mutable` and never used** — `expr.h:157`.
`expr_product` caches into its `_sign`; `expr_sum::sign()` — much the more expensive of
the two — does not.

**5.10 `expr_sum::sign()` recurses without a depth bound** — `expr.h:881-890`. It squares
both partial sums and recurses on the difference; the only guard is "either side is a power
above 1". Nested radicals can grow the expression on the way down.

**5.11 `pow(expr_ptr, rational)` distributes roots over powers and products** —
`expr.h:537-543`. `pow(pow(x,2), 1/2)` folds to `x`, and `sqrt(a*b)` to
`sqrt(a)*sqrt(b)`; both are false for negative operands, and `expr_var` makes negative
operands representable.

**5.12 `MultiPolygon2`/`ArcPolygon2` `operator==` is structural, not geometric** —
`polygon2.h:139`, `polygon2_arc.h:294`. Rotating a ring's vertex list, or reordering the
rings, makes two representations of the same region compare unequal. Nothing says so.

**5.13 `boolean_op` does not simplify its result** — collinear vertices introduced by the
cutting stage survive: the union of `[0,1]²` and `[1,2]x[0,1]` comes back as
`(0,1)(0,0)(1,0)(2,0)(2,1)(1,1)`, with `(1,0)` and `(1,1)` redundant. **[verified]**

**5.14 `__stitch` picks an arbitrary outgoing fragment at a shared vertex** —
`polygon2_boolean.h:148-155`. Where four fragments meet, the "leftmost turn" rule is what
separates two touching rings from one pinched ring. Membership and area are unaffected
(winding is orientation-based), so this is a representation-quality issue, but the result
depends on the sort order of `frags` rather than on geometry.

**5.15 `try_fermat_factorize` overloads its return value** — `integer.h:296`. `0` means
"not found", but `1` (a legitimate return when `n = a² - b²` with `a - b = 1`) means
"prime" and is not a usable factor. No caller in tree; a caller looping on the result
would spin.

**5.16 `iroot(a, 0)` answers two different ways** — `integer.h:667-672`. `a <= 1` returns
`a` (so `iroot(0,0) == 0`) and everything else returns `1`, because the `n == 0` test comes
after the `a <= 1` test.

**5.17 `integer::bit()` and `popcount()` disagree about representation** —
`integer_class.h:281-312`. `popcount()` documents (and implements) two's complement for
negative values; `bit(i)` reads the magnitude. `integer(-1).bit(0)` and
`(~integer(0)).bit(0)` both return true for different reasons. **[verified]**

**5.18 `integer % int128_t` does not compile.** The `%` overload set is
`(integer&, int64_t)`, `(integer&, int)`, `(integer&, unsigned)`, `(integer&, uint64_t)`,
`(integer&, const integer&)`; a 128-bit right operand is an ambiguous integral conversion.
The overloads also return four different types (`int`, `int64_t`, `integer`).

**5.19 `std::format("{:>10}", integer)` is a compile error** — `integer_class.h:1817`
**[verified]**. `parse()` only recognizes an alignment when it is preceded by a fill
character, so the common `{:>10}` / `{:<10}` forms fail the `Check` at the end of
`parse()`; `{:*>10}` works. Since C++23 checks format strings at compile time, this is a
hard error rather than a throw, but it is still a gap versus `std::format` for builtins
(as are `{:+}`, `{:#x}` and `{:_}`).

**5.20 Formatters discard `format_to`'s returned iterator** — `xrational.h:357`,
`expr.h:588`, `real_class.h:254`, `vector.h`, `dual.h`. They call `std::format_to(ctx.out(), ...)`
and then `return ctx.out()`. That works for `std::format`'s back-insert iterator and
breaks for `format_to` into a `char*` or a `format_to_n` counting iterator.

**5.21 `constexpr` on functions that can never be constant-evaluated** —
`operator<<(std::ostream&, ...)` in four headers, and `dcast` (`expr.h:62`), which is
`dynamic_cast`. Harmless, but it dilutes the "full constexpr support" claim.

**5.22 `<ostream>` is never included** although five headers define
`operator<<(std::ostream&, ...)`. It compiles only because `<print>`/`<format>` pull it in
transitively on this toolchain.

---

## 6. Performance

**6.1 `isqrt`, `isqrt2` and `isqrt3` each end with a big multiplication that appears never
to be needed** — `integer.h:583-588`, `618-624`, `656-664`. Newton from an initial guess
`>= sqrt(a)` converges to exactly `floor(sqrt(a))`, so the trailing `(x+1)² <= a` check
should be dead; replaying `isqrt`'s loop over 20k random values of 129..600 bits, the
correction fired 0 times. It costs a full multiply per call (`isqrt` is called from `factorize`,
`exact_sqrt`, `is_power_of_three`, `sqrt_bits` and `xrational::simplify`).

**6.2 `is_likely_prime` rebuilds its 40-entry base array on every call** —
`integer.h:841`. `std::array<int,40> primes = {...}` is a local, not `static constexpr` —
the same fix `REVIEW.md` applied to `is_prime`.

**6.3 `str(const integer&)` is ~19x slower than `integer::str()`** — see 3.3.

**6.4 `pow_mod` / `is_likely_prime` reduce with `a *= b; a %= m`** — `integer.h:767-772`.
The full double-width product is formed before reduction, and `mul_mod` (the routine with
the size-aware fast paths) is right there and unused. See 4.10.

**6.5 `mod63_65` extracts 64 bits to consume 12** — `integer.h:873-898`: `extract_u64` is
called once per 12 bits, so ~5x the necessary work, and `a.num_bits()` is recomputed in the
loop condition.

**6.6 `rational operator+/-` with an integral operand call `simplify()` needlessly** —
`rational_class.h:157`, `343`, `344`. `gcd(num + b*den, den) == gcd(num, den) == 1`, so the
result is already in lowest terms; `normalized()` exists for exactly this.

**6.7 `dilate()` runs one full boolean union per edge** — `polygon2_buffer.h:104-114`. Each
union is quadratic in the accumulated edge count, so buffering an n-gon is roughly cubic.
Building the union of all the per-edge hulls in one pass, or unioning pairwise in a tree,
would help; the header documents the algorithm but not the cost.

**6.8 `simplify(MultiPolygon2&)` is O(n³)** — `polygon2.h:200-212`: the collinear scan
restarts from index 0 after every removal, and each removal is an O(n) `vector::erase`.

**6.9 `expr::make_sum` compares terms with `safe_sign(v[i] - v[j])`** —
`expr.h:335`. That builds and signs a difference expression for every pair (the cheaper
`identical(v[i], v[j])` is right there, commented out on the next line), and `operator-`
can recurse back into `make_sum`.

**6.10 `ArcRegion` deep-copies both operands on every combination** —
`polygon2_arc_boolean.h:74`. `make_shared<const ArcRegion>(a)` copies the whole subtree, so
building an n-node tree costs O(n²). The nodes are already immutable and shared_ptr-held;
the operators could take `shared_ptr` and share instead of copy.

**6.11 `interval::pow` constructs an `integer` per loop iteration** — `expr.h:754`
(`for (int i = 0; i < abs(b); i++)`), and rebuilds `abs(b)` each time.

**6.12 `.bazelrc` passes `-mavx2 -mfma`** for a library whose hot loops are 64-bit integer
multiply/divide and inline `divq`. It buys nothing measurable here and makes the binaries
non-portable to pre-Haswell hardware.

---

## 7. Test gaps and weak tests

**7.1 The random big-integer generator never produces a small top word** —
`integer_magnitude_test.cc:21-41`. `rand_integer` fills *every* word, top included, from
the full `uint64_t` range, so `clz(a.back())` is 0 with overwhelming probability. Every
bit-length-sensitive path — `mul_max_size`, `div_max_size`, `add_max_size`,
`__saturated_div`'s two corrections, `__less_a_bc*`'s bit-count screens — is therefore
exercised on one narrow input class only. This is exactly why 10M iterations of
`mul_karatsuba 4` miss finding 1.2. Generating a random *bit length* first (as
`integer_number_theory_test.cc` does with `uniform_sample_bits`) would cover it.

**7.2 Base 8 and base 2 parsing have one assertion between them** —
`integer_class_test.cc:656` (`integer("777", 8) == 511`). Chunk boundaries (21, 22, 42, 43
octal digits; 64, 65 binary; 16, 17 hex) are untested, which is why 1.1 survived.

**7.3 No test drives `integer` division by a 128-bit divisor** — hence 1.3. Nor is
`div(integer, T, integer&)` tested for any unsigned `T`.

**7.4 The integer→builtin conversion operators are almost entirely untested.**
`operator char`, `short`, `uint8_t`, `uint16_t`, `int128_t`, `is_uint16()`, `low_word()`,
`size_of()`, `unsafe_u128()`, `str_size_upper_bound()` never appear in any test file — the
overflow `Check`s and the negation-in-unsigned logic they document go unverified.

**7.5 `expr` is the least covered public component.** 147 assertions total, and
`make_sum`, `make_product`, `subvec`, `identical` on most node types, `bounds`,
`interval` arithmetic, `expr_matrix`, `expr_var`, `expr_rel` and `needs_parenthesis` have
no direct test. Findings 1.4 and 1.5 are both one assertion away. Three assertions are
commented out at `expr_test.cc:25-27`.

**7.6 Assertion-free tests.** `real_test.cc` has an empty `TEST_CASE("empty")`;
`real_class_test.cc`'s `TEST_CASE("hash")` inserts into an `unordered_map` and asserts
nothing; `integer_class_test.cc` and `integer_magnitude_test.cc` both contain
`TEST_CASE("<<") { integer(1) << 64; }` (a duplicated regression test that would only
catch a crash).

**7.7 ~100 lines of dead code in `integer_number_theory_test.cc`.** `fast_isqrt`
(lines 300-387), `double_to_integer` and `doubleToLongBits` — the last of which is
strict-aliasing UB (`*reinterpret_cast<const ulong*>(&a)`) — are reachable only from an
`#if 0` test case (lines 389-397) that is itself an infinite `while (true)` print loop.
Two more `#if 0` blocks sit at lines 73-97 and 247-284.

**7.8 `rational_test.cc:31` (`round(aa, 100) == 2`) passes only because of which side
Newton converges from.** `round()` truncates, so the assertion holds for
`sqrt(2)²  > 2` and would fail if the iteration ever approached from below. An
`abs(aa - 2) < eps` form would be robust.

**7.9 Untested error paths.** Nothing asserts that `Check`/`throw` fires for: division by
zero on `rational`/`real`/`xrational`, `signed_area`/`bounding_box` of an unbounded region,
`str()` with a base outside 2..36, "buffer too small", `dilate` with an element that
misses the origin, `nth_root` with a zero exponent, `real(rational)` inexact conversion,
or the `integer -> intN` overflow checks.

**7.10 No test for negative or zero inputs to `sqrt(integer, iterations)` /
`sqrt(rational, iterations)`**, which unlike `sqrt_bits` have no `Check` at all
(`rational.h:10`, `rational.h:28`) and return nonsense for a negative argument.

**7.11 `abs_greater` is tested only for `rational`** (`rational_test.cc`), not for
`integer`, and not at all for the generic template in
`segment_segment_intersection.h`, which is what `argmax_abs`/`div_colinear` use for
non-`rational` coordinate types.

**7.12 Stale test logs from the rename.** `bazel-testlogs/algebra/natural_test/` and
`natural_class_test/` are still present from before `af34c48`; harmless, but they make
`bazel-testlogs` misleading when read directly.

---

## 8. Build and tooling

**8.1 No `-Wall -Wextra`.** `.bazelrc` sets only `-Wfatal-errors`, `-DNDEBUG` and `-O3`.
With `-Wall -Wextra` the umbrella header produces 12 `-Wsign-compare` warnings
(`expr.h:319,430,628,648,850`, `integer.h:817,820,843`, `rational.h:127`,
`integer_backend.h:24`), 1.9 (`min(int32_t, uint32_t)`) is exactly a `-Wsign-compare` at the
first call site that instantiates it.

**8.2 `UINT128_MAX` is `static const auto` at namespace scope** — `types.h:13`. One copy per
translation unit and not usable in a constant expression context that requires
`constexpr`; `static constexpr` costs nothing.

**8.3 `library()` in `library.bzl` assumes a `<name>_test.cc` exists** and has no way to
declare extra test-only sources; the four tests that need something different
(`integer_magnitude`, `integer_number_theory`, `integer_backend`, `portable_divq`) are
spelled out by hand in `BUILD`, which is why `integer_class_test.cc` and
`integer_magnitude_test.cc` have drifted into near-duplicates of each other.

---

## Suggested order of work

1. `1.3` (silent wrong quotient), `1.1`, `1.2`, `1.4`, `1.5` — wrong answers and UB in
   public API, each with a one- or two-line fix and an obvious regression test.
2. `1.6`, `1.7`, `1.8`, `1.9`, `1.11`, `1.12`, `1.13` — latent wrong answers; mostly
   `static_assert`/`Check` additions.
3. `8.1` — turn on `-Wall -Wextra` and clear the 12 warnings; this closes the class of bug
   behind 1.8 and 1.9.
4. `7.1` — fix the random generator to vary bit lengths, then re-run the existing stress
   tests; expect them to find more of section 2 on their own.
5. `4.1`-`4.4` — the README describes a library that no longer exists.
6. Sections 3, 5, 6 as ordinary cleanup.
