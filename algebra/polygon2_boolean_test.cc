#include "algebra/polygon2_boolean.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"

using P = MultiPolygon2<rational>;
using R = Ring2<rational>;
using V = Vec2<rational>;

static R box(int x0, int y0, int x1, int y1) {
    return R{V(x0, y0), V(x1, y0), V(x1, y1), V(x0, y1)};
}

// The result is checked by membership rather than by ring shape: a boolean operation has many
// correct representations, but only one correct set of points. Sampling on a half integer grid
// avoids landing on the boundary, where a closed region reports true on both sides.
static void require_same_region(const P& got, auto&& expected, int lo = -4, int hi = 14) {
    for (int i = 2 * lo; i <= 2 * hi; i++)
        for (int j = 2 * lo; j <= 2 * hi; j++) {
            const V p(rational(i, 2) + rational(1, 4), rational(j, 2) + rational(1, 4));
            REQUIRE(contains(got, p) == expected(p));
        }
}

TEST_CASE("union of disjoint boxes") {
    const P a(box(0, 0, 4, 4)), b(box(6, 0, 10, 4));
    const P u = a | b;
    require_same_region(u, [&](const V& p) { return contains(a, p) || contains(b, p); });
    REQUIRE(signed_area(u) == 32);
    REQUIRE(u.rings.size() == 2);
    REQUIRE(!u.complement);
}

TEST_CASE("union of overlapping boxes") {
    const P a(box(0, 0, 6, 6)), b(box(4, 4, 10, 10));
    const P u = a | b;
    require_same_region(u, [&](const V& p) { return contains(a, p) || contains(b, p); });
    REQUIRE(signed_area(u) == 36 + 36 - 4);
    REQUIRE(u.rings.size() == 1);
}

TEST_CASE("intersection") {
    const P a(box(0, 0, 6, 6)), b(box(4, 4, 10, 10));
    const P i = a & b;
    require_same_region(i, [&](const V& p) { return contains(a, p) && contains(b, p); });
    REQUIRE(signed_area(i) == 4);

    // disjoint boxes intersect in nothing
    const P d = P(box(0, 0, 4, 4)) & P(box(6, 0, 10, 4));
    REQUIRE(d.rings.empty());
    REQUIRE(d.is_empty());
}

TEST_CASE("difference") {
    const P a(box(0, 0, 6, 6)), b(box(4, 4, 10, 10));
    const P d = a - b;
    require_same_region(d, [&](const V& p) { return contains(a, p) && !contains(b, p); });
    REQUIRE(signed_area(d) == 36 - 4);

    // subtracting a superset leaves nothing
    REQUIRE((P(box(2, 2, 4, 4)) - P(box(0, 0, 10, 10))).rings.empty());
    // subtracting a disjoint region changes nothing
    const P e = P(box(0, 0, 4, 4)) - P(box(6, 0, 10, 4));
    REQUIRE(signed_area(e) == 16);
}

TEST_CASE("symmetric difference") {
    const P a(box(0, 0, 6, 6)), b(box(4, 4, 10, 10));
    const P s = a ^ b;
    require_same_region(s, [&](const V& p) { return contains(a, p) != contains(b, p); });
    REQUIRE(signed_area(s) == 36 + 36 - 2 * 4);

    // a region with itself cancels
    REQUIRE((a ^ a).rings.empty());
}

TEST_CASE("a hole is punched by difference") {
    const P a(box(0, 0, 10, 10)), b(box(3, 3, 7, 7));
    const P d = a - b;
    require_same_region(d, [&](const V& p) { return contains(a, p) && !contains(b, p); });
    REQUIRE(signed_area(d) == 100 - 16);
    REQUIRE(d.rings.size() == 2);
    // the outer ring winds counter clockwise and the hole the other way
    REQUIRE(signed_area2(d.rings[0]) + signed_area2(d.rings[1]) == 2 * (100 - 16));
    REQUIRE((is_ccw(d.rings[0]) != is_ccw(d.rings[1])));
}

TEST_CASE("operations on a region that already has a hole") {
    R hole = box(3, 3, 7, 7);
    reverse(hole);
    const P ring({box(0, 0, 10, 10), hole});
    REQUIRE(signed_area(ring) == 84);

    // filling the hole back in
    const P filled = ring | P(box(3, 3, 7, 7));
    require_same_region(filled, [&](const V& p) { return contains(ring, p) || contains(P(box(3, 3, 7, 7)), p); });
    REQUIRE(signed_area(filled) == 100);

    // intersecting with a box that reaches into the hole
    const P cut = ring & P(box(5, 5, 15, 15));
    require_same_region(cut, [&](const V& p) { return contains(ring, p) && contains(P(box(5, 5, 15, 15)), p); });
    REQUIRE(signed_area(cut) == 25 - 4);
}

TEST_CASE("touching boxes") {
    // sharing a whole edge
    const P a(box(0, 0, 5, 5)), b(box(5, 0, 10, 5));
    const P u = a | b;
    require_same_region(u, [&](const V& p) { return contains(a, p) || contains(b, p); });
    REQUIRE(signed_area(u) == 50);
    REQUIRE(u.rings.size() == 1);

    // the shared edge is the whole intersection, which encloses no area
    const P i = a & b;
    REQUIRE(signed_area(i) == 0);

    // sharing part of an edge
    const P c(box(0, 0, 5, 5)), d(box(5, 2, 10, 8));
    const P u2 = c | d;
    require_same_region(u2, [&](const V& p) { return contains(c, p) || contains(d, p); });
    REQUIRE(signed_area(u2) == 25 + 30);

    // meeting at a single corner
    const P e(box(0, 0, 5, 5)), f(box(5, 5, 10, 10));
    const P u3 = e | f;
    require_same_region(u3, [&](const V& p) { return contains(e, p) || contains(f, p); });
    REQUIRE(signed_area(u3) == 50);
}

TEST_CASE("identical operands") {
    const P a(box(0, 0, 5, 5));
    REQUIRE(signed_area(a | a) == 25);
    REQUIRE(signed_area(a & a) == 25);
    REQUIRE((a - a).rings.empty());
    REQUIRE((a ^ a).rings.empty());
    require_same_region(a | a, [&](const V& p) { return contains(a, p); });
    require_same_region(a & a, [&](const V& p) { return contains(a, p); });
}

TEST_CASE("empty and whole plane operands") {
    const P a(box(0, 0, 5, 5));
    const P empty;
    const P plane = ~empty;

    REQUIRE(signed_area(a | empty) == 25);
    REQUIRE((a & empty).rings.empty());
    REQUIRE((a & empty).is_empty());
    REQUIRE(signed_area(a - empty) == 25);
    REQUIRE(signed_area(a ^ empty) == 25);

    // anything united with the plane is the plane, and intersected with it is itself
    REQUIRE((a | plane).is_whole_plane());
    require_same_region(a & plane, [&](const V& p) { return contains(a, p); });
    // the plane minus a box is the box inverted
    const P inv = plane - a;
    require_same_region(inv, [&](const V& p) { return !contains(a, p); });
    REQUIRE(inv.complement);
}

TEST_CASE("unbounded operands") {
    const P a(box(0, 0, 10, 10));
    const P b(box(4, 4, 14, 14));
    const P na = ~a;

    // De Morgan, checked by membership
    require_same_region(na | ~b, [&](const V& p) { return !contains(a, p) || !contains(b, p); });
    require_same_region(na & ~b, [&](const V& p) { return !contains(a, p) && !contains(b, p); });
    REQUIRE((na & ~b).complement);
    REQUIRE((na | ~b).complement);

    // an unbounded region intersected with a bounded one is bounded
    const P c = na & b;
    require_same_region(c, [&](const V& p) { return !contains(a, p) && contains(b, p); });
    REQUIRE(!c.complement);
    REQUIRE(signed_area(c) == 100 - 36);
}

TEST_CASE("non convex operands") {
    // an L shape
    const R l{V(0, 0), V(10, 0), V(10, 4), V(4, 4), V(4, 10), V(0, 10)};
    const P a(l);
    REQUIRE(signed_area(a) == 64);
    const P b(box(2, 2, 8, 8));

    const P i = a & b;
    require_same_region(i, [&](const V& p) { return contains(a, p) && contains(b, p); });
    const P u = a | b;
    require_same_region(u, [&](const V& p) { return contains(a, p) || contains(b, p); });
    const P d = a - b;
    require_same_region(d, [&](const V& p) { return contains(a, p) && !contains(b, p); });
    const P s = a ^ b;
    require_same_region(s, [&](const V& p) { return contains(a, p) != contains(b, p); });

    // the four results have to account for the same total area
    REQUIRE(signed_area(i) + signed_area(s) == signed_area(u));
    REQUIRE(signed_area(d) + signed_area(i) == signed_area(a));
}

TEST_CASE("a result that splits into two pieces") {
    // a bar crossed by a taller box, subtracting the middle
    const P bar(box(0, 0, 12, 4));
    const P mid(box(5, -2, 7, 6));
    const P d = bar - mid;
    require_same_region(d, [&](const V& p) { return contains(bar, p) && !contains(mid, p); });
    REQUIRE(d.rings.size() == 2);
    REQUIRE(signed_area(d) == 12 * 4 - 2 * 4);
}

TEST_CASE("non integral coordinates") {
    const P a(R{V(0, 0), V(rational(7, 2), 0), V(rational(7, 2), rational(7, 2)), V(0, rational(7, 2))});
    const P b(R{V(rational(3, 2), rational(3, 2)), V(5, rational(3, 2)), V(5, 5), V(rational(3, 2), 5)});
    const P i = a & b;
    require_same_region(i, [&](const V& p) { return contains(a, p) && contains(b, p); }, -1, 6);
    REQUIRE(signed_area(i) == rational(4, 1)); // a 2x2 square from 3/2 to 7/2
}

TEST_CASE("De Morgan and idempotence over a set of shapes") {
    const std::vector<P> shapes = {
        P(box(0, 0, 6, 6)),
        P(box(3, 3, 9, 9)),
        P(box(-2, 1, 2, 5)),
        P(R{V(0, 0), V(8, 0), V(4, 6)}),
        P({box(0, 0, 10, 10), [] { R h = box(4, 4, 6, 6); reverse(h); return h; }()}),
    };
    for (const P& a : shapes)
        for (const P& b : shapes) {
            // ~(a | b) == ~a & ~b, as sets of points
            const P lhs = ~(a | b), rhs = ~a & ~b;
            for (int i = -6; i <= 22; i++)
                for (int j = -6; j <= 22; j++) {
                    const V p(rational(i, 2) + rational(1, 4), rational(j, 2) + rational(1, 4));
                    REQUIRE(contains(lhs, p) == contains(rhs, p));
                }
            // a - b == a & ~b
            const P d1 = a - b, d2 = a & ~b;
            for (int i = -6; i <= 22; i++)
                for (int j = -6; j <= 22; j++) {
                    const V p(rational(i, 2) + rational(1, 4), rational(j, 2) + rational(1, 4));
                    REQUIRE(contains(d1, p) == contains(d2, p));
                }
        }
}

TEST_CASE("randomized boolean operations against pointwise membership") {
    // Boxes on a small integer grid put many edges of one operand exactly on the line that leaves
    // the midpoint of an edge of the other, which is the configuration __safe_step() has to get
    // right: the ray it casts to find a safe distance ignores edges parallel to itself.
    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<int> coord(0, 6), count(1, 3);
    auto rand_region = [&] {
        std::vector<R> rings;
        for (int n = count(rng); n > 0; n--) {
            int x0 = coord(rng), x1 = coord(rng), y0 = coord(rng), y1 = coord(rng);
            if (x0 == x1 || y0 == y1)
                continue;
            if (x0 > x1) std::swap(x0, x1);
            if (y0 > y1) std::swap(y0, y1);
            rings.push_back(box(x0, y0, x1, y1));
        }
        return P(rings);
    };

    for (int iter = 0; iter < 200; iter++) {
        const P a = rand_region(), b = rand_region();
        const P u = a | b, i = a & b, d = a - b, x = a ^ b;
        for (int gx = -1; gx <= 7; gx++)
            for (int gy = -1; gy <= 7; gy++) {
                // quarter offsets stay off every edge, and off every midpoint normal as well
                const V p(rational(4 * gx + 1, 4), rational(4 * gy + 1, 4));
                const bool ca = contains(a, p), cb = contains(b, p);
                REQUIRE(contains(u, p) == (ca || cb));
                REQUIRE(contains(i, p) == (ca && cb));
                REQUIRE(contains(d, p) == (ca && !cb));
                REQUIRE(contains(x, p) == (ca != cb));
            }
    }
}

TEST_CASE("a vertex shared by two lobes gives two rings") {
    // The two squares meet only at (1,1), where four fragments come together. Which outgoing
    // fragment continues the ring decides whether the result is two rings or one ring that visits
    // the shared vertex twice; both bound the same region, but only the first is a boundary.
    const P a(box(0, 1, 1, 2)), b(box(1, 0, 2, 1));
    const P u = a | b;
    require_same_region(u, [&](const V& p) { return contains(a, p) || contains(b, p); });
    REQUIRE(signed_area(u) == 2);
    REQUIRE(u.rings.size() == 2);
    for (const R& r : u.rings) {
        REQUIRE(r.size() == 4);
        REQUIRE(signed_area2(r) == 2);
    }

    // the same vertex on the other diagonal, which the fragment order happened to get right
    const P c(box(0, 0, 1, 1)), d(box(1, 1, 2, 2));
    REQUIRE((c | d).rings.size() == 2);

    // and a difference that leaves two lobes touching at (2,2)
    const P e = (P(box(0, 0, 4, 4)) - P(box(0, 0, 2, 2))) - P(box(2, 2, 4, 4));
    require_same_region(e, [&](const V& p) {
        return contains(P(box(0, 0, 4, 4)), p) && !contains(P(box(0, 0, 2, 2)), p) && !contains(P(box(2, 2, 4, 4)), p);
    });
    REQUIRE(signed_area(e) == 8);
    REQUIRE(e.rings.size() == 2);
}
