#include "algebra/polygon2_arc_boolean.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"

using AR = ArcRegion<rational>;
using AP = ArcPolygon2<rational>;
using Ring = ArcRing2<rational>;
using AV = ArcVertex<rational>;
using V = Vec2<rational>;

static AR disk(int cx, int cy, int r) { return AR(circle_ring<rational>(V(cx, cy), r)); }
static AR box(int x0, int y0, int x1, int y1) {
    return AR(Ring{AV(V(x0, y0)), AV(V(x1, y0)), AV(V(x1, y1)), AV(V(x0, y1))});
}
static bool in_disk(const V& p, int cx, int cy, int r) {
    return (p.x - cx) * (p.x - cx) + (p.y - cy) * (p.y - cy) <= r * r;
}
static bool in_box(const V& p, int x0, int y0, int x1, int y1) {
    return x0 <= p.x && p.x <= x1 && y0 <= p.y && p.y <= y1;
}

// Sampling on quarters avoids sitting exactly on a boundary in most cases, and where it does the
// expected predicate is written to include the boundary too, since every region here is closed.
static void require_same(const AR& got, auto&& want, int lo = -8, int hi = 8) {
    for (int i = 4 * lo; i <= 4 * hi; i++)
        for (int j = 4 * lo; j <= 4 * hi; j++) {
            const V p(rational(i, 4), rational(j, 4));
            REQUIRE(contains(got, p) == want(p));
        }
}

TEST_CASE("union and intersection of two disks") {
    const AR a = disk(0, 0, 3), b = disk(4, 0, 3);
    require_same(a | b, [](const V& p) { return in_disk(p, 0, 0, 3) || in_disk(p, 4, 0, 3); });
    require_same(a & b, [](const V& p) { return in_disk(p, 0, 0, 3) && in_disk(p, 4, 0, 3); });
    // the lens shaped intersection is non empty, and the disks are not nested
    REQUIRE(contains(a & b, V(2, 0)));
    REQUIRE(!contains(a & b, V(0, 0)));
    REQUIRE(!contains(a & b, V(4, 0)));
    REQUIRE(contains(a | b, V(0, 0)));
    REQUIRE(contains(a | b, V(4, 0)));
}

TEST_CASE("disjoint disks intersect in nothing") {
    const AR a = disk(0, 0, 2), b = disk(6, 0, 2);
    for (int i = -16; i <= 32; i++)
        for (int j = -16; j <= 16; j++)
            REQUIRE(!contains(a & b, V(rational(i, 2), rational(j, 2))));
    require_same(a | b, [](const V& p) { return in_disk(p, 0, 0, 2) || in_disk(p, 6, 0, 2); }, -4, 10);
}

TEST_CASE("difference and symmetric difference") {
    const AR a = disk(0, 0, 4), b = disk(3, 0, 2);
    require_same(a - b, [](const V& p) { return in_disk(p, 0, 0, 4) && !in_disk(p, 3, 0, 2); });
    require_same(b - a, [](const V& p) { return in_disk(p, 3, 0, 2) && !in_disk(p, 0, 0, 4); });
    require_same(a ^ b, [](const V& p) { return in_disk(p, 0, 0, 4) != in_disk(p, 3, 0, 2); });

    // a disk minus a concentric smaller one is an annulus
    const AR ring = disk(0, 0, 4) - disk(0, 0, 2);
    REQUIRE(contains(ring, V(3, 0)));
    REQUIRE(!contains(ring, V(1, 0)));
    REQUIRE(!contains(ring, V(0, 0)));
    REQUIRE(!contains(ring, V(5, 0)));
}

TEST_CASE("mixing arcs with straight edges") {
    const AR d = disk(0, 0, 3), q = box(0, 0, 5, 5);
    require_same(d & q, [](const V& p) { return in_disk(p, 0, 0, 3) && in_box(p, 0, 0, 5, 5); });
    require_same(d | q, [](const V& p) { return in_disk(p, 0, 0, 3) || in_box(p, 0, 0, 5, 5); });
    require_same(d - q, [](const V& p) { return in_disk(p, 0, 0, 3) && !in_box(p, 0, 0, 5, 5); });
    require_same(q - d, [](const V& p) { return in_box(p, 0, 0, 5, 5) && !in_disk(p, 0, 0, 3); });
}

TEST_CASE("complement is a strict negation") {
    const AR a = disk(0, 0, 3);
    // ~ on a region negates contains(), so unlike ArcPolygon2's flag flipping ~ it does not leave
    // the boundary belonging to both sides. That is what makes this a proper Boolean algebra.
    REQUIRE(contains(a, V(3, 0)));
    REQUIRE(!contains(~a, V(3, 0)));
    REQUIRE(contains(a, V(0, 3)));
    REQUIRE(!contains(~a, V(0, 3)));
    require_same(~a, [](const V& p) { return !in_disk(p, 0, 0, 3); });
    REQUIRE((~a).leaf_count() == 1);
    // and it is an involution, without stacking two nodes
    require_same(~~a, [](const V& p) { return in_disk(p, 0, 0, 3); });
    REQUIRE((~~a).kind == AR::Kind::Leaf);

    // complementing a combination wraps it, and unwraps again
    const AR u = disk(0, 0, 2) | disk(4, 0, 2);
    REQUIRE((~u).kind == AR::Kind::Complement);
    REQUIRE((~~u).kind == AR::Kind::Union);
    // ~ of a combination is a strict negation of contains(), unlike ~ of a leaf
    require_same(~u, [&](const V& p) { return !contains(u, p); }, -4, 8);
}

TEST_CASE("De Morgan and the difference identity") {
    const AR a = disk(0, 0, 3), b = disk(3, 2, 3), c = box(-1, -1, 4, 2);
    for (const AR& x : {a, b, c})
        for (const AR& y : {a, b, c}) {
            for (int i = -24; i <= 32; i++)
                for (int j = -24; j <= 32; j++) {
                    const V p(rational(i, 4), rational(j, 4));

                    // These four hold everywhere, since none of them complements a whole
                    // combination -- the shared boundary cancels on both sides.
                    REQUIRE(contains(x ^ y, p) == contains((x - y) | (y - x), p));
                    REQUIRE(contains(x | y, p) == (contains(x, p) || contains(y, p)));
                    REQUIRE(contains(x & y, p) == (contains(x, p) && contains(y, p)));
                    REQUIRE(contains(x - y, p) == (contains(x, p) && !contains(y, p)));

                    // With ~ a strict negation these hold at every point, boundary included.
                    REQUIRE(contains(~(x | y), p) == contains(~x & ~y, p));
                    REQUIRE(contains(~(x & y), p) == contains(~x | ~y, p));
                    REQUIRE(contains(x - y, p) == contains(x & ~y, p));
                }
        }
}

TEST_CASE("nesting to several levels stays exact") {
    const AR a = disk(0, 0, 4), b = disk(3, 0, 3), c = disk(0, 3, 3), d = box(-2, -2, 2, 2);
    const AR e = ((a - b) | (c & d)) ^ disk(1, 1, 1);
    REQUIRE(e.leaf_count() == 5);
    require_same(e, [](const V& p) {
        const bool lhs = (in_disk(p, 0, 0, 4) && !in_disk(p, 3, 0, 3)) ||
                         (in_disk(p, 0, 3, 3) && in_box(p, -2, -2, 2, 2));
        return lhs != in_disk(p, 1, 1, 1);
    });
}

TEST_CASE("idempotence and absorption") {
    const AR a = disk(0, 0, 3), b = box(0, 0, 4, 4);
    for (int i = -20; i <= 24; i++)
        for (int j = -20; j <= 24; j++) {
            const V p(rational(i, 4), rational(j, 4));
            REQUIRE(contains(a | a, p) == contains(a, p));
            REQUIRE(contains(a & a, p) == contains(a, p));
            REQUIRE(!contains(a - a, p));
            REQUIRE(!contains(a ^ a, p));
            REQUIRE(contains(a | (a & b), p) == contains(a, p));   // absorption
            REQUIRE(contains(a & (a | b), p) == contains(a, p));
        }
}

TEST_CASE("an arc region with a hole, combined") {
    // a disk with a smaller disk removed, then intersected with a half plane like box
    const AR annulus = disk(0, 0, 4) - disk(0, 0, 2);
    const AR half = box(0, -5, 5, 5);
    require_same(annulus & half, [](const V& p) {
        return in_disk(p, 0, 0, 4) && !in_disk(p, 0, 0, 2) && in_box(p, 0, -5, 5, 5);
    });
    REQUIRE(contains(annulus & half, V(3, 0)));
    REQUIRE(!contains(annulus & half, V(-3, 0)));
    REQUIRE(!contains(annulus & half, V(1, 0)));
}

TEST_CASE("area_bounds brackets the true area") {
    // a disk of radius 2 has area 4*pi == 12.566...
    const AR a = disk(0, 0, 2);
    rational lower, undecided;
    area_bounds(a, V(-3, -3), V(3, 3), 7, lower, undecided);
    REQUIRE(lower <= rational(4) * rational(314160, 100000));
    REQUIRE(lower + undecided >= rational(4) * rational(314159, 100000));
    REQUIRE(undecided > 0); // a curved boundary always leaves something straddling

    // deeper subdivision cannot make the bracket worse
    rational lower2, undecided2;
    area_bounds(a, V(-3, -3), V(3, 3), 9, lower2, undecided2);
    REQUIRE(lower2 >= lower);
    REQUIRE(undecided2 <= undecided);

    // a straight edged region is bracketed too, and its exact area is known independently
    const AR b = box(0, 0, 4, 2);
    area_bounds(b, V(-1, -1), V(5, 3), 8, lower, undecided);
    REQUIRE(lower <= 8);
    REQUIRE(lower + undecided >= 8);
}

TEST_CASE("combining shares the subtree instead of copying it") {
    // A node holds its two children behind shared_ptr, so combining does not duplicate the tree
    // below the operands: the leaves of the result are the very leaf objects of the inputs.
    const AR leaf = disk(0, 0, 1);
    const AR a = leaf | leaf;
    const AR b = a | a;
    REQUIRE(b.leaf_count() == 4);
    // each operand is copied into a node of its own, and that copy is one node deep:
    REQUIRE(b.a.get() != b.b.get());
    // both copies of `a` point at the very same children, rather than at copies of them
    REQUIRE(b.a->a.get() == b.b->a.get());
    REQUIRE(b.a->b.get() == b.b->b.get());
    // which are the leaves the inputs were built from, so no leaf was duplicated by the second |
    REQUIRE(b.a->a->kind == AR::Kind::Leaf);
}
