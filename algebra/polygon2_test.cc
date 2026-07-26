#include "algebra/polygon2.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"

using P = MultiPolygon2<rational>;
using R = Ring2<rational>;
using V = Vec2<rational>;

// counter clockwise unit-ish square [0,10]^2
static R square(int lo = 0, int hi = 10) {
    return R{V(lo, lo), V(hi, lo), V(hi, hi), V(lo, hi)};
}

TEST_CASE("signed_area2 and orientation") {
    const R ccw = square();
    REQUIRE(signed_area2(ccw) == 200);
    REQUIRE(signed_area(ccw) == 100);
    REQUIRE(is_ccw(ccw));

    R cw = ccw;
    reverse(cw);
    REQUIRE(signed_area2(cw) == -200);
    REQUIRE(!is_ccw(cw));

    // a triangle, and the area is independent of where the ring starts
    R t{V(0, 0), V(4, 0), V(0, 3)};
    REQUIRE(signed_area(t) == 6);
    R t2{V(4, 0), V(0, 3), V(0, 0)};
    REQUIRE(signed_area(t2) == 6);

    // degenerate rings enclose nothing
    REQUIRE(signed_area2(R{}) == 0);
    REQUIRE(signed_area2(R{V(1, 1)}) == 0);
    REQUIRE(signed_area2(R{V(1, 1), V(2, 2)}) == 0);
    // a ring that doubles back has zero area
    REQUIRE(signed_area2(R{V(0, 0), V(5, 0), V(0, 0), V(5, 0)}) == 0);

    // non integral coordinates
    R h{V(0, 0), V(rational(1, 2), 0), V(0, rational(1, 2))};
    REQUIRE(signed_area(h) == rational(1, 8));
}

TEST_CASE("signed_area of a region with a hole") {
    R hole = square(2, 8);
    reverse(hole); // a hole winds the other way
    P a({square(), hole});
    REQUIRE(signed_area(a) == 100 - 36);

    // an unbounded region has no finite area
    REQUIRE_THROWS(signed_area(~a));
}

TEST_CASE("winding_number") {
    const P a(square());
    REQUIRE(winding_number(a, V(5, 5)) == 1);
    REQUIRE(winding_number(a, V(15, 5)) == 0);
    REQUIRE(winding_number(a, V(-1, 5)) == 0);
    REQUIRE(winding_number(a, V(5, 15)) == 0);

    // reversing the ring negates the winding number, it does not change what is enclosed
    R cw = square();
    reverse(cw);
    REQUIRE(winding_number(P(cw), V(5, 5)) == -1);

    // a point level with a vertex must not be counted twice by the two edges meeting there
    R diamond{V(0, 0), V(5, -5), V(10, 0), V(5, 5)};
    REQUIRE(winding_number(P(diamond), V(5, 0)) == 1);
    REQUIRE(winding_number(P(diamond), V(-1, 0)) == 0);
    REQUIRE(winding_number(P(diamond), V(11, 0)) == 0);

    // two nested counter clockwise rings wind twice
    P twice({square(0, 10), square(2, 8)});
    REQUIRE(winding_number(twice, V(5, 5)) == 2);
    REQUIRE(winding_number(twice, V(1, 1)) == 1);
}

TEST_CASE("on_boundary") {
    const P a(square());
    REQUIRE(on_boundary(a, V(0, 0)));        // vertex
    REQUIRE(on_boundary(a, V(10, 10)));
    REQUIRE(on_boundary(a, V(5, 0)));        // edge interior
    REQUIRE(on_boundary(a, V(0, 5)));
    REQUIRE(on_boundary(a, V(10, 5)));
    REQUIRE(on_boundary(a, V(5, 10)));       // the closing edge
    REQUIRE(on_boundary(a, V(rational(1, 2), 0)));
    REQUIRE(!on_boundary(a, V(5, 5)));       // interior
    REQUIRE(!on_boundary(a, V(11, 0)));      // collinear with an edge but outside it
    REQUIRE(!on_boundary(a, V(-1, 0)));

    // a single point ring is its own boundary
    REQUIRE(on_boundary(R{V(1, 2)}, V(1, 2)));
    REQUIRE(!on_boundary(R{V(1, 2)}, V(1, 3)));
    // a two point ring is a segment
    REQUIRE(on_boundary(R{V(0, 0), V(4, 4)}, V(2, 2)));
    REQUIRE(!on_boundary(R{V(0, 0), V(4, 4)}, V(2, 3)));
}

TEST_CASE("contains") {
    const P a(square());
    REQUIRE(contains(a, V(5, 5)));
    REQUIRE(contains(a, V(rational(1, 100), rational(1, 100))));
    REQUIRE(!contains(a, V(11, 5)));
    REQUIRE(!contains(a, V(5, -1)));
    // the boundary belongs to the closed region
    REQUIRE(contains(a, V(0, 0)));
    REQUIRE(contains(a, V(5, 10)));

    // a hole is not contained, but its boundary still is
    R hole = square(2, 8);
    reverse(hole);
    const P withhole({square(), hole});
    REQUIRE(contains(withhole, V(1, 1)));
    REQUIRE(!contains(withhole, V(5, 5)));
    REQUIRE(contains(withhole, V(2, 5)));   // on the hole boundary
    REQUIRE(contains(withhole, V(8, 8)));

    // orientation of the outer ring does not decide membership, the winding rule does
    R cw = square();
    reverse(cw);
    REQUIRE(contains(P(cw), V(5, 5)));
}

TEST_CASE("inversion") {
    const P a(square());
    const P b = ~a;
    REQUIRE(b.complement);
    REQUIRE(b.rings == a.rings);   // inversion shares the boundary, it does not move it

    // the interior and the exterior swap
    REQUIRE(contains(a, V(5, 5)));
    REQUIRE(!contains(b, V(5, 5)));
    REQUIRE(!contains(a, V(15, 5)));
    REQUIRE(contains(b, V(15, 5)));
    // and the boundary stays in both, since both are closed regions
    REQUIRE(contains(a, V(0, 0)));
    REQUIRE(contains(b, V(0, 0)));

    // inversion is an involution
    REQUIRE(~b == a);
    REQUIRE(~~a == a);

    // the empty region and the whole plane are each other's inverse
    const P empty;
    REQUIRE(empty.is_empty());
    REQUIRE(!empty.is_whole_plane());
    REQUIRE(!contains(empty, V(0, 0)));
    const P plane = ~empty;
    REQUIRE(plane.is_whole_plane());
    REQUIRE(!plane.is_empty());
    REQUIRE(contains(plane, V(0, 0)));
    REQUIRE(contains(plane, V(-1000000, rational(1, 3))));
    REQUIRE(~plane == empty);
}

TEST_CASE("bounding_box") {
    V min, max;
    bounding_box(P(square(-3, 7)), min, max);
    REQUIRE(min == V(-3, -3));
    REQUIRE(max == V(7, 7));

    R t{V(0, 0), V(rational(1, 2), 4), V(-2, 1)};
    bounding_box(P(t), min, max);
    REQUIRE(min == V(-2, 0));
    REQUIRE(max == V(rational(1, 2), 4));

    REQUIRE_THROWS(bounding_box(~P(square()), min, max));
    REQUIRE_THROWS(bounding_box(P(), min, max));
}

TEST_CASE("simplify") {
    // a repeated vertex, including across the closing edge
    P a(R{V(0, 0), V(0, 0), V(10, 0), V(10, 10), V(0, 10), V(0, 0)});
    simplify(a);
    REQUIRE(a.rings.size() == 1);
    REQUIRE(a.rings[0] == square());

    // a vertex in the middle of a straight edge
    P b(R{V(0, 0), V(5, 0), V(10, 0), V(10, 10), V(0, 10)});
    simplify(b);
    REQUIRE(b.rings[0].size() == 4);
    REQUIRE(signed_area(b) == 100);

    // several collinear vertices in a row, so removing one exposes the next
    P c(R{V(0, 0), V(2, 0), V(4, 0), V(6, 0), V(10, 0), V(10, 10), V(0, 10)});
    simplify(c);
    REQUIRE(c.rings[0].size() == 4);
    REQUIRE(signed_area(c) == 100);

    // rings that enclose nothing are dropped: too few vertices, or every vertex on one line.
    // the last one doubles back on itself, so no vertex lies between its neighbours and the
    // collinear rule above cannot remove any of them.
    P d({R{V(1, 1)}, R{V(0, 0), V(4, 4)}, square(), R{V(0, 0), V(5, 0), V(0, 0), V(5, 0)},
         R{V(0, 0), V(1, 1), V(3, 3)}});
    simplify(d);
    REQUIRE(d.rings.size() == 1);
    REQUIRE(d.rings[0] == square());

    // simplify does not change membership
    P e(R{V(0, 0), V(5, 0), V(10, 0), V(10, 10), V(0, 10), V(0, 5)});
    P f = e;
    simplify(f);
    for (int x = -1; x <= 11; x++)
        for (int y = -1; y <= 11; y++)
            REQUIRE(contains(e, V(x, y)) == contains(f, V(x, y)));

    // the complement flag survives
    P g = ~P(square());
    simplify(g);
    REQUIRE(g.complement);
}
