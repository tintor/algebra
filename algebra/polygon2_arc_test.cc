#include "algebra/polygon2_arc.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"

using A = ArcPolygon2<rational>;
using AR = ArcRing2<rational>;
using AV = ArcVertex<rational>;
using V = Vec2<rational>;
static const rational half(1, 2);

TEST_CASE("arc centre, radius and midpoint are rational") {
    // a half circle: chord (-1,0)-(1,0) with the arc above it
    const V a(-1, 0), b(1, 0);
    REQUIRE(arc_midpoint(a, b, rational(1)) == V(0, 1));
    REQUIRE(arc_center(a, b, rational(1)) == V(0, 0));
    REQUIRE(arc_radius2(a, b, rational(1)) == 1);
    // the other way round the arc is below
    REQUIRE(arc_midpoint(a, b, rational(-1)) == V(0, -1));
    REQUIRE(arc_center(a, b, rational(-1)) == V(0, 0));

    // a quarter circle has bulge tan(pi/8), which is irrational, so use a major arc with a
    // rational bulge instead: chord (0,0)-(4,0) with bulge 2
    const V c(0, 0), d(4, 0);
    REQUIRE(arc_midpoint(c, d, rational(2)) == V(2, 4));
    REQUIRE(arc_center(c, d, rational(2)) == V(2, rational(3, 2)));
    REQUIRE(arc_radius2(c, d, rational(2)) == rational(25, 4));
    // the centre being on the same side as the arc midpoint is what makes it a major arc
    REQUIRE(cross(d - c, arc_center(c, d, rational(2)) - c) > 0);
    REQUIRE(cross(d - c, arc_midpoint(c, d, rational(2)) - c) > 0);

    // a shallow arc, where the centre is far away on the other side
    REQUIRE(arc_center(a, b, half) == V(0, rational(-3, 4)));
    REQUIRE(arc_radius2(a, b, half) == rational(25, 16));

    REQUIRE_THROWS(arc_center(a, b, rational(0)));
}

TEST_CASE("on_arc") {
    const V a(-1, 0), b(1, 0);
    // upper half circle
    REQUIRE(on_arc(a, b, rational(1), V(0, 1)));
    REQUIRE(on_arc(a, b, rational(1), a));
    REQUIRE(on_arc(a, b, rational(1), b));
    REQUIRE(on_arc(a, b, rational(1), V(rational(3, 5), rational(4, 5))));  // 3-4-5 on the circle
    REQUIRE(!on_arc(a, b, rational(1), V(rational(3, 5), rational(-4, 5)))); // on the circle, lower
    REQUIRE(!on_arc(a, b, rational(1), V(0, -1)));
    REQUIRE(!on_arc(a, b, rational(1), V(0, 0)));      // the centre is not on the arc
    REQUIRE(!on_arc(a, b, rational(1), V(0, 2)));      // off the circle

    // lower half circle
    REQUIRE(on_arc(a, b, rational(-1), V(0, -1)));
    REQUIRE(!on_arc(a, b, rational(-1), V(0, 1)));

    // a straight edge falls back to the segment test
    REQUIRE(on_arc(a, b, rational(0), V(0, 0)));
    REQUIRE(on_arc(a, b, rational(0), V(half, 0)));
    REQUIRE(!on_arc(a, b, rational(0), V(2, 0)));
    REQUIRE(!on_arc(a, b, rational(0), V(0, 1)));

    // a major arc covers points on both sides of its chord
    const V c(0, 0), d(4, 0);
    const rational bulge(2);
    REQUIRE(on_arc(c, d, bulge, V(2, 4)));           // the arc midpoint, above the chord
    REQUIRE(!on_arc(c, d, bulge, V(2, -1)));         // the far side of the circle, below
    REQUIRE(on_arc(c, d, bulge, arc_center(c, d, bulge) + V(rational(-5, 2), 0)));  // left extreme
    REQUIRE(on_arc(c, d, bulge, arc_center(c, d, bulge) + V(rational(5, 2), 0)));   // right extreme
}

TEST_CASE("a disk from two half circle edges") {
    const AR ring = circle_ring<rational>(V(0, 0), 2);
    REQUIRE(ring.size() == 2);
    const A disk(ring);

    // the boundary is the circle of radius 2
    REQUIRE(on_boundary(disk, V(2, 0)));
    REQUIRE(on_boundary(disk, V(-2, 0)));
    REQUIRE(on_boundary(disk, V(0, 2)));
    REQUIRE(on_boundary(disk, V(0, -2)));
    REQUIRE(on_boundary(disk, V(rational(6, 5), rational(8, 5)))); // 3-4-5 scaled
    REQUIRE(!on_boundary(disk, V(0, 0)));
    REQUIRE(!on_boundary(disk, V(0, 3)));

    // inside
    REQUIRE(contains(disk, V(0, 1)));
    REQUIRE(contains(disk, V(0, -1)));
    REQUIRE(contains(disk, V(1, 0)));
    REQUIRE(contains(disk, V(-1, 0)));
    REQUIRE(contains(disk, V(0, 0)));            // the centre lies on both chords
    REQUIRE(contains(disk, V(1, 1)));
    REQUIRE(contains(disk, V(rational(-7, 5), rational(1, 3))));
    // outside
    REQUIRE(!contains(disk, V(0, 3)));
    REQUIRE(!contains(disk, V(3, 0)));
    REQUIRE(!contains(disk, V(2, 2)));
    REQUIRE(!contains(disk, V(-3, -3)));

    // every point of the plane agrees with the exact disk test
    for (int i = -6; i <= 6; i++)
        for (int j = -6; j <= 6; j++) {
            const V p(rational(i, 2), rational(j, 2));
            REQUIRE(contains(disk, p) == (p.x * p.x + p.y * p.y <= 4));
        }

    // the chord polygon of a two edge circle is degenerate, so it has no area of its own
    REQUIRE(signed_area_chords(disk) == 0);
}

TEST_CASE("an arc bulging outwards adds area") {
    // a square with its bottom edge replaced by a half circle below it
    const A a(AR{AV(V(0, 0), rational(-1)), AV(V(4, 0)), AV(V(4, 4)), AV(V(0, 4))});
    REQUIRE(signed_area_chords(a) == 16);

    REQUIRE(contains(a, V(2, 2)));       // inside the square
    REQUIRE(contains(a, V(2, -1)));      // inside the added half disk
    REQUIRE(contains(a, V(2, -2)));      // on its boundary
    REQUIRE(!contains(a, V(2, -3)));     // beyond it
    REQUIRE(!contains(a, V(2, 5)));
    // the half disk has centre (2,0) and radius 2
    REQUIRE(contains(a, V(rational(4, 5), rational(-8, 5))));
    REQUIRE(!contains(a, V(rational(1, 2), -2)));

    for (int i = -2; i <= 12; i++)
        for (int j = -12; j <= 12; j++) {
            const V p(rational(i, 2), rational(j, 2));
            const bool in_square = 0 <= p.x && p.x <= 4 && 0 <= p.y && p.y <= 4;
            const bool in_half = p.y <= 0 && (p.x - 2) * (p.x - 2) + p.y * p.y <= 4;
            REQUIRE(contains(a, p) == (in_square || in_half));
        }
}

TEST_CASE("an arc bulging inwards removes area") {
    // the same square, with the bottom edge bulging up into it
    const A a(AR{AV(V(0, 0), rational(1)), AV(V(4, 0)), AV(V(4, 4)), AV(V(0, 4))});
    REQUIRE(signed_area_chords(a) == 16);

    REQUIRE(!contains(a, V(2, 1)));      // inside the removed half disk
    REQUIRE(contains(a, V(2, 3)));       // above it
    REQUIRE(contains(a, V(2, 2)));       // on its boundary
    REQUIRE(!contains(a, V(rational(1, 4), rational(1, 4)))); // also inside the removed half disk

    for (int i = 0; i <= 8; i++)
        for (int j = 0; j <= 8; j++) {
            const V p(rational(i, 2), rational(j, 2));
            const bool in_square = 0 <= p.x && p.x <= 4 && 0 <= p.y && p.y <= 4;
            const bool in_half = p.y >= 0 && (p.x - 2) * (p.x - 2) + p.y * p.y < 4;
            REQUIRE(contains(a, p) == (in_square && !in_half));
        }
}

TEST_CASE("an arc hole") {
    // a square with a circular hole, the hole wound the other way
    AR hole = circle_ring<rational>(V(5, 5), 2);
    std::reverse(hole.begin(), hole.end());
    // reversing the vertex order also has to move each bulge to the edge that now starts there
    hole = AR{AV(hole[0].p, rational(1)), AV(hole[1].p, rational(1))};

    const A a({AR{AV(V(0, 0)), AV(V(10, 0)), AV(V(10, 10)), AV(V(0, 10))}, hole});
    REQUIRE(contains(a, V(1, 1)));
    REQUIRE(!contains(a, V(5, 5)));        // in the hole
    REQUIRE(!contains(a, V(5, rational(13, 2))));
    REQUIRE(contains(a, V(5, 8)));
    REQUIRE(contains(a, V(5, 7)));         // on the hole boundary
    REQUIRE(!contains(a, V(11, 5)));

    for (int i = -2; i <= 24; i++)
        for (int j = -2; j <= 24; j++) {
            const V p(rational(i, 2), rational(j, 2));
            const bool in_box = 0 <= p.x && p.x <= 10 && 0 <= p.y && p.y <= 10;
            const rational d = (p.x - 5) * (p.x - 5) + (p.y - 5) * (p.y - 5);
            REQUIRE(contains(a, p) == (in_box && d >= 4));
        }
}

TEST_CASE("inversion of an arc region") {
    const A disk(circle_ring<rational>(V(0, 0), 2));
    const A out = ~disk;
    REQUIRE(out.complement);
    REQUIRE(out.rings == disk.rings);
    REQUIRE(contains(disk, V(0, 1)));
    REQUIRE(!contains(out, V(0, 1)));
    REQUIRE(!contains(disk, V(3, 0)));
    REQUIRE(contains(out, V(3, 0)));
    // the boundary belongs to both closed regions
    REQUIRE(contains(disk, V(2, 0)));
    REQUIRE(contains(out, V(2, 0)));
    REQUIRE(~out == disk);
    REQUIRE_THROWS(signed_area_chords(out));
}

TEST_CASE("a straight edged arc ring matches the straight edged type") {
    const AR ring{AV(V(0, 0)), AV(V(6, 0)), AV(V(6, 6)), AV(V(0, 6))};
    const A a(ring);
    const MultiPolygon2<rational> b(Ring2<rational>{V(0, 0), V(6, 0), V(6, 6), V(0, 6)});
    REQUIRE(signed_area_chords(a) == signed_area(b));
    REQUIRE(chord_polygon(a) == b);
    for (int i = -2; i <= 16; i++)
        for (int j = -2; j <= 16; j++) {
            const V p(rational(i, 2), rational(j, 2));
            REQUIRE(contains(a, p) == contains(b, p));
        }
}

TEST_CASE("contains a point where two chords cross") {
    // A chord is not part of the boundary, so contains() has to step off one before applying the
    // winding formula. Stepping along a chord's own direction keeps the point on it, so a point
    // lying on a horizontal chord and a vertical one at once needs a direction that is neither.

    // a disk of radius 2 whose chord runs vertically through the origin
    const AR disk{AV{V(0, 2), -1}, AV{V(0, -2), -1}};
    // a lens whose two arcs both bulge below the x axis, sharing the endpoints (-1,0) and (1,0):
    // its chords pass through the origin, its body does not
    const AR lens{AV{V(-1, 0), -half}, AV{V(1, 0), half}};
    const A a{std::vector<AR>{disk, lens}};

    const V o(0, 0);
    REQUIRE(!on_boundary(a, o));
    REQUIRE(__on_any_chord(a, o));

    // the origin is well inside the disk and above the lens, like everything around it
    REQUIRE(contains(a, V(rational(1, 8), rational(1, 8))));
    REQUIRE(contains(a, V(rational(-1, 8), rational(1, 8))));
    REQUIRE(contains(a, V(rational(1, 8), rational(-1, 8))));
    REQUIRE(contains(a, V(rational(-1, 8), rational(-1, 8))));
    REQUIRE(contains(a, o));

    // the same for a point on the two chords that is not their crossing point
    REQUIRE(contains(a, V(half, 0)));
    REQUIRE(contains(a, V(0, half)));

    // and a point outside is still outside, chord or no chord
    REQUIRE(!contains(a, V(0, 3)));
    REQUIRE(!contains(a, V(3, 0)));
}
