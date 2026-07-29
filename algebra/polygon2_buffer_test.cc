#include "algebra/polygon2_buffer.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"

using P = MultiPolygon2<rational>;
using R = Ring2<rational>;
using V = Vec2<rational>;

static R box(int x0, int y0, int x1, int y1) {
    return R{V(x0, y0), V(x1, y0), V(x1, y1), V(x0, y1)};
}

TEST_CASE("convex_hull") {
    // a square with interior and edge points thrown in
    const R h = convex_hull(std::vector<V>{V(0, 0), V(4, 0), V(4, 4), V(0, 4), V(2, 2), V(2, 0), V(1, 3)});
    REQUIRE(h.size() == 4);
    REQUIRE(is_ccw(h));
    REQUIRE(signed_area(h) == 16);

    // collinear points collapse to a segment, and duplicates are dropped
    REQUIRE(convex_hull(std::vector<V>{V(0, 0), V(1, 1), V(2, 2)}).size() == 2);
    REQUIRE(convex_hull(std::vector<V>{V(1, 1), V(1, 1)}).size() == 1);
    REQUIRE(convex_hull(std::vector<V>{}).empty());

    // a triangle keeps its three vertices and comes out counter clockwise
    const R t = convex_hull(std::vector<V>{V(0, 3), V(0, 0), V(4, 0)});
    REQUIRE(t.size() == 3);
    REQUIRE(is_ccw(t));
    REQUIRE(signed_area(t) == 6);
}

TEST_CASE("structuring elements") {
    const R s = square_element<rational>(2);
    REQUIRE(signed_area(s) == 16);
    REQUIRE(contains(P(s), V(0, 0)));
    REQUIRE(contains(P(s), V(2, 2)));
    REQUIRE(!contains(P(s), V(rational(21, 10), 0)));

    const R d = diamond_element<rational>(2);
    REQUIRE(signed_area(d) == 8);
    REQUIRE(contains(P(d), V(0, 0)));
    REQUIRE(contains(P(d), V(1, 1)));            // |1|+|1| == 2
    REQUIRE(!contains(P(d), V(rational(11, 10), 1)));

    REQUIRE_THROWS(square_element<rational>(0));
    REQUIRE_THROWS(diamond_element<rational>(-1));

    // an inscribed polygon stays inside the disk and contains the origin
    const R e = polygon_element<rational>(2, 6);
    REQUIRE(e.size() == 12);
    REQUIRE(is_ccw(e));
    REQUIRE(contains(P(e), V(0, 0)));
    for (const V& p : e)
        REQUIRE(p.x * p.x + p.y * p.y == 4);     // every vertex is exactly on the circle
    // inscribed, so strictly inside the disk, and closer with more sides
    REQUIRE(signed_area(e) < 4 * rational(314159, 100000));
    REQUIRE(signed_area(e) > rational(11));
    REQUIRE(signed_area(polygon_element<rational>(2, 16)) > signed_area(e));
    REQUIRE(signed_area(polygon_element<rational>(2, 16)) < 4 * rational(314159, 100000));
    REQUIRE_THROWS(polygon_element<rational>(2, 2));
}

TEST_CASE("reflect") {
    // negating every vertex is a rotation by half a turn, which keeps the orientation
    const R a{V(0, 0), V(3, 0), V(0, 4)};
    REQUIRE(is_ccw(a));
    const R r = reflect(a);
    REQUIRE(is_ccw(r));
    REQUIRE(signed_area(r) == signed_area(a));

    // and it is the same ring, negated, starting from the same vertex
    REQUIRE(r.size() == a.size());
    for (size_t i = 0; i < a.size(); i++) {
        REQUIRE(r[i].x == -a[i].x);
        REQUIRE(r[i].y == -a[i].y);
    }

    // an involution
    const R rr = reflect(r);
    REQUIRE(rr == a);

    // a clockwise ring stays clockwise
    R cw = a;
    reverse(cw);
    REQUIRE(!is_ccw(cw));
    REQUIRE(!is_ccw(reflect(cw)));
}

TEST_CASE("dilate a square by a square") {
    const P a(box(0, 0, 10, 10));
    const P g = dilate(a, square_element<rational>(2));
    // Chebyshev growth turns a box into a bigger box
    REQUIRE(signed_area(g) == 14 * 14);
    for (int i = -8; i <= 36; i++)
        for (int j = -8; j <= 36; j++) {
            const V p(rational(i, 2), rational(j, 2));
            const bool want = -2 <= p.x && p.x <= 12 && -2 <= p.y && p.y <= 12;
            REQUIRE(contains(g, p) == want);
        }
}

TEST_CASE("dilate a square by a diamond") {
    const P a(box(0, 0, 10, 10));
    const P g = dilate(a, diamond_element<rational>(2));
    // an octagon: the square grown by 2 on each side with the corners cut back
    REQUIRE(signed_area(g) == 100 + 4 * 20 + 4 * 2);
    for (int i = -8; i <= 36; i++)
        for (int j = -8; j <= 36; j++) {
            const V p(rational(i, 2), rational(j, 2));
            // the Manhattan distance from p to the box has to be at most 2
            rational dx = 0, dy = 0;
            if (p.x < 0) dx = -p.x; else if (p.x > 10) dx = p.x - 10;
            if (p.y < 0) dy = -p.y; else if (p.y > 10) dy = p.y - 10;
            REQUIRE(contains(g, p) == (dx + dy <= 2));
        }
}

TEST_CASE("erode is the dual of dilate") {
    const P a(box(0, 0, 10, 10));
    const P e = erode(a, square_element<rational>(2));
    REQUIRE(signed_area(e) == 6 * 6);
    for (int i = -4; i <= 28; i++)
        for (int j = -4; j <= 28; j++) {
            const V p(rational(i, 2), rational(j, 2));
            const bool want = 2 <= p.x && p.x <= 8 && 2 <= p.y && p.y <= 8;
            REQUIRE(contains(e, p) == want);
        }

    // eroding by more than half the width leaves nothing
    const P gone = erode(a, square_element<rational>(6));
    REQUIRE(gone.rings.empty());
    REQUIRE(gone.is_empty());
}

TEST_CASE("buffer dispatches on the sign of the size") {
    const P a(box(0, 0, 10, 10));
    REQUIRE(signed_area(buffer(a, rational(2))) == 196);
    REQUIRE(signed_area(buffer(a, rational(-2))) == 36);
    REQUIRE(buffer(a, rational(0)) == a);
    REQUIRE(signed_area(buffer(a, rational(2), diamond_element<rational>)) == 100 + 80 + 8);

    // a fractional size
    REQUIRE(signed_area(buffer(a, rational(1, 2))) == 11 * 11);
}

TEST_CASE("dilating closes a gap and eroding splits a bar") {
    // two boxes 2 apart, dilated by 1 in Chebyshev, become one
    const P two({box(0, 0, 4, 4), box(6, 0, 10, 4)});
    const P joined = dilate(two, square_element<rational>(1));
    REQUIRE(joined.rings.size() == 1);
    REQUIRE(signed_area(joined) == 12 * 6);

    // an H with two bars 6 wide joined by a crossbar only 2 tall. Eroding by 1 takes the crossbar
    // down to nothing, so the result is the two bars, each shrunk from 6x16 to 4x14.
    const P h(R{V(0, 0), V(6, 0), V(6, 7), V(10, 7), V(10, 0), V(16, 0),
                V(16, 16), V(10, 16), V(10, 9), V(6, 9), V(6, 16), V(0, 16)});
    REQUIRE(signed_area(h) == 6 * 16 + 6 * 16 + 4 * 2);
    const P thin = erode(h, square_element<rational>(1));
    REQUIRE(thin.rings.size() == 2);
    REQUIRE(signed_area(thin) == 2 * (4 * 14));

    // the sound direction of the erosion property: anything left in the result has its whole
    // square inside the original. Sampling the square only *disproves* membership, so the check
    // runs this way round and not the other.
    for (int i = -2; i <= 34; i++)
        for (int j = -2; j <= 34; j++) {
            const V p(rational(i, 2), rational(j, 2));
            if (!contains(thin, p))
                continue;
            for (int dx = -1; dx <= 1; dx++)
                for (int dy = -1; dy <= 1; dy++)
                    REQUIRE(contains(h, V(p.x + dx, p.y + dy)));
        }
}

TEST_CASE("a region with a hole") {
    R hole = box(3, 3, 7, 7);
    reverse(hole);
    const P a({box(0, 0, 10, 10), hole});
    REQUIRE(signed_area(a) == 84);

    // dilating grows the outside and shrinks the hole
    const P g = dilate(a, square_element<rational>(1));
    REQUIRE(signed_area(g) == 144 - 4);
    REQUIRE(g.rings.size() == 2);

    // dilating by enough closes the hole completely
    const P closed = dilate(a, square_element<rational>(2));
    REQUIRE(closed.rings.size() == 1);
    REQUIRE(signed_area(closed) == 196);
}

TEST_CASE("dilate and erode do not have to be inverses") {
    // a thin spike is lost by eroding, so growing it back does not restore it
    const P a({R{V(0, 0), V(10, 0), V(10, 6), V(0, 6)}, R{V(4, 6), V(6, 6), V(5, 20)}});
    const P once = erode(a, square_element<rational>(1));
    const P back = dilate(once, square_element<rational>(1));
    // the opening is contained in the original, which is the property that always holds
    for (int i = -2; i <= 24; i++)
        for (int j = -2; j <= 44; j++) {
            const V p(rational(i, 2), rational(j, 2));
            if (contains(back, p))
                REQUIRE(contains(a, p));
        }
}

TEST_CASE("the element has to contain the origin") {
    const P a(box(0, 0, 10, 10));
    REQUIRE_THROWS(dilate(a, box(5, 5, 7, 7)));
    REQUIRE_THROWS(erode(a, box(5, 5, 7, 7)));
    // one that only touches the origin is allowed, since the region is closed
    REQUIRE_NOTHROW(dilate(a, box(0, 0, 2, 2)));
}

TEST_CASE("buffer takes any callable as the element") {
    const P a(box(0, 0, 4, 4));

    // the two elements that are plain functions, as before
    REQUIRE(buffer(a, rational(1)) == dilate(a, square_element(rational(1))));
    REQUIRE(buffer(a, rational(1), diamond_element<rational>) == dilate(a, diamond_element(rational(1))));
    REQUIRE(buffer(a, rational(-1), diamond_element<rational>) == erode(a, diamond_element(rational(1))));

    // and one that needs more than a size, which a function pointer parameter could not express
    auto disk = [](const rational& r) { return polygon_element(r, 4); };
    const P grown = buffer(a, rational(1), disk);
    REQUIRE(grown == dilate(a, polygon_element(rational(1), 4)));
    REQUIRE(signed_area(grown) > signed_area(a));

    // a capturing lambda works as well
    int sides = 6;
    auto captured = [sides](const rational& r) { return polygon_element(r, sides); };
    REQUIRE(buffer(a, rational(1), captured) == dilate(a, polygon_element(rational(1), 6)));

    REQUIRE(buffer(a, rational(0), disk) == a); // zero is still the identity
}
