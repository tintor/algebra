#include "algebra/geometry.h"
#include "algebra/rational_vector.h"
#include "algebra/__test.h"

using V = qvec3;
using P = Plane3<rational>;
using L = Line3<rational>;

static bool on_plane(const P& p, const V& x) { return dot(p.n, x) + p.d == 0; }

TEST_CASE("plane_plane_intersection") {
    const P x5{{1, 0, 0}, -5, 1}; // x == 5
    const P y7{{0, 1, 0}, -7, 1}; // y == 7

    auto res = plane_intersection(x5, y7);
    REQUIRE(std::holds_alternative<L>(res));
    const L l = std::get<L>(res);
    REQUIRE(!is_zero(l.dir));
    for (rational t : {rational(0), rational(1), rational(-3, 2), rational(11, 7)}) {
        const V p = l.orig + l.dir * t;
        REQUIRE(on_plane(x5, p));
        REQUIRE(on_plane(y7, p));
    }

    // slanted planes
    const P a{{1, 1, 0}, -2, 2};
    const P b{{0, 1, 1}, -3, 2};
    res = plane_intersection(a, b);
    REQUIRE(std::holds_alternative<L>(res));
    const L m = std::get<L>(res);
    for (rational t : {rational(0), rational(1), rational(-5, 3)}) {
        const V p = m.orig + m.dir * t;
        REQUIRE(on_plane(a, p));
        REQUIRE(on_plane(b, p));
    }

    // same plane
    res = plane_intersection(x5, x5);
    REQUIRE(std::holds_alternative<P>(res));
    // same plane, opposite orientation
    res = plane_intersection(x5, -x5);
    REQUIRE(std::holds_alternative<P>(res));
    // parallel, but distinct
    const P x6{{1, 0, 0}, -6, 1};
    res = plane_intersection(x5, x6);
    REQUIRE(std::holds_alternative<None>(res));
}

TEST_CASE("line_plane_intersection") {
    const P z5{{0, 0, 1}, -5, 1}; // z == 5

    // crossing
    auto res = line_plane_intersection(L{{0, 0, 0}, {0, 0, 1}}, z5);
    REQUIRE(std::holds_alternative<V>(res));
    REQUIRE(std::get<V>(res) == V{0, 0, 5});

    res = line_plane_intersection(L{{1, 2, 3}, {1, 1, 1}}, z5);
    REQUIRE(std::holds_alternative<V>(res));
    REQUIRE(std::get<V>(res) == V{3, 4, 5});

    // parallel, outside the plane
    res = line_plane_intersection(L{{0, 0, 0}, {1, 0, 0}}, z5);
    REQUIRE(std::holds_alternative<None>(res));

    // inside the plane
    res = line_plane_intersection(L{{0, 0, 5}, {1, 0, 0}}, z5);
    REQUIRE(std::holds_alternative<L>(res));
}

TEST_CASE("plane_plane_plane_intersection") {
    const P x1{{1, 0, 0}, -1, 1};
    const P y2{{0, 1, 0}, -2, 1};
    const P z3{{0, 0, 1}, -3, 1};

    // single point
    auto res = plane_intersection(x1, y2, z3);
    REQUIRE(std::holds_alternative<V>(res));
    REQUIRE(std::get<V>(res) == V{1, 2, 3});

    // two identical planes and a crossing one -> line
    res = plane_intersection(x1, x1, z3);
    REQUIRE(std::holds_alternative<L>(res));
    const L l = std::get<L>(res);
    for (rational t : {rational(0), rational(2)}) {
        const V p = l.orig + l.dir * t;
        REQUIRE(on_plane(x1, p));
        REQUIRE(on_plane(z3, p));
    }

    // three identical planes -> plane
    res = plane_intersection(x1, x1, x1);
    REQUIRE(std::holds_alternative<P>(res));

    // parallel and distinct -> empty
    const P x2{{1, 0, 0}, -2, 1};
    res = plane_intersection(x1, x2, z3);
    REQUIRE(std::holds_alternative<None>(res));

    // triangular prism: normals are linearly dependent, no two planes parallel
    const P p1{{1, 0, 0}, 0, 1};
    const P p2{{0, 1, 0}, 0, 1};
    const P p3{{1, 1, 0}, -1, 2};
    res = plane_intersection(p1, p2, p3);
    REQUIRE(std::holds_alternative<None>(res));

    // three planes sharing a line
    const P q3{{1, 1, 0}, 0, 2};
    res = plane_intersection(p1, p2, q3);
    REQUIRE(std::holds_alternative<L>(res));
}
