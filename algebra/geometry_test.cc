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

TEST_CASE("are_parallel") {
    REQUIRE(are_parallel(V{1, 2, 3}, V{2, 4, 6}));
    REQUIRE(are_parallel(V{1, 2, 3}, V{-1, -2, -3})); // anti-parallel counts
    REQUIRE(are_parallel(V{1, 2, 3}, V{1/2_q, 1, 3/2_q}));
    REQUIRE(!are_parallel(V{1, 2, 3}, V{2, 4, 7}));
    REQUIRE(!are_parallel(V{1, 0, 0}, V{0, 1, 0}));
    // a zero vector is parallel to everything
    REQUIRE(are_parallel(V{0, 0, 0}, V{1, 2, 3}));
    REQUIRE(are_parallel(V{1, 2, 3}, V{0, 0, 0}));
    REQUIRE(are_parallel(V{0, 0, 0}, V{0, 0, 0}));
    // one zero component must not short circuit the other two
    REQUIRE(are_parallel(V{0, 2, 4}, V{0, 1, 2}));
    REQUIRE(!are_parallel(V{0, 2, 4}, V{0, 1, 3}));
}

TEST_CASE("Plane3 operator-") {
    const P x5{{1, 0, 0}, -5, 1};
    const P n = -x5;
    REQUIRE(n.n == V{-1, 0, 0});
    REQUIRE(n.d == 5);
    REQUIRE(n.den == 1); // den is not negated
    REQUIRE(-n == x5);
    // same set of points, opposite orientation
    REQUIRE(on_plane(x5, V{5, 1, 2}));
    REQUIRE(on_plane(n, V{5, 1, 2}));
}

TEST_CASE("Plane3 operator==") {
    const P x5{{1, 0, 0}, -5, 1};
    REQUIRE(x5 == x5);
    REQUIRE(!(x5 == -x5)); // equality includes orientation
    REQUIRE(!(x5 == P{{1, 0, 0}, -6, 1}));
    REQUIRE(!(x5 == P{{0, 1, 0}, -5, 1}));

    // den != 1 takes the scaling aware path: (2x - 10)/sqrt(4) is the same equation as x - 5
    REQUIRE(P{{2, 0, 0}, -10, 4} == x5);
    REQUIRE(x5 == P{{2, 0, 0}, -10, 4});
    REQUIRE(!(P{{2, 0, 0}, 10, 4} == x5)); // x == -5
    REQUIRE(!(P{{2, 0, 0}, -10, 4} == -x5));
    REQUIRE(P{{-2, 0, 0}, 10, 4} == -x5);

    // through the origin, where d == 0 on both sides
    const P x0{{1, 0, 0}, 0, 1};
    REQUIRE(P{{3, 0, 0}, 0, 9} == x0);
    REQUIRE(!(P{{-3, 0, 0}, 0, 9} == x0));
    REQUIRE(P{{-3, 0, 0}, 0, 9} == -x0);
}

TEST_CASE("plane_plane_plane_intersection - remaining parallel branches") {
    const P x1{{1, 0, 0}, -1, 1};
    const P y2{{0, 1, 0}, -2, 1};
    const P z3{{0, 0, 1}, -3, 1};

    // b parallel to c, and equal: result is a into b
    auto res = plane_intersection(x1, y2, y2);
    REQUIRE(std::holds_alternative<L>(res));
    for (rational t : {rational(0), rational(3)}) {
        const V p = std::get<L>(res).orig + std::get<L>(res).dir * t;
        REQUIRE(on_plane(x1, p));
        REQUIRE(on_plane(y2, p));
    }

    // b parallel to c, and distinct: empty
    REQUIRE(std::holds_alternative<None>(plane_intersection(x1, y2, P{{0, 1, 0}, -5, 1})));

    // a parallel to c, and equal
    res = plane_intersection(x1, y2, x1);
    REQUIRE(std::holds_alternative<L>(res));

    // a parallel to c, and distinct: empty
    REQUIRE(std::holds_alternative<None>(plane_intersection(x1, y2, P{{1, 0, 0}, -7, 1})));

    // an opposingly oriented duplicate is still the same set of points
    res = plane_intersection(x1, -x1, z3);
    REQUIRE(std::holds_alternative<L>(res));
    REQUIRE(std::holds_alternative<P>(plane_intersection(x1, -x1, x1)));
}
