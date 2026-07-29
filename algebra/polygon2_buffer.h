#pragma once
#include "algebra/polygon2_boolean.h"

namespace algebra {

// Buffering (offsetting) a 2d region, exact for an exact T.
//
// The operation is a Minkowski sum with a convex structuring element that contains the origin:
// dilating grows the region by that element, eroding shrinks it by the same amount. Note what is
// *not* offered: buffering by a Euclidean distance r. Moving an edge out by r needs the unit normal,
// i.e. sqrt(dx*dx + dy*dy), which is irrational for rational input, and a round join needs a
// circular arc whose intersections with its neighbours are irrational too. Neither is representable
// in rational, so the element is given explicitly instead and the result stays exact.
//
// square_element(r) buffers in the Chebyshev metric (max(|dx|, |dy|) <= r) and diamond_element(r) in
// the Manhattan one (|dx| + |dy| <= r). A finer convex polygon approximates a disk as closely as
// wanted while staying exact, at the cost of more edges.
//
// Dilation uses A + B == A union (boundary of A) + B, which holds because B is convex and contains
// the origin: for q = p + b with p in A, the whole segment from p to q lies in p + B, so if q is
// outside A that segment leaves A at some p', and q - p' is still a multiple of b in B. The sum of a
// segment with a convex B is the convex hull of the two translated copies, so each edge contributes
// one convex polygon and the result is their union with A.
//
// Erosion is the dual, A - B == complement of (complement of A) + (-B), which the complement flag
// makes exact and cheap.

// The counter clockwise convex hull of a set of points, by monotone chain. Exact for an exact T.
template<typename T>
constexpr Ring2<T> convex_hull(std::vector<Vec2<T>> p) {
    std::sort(p.begin(), p.end(), [](const Vec2<T>& a, const Vec2<T>& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    p.erase(std::unique(p.begin(), p.end(), [](const Vec2<T>& a, const Vec2<T>& b) { return a == b; }), p.end());
    if (p.size() < 3)
        return Ring2<T>(p.begin(), p.end());

    Ring2<T> h;
    h.reserve(p.size() + 1);
    // the lower hull, then the upper one; a turn that is not strictly left pops the middle point
    for (int pass = 0; pass < 2; pass++) {
        const size_t base = h.size();
        for (size_t i = 0; i < p.size(); i++) {
            const Vec2<T>& q = pass == 0 ? p[i] : p[p.size() - 1 - i];
            while (h.size() >= base + 2 && cross(h[h.size() - 1] - h[h.size() - 2], q - h[h.size() - 2]) <= 0)
                h.pop_back();
            h.push_back(q);
        }
        h.pop_back(); // the last point of each pass starts the next one
    }
    return h;
}

// max(|dx|, |dy|) <= r, i.e. buffering in the Chebyshev metric.
template<typename T>
constexpr Ring2<T> square_element(const T& r) {
    Check(r > 0, "square_element() needs a positive size");
    return Ring2<T>{Vec2<T>{-r, -r}, Vec2<T>{r, -r}, Vec2<T>{r, r}, Vec2<T>{-r, r}};
}

// |dx| + |dy| <= r, i.e. buffering in the Manhattan metric.
template<typename T>
constexpr Ring2<T> diamond_element(const T& r) {
    Check(r > 0, "diamond_element() needs a positive size");
    return Ring2<T>{Vec2<T>{-r, 0}, Vec2<T>{0, -r}, Vec2<T>{r, 0}, Vec2<T>{0, r}};
}

// A convex polygon with 2*sides vertices inscribed in the circle of radius r, with rational
// coordinates throughout. It is a subset of the disk, so it under-approximates a round buffer, and
// more sides brings it closer while staying exact.
template<typename T>
constexpr Ring2<T> polygon_element(const T& r, int sides) {
    Check(r > 0, "polygon_element() needs a positive size");
    Check(sides >= 3, "polygon_element() needs at least 3 sides");
    // rational points on the circle come from the Pythagorean parametrisation
    // (r*(1-t*t)/(1+t*t), r*2*t/(1+t*t)) for rational t, which sweeps the whole circle
    std::vector<Vec2<T>> pts;
    for (int i = 0; i < sides; i++) {
        const T t = T(2 * i - sides) / T(sides); // t over (-1, 1] covers half the circle
        const T d = T(1) + t * t;
        const Vec2<T> q{r * (T(1) - t * t) / d, r * (T(2) * t) / d};
        pts.push_back(q);
        pts.push_back(Vec2<T>{-q.x, -q.y}); // the opposite point, to cover the other half
    }
    return convex_hull(pts);
}

// -B, which is what eroding by B needs. Negating every vertex is a rotation by half a turn, so the
// orientation carries over on its own; reversing the ring as well used to flip it.
template<typename T>
constexpr Ring2<T> reflect(Ring2<T> a) {
    for (Vec2<T>& p : a)
        p = Vec2<T>{-p.x, -p.y};
    return a;
}

// Grows `a` by the convex element `b`, which has to contain the origin.
template<typename T>
constexpr MultiPolygon2<T> dilate(const MultiPolygon2<T>& a, const Ring2<T>& b) {
    Check(contains(MultiPolygon2<T>(b), Vec2<T>{T(0), T(0)}), "the structuring element must contain the origin");

    // One boolean operation, not one per edge. Every hull below comes out of convex_hull() counter
    // clockwise, so their winding numbers only add: a point is inside at least one of them exactly
    // when the winding number of the whole set is non-zero, which is what the nonzero rule already
    // says. So the hulls can be handed over as one region and unioned with `a` in a single pass.
    std::vector<Ring2<T>> hulls;
    for (const Ring2<T>& ring : a.rings) {
        if (ring.size() < 2)
            continue;
        for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
            // the sum of the edge with b is the hull of b translated to each end
            std::vector<Vec2<T>> pts;
            pts.reserve(b.size() * 2);
            for (const Vec2<T>& q : b) {
                pts.push_back(ring[j] + q);
                pts.push_back(ring[i] + q);
            }
            hulls.push_back(convex_hull(pts));
        }
    }
    if (hulls.empty())
        return a;
    // Union the hulls in a balanced tree rather than one at a time into an accumulator: each
    // boolean_op is quadratic in the edges it is given, so pairing small operands first keeps the
    // large operands to the last few rounds.
    std::vector<MultiPolygon2<T>> level;
    level.reserve(hulls.size());
    for (Ring2<T>& h : hulls)
        level.push_back(MultiPolygon2<T>(std::move(h)));
    while (level.size() > 1) {
        std::vector<MultiPolygon2<T>> next;
        next.reserve((level.size() + 1) / 2);
        for (size_t i = 0; i + 1 < level.size(); i += 2)
            next.push_back(level[i] | level[i + 1]);
        if (level.size() % 2)
            next.push_back(std::move(level.back()));
        level = std::move(next);
    }
    return a | level.front();
}

// Shrinks `a` by the convex element `b`, the dual of dilate through the complement.
template<typename T>
constexpr MultiPolygon2<T> erode(const MultiPolygon2<T>& a, const Ring2<T>& b) {
    MultiPolygon2<T> out = ~dilate(~a, reflect(b));
    // dilate() leaves each ring with the interior of the *complemented* region on its left, and
    // flipping the flag swaps which side that is, so the rings have to be reversed to restore the
    // usual orientation. Membership does not depend on this, since the winding rule only tests
    // against zero, but the sign of signed_area() does.
    for (Ring2<T>& ring : out.rings)
        reverse(ring);
    return out;
}

// Positive size grows, negative shrinks, zero is the identity. `element` maps a size to a convex
// structuring element; the default buffers in the Chebyshev metric.
template<typename T>
constexpr MultiPolygon2<T> buffer(const MultiPolygon2<T>& a, const T& size,
                                  Ring2<T> (*element)(const T&) = square_element<T>) {
    if (size == 0)
        return a;
    if (size > 0)
        return dilate(a, element(size));
    return erode(a, element(-size));
}

}
