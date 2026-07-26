#pragma once
#include "algebra/polygon2.h"
#include "algebra/solve_linear.h"

namespace algebra {

// A 2d region whose edges are line segments or circular arcs, with everything exact in T.
//
// An edge is stored as its start vertex plus a bulge, the DXF convention: bulge == tan(theta/4)
// where theta is the arc's included angle, and 0 means a straight edge. The point of that choice is
// that a rational bulge and rational endpoints give a *rational* centre and squared radius, so no
// irrational coordinate is ever needed:
//
//     midpoint of the arc  M = mid + perp(B-A) * bulge/2
//     centre               C = mid - perp(B-A) * (1 - bulge^2)/(4*bulge)
//     squared radius       R2 = dot2(A - C)
//
// A positive bulge puts the arc on the left of A->B, a negative one on the right. A bulge cannot
// describe a full circle, since theta == 2*pi needs tan(pi/2); use two half circle edges instead.
//
// Membership works by chord decomposition. Replacing the chord A->B by the arc A->B changes the
// enclosed region by exactly the circular segment between them, added or removed depending on which
// side the arc bulges. Winding numbers add over closed curves, so
//
//     winding(arc ring, p) == winding(chord polygon, p) - sum over arcs of sign(bulge) [p in S]
//
// where S is the circular segment. Both terms are rational: S is the intersection of the open disk
// with the open half plane on the arc's side of the chord, which also handles a major arc, where S
// is the larger piece and contains the centre. Nothing here needs a square root.

template<typename T>
struct ArcVertex {
    Vec2<T> p;
    T bulge = 0; // for the edge from this vertex to the next one, 0 for a straight edge

    constexpr ArcVertex() {}
    constexpr ArcVertex(Vec2<T> p, T bulge = T(0)) : p(std::move(p)), bulge(std::move(bulge)) {}
    constexpr bool operator==(const ArcVertex& o) const { return p == o.p && bulge == o.bulge; }
};

template<typename T>
using ArcRing2 = std::vector<ArcVertex<T>>;

template<typename T = rational>
struct ArcPolygon2 {
    std::vector<ArcRing2<T>> rings;
    bool complement = false;

    constexpr ArcPolygon2() {}
    constexpr ArcPolygon2(ArcRing2<T> ring) { rings.push_back(std::move(ring)); }
    constexpr ArcPolygon2(std::vector<ArcRing2<T>> r, bool c = false) : rings(std::move(r)), complement(c) {}

    constexpr bool is_empty() const { return rings.empty() && !complement; }
    constexpr bool is_whole_plane() const { return rings.empty() && complement; }
};

template<typename T>
constexpr Vec2<T> __perp(const Vec2<T>& v) { return Vec2<T>{-v.y, v.x}; }

// The point halfway along the arc. Rational, unlike the radius.
template<typename T>
constexpr Vec2<T> arc_midpoint(const Vec2<T>& a, const Vec2<T>& b, const T& bulge) {
    return (a + b) / T(2) + __perp(b - a) * (bulge / T(2));
}

template<typename T>
constexpr Vec2<T> arc_center(const Vec2<T>& a, const Vec2<T>& b, const T& bulge) {
    Check(bulge != 0, "arc_center() of a straight edge");
    return (a + b) / T(2) - __perp(b - a) * ((T(1) - bulge * bulge) / (T(4) * bulge));
}

template<typename T>
constexpr T arc_radius2(const Vec2<T>& a, const Vec2<T>& b, const T& bulge) {
    return dot2(a - arc_center(a, b, bulge));
}

// True when p is strictly inside the circular segment cut off by the chord a->b: strictly inside
// the circle, and strictly on the side the arc bulges towards.
template<typename T>
constexpr bool __in_circular_segment(const Vec2<T>& a, const Vec2<T>& b, const T& bulge, const Vec2<T>& p) {
    if (bulge == 0)
        return false;
    const Vec2<T> c = arc_center(a, b, bulge);
    if (!(dot2(p - c) < dot2(a - c)))
        return false;
    // the arc's side of the chord is the side its midpoint is on, which is sign(bulge)
    const T side = cross(b - a, p - a);
    return (bulge > 0) ? (side > 0) : (side < 0);
}

// The chord polygon: the same vertices with every arc replaced by its chord.
template<typename T>
constexpr Ring2<T> chord_ring(const ArcRing2<T>& ring) {
    Ring2<T> out;
    out.reserve(ring.size());
    for (const ArcVertex<T>& v : ring)
        out.push_back(v.p);
    return out;
}

template<typename T>
constexpr MultiPolygon2<T> chord_polygon(const ArcPolygon2<T>& a) {
    MultiPolygon2<T> out;
    out.complement = a.complement;
    for (const ArcRing2<T>& ring : a.rings)
        out.rings.push_back(chord_ring(ring));
    return out;
}

// Is w within the counter clockwise sweep from u to v, all measured from the centre?
template<typename T>
constexpr bool __in_ccw_sweep(const Vec2<T>& u, const Vec2<T>& v, const Vec2<T>& w) {
    const T uv = cross(u, v);
    if (uv > 0) // the sweep is less than half a turn
        return cross(u, w) >= 0 && cross(w, v) >= 0;
    if (uv < 0) // more than half a turn, so take the complement of the short way back
        return !(cross(v, w) > 0 && cross(w, u) > 0);
    // exactly half a turn, or none at all
    if (dot(u, v) < 0)
        return cross(u, w) >= 0;
    return cross(u, w) == 0 && dot(u, w) > 0;
}

// True when p lies on the arc (or the straight edge) from a to b, endpoints included.
template<typename T>
constexpr bool on_arc(const Vec2<T>& a, const Vec2<T>& b, const T& bulge, const Vec2<T>& p) {
    if (bulge == 0)
        return cross(b - a, p - a) == 0 && loose_order(a, p, b);
    const Vec2<T> c = arc_center(a, b, bulge);
    if (dot2(p - c) != dot2(a - c))
        return false;
    if (p == a || p == b)
        return true;
    // a positive bulge leaves the arc on the left of a->b, which is clockwise about the centre,
    // so the counter clockwise sweep runs from b to a in that case
    return (bulge > 0) ? __in_ccw_sweep(b - c, a - c, p - c) : __in_ccw_sweep(a - c, b - c, p - c);
}

template<typename T>
constexpr bool on_boundary(const ArcRing2<T>& ring, const Vec2<T>& p) {
    if (ring.size() == 1)
        return ring[0].p == p;
    for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++)
        if (on_arc(ring[j].p, ring[i].p, ring[j].bulge, p))
            return true;
    return false;
}

template<typename T>
constexpr bool on_boundary(const ArcPolygon2<T>& a, const Vec2<T>& p) {
    for (const ArcRing2<T>& ring : a.rings)
        if (on_boundary(ring, p))
            return true;
    return false;
}

// True when p is on a chord of an arc, where the winding formula below does not apply. Such a point
// is not on the region's boundary; contains() steps off it instead.
template<typename T>
constexpr bool __on_any_chord(const ArcPolygon2<T>& a, const Vec2<T>& p) {
    for (const ArcRing2<T>& ring : a.rings) {
        if (ring.size() < 2)
            continue;
        for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
            if (ring[j].bulge == 0)
                continue;
            if (cross(ring[i].p - ring[j].p, p - ring[j].p) == 0 && loose_order(ring[j].p, p, ring[i].p))
                return true;
        }
    }
    return false;
}

// Undefined when p is on the boundary or on a chord of an arc, in the same way that the straight
// edged winding_number is undefined on the boundary.
template<typename T>
constexpr int winding_number(const ArcPolygon2<T>& a, const Vec2<T>& p) {
    int wn = 0;
    for (const ArcRing2<T>& ring : a.rings) {
        wn += winding_number(chord_ring(ring), p);
        if (ring.size() < 2)
            continue;
        for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
            const T& bulge = ring[j].bulge;
            if (bulge == 0)
                continue;
            if (__in_circular_segment(ring[j].p, ring[i].p, bulge, p))
                wn -= (bulge > 0) ? 1 : -1;
        }
    }
    return wn;
}

// Does p + dir*t meet the segment u..v for some t in (0, tmax]?
template<typename T>
constexpr bool __ray_meets_segment(const Vec2<T>& p, const Vec2<T>& dir, const T& tmax,
                                   const Vec2<T>& u, const Vec2<T>& v) {
    T t, s, det;
    if (!__solve_linear<T>(p - u, dir, u - v, t, s, det))
        return false; // parallel: it either misses, or overlaps and an endpoint is caught elsewhere
    if (det < 0) {
        negate(t);
        negate(s);
        negate(det);
    }
    return t > 0 && t <= tmax * det && s >= 0 && s <= det;
}

// Does p + dir*t meet the circle for some t in (0, tmax]? Solved by sign analysis of the quadratic
// rather than by taking a square root, so it stays in T. Conservative on purpose: it answers for the
// whole circle, not just the arc, which is all that the step below needs.
template<typename T>
constexpr bool __ray_meets_circle(const Vec2<T>& p, const Vec2<T>& dir, const T& tmax,
                                  const Vec2<T>& c, const T& r2) {
    // f(t) = A t^2 + 2 B t + C, with A > 0
    const T A = dot2(dir);
    if (A == 0)
        return false;
    const T B = dot(dir, p - c);
    const T C = dot2(p - c) - r2;
    const T f1 = A * tmax * tmax + 2 * B * tmax + C;
    if (C == 0) {
        // p is on the circle, so t == 0 is a root; the other one is at -2B/A
        const T other = -2 * B;
        return other > 0 && other <= tmax * A;
    }
    if ((C > 0) != (f1 > 0) || f1 == 0)
        return true; // a sign change, so a root in (0, tmax]
    // both ends on the same side, so a root needs the extremum inside the interval and past zero
    const T tv = -B; // the extremum is at -B/A
    if (!(tv > 0 && tv < tmax * A))
        return false;
    return B * B >= A * C; // discriminant of f, which touches or crosses zero
}

// A step along dir that provably crosses nothing: no chord, and no circle of any arc. Found by
// halving, which terminates because p is not on the boundary and any point off a curve has a
// positive distance to it.
template<typename T>
constexpr T __arc_safe_step(const ArcPolygon2<T>& a, const Vec2<T>& p, const Vec2<T>& dir) {
    T t = 1;
    for (int i = 0; i < 200; i++) {
        bool clear = true;
        for (const ArcRing2<T>& ring : a.rings) {
            if (ring.size() < 2)
                continue;
            for (size_t k = 0, j = ring.size() - 1; k < ring.size() && clear; j = k++) {
                const Vec2<T>& u = ring[j].p;
                const Vec2<T>& v = ring[k].p;
                if (__ray_meets_segment(p, dir, t, u, v))
                    clear = false;
                else if (ring[j].bulge != 0 &&
                         __ray_meets_circle(p, dir, t, arc_center(u, v, ring[j].bulge),
                                            arc_radius2(u, v, ring[j].bulge)))
                    clear = false;
            }
            if (!clear)
                break;
        }
        if (clear)
            return t;
        t /= 2;
    }
    return t;
}

// Closed region membership. The boundary belongs to the region on either side of a complement, the
// same as for the straight edged type.
template<typename T>
constexpr bool contains(const ArcPolygon2<T>& a, const Vec2<T>& p) {
    if (on_boundary(a, p))
        return true;
    Vec2<T> q = p;
    if (__on_any_chord(a, p)) {
        // A chord is not part of the boundary, but the winding formula is undefined on one, so step
        // off it. The step provably crosses nothing, so it stays in the same face. Trying two
        // directions is enough, since a point cannot lie on chords running both ways at once.
        for (const Vec2<T>& dir : {Vec2<T>{T(1), T(0)}, Vec2<T>{T(0), T(1)}}) {
            q = p + dir * __arc_safe_step(a, p, dir);
            if (!__on_any_chord(a, q) && !on_boundary(a, q))
                break;
        }
    }
    return (winding_number(a, q) != 0) != a.complement;
}

template<typename T>
constexpr ArcPolygon2<T> operator~(ArcPolygon2<T> a) {
    a.complement = !a.complement;
    return a;
}

template<typename T>
constexpr bool operator==(const ArcPolygon2<T>& a, const ArcPolygon2<T>& b) {
    return a.complement == b.complement && a.rings == b.rings;
}

// The area of the chord polygon. The true area differs by the circular segments, whose area is
// r^2*(theta - sin theta)/2 and so is transcendental -- not representable in rational, which is why
// there is no exact signed_area() for an arc region here.
template<typename T>
constexpr T signed_area_chords(const ArcPolygon2<T>& a) {
    Check(!a.complement, "signed_area_chords() of an unbounded region");
    T sum = 0;
    for (const ArcRing2<T>& ring : a.rings)
        sum += signed_area2(chord_ring(ring));
    return sum / 2;
}

// A closed circle, as the two half circle edges that a single bulge cannot express.
template<typename T>
constexpr ArcRing2<T> circle_ring(const Vec2<T>& center, const T& radius) {
    const Vec2<T> l{center.x - radius, center.y};
    const Vec2<T> r{center.x + radius, center.y};
    // counter clockwise: right vertex bulging below, left vertex bulging above
    return ArcRing2<T>{ArcVertex<T>{r, T(-1)}, ArcVertex<T>{l, T(-1)}};
}

}
