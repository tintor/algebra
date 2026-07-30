#pragma once
#include "algebra/rational.h"
#include "algebra/vector.h"
#include <vector>

namespace algebra {

// A closed ring of vertices. The closing edge from back() to front() is implicit, so a ring is
// never expected to repeat its first vertex at the end. Fewer than 3 vertices encloses no area.
template<typename T>
using Ring2 = std::vector<Vec2<T>>;

// A region of the plane, as a set of rings plus a complement flag.
//
// Membership uses the nonzero winding rule: a point is inside when its winding number around the
// rings is nonzero, which is then flipped by `complement`. An outer boundary winds counter
// clockwise (positive signed area) and a hole winds clockwise, so a hole cancels the shell that
// contains it. Nesting to any depth works without tracking which ring is whose hole.
//
// `complement` is what makes inversion exact and free: the complement of a bounded region is
// unbounded and has no finite ring representation, so it is carried as a flag instead. The empty
// region is no rings, and the whole plane is no rings with complement set.
//
// T needs exact arithmetic for the predicates to be reliable; rational is the default.
// True when halving is exact in T, which the boolean and buffer operations need: they take the
// midpoint of a fragment and step off it by a fraction of the normal. With a truncating T both land
// back on the boundary, where the winding rule says nothing.
template<typename T>
constexpr bool __exact_halving() { return T(1) / T(2) * T(2) == T(1); }

template<typename T = rational>
struct MultiPolygon2 {
    std::vector<Ring2<T>> rings;
    bool complement = false;

    constexpr MultiPolygon2() {}
    constexpr MultiPolygon2(Ring2<T> ring) { rings.push_back(std::move(ring)); }
    constexpr MultiPolygon2(std::vector<Ring2<T>> r, bool c = false) : rings(std::move(r)), complement(c) {}

    constexpr bool is_empty() const { return rings.empty() && !complement; }
    constexpr bool is_whole_plane() const { return rings.empty() && complement; }
};

// Twice the signed area, which stays integral when T is. Positive is counter clockwise.
template<typename T>
constexpr T signed_area2(const Ring2<T>& ring) {
    T sum = 0;
    if (ring.size() < 3)
        return sum;
    // the shoelace formula, taken around the closing edge as well
    for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++)
        sum += cross(ring[j], ring[i]);
    return sum;
}

template<typename T>
constexpr T signed_area(const Ring2<T>& ring) { return signed_area2(ring) / 2; }

template<typename T>
constexpr bool is_ccw(const Ring2<T>& ring) { return signed_area2(ring) > 0; }

template<typename T>
constexpr void reverse(Ring2<T>& ring) { std::reverse(ring.begin(), ring.end()); }

// The signed area of the whole region. Throws for an unbounded one, which has no finite area.
template<typename T>
constexpr T signed_area(const MultiPolygon2<T>& a) {
    Check(!a.complement, "signed_area() of an unbounded region");
    T sum = 0;
    for (const Ring2<T>& ring : a.rings)
        sum += signed_area2(ring);
    return sum / 2;
}

// True when p lies exactly on the closing edge or on any edge of the ring.
template<typename T>
constexpr bool on_boundary(const Ring2<T>& ring, const Vec2<T>& p) {
    if (ring.size() == 1)
        return ring[0] == p;
    for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
        const Vec2<T>& a = ring[j];
        const Vec2<T>& b = ring[i];
        // collinear with the edge, and inside its bounding box
        if (cross(b - a, p - a) == 0 && loose_order(a, p, b))
            return true;
    }
    return false;
}

template<typename T>
constexpr bool on_boundary(const MultiPolygon2<T>& a, const Vec2<T>& p) {
    for (const Ring2<T>& ring : a.rings)
        if (on_boundary(ring, p))
            return true;
    return false;
}

// How many times the ring wraps counter clockwise around p. Undefined when p is on the boundary,
// so callers test on_boundary() first. Counts upward crossings of the ray going right from p.
template<typename T>
constexpr int winding_number(const Ring2<T>& ring, const Vec2<T>& p) {
    int wn = 0;
    if (ring.size() < 3)
        return wn;
    for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
        const Vec2<T>& a = ring[j];
        const Vec2<T>& b = ring[i];
        // a half open rule on y, so a vertex exactly at p.y is counted by one edge and not two
        if (a.y <= p.y) {
            if (b.y > p.y && cross(b - a, p - a) > 0)
                wn += 1;
        } else {
            if (b.y <= p.y && cross(b - a, p - a) < 0)
                wn -= 1;
        }
    }
    return wn;
}

template<typename T>
constexpr int winding_number(const MultiPolygon2<T>& a, const Vec2<T>& p) {
    int wn = 0;
    for (const Ring2<T>& ring : a.rings)
        wn += winding_number(ring, p);
    return wn;
}

// Closed region membership: the boundary belongs to the region, on either side of a complement.
template<typename T>
constexpr bool contains(const MultiPolygon2<T>& a, const Vec2<T>& p) {
    if (on_boundary(a, p))
        return true;
    return (winding_number(a, p) != 0) != a.complement;
}

// Inversion. Exact and free, since it only flips the flag; the boundary stays shared, which is
// what makes `contains` true for it on both sides.
template<typename T>
constexpr MultiPolygon2<T> operator~(MultiPolygon2<T> a) {
    a.complement = !a.complement;
    return a;
}

// Equality of the representation, not of the region: the same set of points has many
// representations, and rotating a ring's vertex list or reordering the rings gives an unequal one.
// Two regions are the same region when their symmetric difference is empty, which costs a boolean
// operation -- (a ^ b).is_empty() spells that out at the call site.
template<typename T>
constexpr bool operator==(const MultiPolygon2<T>& a, const MultiPolygon2<T>& b) {
    return a.complement == b.complement && a.rings == b.rings;
}

// The smallest axis aligned box containing every vertex. Throws for an unbounded region.
template<typename T>
constexpr void bounding_box(const MultiPolygon2<T>& a, Vec2<T>& min, Vec2<T>& max) {
    Check(!a.complement, "bounding_box() of an unbounded region");
    bool first = true;
    for (const Ring2<T>& ring : a.rings)
        for (const Vec2<T>& p : ring) {
            if (first) {
                min = max = p;
                first = false;
                continue;
            }
            if (p.x < min.x) min.x = p.x;
            if (p.y < min.y) min.y = p.y;
            if (max.x < p.x) max.x = p.x;
            if (max.y < p.y) max.y = p.y;
        }
    Check(!first, "bounding_box() of an empty region");
}

// True when every vertex lies on one line, so the ring encloses no area whatever its shape.
template<typename T>
constexpr bool __all_collinear(const Ring2<T>& ring) {
    if (ring.size() < 3)
        return true;
    size_t i = 1;
    while (i < ring.size() && ring[i] == ring[0])
        i += 1;
    if (i == ring.size())
        return true; // every vertex is the same point
    const Vec2<T> dir = ring[i] - ring[0];
    for (const Vec2<T>& p : ring)
        if (cross(p - ring[0], dir) != 0)
            return false;
    return true;
}

// Drops repeated consecutive vertices, including across the closing edge, and then vertices that
// lie exactly on the segment between their neighbours. Rings that enclose no area at all are
// removed: fewer than 3 vertices, or every vertex on one line.
//
// The interior is unchanged. Note that removing a degenerate ring does drop its boundary, so a
// point that was only on such a sliver stops being contained -- that is the point of the call, and
// it is why zero *signed* area is not the test used here: a ring shaped like a figure eight with
// equal lobes also has zero signed area, yet each lobe has a real interior.
template<typename T>
constexpr void simplify(MultiPolygon2<T>& a) {
    std::vector<Ring2<T>> out;
    for (Ring2<T>& ring : a.rings) {
        Ring2<T> r;
        for (const Vec2<T>& p : ring)
            if (r.empty() || !(r.back() == p))
                r.push_back(p);
        while (r.size() > 1 && r.front() == r.back())
            r.pop_back();

        // Remove the vertices that lie on the segment between their neighbours. Removing one can
        // expose another, so the scan carries a stack: after dropping the middle of the last three,
        // the new last three are checked again before moving on. That is one pass instead of a
        // restart from the beginning per removal.
        auto on_segment = [](const Vec2<T>& x, const Vec2<T>& y, const Vec2<T>& z) {
            return cross(z - x, y - x) == 0 && loose_order(x, y, z);
        };
        Ring2<T> s;
        for (const Vec2<T>& p : r) {
            s.push_back(p);
            while (s.size() >= 3 && on_segment(s[s.size() - 3], s[s.size() - 2], s.back()))
                s.erase(s.end() - 2);
        }
        r = std::move(s);

        // The closing edge, which the pass above cannot see: the last vertex has the first as its
        // neighbour and the other way round, and each removal can expose the other end again.
        for (bool again = true; again && r.size() >= 3;) {
            again = false;
            while (r.size() >= 3 && on_segment(r[r.size() - 2], r.back(), r.front())) {
                r.pop_back();
                again = true;
            }
            while (r.size() >= 3 && on_segment(r.back(), r.front(), r[1])) {
                r.erase(r.begin());
                again = true;
            }
        }
        if (r.size() >= 3 && !__all_collinear(r))
            out.push_back(std::move(r));
    }
    a.rings = std::move(out);
}

}
