#pragma once
#include "algebra/polygon2.h"
#include "algebra/segment_segment_intersection.h"
#include <algorithm>
#include <map>

namespace algebra {

// Boolean operations on MultiPolygon2, exact for an exact T.
//
// The method is: cut every edge of both regions at every crossing with every other edge, so that
// the resulting fragments have pairwise disjoint interiors; decide for each fragment whether the
// operation's result is inside on one side and outside on the other, which is exactly the condition
// for the fragment to lie on the result's boundary; then stitch the surviving fragments into rings.
//
// The interesting part is classifying the two sides of a fragment without an epsilon. For the
// midpoint M of a fragment and the normal n, the sample point is M + n*d where d is half the
// smallest positive t with M + n*t on any edge. No edge lies strictly between M and the sample, and
// the sample is not on an edge, so the ordinary winding rule applies there and stays exact.
//
// Cost is quadratic in the number of edges for the cutting and for the classification. Correctness
// first: a sweep line would bring it down to O(n log n) and can replace the internals later without
// changing this interface.

enum class BoolOp { Union, Intersection, Difference, SymmetricDifference };

constexpr bool __apply_bool_op(BoolOp op, bool a, bool b) {
    switch (op) {
    case BoolOp::Union: return a || b;
    case BoolOp::Intersection: return a && b;
    case BoolOp::Difference: return a && !b;
    case BoolOp::SymmetricDifference: return a != b;
    }
    return false;
}

template<typename T>
struct __Edge {
    Vec2<T> a, b;
};

template<typename T>
constexpr void __collect_edges(const MultiPolygon2<T>& p, std::vector<__Edge<T>>& out) {
    for (const Ring2<T>& ring : p.rings) {
        if (ring.size() < 2)
            continue;
        for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++)
            if (!(ring[j] == ring[i]))
                out.push_back({ring[j], ring[i]});
    }
}

// Cuts `edges` at every point where they meet, so no two fragments cross or overlap except at
// their endpoints. Reuses segment_segment_intersection, which reports a single point or, for a
// collinear overlap, the two endpoints of the shared part.
template<typename T>
constexpr std::vector<__Edge<T>> __cut_edges(const std::vector<__Edge<T>>& edges) {
    std::vector<__Edge<T>> out;
    for (size_t i = 0; i < edges.size(); i++) {
        const Vec2<T>& a = edges[i].a;
        const Vec2<T>& b = edges[i].b;
        const Vec2<T> dir = b - a;

        // parameters along a->b where this edge has to be cut, as t in [0, 1]
        std::vector<T> cuts{T(0), T(1)};
        for (size_t j = 0; j < edges.size(); j++) {
            if (i == j)
                continue;
            const auto r = segment_segment_intersection<T>(a, b, edges[j].a, edges[j].b);
            if (const Vec2<T>* p = std::get_if<Vec2<T>>(&r)) {
                cuts.push_back(div_colinear(*p - a, dir));
                continue;
            }
            if (const auto* s = std::get_if<std::pair<Vec2<T>, Vec2<T>>>(&r)) {
                cuts.push_back(div_colinear(s->first - a, dir));
                cuts.push_back(div_colinear(s->second - a, dir));
            }
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

        for (size_t k = 0; k + 1 < cuts.size(); k++) {
            if (cuts[k] < 0 || cuts[k + 1] > 1)
                continue; // a crossing outside the segment, which div_colinear can report
            out.push_back({a + dir * cuts[k], a + dir * cuts[k + 1]});
        }
    }
    return out;
}

// Half the smallest positive t with `from + dir*t` on any of `edges`, or 1 when there is none.
// Walking that far from `from` along `dir` cannot leave the face that `from` borders.
template<typename T>
constexpr T __safe_step(const std::vector<__Edge<T>>& edges, const Vec2<T>& from, const Vec2<T>& dir) {
    bool found = false;
    T best = 1;
    for (const __Edge<T>& e : edges) {
        // from + dir*t == e.a + (e.b - e.a)*s, for some s in [0, 1] and t > 0
        T t, s, det;
        if (!__solve_linear<T>(from - e.a, dir, e.a - e.b, t, s, det)) {
            if (det == 0)
                continue; // parallel, so it either misses or overlaps and hits an endpoint anyway
            continue;
        }
        if (det < 0) {
            negate(t);
            negate(s);
            negate(det);
        }
        if (t <= 0 || s < 0 || s > det)
            continue;
        const T tt = t / det;
        if (!found || tt < best) {
            found = true;
            best = tt;
        }
    }
    return found ? best / 2 : T(1);
}

template<typename T>
constexpr bool __inside(const MultiPolygon2<T>& p, const Vec2<T>& sample) {
    return (winding_number(p, sample) != 0) != p.complement;
}

// Joins fragments into closed rings. Each fragment is oriented so the result's interior is on its
// left, which makes the out degree of every vertex equal its in degree, so following unused
// outgoing fragments always closes a ring.
template<typename T>
constexpr std::vector<Ring2<T>> __stitch(std::vector<__Edge<T>> frags) {
    std::map<std::pair<T, T>, std::vector<size_t>> outgoing;
    for (size_t i = 0; i < frags.size(); i++)
        outgoing[{frags[i].a.x, frags[i].a.y}].push_back(i);

    std::vector<bool> used(frags.size(), false);
    std::vector<Ring2<T>> rings;
    for (size_t start = 0; start < frags.size(); start++) {
        if (used[start])
            continue;
        Ring2<T> ring;
        size_t cur = start;
        while (true) {
            used[cur] = true;
            ring.push_back(frags[cur].a);
            const Vec2<T> to = frags[cur].b;
            if (to == frags[start].a)
                break;
            auto it = outgoing.find({to.x, to.y});
            size_t next = frags.size();
            if (it != outgoing.end())
                for (size_t c : it->second)
                    if (!used[c]) {
                        next = c;
                        break;
                    }
            if (next == frags.size())
                break; // no continuation; the input was not a closed boundary
            cur = next;
        }
        if (ring.size() >= 3)
            rings.push_back(std::move(ring));
    }
    return rings;
}

template<typename T>
constexpr MultiPolygon2<T> boolean_op(BoolOp op, const MultiPolygon2<T>& a, const MultiPolygon2<T>& b) {
    MultiPolygon2<T> out;
    // Far from every edge only the complement flags matter, and that fixes whether the result is
    // bounded. The rings below then bound the finite structure inside it.
    out.complement = __apply_bool_op(op, a.complement, b.complement);

    std::vector<__Edge<T>> edges;
    __collect_edges(a, edges);
    __collect_edges(b, edges);
    if (edges.empty())
        return out;

    std::vector<__Edge<T>> frags = __cut_edges(edges);

    // Two coincident fragments, from a shared edge, must contribute one boundary piece and not two.
    auto key = [](const __Edge<T>& e) {
        const bool swap = e.b.x < e.a.x || (e.b.x == e.a.x && e.b.y < e.a.y);
        const Vec2<T>& p = swap ? e.b : e.a;
        const Vec2<T>& q = swap ? e.a : e.b;
        return std::array<T, 4>{p.x, p.y, q.x, q.y};
    };
    std::sort(frags.begin(), frags.end(), [&](const __Edge<T>& x, const __Edge<T>& y) { return key(x) < key(y); });
    frags.erase(std::unique(frags.begin(), frags.end(),
                            [&](const __Edge<T>& x, const __Edge<T>& y) { return key(x) == key(y); }),
                frags.end());

    std::vector<__Edge<T>> kept;
    for (const __Edge<T>& f : frags) {
        const Vec2<T> mid = (f.a + f.b) / T(2);
        const Vec2<T> dir = f.b - f.a;
        const Vec2<T> normal{-dir.y, dir.x}; // to the left of a->b

        const Vec2<T> left = mid + normal * __safe_step(edges, mid, normal);
        const Vec2<T> right = mid - normal * __safe_step(edges, mid, -normal);

        const bool in_left = __apply_bool_op(op, __inside(a, left), __inside(b, left));
        const bool in_right = __apply_bool_op(op, __inside(a, right), __inside(b, right));
        if (in_left == in_right)
            continue; // the same on both sides, so this is not on the result's boundary

        // the interior goes on the left, which orients the rings counter clockwise
        kept.push_back(in_left ? f : __Edge<T>{f.b, f.a});
    }

    out.rings = __stitch(kept);
    return out;
}

template<typename T>
constexpr MultiPolygon2<T> operator|(const MultiPolygon2<T>& a, const MultiPolygon2<T>& b) {
    return boolean_op(BoolOp::Union, a, b);
}

template<typename T>
constexpr MultiPolygon2<T> operator&(const MultiPolygon2<T>& a, const MultiPolygon2<T>& b) {
    return boolean_op(BoolOp::Intersection, a, b);
}

template<typename T>
constexpr MultiPolygon2<T> operator-(const MultiPolygon2<T>& a, const MultiPolygon2<T>& b) {
    return boolean_op(BoolOp::Difference, a, b);
}

template<typename T>
constexpr MultiPolygon2<T> operator^(const MultiPolygon2<T>& a, const MultiPolygon2<T>& b) {
    return boolean_op(BoolOp::SymmetricDifference, a, b);
}

}
