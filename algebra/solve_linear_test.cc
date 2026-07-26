#include "algebra/solve_linear.h"
#include "algebra/rational.h"
#include "algebra/__test.h"

using V2 = Vec2<rational>;
using V3 = Vec3<rational>;

TEST_CASE("solve_linear 2d") {
    rational s, t;

    // (1, 0) + s * (-1, 0) + t * (0, 1) == 0  =>  s == 1, t == 0
    V2 a{1, 0}, b{-1, 0}, c{0, 1};
    REQUIRE(solve_linear(a, b, c, &s, &t));
    REQUIRE(s == 1);
    REQUIRE(t == 0);
    REQUIRE(is_zero(a + b * s + c * t));

    // (5, 7) + s * (1, 0) + t * (0, 1) == 0  =>  s == -5, t == -7
    V2 a2{5, 7}, b2{1, 0}, c2{0, 1};
    REQUIRE(solve_linear(a2, b2, c2, &s, &t));
    REQUIRE(s == -5);
    REQUIRE(t == -7);

    // no unique solution
    V2 b3{1, 2}, c3{2, 4};
    REQUIRE(!solve_linear(a2, b3, c3, &s, &t));

    Random rng(3);
    for (int i = 0; i < 50; i++) {
        auto r = [&] { return rational(rng.Uniform<int>(-9, 9), rng.Uniform<int>(1, 5)); };
        V2 A{r(), r()}, B{r(), r()}, C{r(), r()};
        if (determinant(B, C) == 0)
            continue;
        REQUIRE(solve_linear(A, B, C, &s, &t));
        REQUIRE(is_zero(A + B * s + C * t));

        rational ss, tt, det;
        REQUIRE(__solve_linear(A, B, C, ss, tt, det));
        REQUIRE(ss / det == s);
        REQUIRE(tt / det == t);
    }
}

TEST_CASE("solve_linear 3d") {
    rational s, t, r;

    // (5, 7, 9) + s * (1,0,0) + t * (0,1,0) + r * (0,0,1) == 0
    V3 a{5, 7, 9}, b{1, 0, 0}, c{0, 1, 0}, d{0, 0, 1};
    REQUIRE(solve_linear(a, b, c, d, &s, &t, &r));
    REQUIRE(s == -5);
    REQUIRE(t == -7);
    REQUIRE(r == -9);
    REQUIRE(is_zero(a + b * s + c * t + d * r));

    Random rng(5);
    for (int i = 0; i < 50; i++) {
        auto q = [&] { return rational(rng.Uniform<int>(-9, 9), rng.Uniform<int>(1, 5)); };
        V3 A{q(), q(), q()}, B{q(), q(), q()}, C{q(), q(), q()}, D{q(), q(), q()};
        if (determinant(B, C, D) == 0)
            continue;
        REQUIRE(solve_linear(A, B, C, D, &s, &t, &r));
        REQUIRE(is_zero(A + B * s + C * t + D * r));

        rational ss, tt, rr, det;
        REQUIRE(__solve_linear(A, B, C, D, ss, tt, rr, det));
        REQUIRE(ss / det == s);
        REQUIRE(tt / det == t);
        REQUIRE(rr / det == r);
    }
}

TEST_CASE("solve_linear 1 unknown") {
    // A + B*x == 0
    Vec2<rational> a{6, 9}, b{-2, -3};
    auto res = solve_linear(a, b);
    REQUIRE(std::holds_alternative<rational>(res));
    REQUIRE(std::get<rational>(res) == 3);
}

TEST_CASE("solve_linear 1 unknown - inconsistent systems have no solution") {
    // A + B*x == 0 has to hold in every component, so a system that cannot be satisfied
    // must report None rather than the value that happens to solve the first component
    auto res = solve_linear(Vec2<rational>{1, 1}, Vec2<rational>{1, 2});
    REQUIRE(std::holds_alternative<None>(res));   // x = -1 solves x, leaves y at -1
    res = solve_linear(Vec2<rational>{1, 5}, Vec2<rational>{1, 1});
    REQUIRE(std::holds_alternative<None>(res));
    // the first component being degenerate must not change that
    res = solve_linear(Vec2<rational>{1, 1}, Vec2<rational>{0, 2});
    REQUIRE(std::holds_alternative<None>(res));
    // three components, inconsistent only in the last
    res = solve_linear(Vec3<rational>{2, 4, 7}, Vec3<rational>{1, 2, 1});
    REQUIRE(std::holds_alternative<None>(res));

    // and the consistent cases still solve
    res = solve_linear(Vec3<rational>{2, 4, 6}, Vec3<rational>{1, 2, 3});
    REQUIRE(std::holds_alternative<rational>(res));
    REQUIRE(std::get<rational>(res) == -2);
    res = solve_linear(Vec2<rational>{0, 1}, Vec2<rational>{0, 2});
    REQUIRE(std::holds_alternative<rational>(res));
    REQUIRE(std::get<rational>(res) == rational(-1, 2));

    // A == 0 with B != 0 is solved by x == 0
    res = solve_linear(Vec2<rational>{0, 0}, Vec2<rational>{1, 2});
    REQUIRE(std::holds_alternative<rational>(res));
    REQUIRE(std::get<rational>(res) == 0);
    // B == 0 with A != 0 has no solution, both zero is solved by anything
    REQUIRE(std::holds_alternative<None>(solve_linear(Vec2<rational>{1, 0}, Vec2<rational>{0, 0})));
    REQUIRE(std::holds_alternative<Any>(solve_linear(Vec2<rational>{0, 0}, Vec2<rational>{0, 0})));
}

TEST_CASE("determinant") {
    using V2 = Vec2<rational>;
    using V3 = Vec3<rational>;
    REQUIRE(determinant(V2{1, 0}, V2{0, 1}) == 1);
    REQUIRE(determinant(V2{0, 1}, V2{1, 0}) == -1);   // swapping negates
    REQUIRE(determinant(V2{2, 3}, V2{4, 5}) == -2);
    REQUIRE(determinant(V2{1, 2}, V2{2, 4}) == 0);    // parallel

    REQUIRE(determinant(V3{1, 0, 0}, V3{0, 1, 0}, V3{0, 0, 1}) == 1);
    REQUIRE(determinant(V3{0, 1, 0}, V3{1, 0, 0}, V3{0, 0, 1}) == -1);
    // a cyclic rotation preserves it, a single swap negates it
    const V3 a{1, 2, 3}, b{4, 5, 6}, c{7, 8, 10};
    const rational d = determinant(a, b, c);
    REQUIRE(d != 0);
    REQUIRE(determinant(b, c, a) == d);
    REQUIRE(determinant(c, a, b) == d);
    REQUIRE(determinant(b, a, c) == -d);
    REQUIRE(determinant(a, c, b) == -d);
    // linearly dependent rows
    REQUIRE(determinant(V3{1, 2, 3}, V3{2, 4, 6}, V3{7, 8, 9}) == 0);
}

TEST_CASE("__solve_linear returns the solution scaled by the determinant") {
    using V2 = Vec2<rational>;
    using V3 = Vec3<rational>;
    // the scaled variants must agree with the dividing ones
    const V2 a{-5, -6}, b{1, 0}, c{0, 1};
    rational s, t, det;
    REQUIRE(__solve_linear(a, b, c, s, t, det));
    rational s2, t2;
    REQUIRE(solve_linear(a, b, c, &s2, &t2));
    REQUIRE(s / det == s2);
    REQUIRE(t / det == t2);
    REQUIRE(a + b * s2 + c * t2 == V2{0, 0});

    // degenerate: B and C parallel
    REQUIRE(!__solve_linear(V2{1, 1}, V2{1, 2}, V2{2, 4}, s, t, det));
    REQUIRE(!solve_linear(V2{1, 1}, V2{1, 2}, V2{2, 4}, &s, &t));

    const V3 a3{-1, -2, -3}, b3{1, 0, 0}, c3{0, 1, 0}, d3{0, 0, 1};
    rational r3, s3, t3, det3;
    REQUIRE(__solve_linear(a3, b3, c3, d3, s3, t3, r3, det3));
    rational s4, t4, r4;
    REQUIRE(solve_linear(a3, b3, c3, d3, &s4, &t4, &r4));
    REQUIRE(s3 / det3 == s4);
    REQUIRE(t3 / det3 == t4);
    REQUIRE(r3 / det3 == r4);
    REQUIRE(a3 + b3 * s4 + c3 * t4 + d3 * r4 == V3{0, 0, 0});

    // the output pointers are optional
    REQUIRE(solve_linear(a, b, c, static_cast<rational*>(nullptr), static_cast<rational*>(nullptr)));
    REQUIRE(solve_linear(a3, b3, c3, d3, static_cast<rational*>(nullptr), static_cast<rational*>(nullptr),
                         static_cast<rational*>(nullptr)));
}
