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
