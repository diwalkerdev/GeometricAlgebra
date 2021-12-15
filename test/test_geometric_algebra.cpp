#include "GeometricAlgebra/geometric_algebra.h"

#include <cassert>
#include <math.h>
#include <stdio.h>


void
Test_SetGetElements()
{
    printf(__func__);
    printf("\n");

    Vec v1 { 1.f, 2.f, 3.f };
    assert(v1.x == 1.f);
    assert(v1[0] == 1.f);

    v1.x = 4.f;
    assert(v1.x == 4.f);

    v1[0] = 5.f;
    assert(v1[0] == 5.f);
}


void
Test_VectorAddition()
{
    printf(__func__);
    printf("\n");

    Vec v1 { 1.f, 2.f, 3.f };
    Vec v2 { 1.f, 1.f, 1.f };

    auto v3 = v1 + v2;
    assert(v3.x == 2.f);
    assert(v3.y == 3.f);
    assert(v3.z == 4.f);
}


void
Test_VectorSubtraction()
{
    printf(__func__);
    printf("\n");

    Vec v1 { 1.f, 2.f, 3.f };
    Vec v2 { 1.f, 1.f, 1.f };

    auto v3 = v1 - v2;
    assert(v3.x == 0.f);
    assert(v3.y == 1.f);
    assert(v3.z == 2.f);
}


void
Test_DotProduct()
{
    printf(__func__);
    printf("\n");

    Vec v1 { 1.0, 0.0, 0.0 };
    Vec v2 { 0.0, 1.0, 0.0 };

    assert(Vec_Dot(v1, v1) == 1.0); // Dot of parallel vectors == 1.
    assert(Vec_Dot(v1, v2) == 0.0); // Dot of orthogonal vectors == 0.
}


// void
// Test_WedgeProduct()
// {
//     printf(__func__);
//     printf("\n");

//     Vec v1 { 1.0, 0.0, 0.0 };
//     Vec v2 { 0.0, 1.0, 0.0 };

//     assert(Vec_Wedge(v1, v1) == 0.0); // Wedge (or area) of parallel vectors == 0.
//     assert(Vec_Wedge(v1, v2) == 1.0); // Wedge (or area) of orthogonal vectors == 1.
// }


// void
// Test_RotateByMultiVector()
// {
//     printf(__func__);
//     printf("\n");

//     // Define a multivector using two orthoginal vectors.
//     Vec  u { 1.0, 0.0, 0.0 };
//     Vec  v { 0.0, 1.0, 0.0 };
//     auto M = Vec_Mul(u, v);

//     // Then rotate u by the multivector M. Which is equivalent to rotating
//     // by 90 degrees.
//     Vec u_prime = Vec_Mul(M, u);
//     assert(u_prime.x < +0.01);
//     assert(u_prime.x > -0.01);
//     assert(u_prime.y < +1.01);
//     assert(u_prime.y > +0.99);

//     // Test rotor constructed with non-unit vector doesn't scale rotated vector.
//     Vec  w { 100.0, 0.0, 0.0 };
//     auto W  = Vec_Mul(w, v);
//     u_prime = Vec_Mul(u, W);
//     assert(u_prime.x < +0.01);
//     assert(u_prime.x > -0.01);
//     assert(u_prime.y < +1.01);
//     assert(u_prime.y > +0.99);

//     // Swapping u and v rotates the other direction.
//     auto N  = Vec_Mul(v, u);
//     u_prime = Vec_Mul(u, N);
//     assert(u_prime.x < +0.01);
//     assert(u_prime.x > -0.01);
//     assert(u_prime.y > -1.01);
//     assert(u_prime.y < -0.99);
// }


void
Test_Rotate3D()
{
    printf(__func__);
    printf("\n");

    // Define two orthoginal vectors in the XY plane.
    {
        Vec u = { 1.0, 0.0, 0.0 };
        Vec v = { 0.0, 1.0, 0.0 };

        float    s;
        BiVector B;
        std::tie(s, B) = Vec_Mul(u, v);

        assert(s == 0.0);
        assert(B.e12 == 1.0);
        assert(B.e13 == 0.0);
        assert(B.e23 == 0.0);
    }

    // Define two orthoginal vectors in the XZ plane.
    {
        Vec u = { 1.0, 0.0, 0.0 };
        Vec v = { 0.0, 0.0, 1.0 };

        float    s;
        BiVector B;
        std::tie(s, B) = Vec_Mul(u, v);

        assert(s == 0.0);
        assert(B.e12 == 0.0);
        assert(B.e13 == 1.0);
        assert(B.e23 == 0.0);
    }

    {
        // Define two orthoginal vectors in the ZY plane. (Yes, these are reversed).
        Vec u = { 0.0, 0.0, 1.0 };
        Vec v = { 0.0, 1.0, 0.0 };

        float    s;
        BiVector B;
        std::tie(s, B) = Vec_Mul(u, v);

        assert(s == 0.0);
        assert(B.e12 == 0.0);
        assert(B.e13 == 0.0);
        assert(B.e23 == -1.0);
    }

    // Define a pair of parallel vectors.
    {
        Vec u = { 0.0, 0.0, 1.0 };
        Vec v = { 0.0, 0.0, 1.0 };

        float    s;
        BiVector B;
        std::tie(s, B) = Vec_Mul(u, v);

        assert(s == 1.0);
        assert(B.e12 == 0.0);
        assert(B.e13 == 0.0);
        assert(B.e23 == 0.0);
    }

    // BiVector-Vector product.
    // Gives a Vec and a TriVector. However, since u is already in the plane
    // (uvu), the TriVector will always be 0.
    //
    {
        Vec u { 1.0, 0.0, 0.0 };
        // Vec v { 0.0, 1.0, 0.0 };
        // Vec v { 0.707, 0.707, 0.0 };
        Vec v { 0.707, 0.707, 0.707 };

        // uv
        float    s;
        BiVector B;
        std::tie(s, B) = Vec_Mul(u, v);

        // Reflect = -uvu
        // uvu
        Vec       w;
        TriVector T;
        std::tie(T, w) = Vec_Mul(Rotor { s, B }, u);

        // Remember the reflection formula.
        w *= -1;
        Print("w (uvu): ", w);
        Print("T: ", T);
        // TODO(DW): Check reflection result.
        // assert(w.x <= -0.707 && w.x >= -0.708);
        // assert(w.y >= +0.707 && w.y <= +0.708);
        // assert(w.z >= +0.707 && w.z <= +0.708);
    }


    // BiVector-Vector product.
    // Trivector result should be non-zero as w is not in the plane of uv.
    //
    {
        Vec const u { 1.0, 0.0, 0.0 };
        Vec const v { 0.0, 1.0, 0.0 };

        // uv
        auto R = Rotor(Vec_Mul(u, v));

        // Rw
        Vec       x;
        TriVector T;
        std::tie(T, x) = Vec_Mul(R, x);

        Print("x: ", x);
        Print("T: ", T);

        assert(T.e123 != 0.0);
    }

    // Rotate 3D manually => uvwvu => Ruv w Rvu
    // This requires that you calculate both rotors Ruv and Rvu.
    //
    {
        Vec u { 1.0, 0.0, 0.0 };
        Vec v { 0.0, 1.0, 0.0 };
        Vec w { 0.707, 0.707, 0.707 };

        Rotor Ruv { Vec_Mul(u, v) };
        Rotor Rvu { Vec_Mul(v, u) };
        Print("Ruv: ", Ruv);
        Print("Rvu: ", Rvu);

        Vec       w_hat;
        TriVector T;
        std::tie(T, w_hat) = Vec_Mul(Ruv, w);

        Vec       w_dot;
        TriVector Q;
        std::tie(w_dot, Q) = Vec_Mul({ T, w_hat }, Rvu);

        Print("w (manual): ", w_dot);
        Print("Q: ", Q);

        assert(w_dot.x <= -0.707 && w_dot.x >= -0.708);
        assert(w_dot.y <= -0.707 && w_dot.y >= -0.708);
        assert(w_dot.z >= +0.707 && w_dot.z <= +0.708);
    }

    // Rotate 3D uvwvu using shortcut => Ruv w Rvu => Ruv w R*uv.
    // This requires only 1 rotor.
    {
        Vec u { 1.0, 0.0, 0.0 };
        Vec v { 0.0, 1.0, 0.0 };
        Vec w { 0.707, 0.707, 0.707 };

        Rotor Ruv { Vec_Mul(u, v) };
        Vec   w_dot = Vec_Rotate(Ruv, w);
        Print("w (rotate): ", w_dot);

        assert(w_dot.x <= -0.707 && w_dot.x >= -0.708);
        assert(w_dot.y <= -0.707 && w_dot.y >= -0.708);
        assert(w_dot.z >= +0.707 && w_dot.z <= +0.708);
    }

    // Check that composing rotations works as expected. There are two
    // ways to do this.
    // 1) Sandwitch a rotation in another rotation R2R1wR1'R2'.
    // 2) Multiply the two rotors together and do one roation:
    //        R3 = R1*R2
    //        w_dot = R3wR3'
    {
        // Note the vectors chosen causes the starting vector u to be rotated
        // back to its starting position.
        auto const u = Vec { 1.0f, 0.0f, 0.0f };
        auto const v = Vec { 0.0f, 1.0f, 0.0f };
        auto const w = Vec { 0.0f, 0.0f, 1.0f };

        auto Ruv = Rotor { Vec_Mul(u, v) };
        auto Ruw = Rotor { Vec_Mul(u, w) };

        // 1
        {
            auto w_dot = Vec_Rotate(Ruv, u);
            w_dot      = Vec_Rotate(Ruw, w_dot);

            assert(w_dot.x == 1.0f);
            assert(w_dot.y == 0.0f);
            assert(w_dot.z == 0.0f);
        }

        // 2
        {
            auto Ruwuv   = Geo_Mul(Ruw, Ruv);
            auto w_Ruwuv = Vec_Rotate(Ruwuv, u);

            assert(w_Ruwuv.x == 1.0f);
            assert(w_Ruwuv.y == 0.0f);
            assert(w_Ruwuv.z == 0.0f);
        }
    }
}


int
main(void)
{
    Test_SetGetElements();
    Test_VectorAddition();
    Test_VectorSubtraction();
    Test_DotProduct();
    // Test_WedgeProduct();
    // Test_RotateByTheta();
    // Test_RotateByMultiVector();
    Test_Rotate3D();

    printf("%s PASSED\n", "test_basic_operators.cpp");
}