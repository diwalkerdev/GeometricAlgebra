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

    assert(Vec_Dot2D(v1, v1) == 1.0); // Dot of parallel vectors == 1.
    assert(Vec_Dot2D(v1, v2) == 0.0); // Dot of orthogonal vectors == 0.
}


void
Test_WedgeProduct()
{
    printf(__func__);
    printf("\n");

    Vec v1 { 1.0, 0.0, 0.0 };
    Vec v2 { 0.0, 1.0, 0.0 };

    assert(Vec_Wedge2D(v1, v1) == 0.0); // Wedge (or area) of parallel vectors == 0.
    assert(Vec_Wedge2D(v1, v2) == 1.0); // Wedge (or area) of orthogonal vectors == 1.
}

void
Test_RotateByTheta()
{
    printf(__func__);
    printf("\n");

    Vec u { 1.0, 0.0, 0.0 };

    Vec u_prime = Vec_Rotate(u, M_PI / 2.0);
    assert(u_prime.x < +0.01);
    assert(u_prime.x > -0.01);
    assert(u_prime.y < +1.01);
    assert(u_prime.y > -1.01);
}


void
Test_RotateByMultiVector()
{
    printf(__func__);
    printf("\n");

    // Define a multivector using two orthoginal vectors.
    Vec  u { 1.0, 0.0, 0.0 };
    Vec  v { 0.0, 1.0, 0.0 };
    auto M = Vec_Mul2D(u, v);

    // Then rotate u by the multivector M. Which is equivalent to rotating
    // by 90 degrees.
    Vec u_prime = Vec_Rotate(u, M);
    assert(u_prime.x < +0.01);
    assert(u_prime.x > -0.01);
    assert(u_prime.y < +1.01);
    assert(u_prime.y > +0.99);

    // Test rotor constructed with non-unit vector doesn't scale rotated vector.
    Vec  w { 100.0, 0.0, 0.0 };
    auto W  = Vec_Mul2D(w, v);
    u_prime = Vec_Rotate(u, W);
    assert(u_prime.x < +0.01);
    assert(u_prime.x > -0.01);
    assert(u_prime.y < +1.01);
    assert(u_prime.y > +0.99);

    // Swapping u and v rotates the other direction.
    auto N  = Vec_Mul2D(v, u);
    u_prime = Vec_Rotate(u, N);
    assert(u_prime.x < +0.01);
    assert(u_prime.x > -0.01);
    assert(u_prime.y > -1.01);
    assert(u_prime.y < -0.99);
}


int
main(void)
{
    Test_SetGetElements();
    Test_VectorAddition();
    Test_VectorSubtraction();
    Test_DotProduct();
    Test_WedgeProduct();
    Test_RotateByTheta();
    Test_RotateByMultiVector();

    printf("%s PASSED\n", "test_basic_operators.cpp");
}