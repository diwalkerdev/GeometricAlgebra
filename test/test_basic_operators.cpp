#include "GeometricAlgebra/geometric_algebra.h"

#include <cassert>
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


int
main(void)
{
    Test_SetGetElements();
    Test_VectorAddition();
    printf("%s PASSED\n", "test_basic_operators.cpp");
}