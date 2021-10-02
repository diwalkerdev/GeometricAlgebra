#include "GeometricAlgebra/geometric_algebra.h"
#include <stdio.h>


void
Vec_Print(Vec const& v1)
{
    printf("%0.3f<x> %0.3f<y> %0.3f<z>\n", v1.x, v1.y, v1.z);
}


Vec
Vec_Zero()
{
    return { 0.f, 0.f, 0.f };
}


void
Vec_Add_(Vec& v1, Vec const& v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
}
