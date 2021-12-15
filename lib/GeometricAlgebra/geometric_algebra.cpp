#include "GeometricAlgebra/geometric_algebra.h"
#include <stdio.h>


void
Print(char const* text, Vec const& v1)
{
    printf("%s %0.3f<x> %0.3f<y> %0.3f<z>\n", text, v1.x, v1.y, v1.z);
}


void
Print(char const* text, Rotor const& R)
{
    printf("%s %0.3f<r> %0.3f<e12> %0.3f<e13> %0.3f<e23>\n", text, R.s, R.B.e12, R.B.e13, R.B.e23);
}


void
Print(char const* text, TriVector const& T)
{
    printf("%s %0.3f<e123>\n", text, T.e123);
}


Vec
Vec_Zero()
{
    return { 0.f, 0.f, 0.f };
}


float
Vec_Magnitude(Vec const& v1)
{
    return sqrtf(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
}


Vec
Vec_Normalise(Vec const& u)
{
    auto m = Vec_Magnitude(u);
    return {
        u.x / m,
        u.y / m,
        u.z / m
    };
}
