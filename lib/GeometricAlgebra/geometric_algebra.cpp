#include "GeometricAlgebra/geometric_algebra.h"
#include <math.h>
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


float
Vec_Magnitude(Vec const& v1)
{
    return sqrtf(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
}


float
Vec_Dot2D(Vec const& a, Vec const& b)
{
    return a.x * b.x + a.y * b.y;
}


float
Vec_Wedge2D(Vec const& a, Vec const& b)
{
    return a.x * b.y - a.y * b.x;
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


Rotor2D
Vec_Mul2D(Vec a, Vec b)
{
    // Make sure the scalars are normalised, otherwise the rotor will also
    // scale the vectors it is applied to.
    a = Vec_Normalise(a);
    b = Vec_Normalise(b);

    auto s = Vec_Dot2D(a, b);
    auto I = Vec_Wedge2D(a, b);
    return { s, I };
}


static inline Vec
_Vec_Rotate(Vec const& u, float r, float i)
{
    auto x = (u.x * r) - (u.y * i);
    auto y = (u.x * i) + (u.y * r);

    return { x, y, 0.0 };
}


Vec
Vec_Rotate(Vec const& u, float theta)
{
    auto r = cosf(theta);
    auto i = sinf(theta);

    return _Vec_Rotate(u, r, i);
}


Vec
Vec_Rotate(Vec const& u, Rotor2D const& M)
{
    auto r = M.s;
    auto i = M.I;

    return _Vec_Rotate(u, r, i);
}


Vec
Vec_Rotate(Rotor2D const& M, Vec const& u)
{
    auto r = M.s;
    auto i = M.I;

    return _Vec_Rotate(u, r, i);
}