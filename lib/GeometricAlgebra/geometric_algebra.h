#pragma once
#include <cstddef>

struct Vec
{
    // The union, in conjunction with the [] operator, allows Vec to be
    // accessed as if it was a struct, or as an array of elements.
    // Both of the following are valid and equivalent:
    //   v.x = 0;
    //   v[0] = 0;
    union
    {
        struct
        {
            float x, y, z;
        };
        float data[3];
    };


    float&
    operator[](size_t index)
    {
        return this->data[index];
    }

    float
    operator[](size_t index) const
    {
        return this->data[index];
    }


    Vec
    operator+(Vec const& other)
    {
        return {
            this->x + other.x,
            this->y + other.y,
            this->z + other.z
        };
    }


    Vec
    operator+=(Vec const& other)
    {
        this->x += other.x;
        this->y += other.y;
        this->z += other.z;
        return *this;
    }

    Vec
    operator-(Vec const& other)
    {
        return {
            this->x - other.x,
            this->y - other.y,
            this->z - other.z
        };
    }


    Vec
    operator-=(Vec const& other)
    {
        this->x -= other.x;
        this->y -= other.y;
        this->z -= other.z;
        return *this;
    }


    Vec
    operator*(float scalar)
    {
        return {
            this->x * scalar,
            this->y * scalar,
            this->z * scalar
        };
    }


    Vec
    operator*=(float scalar)
    {
        this->x *= scalar;
        this->y *= scalar;
        this->z *= scalar;
        return *this;
    }
};


// Rotor2D is really a Multivector without the basis components.
//
// We've removed the basis components, as _geneally_ it is clearer if we keep
// vectors seperate from the 'geometric parts'.
//
// To make this more concrete, we _could_ represent two othoginal vectors as:
// u = 0s, 1e1, 0e2, 0I and v = 0s, 0e1, 1e1, 0I. Multiply them together to get
// w = 0s, 0e1, 0e2, 1I. Finally, we could rotate u by w giving u_prime.
//
// The problem with this approach, is that u, v, w and u_prime are all
// multivectors. This introduces ambiguity in that u and v are likely to
// reprensent things we want to draw on the screen, like position and veloctity,
// where w is really just an operator.
//
// Having a different type, which represents the result of a geometric product,
// keeps a clear and distinct different between vectors and a geometric operation.
struct Rotor2D
{
    float s;
    float I; // The bivector e1e2.
};


// Utility function for printing a Vec.
// The format is "a<x> b<y> c<z>" where a,b and c are floating point numbers
// formatted to 3dp.
void
Vec_Print(Vec const& v1);


// Creates a zero vector (all elements set to zero).
Vec
Vec_Zero();


float
Vec_Magnitude(Vec const& v1);


float
Vec_Dot2D(Vec const& a, Vec const& b);


float
Vec_Wedge2D(Vec const& a, Vec const& b);


Rotor2D
Vec_Mul2D(Vec a, Vec b);


// Rotate the vector u by the angle theta.
//
// Internally, this function uses the standard form of eulers formula to
// perform the rotation.
Vec
Vec_Rotate(Vec const& u, float theta);


// Rotate the vector u by the multivector M.
//
// This function is useful if you to need to rotate a vector u by an
// angle defined by two different vectors a and b.
// a * b generates a multivector with a scalar and bivector part. The
// scalar and bivector are closely related to the real and imaginary
// parts of a complex number. Therefore the multivector can be used to
// rotated the vector u. See below for more details.
//
// Rotation by a multivector is closely related to rotating by a complex
// number - a multivector can be thought of having a 'real' part (s) and
// an imaginary part (I).
//
// Internally, this function still uses eulers formula to perform the
// rotation, however, we now no longer need to use sin and cosine to
// calculate the real and imaginary parts, since they are already
// contained within the multivector.
//
Vec
Vec_Rotate(Vec const& u, Rotor2D const& M);


Vec
Vec_Rotate(Rotor2D const& M, Vec const& u);


Vec
Vec_Normalise(Vec const& u);