#pragma once
#include <cstddef>
#include <math.h>
#include <tuple>

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


struct BiVector
{
    union
    {
        struct
        {
            float e12, e13, e23;
        };
        float data[3];
    };
};


struct TriVector
{
    float e123;
};


struct Rotor
{
    float    s;
    BiVector B;

    Rotor(float s, BiVector B)
    {
        this->s = s;
        this->B = B;
    }
    Rotor(std::tuple<float, BiVector> s_B)
    {
        this->s = std::get<0>(s_B);
        this->B = std::get<1>(s_B);
    }
};


void
Print(char const* text, TriVector const& T);


void
Print(char const* text, Rotor const& R);


// Utility function for printing a Vec.
// The format is "a<x> b<y> c<z>" where a,b and c are floating point numbers
// formatted to 3dp.
void
Print(const char* text, Vec const& v1);

template <typename Tp>
Tp
Square(Tp const& x)
{
    return x * x;
}

inline float
Geo_Length(Rotor const& R)
{
    return sqrtf(Square(R.s) + Square(R.B.e12) + Square(R.B.e13) + Square(R.B.e23));
}


inline void
Geo_Normalise(Rotor& R)
{
    auto l = Geo_Length(R);
    R.s /= l;
    R.B.e12 /= l;
    R.B.e13 /= l;
    R.B.e23 /= l;
}

// Creates a zero vector (all elements set to zero).
Vec
Vec_Zero();


float
Vec_Magnitude(Vec const& v1);


// Rotate the vector u by the angle theta.
//
// Internally, this function uses the standard form of eulers formula to
// perform the rotation.
Vec
Vec_Rotate(Vec const& u, float theta);


Vec
Vec_Normalise(Vec const& u);


inline float
Vec_Dot(Vec const& a, Vec const& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


inline BiVector
Vec_Wedge(Vec const& a, Vec const& b)
{
    return {
        a.x * b.y - b.x * a.y,
        a.x * b.z - b.x * a.z,
        a.y * b.z - b.y * a.z
    };
}


inline std::tuple<float, BiVector>
Vec_Mul(Vec const& a, Vec const& b)
{
    Rotor R = {
        Vec_Dot(a, b),
        Vec_Wedge(a, b)
    };
    Geo_Normalise(R);

    // TODO(DW): We should just return a rotor here.
    return { R.s, R.B };
}

inline float
Vec_Times(Vec const& u, Vec&& v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

inline std::tuple<TriVector, Vec>
Vec_Mul(std::tuple<float, BiVector> const& M, Vec const& v)
{
    auto& s = std::get<0>(M);
    auto& B = std::get<1>(M);

    auto x = Vec_Times(v, { s, +B.e12, +B.e13 });
    auto y = Vec_Times(v, { -B.e12, s, +B.e23 });
    auto z = Vec_Times(v, { -B.e13, -B.e23, s });

    auto T = Vec_Times(v, { B.e23, -B.e13, B.e12 });

    return {
        TriVector { T },
        { x, y, z },
    };
}


inline std::tuple<TriVector, Vec>
Vec_Mul(Rotor const& R, Vec const& v)
{
    return Vec_Mul(std::tuple<float, BiVector> { R.s, R.B }, v);
}


inline std::tuple<Vec, TriVector>
Vec_Mul(std::tuple<TriVector, Vec> const& M, std::tuple<float, BiVector> R)
{
    auto& T = std::get<0>(M);
    auto& w = std::get<1>(M);
    auto& s = std::get<0>(R);
    auto& B = std::get<1>(R);

    // auto x = s * w[0] + w[1] * B.b1 + w[2] * B.b2 + B.b3 * T.e123;
    // auto y = s * w[1] - w[0] * B.b1 + w[2] * B.b3 - B.b2 * T.e123;
    // auto z = s * w[2] - w[0] * B.b2 - w[1] * B.b3 + B.b1 * T.e123;
    // auto Q = w[0] * B.b3 - w[1] * B.b2 + w[2] * B.b1;
    auto x = Vec_Times(w, { s, -B.e12, -B.e13 }) - T.e123 * B.e23;
    auto y = Vec_Times(w, { B.e12, s, -B.e23 }) + T.e123 * B.e13;
    auto z = Vec_Times(w, { B.e13, B.e23, s }) - T.e123 * B.e12;
    auto Q = w[0] * B.e23 - w[1] * B.e13 + w[2] * B.e12;

    return {
        Vec { x, y, z },
        TriVector { (T.e123 * s) + Q },
    };
}


inline Vec
Vec_Rotate(Rotor const& Ruv, Vec const& q)
{
    auto& s = Ruv.s;
    auto& B = Ruv.B;

    Vec       w;
    TriVector T;
    std::tie(T, w) = Vec_Mul(Ruv, q);

    auto x = s * w[0] + w[1] * B.e12 + w[2] * B.e13 + B.e23 * T.e123;
    auto y = s * w[1] - w[0] * B.e12 + w[2] * B.e23 - B.e13 * T.e123;
    auto z = s * w[2] - w[0] * B.e13 - w[1] * B.e23 + B.e12 * T.e123;

    return { x, y, z };

    // Trivector should always be 0 for a rotation.
    // auto Q = w[0] * B.b3 - w[1] * B.b2 + w[2] * B.b1;
    // return {
    //     Vec { x, y, z },
    //     TriVector { (T.e123 * s) + Q },
    // };
}


inline std::tuple<Vec, TriVector>
Vec_Mul(std::tuple<TriVector, Vec> const& M, Rotor const& R)
{
    return Vec_Mul(M, std::tuple<float, BiVector> { R.s, R.B });
}


inline Rotor
Geo_Mul(Rotor const& X, Rotor const& Y)
{
    auto const& a  = X.s;
    auto const& B1 = X.B;

    auto const& b  = Y.s;
    auto const& B2 = Y.B;

    Rotor R = {
        a * b + B1.e12 * B2.e12 + B1.e13 * B2.e13 + B1.e23 * B2.e23, // scalar
        { a * B2.e12 + b * B1.e12 + B1.e13 * B2.e23 - B1.e23 * B2.e13, // e12
          a * B2.e13 + B1.e12 * B2.e23 + b * B1.e13 - B1.e23 * B2.e12, // e13
          a * B2.e23 - B1.e12 * B2.e23 + B1.e13 * B2.e12 + b * B1.e23 } // e23
    };
    Geo_Normalise(R);
    return R;
}


struct Matrix4
{
    float data[4 * 4];

    float&
    operator[](int x)
    {
        return data[x];
    }
};


struct Quat
{
    float data[4];
};


inline Matrix4
ToMatrix4(Rotor const& R)
{
    auto v0 = Vec_Rotate(R, { 1, 0, 0 });
    auto v1 = Vec_Rotate(R, { 0, 1, 0 });
    auto v2 = Vec_Rotate(R, { 0, 0, 1 });

    Matrix4 mat;

    mat[0] = v0[0];
    mat[1] = v0[1];
    mat[2] = v0[2];
    mat[3] = 0;

    mat[4] = v1[0];
    mat[5] = v1[1];
    mat[6] = v1[2];
    mat[7] = 0;

    mat[8]  = v2[0];
    mat[9]  = v2[1];
    mat[10] = v2[2];
    mat[11] = 0;

    mat[12] = 0;
    mat[13] = 0;
    mat[14] = 0;
    mat[15] = 0;

    return mat;
}