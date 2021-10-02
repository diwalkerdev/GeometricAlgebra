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
    operator-(Vec const& other)
    {
        return {
            this->x - other.x,
            this->y - other.y,
            this->z - other.z
        };
    }
};


// Utility function for printing a Vec.
// The format is "a<x> b<y> c<z>" where a,b and c are floating point numbers
// formatted to 3dp.
void
Vec_Print(Vec const& v1);


// For performance, we could also provide "inplace" operators.
// Inplace operators prevent the allocation of a new vector by modifying v1.
void
Vec_Add_(Vec& v1, Vec const& v2);


// Creates a zero vector (all elements set to zero).
Vec
Vec_Zero();