#ifndef STDINCLUDE
#define STDINCLUDE

#ifdef NVCC
#define HOSTDEVICE __host__ __device__ inline
#else
#define HOSTDEVICE inline __attribute__((always_inline))
#endif

#include "vector_types.h"
#include "vector_functions.h"

#define PI 3.14159265358979323846

#define Dscalar double
#define Dscalar2 double2
#define ncDscalar ncDouble

#define cur_norm curand_normal_double

HOSTDEVICE bool operator<(const Dscalar2 &a, const Dscalar2 &b)
    {
    return a.x<b.x;
    }

HOSTDEVICE Dscalar2 make_Dscalar2(Dscalar x, Dscalar y)
    {
    Dscalar2 ans;
    ans.x =x;
    ans.y=y;
    return ans;
    }




#undef HOSTDEVICE

#endif

