#define __SINGLE_PRECISION 32
#define __DOUBLE_PRECISION 64

#define __KIND __SINGLE_PRECISION
#include "m2p_interpolation_s.c"
#undef   __KIND
#define __KIND __DOUBLE_PRECISION
#include "m2p_interpolation_d.c"
#undef   __KIND

#undef __SINGLE_PRECISION
#undef __DOUBLE_PRECISION
