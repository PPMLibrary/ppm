#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

// C/Fortran wrapper code for Joachim Stadel's Hilbert key sorting code
// This is just prototypic code and probably needs to be extended to work
// properly for doubles.
//
// Author: Omar Awile
// Date: 15.11.2012

typedef struct idx {
  uint64_t d;
  int i;
  double x[3];
} idx_t;

int cmp (const void * elem1, const void * elem2) {
  idx_t f = *((idx_t*)elem1);
  idx_t s = *((idx_t*)elem2);
  if (f.d > s.d) return  1;
  if (f.d < s.d) return -1;
  return 0;
}

/**
 * The following code was kindly provided by Joachim Stadel (University of Zurich)
 * after a discussion at the HP2C Workshop in Lugano 20.3.12
 **/
/*
 ** x and y must have range [1,2) !
 */
uint64_t hilbert2d(float x,float y) {
  uint64_t s = 0;
  uint32_t m,ux,uy,ut;

  union {
    float    f;
    uint32_t u;
  } punner;

  punner.f = x; ux = punner.u >> 2;
  punner.f = y; uy = punner.u >> 2;

  m = 0x00100000;

  while (m) {
    s = s << 2;
    if (ux&m) {
      if (uy&m) {
        s |= 2;
      }
      else {
        ut = ux;
        ux = ~uy;
        uy = ~ut;
        s |= 3;
      }
    }
    else {
      if (uy&m) {
        s |= 1;
      }
      else {
        ut = ux;
        ux = uy;
        uy = ut;
      }
    }
    m = m >> 1;
  }
  return s;
}

/*
 ** x, y and z must have range [1,2) !
 */
uint64_t hilbert3d(float x,float y,float z) {
  uint64_t s = 0;
  uint32_t m,ux,uy,uz,ut;

  ux = (*(uint32_t *)&x)>>2;
  uy = (*(uint32_t *)&y)>>2;
  uz = (*(uint32_t *)&z)>>2;

  m = 0x00100000;

  while (m) {
    s = s << 3;

    if (ux&m) {
      if (uy&m) {
        if (uz&m) {
          ut = ux;
          ux = uy;
          uy = ~uz;
          uz = ~ut;
          s |= 5;
        }
        else {
          ut = uz;
          uz = ux;
          ux = uy;
          uy = ut;
          s |= 2;
        }
      }
      else {
        ux = ~ux;
        uy = ~uy;
        if (uz&m) {
          s |= 4;
        }
        else {
          s |= 3;
        }
      }
    }
    else {
      if (uy&m) {
        if (uz&m) {
          ut = ux;
          ux = uy;
          uy = ~uz;
          uz = ~ut;
          s |= 6;
        }
        else {
          ut = uz;
          uz = ux;
          ux = uy;
          uy = ut;
          s |= 1;
        }
      }
      else {
        if (uz&m) {
          ut = uy;
          uy = ux;
          ux = ~uz;
          uz = ~ut;
          s |= 7;
        }
        else {
          ut = uy;
          uy = ux;
          ux = uz;
          uz = ut;
          s |= 0;
        }
      }
    }
    m = m >> 1;
  }
  return s;
}

/**************************/

// variable transformation from [min,max] to [1,2) and casting
// from double to float.
// This prepares our data for the morton key sorting routines
float vtf(double x, double min, double max) {
  return (float)((x - min) / (max - min + FLT_EPSILON*(max-min)) + 1.0);

}


void hilbert_sort_s(float* xp, int np, float* mind, float* maxd, int dim, int* info) {

  int i = 0;
  // allocate indices array 
  idx_t* indices = (idx_t*)malloc(sizeof(idx_t)*np);

  // iterate over all particle positions and compute their 
  // morton curve indices
  if (dim == 2) {
    // 2D case
    for(i = 0; i < np*2; i+=2) {
      indices[i/2].d = hilbert2d(
          vtf(xp[i],mind[0],maxd[0]),
          vtf(xp[i+1],mind[1],maxd[1]));
      indices[i/2].i = i;
      indices[i/2].x[0] = xp[i];
      indices[i/2].x[1] = xp[i+1];
    }
    // sort the indices according to the morton key
    qsort(indices,np,sizeof(idx_t),cmp);

    // now iterate through the indices and swap the xp positions
    for(i = 0; i < np; i++) {
      xp[2*i] = indices[i].x[0];
      xp[2*i+1] = indices[i].x[1];
    }

  } else {
    // 3D case
    for(i = 0; i < np*3; i+=3) {
      indices[i/3].d = hilbert3d(
          vtf(xp[i],mind[0],maxd[0]),
          vtf(xp[i+1],mind[1],maxd[1]),
          vtf(xp[i+2],mind[2],maxd[2]));
      indices[i/3].i = i;
      indices[i/3].x[0] = xp[i];
      indices[i/3].x[1] = xp[i+1];
      indices[i/3].x[2] = xp[i+2];
    }
    
    // sort the indices according to the morton key
    qsort(indices,np,sizeof(idx_t),cmp);
    
    // now iterate through the indices and swap the xp positions
    for(i = 0; i < np; i++) {
      xp[3*i] = indices[i].x[0];
      xp[3*i+1] = indices[i].x[1];
      xp[3*i+2] = indices[i].x[2];
    }
  }

  free(indices);

}

void hilbert_sort_d(double* xp, int np, double* mind, double* maxd, int dim, int* info) {

  int i = 0;
  // allocate indices array 
  idx_t* indices = (idx_t*)malloc(sizeof(idx_t)*np);

  // iterate over all particle positions and compute their 
  // morton curve indices
  if (dim == 2) {
    // 2D case
    for(i = 0; i < np*2; i+=2) {
      indices[i/2].d = hilbert2d(
          vtf(xp[i],mind[0],maxd[0]),
          vtf(xp[i+1],mind[1],maxd[1]));
      indices[i/2].i = i;
      indices[i/2].x[0] = xp[i];
      indices[i/2].x[1] = xp[i+1];
    }
    // sort the indices according to the morton key
    qsort(indices,np,sizeof(idx_t),cmp);

    // now iterate through the indices and swap the xp positions
    for(i = 0; i < np; i++) {
      xp[2*i] = indices[i].x[0];
      xp[2*i+1] = indices[i].x[1];
    }

  } else {
    // 3D case
    for(i = 0; i < np*3; i+=3) {
      indices[i/3].d = hilbert3d(
          vtf(xp[i],mind[0],maxd[0]),
          vtf(xp[i+1],mind[1],maxd[1]),
          vtf(xp[i+2],mind[2],maxd[2]));
      indices[i/3].i = i;
      indices[i/3].x[0] = xp[i];
      indices[i/3].x[1] = xp[i+1];
      indices[i/3].x[2] = xp[i+2];
    }
    
    // sort the indices according to the morton key
    qsort(indices,np,sizeof(idx_t),cmp);
    
    // now iterate through the indices and swap the xp positions
    for(i = 0; i < np; i++) {
      xp[3*i] = indices[i].x[0];
      xp[3*i+1] = indices[i].x[1];
      xp[3*i+2] = indices[i].x[2];
    }
  }

  free(indices);

}
