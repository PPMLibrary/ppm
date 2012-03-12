/////////////////////////////////////////////////////
// File: checkers.h                                //
// Test driver functions for GPU implementations   //
// Author, date: Ferit Buyukkececi, 26.03.2011     //
// todo:                                           //
//      1. Preprocessor directives for timings     //
//      2. Releasing memory objects in the end     //
//                                                 //
/////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../opencl_utils.h"
#include "m2p_common.h"

void createRandomParticlePosition_2d(
			        int   p_per_cell, 
			        float *physmin, 
	       			float *physmax,
					int   *mesh_size,
		    		float *xp,
					int    ndim
);

void createRandomParticlePosition_3d(
			        int   p_per_cell, 
			        float *physmin, 
	       			float *physmax,
					int   *mesh_size,
		    		float *xp,
					int    ndim
);


void createRandomParticleMass_2d(
			        int   p_per_cell, 
					int   *mesh_size,
                    float *mass
);
void createRandomParticleMass_3d(
			        int   p_per_cell, 
					int   *mesh_size,
                    float *mass
);

void createRandomMesh_2d_m2p(
					int   *mesh_size,
					int    nmass,
					float *mesh
);

void createRandomMesh_3d_m2p(
					int   *mesh_size,
					int    nmass,
					float *mesh
);

void check_interpolation_2d_m2p_mp4(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*   nmass_in
);

void check_interpolation_3d_m2p_mp4(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*   nmass_in
);

void check_interpolation_2d_m2p_bsp2(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
                             int*   ndim_in,
                             int*   nmass_in
);

void check_interpolation_3d_m2p_bsp2(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
                             int*   ndim_in,
                             int*   nmass_in
);
