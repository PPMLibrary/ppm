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
#include "p2m_common.h"

void createRandomParticlePosition_2d(
                    int   p_per_cell,
                    float *physmin,
                    float *physmax,
                    int   *mesh_size,
                    float *xp,
					int   ndim
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
                    float *mass,
					int    nmass
);

void createRandomParticleMass_3d(
			        int   p_per_cell, 
					int   *mesh_size,
                    float *mass,
					int    nmass
);

void createRandomMesh_2d_m2p(
					int   *mesh_size,
					float *mesh,
					int	   nmass
);

void createRandomMesh_3d_m2p(
					int   *mesh_size,
					float *mesh,
					int	   nmass
);

void check_interpolation_2d_p2m_mp4_
(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   grid_dim,
                             int*   np,
							 int*   ndim_in,
							 int*	nmass_in
);

void check_interpolation_3d_p2m_mp4_
(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   grid_dim,
                             int*   np,
							 int*   ndim_in,
							 int*	nmass_in
);

void check_interpolation_2d_m2p_mp4
(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*	nmass_in
);

void check_interpolation_3d_m2p_mp4
(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*	nmass_in
);

void check_interpolation_2d_p2m_bsp2_
(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin_wo_ghost,
                             float* physmax_wo_ghost,
                             int*   mesh_size,
                             int*   np,
                             int*   ndim_in,
                             int*   nmass_in
);

void check_interpolation_3d_p2m_bsp2_
(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin_wo_ghost,
                             float* physmax_wo_ghost,
                             int*   mesh_size,
                             int*   np,
                             int*   ndim_in,
                             int*   nmass_in
							
);
