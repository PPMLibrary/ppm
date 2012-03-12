////////////////////////////////////////////////////
// File: p2m_interpolation.h                       //
// p2m interpolation on the GPU                    //
// Author, date: Ferit Buyukkececi, 26.03.2011     //
// todo:                                           //
//      1. Preprocessor directives for timings     //
//      2. Releasing memory objects in the end     //

/////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../common_vars/interp_common.h"
#include "../common_vars/interp_common_ocl.h"

int ppm_gpu_m2p_initialize_(
					int*   mesh_size,
                    float* minphys_with_ghost,
                    float* maxphys_with_ghost,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
);

int ppm_gpu_m2p_interpolate_(
                    float* xp,
                    float* mass,
                    float* mesh,
                    float* minphys_with_ghost,
                    float* maxphys_with_ghost,
					int*   mesh_size,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
); 

int m2p_initialize_np_cell();

int initialize_padded_mesh(int nmass);

int m2p_get_particle_positions
(
						int nparticle,
						int ndim,
						int nmass
);

int m2p_reduce_np_cell();

int m2p_allocate_data_buffers
(
					int ndim,
					int nmass
);

int duplicate_mesh
(
				int ndim,
				int nmass
);

int place_particle_pos
(
					int nparticle,
					int ndim,
					int nmass
);

int m2p_interpolate
(
				int ndim,
				int nmass
);

int collect_particles
(
					int nparticle,
					int nmass
);

int ppm_gpu_m2p_clean_up();
