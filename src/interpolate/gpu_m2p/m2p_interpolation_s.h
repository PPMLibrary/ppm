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

int ppm_gpu_m2p_initialize_s_(
					int*   mesh_size,
                    float* minphys_with_ghost,
                    float* maxphys_with_ghost,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
);

int ppm_gpu_m2p_interpolate_s_(
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

int m2p_initialize_np_cell_s();

int initialize_padded_mesh_s(int nmass);

int m2p_get_particle_positions_s
(
						int nparticle,
						int ndim,
						int nmass
);

int m2p_reduce_np_cell_s();

int m2p_allocate_data_buffers_s
(
					int ndim,
					int nmass
);

int duplicate_mesh_s
(
				int ndim,
				int nmass
);

int place_particle_pos_s
(
					int nparticle,
					int ndim,
					int nmass
);

int m2p_interpolate_s
(
				int ndim,
				int nmass
);

int collect_particles_s
(
					int nparticle,
					int nmass
);

int ppm_gpu_m2p_clean_up_s();
