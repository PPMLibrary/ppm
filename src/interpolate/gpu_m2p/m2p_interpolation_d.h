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

int ppm_gpu_m2p_initialize_d_(
					int*   mesh_size,
                    double* minphys_with_ghost,
                    double* maxphys_with_ghost,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
);

int ppm_gpu_m2p_interpolate_d_(
                    double* xp,
                    double* mass,
                    double* mesh,
                    double* minphys_with_ghost,
                    double* maxphys_with_ghost,
					int*   mesh_size,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
); 

int m2p_initialize_np_cell_d();

int initialize_padded_mesh_d(int nmass);

int m2p_get_particle_positions_d
(
						int nparticle,
						int ndim,
						int nmass
);

int m2p_reduce_np_cell_d();

int m2p_allocate_data_buffers_d
(
					int ndim,
					int nmass
);

int duplicate_mesh_d
(
				int ndim,
				int nmass
);

int place_particle_pos_d
(
					int nparticle,
					int ndim,
					int nmass
);

int m2p_interpolate_d
(
				int ndim,
				int nmass
);

int collect_particles_d
(
					int nparticle,
					int nmass
);

int ppm_gpu_m2p_clean_up_d();
