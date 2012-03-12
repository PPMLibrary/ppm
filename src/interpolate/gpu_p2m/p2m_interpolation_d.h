/////////////////////////////////////////////////////
// File: p2m_interpolation.h                       //
// p2m interpolation on the GPU                    //
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
#include "../common_vars/interp_common.h"
#include "../common_vars/interp_common_ocl.h"

int ppm_gpu_p2m_initialize_d_(
                    int*   mesh_size,
                    double* minphys_wo_ghost,
                    double* maxphys_wo_ghost,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
);

int ppm_gpu_p2m_interpolation_d_(
                    double* xp,
                    double* mass,
                    double* mesh,
                    double* minphys_wo_ghost,
                    double* maxphys_wo_ghost,
                    int*   mesh_size,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel);

int p2m_initialize_np_cell_d();

int	p2m_get_particle_positions_d(
						int nparticle, 
						int ndim, 
						int nmass 
);

int p2m_reduce_np_cell_d();

int	p2m_allocate_data_buffers_d(
				int ndim,
				int nmass
);

int	init_sorted_particle_mass_d();

int init_padded_mesh_d(int nmass);

int place_particle_positions_d(
						int nparticle,
						int ndim,
						int nmass
);

int	place_particle_mass_d(
						int nparticle,
						int ndim,
						int nmass
);

int	p2m_interpolate_d(
					int ndim,
					int nmass
);

int	sew_padded_mesh_d(
					int *mesh_size,
					int  ndim
);

int ppm_gpu_p2m_clean_up_d();
