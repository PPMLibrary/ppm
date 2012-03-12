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

cl_platform_id   platform;
cl_device_id     device;
cl_context       context;
cl_command_queue queue;
cl_program       program;

int ppm_gpu_p2m_initialize_(
                    int*   mesh_size,
                    float* minphys_wo_ghost,
                    float* maxphys_wo_ghost,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
);

int ppm_gpu_p2m_interpolation_(
                    float* xp,
                    float* mass,
                    float* mesh,
                    float* minphys_wo_ghost,
                    float* maxphys_wo_ghost,
                    int*   mesh_size,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel);

int p2m_initialize_np_cell();

int	p2m_get_particle_positions(
						int nparticle, 
						int ndim, 
						int nmass 
);

int p2m_reduce_np_cell();

int	p2m_allocate_data_buffers
(
				int ndim,
				int nmass
);

int	init_sorted_particle_mass();

int init_padded_mesh(int nmass);

int place_particle_positions(
						int nparticle,
						int ndim,
						int nmass
);

int	place_particle_mass(
						int nparticle,
						int ndim,
						int nmass
);

int	p2m_interpolate
(
					int ndim,
					int nmass
);

int	sew_padded_mesh(
					int *mesh_size,
					int  ndim
);

int ppm_gpu_p2m_clean_up();
