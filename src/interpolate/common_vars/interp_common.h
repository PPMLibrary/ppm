#ifndef INTERP_COMMON_H
#define INTERP_COMMON_H

#include "../../gpu/opencl_utils.h"
#include "../../gpu/state.h"
#include <stdio.h>

extern int         last_kernel;
extern int			MP4_KERNEL;
extern int			BSP2_KERNEL;

extern int			ndim_prev;
extern int			nmass_prev;
extern int         mesh_size_prev[3];

extern int         ncell_int;
extern int         ncell_ext;
extern int         ncell_ext_padded;
extern int         ncell_int_padded;

extern int         cell_size_int[3];
extern int         cell_size_ext[3];
extern int         cell_size_int_padded[3];
extern int         cell_size_ext_padded[3];

extern int         num_groups[3];

extern int	 		nmesh;
extern int			nmesh_int;
extern int			nmesh_ext;
extern int         nmesh_padded;
extern int         nmesh_int_padded;
extern int         nmesh_ext_padded;

extern int         mesh_size_int[3];
extern int         mesh_size_ext[3];
extern int         mesh_size_padded[3];
extern int         mesh_size_int_padded[3];
extern int         mesh_size_ext_padded[3];

extern int         n_sorted_particle_position;
extern int         n_sorted_particle_mass;
extern size_t      size_sorted_particle_position;
extern size_t      size_sorted_particle_mass;

extern int         kernel_size_m2p_mp4[3];
extern int         kernel_size_m2p_bsp2[3];

extern int         max_np;

extern int			kernel_size[3];
extern int         kernel_size_p2m_mp4[3];
extern int         kernel_size_p2m_bsp2[3];
extern int         kernel_size_m2p_mp4[3];
extern int         kernel_size_m2p_bsp2[3];

extern int     blocksize[2];

extern int 		BLOCK_SIZE_interior[3];
extern int         BLOCK_SIZE_exterior[3]; 
extern int         BLOCK_SIZE_interior_2d[3];
extern int         BLOCK_SIZE_interior_3d[3];

extern size_t		globalWorkSize;
extern size_t    	localWorkSize;
extern size_t 		num_wg;

extern size_t		globalWorkSizeND[3];
extern size_t    	localWorkSizeND[3];
extern size_t 		num_wgND[3];

extern float		h_s[3];
extern float        minphys_s[3];
extern float        maxphys_s[3];

extern double		h_d[3];
extern double       minphys_d[3];
extern double       maxphys_d[3];

extern int         size_max_np;
extern int         size_max_np2;
extern int         size_max_np3;

extern int			interp_lmem_size;

#endif
