#include "interp_common.h"

int         last_kernel = -99;
int			MP4_KERNEL  = -1;
int			BSP2_KERNEL = -2;

int			ndim_prev  = 0;
int			nmass_prev = 0;
int         mesh_size_prev[3] = {0, 0, 0};

int         ncell_int;
int         ncell_ext;
int         ncell_ext_padded;
int         ncell_int_padded;

int         cell_size_int[3];
int         cell_size_ext[3];
int         cell_size_int_padded[3];
int         cell_size_ext_padded[3];

int         num_groups[3];

int	 		nmesh;
int			nmesh_int;
int			nmesh_ext;
int         nmesh_padded;
int         nmesh_int_padded;
int         nmesh_ext_padded;

int         mesh_size_int[3];
int         mesh_size_ext[3];
int         mesh_size_padded[3];
int         mesh_size_int_padded[3];
int         mesh_size_ext_padded[3];

int         n_sorted_particle_position;
int         n_sorted_particle_mass;
size_t      size_sorted_particle_position;
size_t      size_sorted_particle_mass;

int         max_np;

int			kernel_size[3];
int         kernel_size_p2m_mp4[3]  = {2, 2, 2};
int         kernel_size_p2m_bsp2[3] = {1, 1, 1};
int         kernel_size_m2p_mp4[3]  = {1, 1 ,1};
int         kernel_size_m2p_bsp2[3] = {0, 0, 0};

int     blocksize[2];


int 		BLOCK_SIZE_interior[3];
int         BLOCK_SIZE_exterior[3]; 

#if __GPU == NV_TESLA
int         BLOCK_SIZE_interior_2d[3] = {32, 16, 0};
int         BLOCK_SIZE_interior_3d[3] = {32, 4, 4};
#elif __GPU == ATI_CAYMAN
int         BLOCK_SIZE_interior_2d[3] = {16, 16, 0};
int         BLOCK_SIZE_interior_3d[3] = {8, 8, 4};
#elif __GPU == AMD_CPU
int         BLOCK_SIZE_interior_2d[3] = {32, 32, 0};
int         BLOCK_SIZE_interior_3d[3] = {16, 8, 8};
#else
int         BLOCK_SIZE_interior_2d[3] = {16, 16, 0};
int         BLOCK_SIZE_interior_3d[3] = {8, 8, 4};
#endif

size_t		globalWorkSize;
size_t    	localWorkSize;
size_t 		num_wg;

size_t		globalWorkSizeND[3];
size_t    	localWorkSizeND[3];
size_t 		num_wgND[3];

float	  	h_s[3];
float       minphys_s[3];
float       maxphys_s[3];

double	  	h_d[3];
double      minphys_d[3];
double      maxphys_d[3];

int         size_max_np  = 0;
int         size_max_np2 = 0;
int         size_max_np3 = 0;

int			interp_lmem_size;

