#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable

#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 8
#define BLOCK_SIZE_Z 8

#ifdef __SINGLE_PRECISION
typedef float  t_real;
#elif  __DOUBLE_PRECISION
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double t_real;
#endif

__kernel void kernel_clist_init_to_zero_int(
				__global int* buffer,
				int 		  nelements
) {
	int id = get_global_id(0);

	if(id < nelements){
		buffer[id] = 0;	
	}
}

__kernel void kernel_clist_init_to_zero_real(
				__global t_real* buffer,
				int 		    ncell_re
) {
	int id = get_global_id(0);

    if(id < ncell_re){
		buffer[id] = 0.0;	
	}
}

__kernel void kernel_clist_get_particle_positions
(
				const __global t_real* xp,
				int 				  nparticle,
				int					  ndim,
				int					  nval,
				int					  ncell_re,
				const __global t_real* h,
				const __global t_real* physmin,
				const __global int*   num_groups,
				__global int* 		  np_cell,
				__global int* 		  p_cell
) {
	int p_id, cell_ID, image_group, index, i;
	int cell_pos_x, cell_pos_y, cell_pos_z, val;
	int wg_x, wg_y, wg_z;

	p_id = get_global_id(0); 
	
	if(p_id < nparticle) {
        cell_pos_x = (xp[3*p_id]     - physmin[0])/h[0]; 
        cell_pos_y = (xp[3*p_id + 1] - physmin[1])/h[1]; 
        cell_pos_z = (xp[3*p_id + 2] - physmin[2])/h[2]; 

		wg_x = cell_pos_x/BLOCK_SIZE_X;
		wg_y = cell_pos_y/BLOCK_SIZE_Y;
		wg_z = cell_pos_z/BLOCK_SIZE_Z;

		cell_pos_x = cell_pos_x - wg_x*BLOCK_SIZE_X;
		cell_pos_y = cell_pos_y - wg_y*BLOCK_SIZE_Y;
		cell_pos_z = cell_pos_z - wg_z*BLOCK_SIZE_Z;
		cell_ID    = cell_pos_x + BLOCK_SIZE_X*cell_pos_y + BLOCK_SIZE_X*BLOCK_SIZE_Y*cell_pos_z + (wg_x + num_groups[0]*wg_y + num_groups[0]*num_groups[1]*wg_z)*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z; 

		image_group = atom_inc(&np_cell[cell_ID]);    //number of particles in cell_i is incremented
		index = image_group*ncell_re + cell_ID;

		p_cell[p_id] = index;
	}
}

__kernel void kernel_clist_place_particle_pos
(
				const __global t_real* xp,
				const __global t_real* mass,
				int 				  nparticle,
				int					  ndim,
				int					  nval,
				int					  rank,
				int					  ncell_re,
				int					  max_np,
				__global int* 		  p_cell,
				__global t_real*       sorted_particle_position
) {
	int p_id, cell_ID, image_group, index, i;
	t_real pos;

	p_id = get_global_id(0); 
	
	if(p_id < nparticle) {
		index = p_cell[p_id];

		pos   = xp[ndim*p_id + rank];
		sorted_particle_position[index + rank*max_np*ncell_re] = pos;
	}
}

__kernel void kernel_clist_place_particle_mass
(
				const __global t_real* xp,
				const __global t_real* mass,
				int 				  nparticle,
				int					  ndim,
				int					  nval,
				int					  rank,
				int					  ncell_re,
				int					  max_np,
				__global int* 		  p_cell,
				__global t_real*       sorted_particle_mass
) {
	int p_id, cell_ID, image_group, index, i;
	t_real pos_x, pos_y, val;

	p_id = get_global_id(0); 
	
	if(p_id < nparticle) {
		index = p_cell[p_id];

		val   = mass[nval*p_id + rank];
		sorted_particle_mass[index + rank*max_np*ncell_re] = val;
	}
}

__kernel void clist_get_max_np_blockwise(
										__global int* np_cell,
										int			  ncell_ext,
										__local  int* max_np,
										__global int* g_max_np
){
	int l_id, g_id, wg_id, s;

	l_id  = get_local_id(0); 
	g_id  = get_group_id(0)*get_local_size(0)*2 + l_id;

	// Initialize local memory to 0
	max_np[l_id] = 0;
	barrier(CLK_LOCAL_MEM_FENCE);	

	// Reading into local memory 
    max_np[l_id] = (g_id < ncell_ext) ? np_cell[g_id] : 0;
	barrier(CLK_LOCAL_MEM_FENCE);	

	// First reduction in global mem level
    if (g_id + get_local_size(0) < ncell_ext){
        max_np[l_id] = (max_np[l_id] < np_cell[g_id + get_local_size(0)]) ? np_cell[g_id + get_local_size(0)] : max_np[l_id];  
	}
    barrier(CLK_LOCAL_MEM_FENCE);

    // Remaining reductions done in shared memory
    for(s=get_local_size(0)/2; s>0; s>>=1){
        if (l_id < s)
        {
            max_np[l_id] = (max_np[l_id] < max_np[l_id + s]) ? max_np[l_id + s] : max_np[l_id];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write result back to global mem in the index of workgroup
    if (l_id == 0) 
		g_max_np[get_group_id(0)] = max_np[0];
}

__kernel void kernel_interpolate(
		  __global t_real* particle_pos_data,
		  __global t_real* particle_mass_data,
	      __local  t_real* local_mesh,
	               int    ndim,
	               int    nmass,
	               int    max_np,
	               int    ncell_ext_padded,
	               int    nmesh_ext_padded,
	const __global t_real* physmin,
	const __global t_real* h,
	const __global int*   mesh_size,
		  __global t_real* mesh_re
){
	int   i, rank, mass;
	int  m_id_x, m_id_y, m_id_z;
	t_real x, x0, x1, x2, wx_x, wx_y, wx_z;
	t_real pos_x, pos_y, pos_z;
	t_real p_mass[4];

	// Calculate local ID and workgroup ID
	int l_ID   = get_local_id(0) + get_local_id(1)*BLOCK_SIZE_X + get_local_id(2)*BLOCK_SIZE_X*BLOCK_SIZE_Y;
	int wg_ID  = get_group_id(0) + get_group_id(1)*get_num_groups(0) + get_group_id(2)*get_num_groups(0)*get_num_groups(1);

	// Initialize local memory to 0. Barrier is just before the for-loop.
	for(i = l_ID; i < 4*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3); i += BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z){
		local_mesh[i] = 0.0;
	}

	for(rank = 0; rank < max_np; rank++){
		// read particle coordinates and its mass from global memory in a coalesced fashion.
		pos_x  = particle_pos_data[             rank*ncell_ext_padded + wg_ID*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z + l_ID];
		pos_y  = particle_pos_data[(rank +   max_np)*ncell_ext_padded + wg_ID*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z + l_ID];
		pos_z  = particle_pos_data[(rank + 2*max_np)*ncell_ext_padded + wg_ID*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z + l_ID];

		p_mass[0] = particle_mass_data[		        rank*ncell_ext_padded + wg_ID*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z + l_ID];
		p_mass[1] = particle_mass_data[(rank +   max_np)*ncell_ext_padded + wg_ID*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z + l_ID];
		p_mass[2] = particle_mass_data[(rank + 2*max_np)*ncell_ext_padded + wg_ID*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z + l_ID];
		p_mass[3] = particle_mass_data[(rank + 3*max_np)*ncell_ext_padded + wg_ID*BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z + l_ID];

		// Calculate local ID of the mesh point.
		l_ID = get_local_id(0) + get_local_id(1)*(BLOCK_SIZE_X+3) + get_local_id(2)*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3);

		// Make sure local memory is initialized to 0
		barrier(CLK_LOCAL_MEM_FENCE);

		// In for loops, same particle is used to update 4x4x4 mesh points around it. 
		x2 = (pos_z - physmin[2])/h[2];
		x2 = x2 - floor(x2);
		x1 = (pos_y - physmin[1])/h[1];
		x1 = x1 - floor(x1);
		x0 = (pos_x - physmin[0])/h[0];
		x0 = x0 - floor(x0);
		for(m_id_z = -1; m_id_z < 3; m_id_z++){
			x = fabs(x2 - m_id_z);

			if(x < 1.0)
				wx_z = 1.0 - x*x*(2.5 - 1.5*x);
			else
				wx_z = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;
			
// m_id_y = -1
				x = 1.0 + x1;

				if(x < 1.0)
					wx_y = 1.0 - x*x*(2.5 - 1.5*x);
				else 
					wx_y = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

// m_id_x = -1
					x = 1.0 + x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 0
					x = x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 1
					x = 1.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 2
					x = 2.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
				l_ID += BLOCK_SIZE_X-1; // Move on to next mesh point in y-dimension
			
// m_id_y = 0
				x = x1;

				if(x < 1.0)
					wx_y = 1.0 - x*x*(2.5 - 1.5*x);
				else 
					wx_y = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

// m_id_x = -1
					x = 1.0 + x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 0
					x = x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 1
					x = 1.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 2
					x = 2.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
				l_ID += BLOCK_SIZE_X-1; // Move on to next mesh point in y-dimension
			
// m_id_y = 1
				x = 1.0 - x1;

				if(x < 1.0)
					wx_y = 1.0 - x*x*(2.5 - 1.5*x);
				else 
					wx_y = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

// m_id_x = -1
					x = 1.0 + x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 0
					x = x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 1
					x = 1.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 2
					x = 2.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
				l_ID += BLOCK_SIZE_X-1; // Move on to next mesh point in y-dimension
			
// m_id_y = 2
				x = 2.0 - x1;

				if(x < 1.0)
					wx_y = 1.0 - x*x*(2.5 - 1.5*x);
				else 
					wx_y = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

// m_id_x = -1
					x = 1.0 + x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 0
					x = x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 1
					x = 1.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
// m_id_x = 2
					x = 2.0 - x0;

					if(x < 1.0)
						wx_x = 1.0 - x*x*(2.5 - 1.5*x);
					else 
						wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

					// Update mesh value in the local memory
					local_mesh[l_ID] 		+= wx_x*wx_y*wx_z*p_mass[0];
					local_mesh[l_ID + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[1];
					local_mesh[l_ID + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[2];
					local_mesh[l_ID + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)] += wx_x*wx_y*wx_z*p_mass[3];
					barrier(CLK_LOCAL_MEM_FENCE);

					l_ID += 1; // Move on to next mesh point in x-dimension
				
				l_ID += BLOCK_SIZE_X-1; // Move on to next mesh point in y-dimension
			
			l_ID += (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y-1); // Move on to next mesh point in z-dimension
		}
			
		// Write mesh point values in the local memory onto mesh points in the 
		// padded array in the global memory in a coalesced fashion.
		l_ID   = get_local_id(0) + get_local_id(1)*BLOCK_SIZE_X + get_local_id(2)*BLOCK_SIZE_X*BLOCK_SIZE_Y;
	}

	// Coalesced write
	for(i = l_ID; i < (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3); i+=BLOCK_SIZE_X*BLOCK_SIZE_Y*BLOCK_SIZE_Z){
		mesh_re[				     wg_ID*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3) + i] = local_mesh[i];
		mesh_re[nmesh_ext_padded   + wg_ID*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3) + i] = local_mesh[i + (BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)];
		mesh_re[2*nmesh_ext_padded + wg_ID*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3) + i] = local_mesh[i + 2*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)];
		mesh_re[3*nmesh_ext_padded + wg_ID*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3) + i] = local_mesh[i + 3*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3)];
	}
}

void get_mesh_ID_wg_padded(
	__global int* mesh_size_padded,
	__global int* num_groups,
	int  global_id_x,
	int  global_id_y,
	int  global_id_z,
	int *mesh_ID_padded,
	int *wg_id_x,
	int *wg_id_y,
	int *wg_id_z
) {
	int wg_ID, local_id_x, local_id_y, local_id_z, local_ID;

	*wg_id_x = global_id_x/BLOCK_SIZE_X; 
	*wg_id_y = global_id_y/BLOCK_SIZE_Y;
	*wg_id_z = global_id_z/BLOCK_SIZE_Z;
	wg_ID    = (*wg_id_x) + (*wg_id_y)*num_groups[0] + (*wg_id_z)*num_groups[0]*num_groups[1];
	
	local_id_x = global_id_x  - (*wg_id_x)*BLOCK_SIZE_X; 
	local_id_y = global_id_y  - (*wg_id_y)*BLOCK_SIZE_Y; 
	local_id_z = global_id_z  - (*wg_id_z)*BLOCK_SIZE_Z; 

	global_id_x = 3 + local_id_x + (*wg_id_x)*(BLOCK_SIZE_X+3);
	global_id_y = 3 + local_id_y + (*wg_id_y)*(BLOCK_SIZE_Y+3);
	global_id_z = 3 + local_id_z + (*wg_id_z)*(BLOCK_SIZE_Z+3);

	local_id_x = global_id_x  - (*wg_id_x)*(BLOCK_SIZE_X+3); 
	local_id_y = global_id_y  - (*wg_id_y)*(BLOCK_SIZE_Y+3); 
	local_id_z = global_id_z  - (*wg_id_z)*(BLOCK_SIZE_Z+3); 
	local_ID   = local_id_x + local_id_y*(BLOCK_SIZE_X+3) + local_id_z*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3);
	*mesh_ID_padded = wg_ID*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3) + local_ID;
}

void get_mesh_ID_to_sew(
	__global int* mesh_size_padded,
	__global int* num_groups,
	int  global_id_x,
	int  global_id_y,
	int  global_id_z,
	int  wg_id_x_prev,
	int  wg_id_y_prev,
	int  wg_id_z_prev,
	int *mesh_ID_padded,
	int *wg_id_x,
	int *wg_id_y,
	int *wg_id_z
) {
	int wg_ID, local_id_x, local_id_y, local_id_z, local_ID;

	*wg_id_x = global_id_x/BLOCK_SIZE_X; 
	*wg_id_y = global_id_y/BLOCK_SIZE_Y;
	*wg_id_z = global_id_z/BLOCK_SIZE_Z;
	wg_ID    = (*wg_id_x) + (*wg_id_y)*num_groups[0] + (*wg_id_z)*num_groups[0]*num_groups[1];

	global_id_x += (wg_id_x_prev - *wg_id_x + 1)*3;
	global_id_y += (wg_id_y_prev - *wg_id_y + 1)*3;
	global_id_z += (wg_id_z_prev - *wg_id_z + 1)*3;
	
	local_id_x = global_id_x  - (*wg_id_x)*BLOCK_SIZE_X; 
	local_id_y = global_id_y  - (*wg_id_y)*BLOCK_SIZE_Y; 
	local_id_z = global_id_z  - (*wg_id_z)*BLOCK_SIZE_Z; 
	local_ID   = local_id_x + local_id_y*(BLOCK_SIZE_X+3) + local_id_z*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3);

	*mesh_ID_padded = wg_ID*(BLOCK_SIZE_X+3)*(BLOCK_SIZE_Y+3)*(BLOCK_SIZE_Z+3) + local_ID;
}
__kernel void kernel_sew_padded_mesh
(
		const __global t_real* padded_mesh,
			  __global t_real* mesh,
			  __global int*   mesh_size,
			  __global int*   mesh_size_padded,
					   int	  nmesh_ext_padded,
			  __global int*   num_groups
) {
	int mesh_id_x, mesh_id_y, mesh_id_z, mesh_ID, mesh_ID_padded;
	int wg_id_x, wg_id_y, wg_id_z, wg_ID;
	int wg_id_x_padded, wg_id_y_padded, wg_id_z_padded;
	int local_id_x, local_id_y, local_id_z, local_ID;
	int mesh_ID_padded_temp;

	t_real total_value[4];

	mesh_id_x = get_global_id(0);
	mesh_id_y = get_global_id(1);
	mesh_id_z = get_global_id(2);
	mesh_ID   = mesh_id_x + mesh_id_y*mesh_size[0] + mesh_id_z*mesh_size[0]*mesh_size[1];

	if((mesh_id_x < mesh_size[0]) && (mesh_id_y < mesh_size[1]) && (mesh_id_z < mesh_size[2])) { 
		get_mesh_ID_wg_padded(mesh_size_padded, num_groups, mesh_id_x, mesh_id_y, mesh_id_z, &mesh_ID_padded, &wg_id_x, &wg_id_y, &wg_id_z);
		total_value[0] = padded_mesh[mesh_ID_padded];
		total_value[1] = padded_mesh[mesh_ID_padded +   nmesh_ext_padded];
		total_value[2] = padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
		total_value[3] = padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];

		// (1,0,0) direction
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x + 3, mesh_id_y, mesh_id_z, wg_id_x, wg_id_y, wg_id_z, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded, &wg_id_z_padded);
		if(wg_id_x != wg_id_x_padded){	
			total_value[0] += padded_mesh[mesh_ID_padded];
			total_value[1] += padded_mesh[mesh_ID_padded +   nmesh_ext_padded];
			total_value[2] += padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
			total_value[3] += padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];
		}

		// (0,1,0) direction
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x, mesh_id_y + 3, mesh_id_z, wg_id_x, wg_id_y, wg_id_z, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded, &wg_id_z_padded);
		if(wg_id_y != wg_id_y_padded){
			total_value[0] += padded_mesh[mesh_ID_padded];
			total_value[1] += padded_mesh[mesh_ID_padded + nmesh_ext_padded];
			total_value[2] += padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
			total_value[3] += padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];
		}

		// (0,0,1) direction
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x, mesh_id_y, mesh_id_z + 3, wg_id_x, wg_id_y, wg_id_z, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded, &wg_id_z_padded);
		if(wg_id_z != wg_id_z_padded){	
			total_value[0] += padded_mesh[mesh_ID_padded];
			total_value[1] += padded_mesh[mesh_ID_padded + nmesh_ext_padded];
			total_value[2] += padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
			total_value[3] += padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];
		}

		// (1,1,0) direction
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x + 3, mesh_id_y + 3, mesh_id_z, wg_id_x, wg_id_y, wg_id_z, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded, &wg_id_z_padded);
		if((wg_id_x != wg_id_x_padded) && (wg_id_y != wg_id_y_padded)){
			total_value[0] += padded_mesh[mesh_ID_padded];
			total_value[1] += padded_mesh[mesh_ID_padded + nmesh_ext_padded];
			total_value[2] += padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
			total_value[3] += padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];
		}

		// (1,0,1) direction
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x + 3, mesh_id_y, mesh_id_z + 3, wg_id_x, wg_id_y, wg_id_z, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded, &wg_id_z_padded);
		if((wg_id_x != wg_id_x_padded) && (wg_id_z != wg_id_z_padded)){
			total_value[0] += padded_mesh[mesh_ID_padded];
			total_value[1] += padded_mesh[mesh_ID_padded + nmesh_ext_padded];
			total_value[2] += padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
			total_value[3] += padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];
		}

		// (0,1,1) direction
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x, mesh_id_y + 3, mesh_id_z + 3, wg_id_x, wg_id_y, wg_id_z, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded, &wg_id_z_padded);
		if((wg_id_y != wg_id_y_padded) && (wg_id_z != wg_id_z_padded)){
			total_value[0] += padded_mesh[mesh_ID_padded];
			total_value[1] += padded_mesh[mesh_ID_padded + nmesh_ext_padded];
			total_value[2] += padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
			total_value[3] += padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];
		}

		// (1,1,1) direction
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x + 3, mesh_id_y + 3, mesh_id_z + 3, wg_id_x, wg_id_y, wg_id_z, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded, &wg_id_z_padded);
		if((wg_id_x != wg_id_x_padded) && (wg_id_y != wg_id_y_padded) && (wg_id_z != wg_id_z_padded)){
			total_value[0] += padded_mesh[mesh_ID_padded];
			total_value[1] += padded_mesh[mesh_ID_padded + nmesh_ext_padded];
			total_value[2] += padded_mesh[mesh_ID_padded + 2*nmesh_ext_padded];
			total_value[3] += padded_mesh[mesh_ID_padded + 3*nmesh_ext_padded];
		}
			
		mesh[4*mesh_ID]     = total_value[0];
		mesh[4*mesh_ID + 1] = total_value[1];
		mesh[4*mesh_ID + 2] = total_value[2];
		mesh[4*mesh_ID + 3] = total_value[3];
	
	}
}
