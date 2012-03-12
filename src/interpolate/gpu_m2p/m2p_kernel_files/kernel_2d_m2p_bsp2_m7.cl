
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable

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
				int 		    ncell
) {
	int id = get_global_id(0);

	if(id < ncell)
		buffer[id] = 0.0;	
}

__kernel void kernel_clist_duplicate_mesh
(
				__global t_real* mesh,
			    __global t_real* mesh_padded,
				__global int*   mesh_size,
				__global int*   num_groups,
							   int    nmesh_padded,
							   int    nmass
) {
	int t_id[2], wg_id[2], local_id[2];
	int mesh_ID, mesh_ID_padded;

	int g_ID, l_ID, wg_ID;	

	g_ID  = get_global_id(0);
	wg_ID = g_ID/561;	
	l_ID  = g_ID - 561*wg_ID;

	wg_id[1] = wg_ID/num_groups[0];
	wg_ID -= wg_id[1]*num_groups[0];
	wg_id[0] = wg_ID;

	local_id[1] = l_ID/33;
	l_ID -= local_id[1]*33;
	local_id[0] = l_ID;

	t_id[0] = wg_id[0]*32 + local_id[0];
	t_id[1] = wg_id[1]*16 + local_id[1];

	if((t_id[0] < mesh_size[0]) && (t_id[1] < mesh_size[1])){
		mesh_ID        = t_id[0] + mesh_size[0]*t_id[1];
		mesh_ID_padded = (wg_id[0] + wg_id[1]*num_groups[0])*561 + local_id[0] + local_id[1]*33;

		mesh_padded[mesh_ID_padded]                  = mesh[7*mesh_ID];
		mesh_padded[mesh_ID_padded +   nmesh_padded] = mesh[7*mesh_ID + 1];
		mesh_padded[mesh_ID_padded + 2*nmesh_padded] = mesh[7*mesh_ID + 2];
		mesh_padded[mesh_ID_padded + 3*nmesh_padded] = mesh[7*mesh_ID + 3];
		mesh_padded[mesh_ID_padded + 4*nmesh_padded] = mesh[7*mesh_ID + 4];
		mesh_padded[mesh_ID_padded + 5*nmesh_padded] = mesh[7*mesh_ID + 5];
		mesh_padded[mesh_ID_padded + 6*nmesh_padded] = mesh[7*mesh_ID + 6];
	}
}

__kernel void kernel_clist_get_particle_positions
(
				const __global t_real* xp,
				int 				  nparticle,
				int					  ndim,
				int					  nval,
				int					  ncell,
				const __global t_real* h,
				const __global t_real* physmin,
				const __global int*   num_groups,
				__global int* 		  np_cell,
				__global int* 		  p_cell
) {
	int p_id, cell_ID, image_group, index, i;
	int cell_pos_x, cell_pos_y, val;
	int wg_x, wg_y;

	p_id = get_global_id(0); 
	
	if(p_id < nparticle) {
        cell_pos_x = (int)((xp[2*p_id] - physmin[0])/h[0]); 
        cell_pos_y = (int)((xp[2*p_id + 1] - physmin[1])/h[1]); 

		wg_x = (int)(cell_pos_x/32);
		wg_y = (int)(cell_pos_y/16);

		cell_pos_x = cell_pos_x - wg_x*32;
		cell_pos_y = cell_pos_y - wg_y*16;
		cell_ID    = cell_pos_x + 32*cell_pos_y + (wg_x + num_groups[0]*wg_y)*512; 

		image_group = atom_inc(&np_cell[cell_ID]);    //number of particles in cell_i is incremented
		index = image_group*ncell + cell_ID;

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
				int					  ncell,
				int					  max_np,
				__global int* 		  p_cell,
				__global t_real*       sorted_particle_data
) {
	int p_id, cell_ID, image_group, index, i;
	t_real pos;

	p_id = get_global_id(0); 
	
	if(p_id < nparticle) {
		index = p_cell[p_id];

		pos   = xp[ndim*p_id + rank];
		sorted_particle_data[index + rank*max_np*ncell] = pos;
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

__kernel void kernel_interpolate
(
		  __global t_real* particle_pos_data,
		  __global t_real* particle_mass_data,
		  __global t_real* mesh_padded,
	      __local  t_real* local_mesh,
	               int    ndim,
	               int    nmass,
	               int    max_np,
	               int    ncell_int_padded,
	               int    nmesh_ext_padded,
	const __global t_real* physmin,
	const __global t_real* h
){
	int   		i, rank;
	int2  		m_id;
	t_real 		x, wx_x, wx_y;
	t_real 		pos_x, pos_y;
	t_real       particle_mass[7];

	int l_ID   = get_local_id(0) + get_local_id(1)*32;
	int wg_ID  = get_group_id(0) + get_group_id(1)*get_num_groups(0);

	for(i = l_ID; i < 561; i+= 512){
		local_mesh[i]        = mesh_padded[wg_ID*561 + i];
		local_mesh[i + 561]  = mesh_padded[wg_ID*561 + i +   nmesh_ext_padded];
		local_mesh[i + 1122] = mesh_padded[wg_ID*561 + i + 2*nmesh_ext_padded];
		local_mesh[i + 1683] = mesh_padded[wg_ID*561 + i + 3*nmesh_ext_padded];
		local_mesh[i + 2244] = mesh_padded[wg_ID*561 + i + 4*nmesh_ext_padded];
		local_mesh[i + 2805] = mesh_padded[wg_ID*561 + i + 5*nmesh_ext_padded];
		local_mesh[i + 3366] = mesh_padded[wg_ID*561 + i + 6*nmesh_ext_padded];
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for(rank = 0; rank < max_np; rank++){
		pos_x  = particle_pos_data[			rank*ncell_int_padded + wg_ID*512 + l_ID];
		pos_y  = particle_pos_data[(rank+max_np)*ncell_int_padded + wg_ID*512 + l_ID];

		particle_mass[0] = 0.0;
		particle_mass[1] = 0.0;
		particle_mass[2] = 0.0;
		particle_mass[3] = 0.0;
		particle_mass[4] = 0.0;
		particle_mass[5] = 0.0;
		particle_mass[6] = 0.0;

		l_ID = get_local_id(0) + get_local_id(1)*33;

		for(m_id.y = 0; m_id.y < 2; m_id.y++){
			x = (pos_y - physmin[1])/h[1];
			x = x - (int)x;
       	    wx_y = 1.0 - fabs(x - (t_real)m_id.y);
		
			for(m_id.x = 0; m_id.x < 2; m_id.x++){	
				x = (pos_x - physmin[0])/h[0];
				x = x - (int)x;
       	    	wx_x = 1.0 - fabs(x - (t_real)m_id.x);

				particle_mass[0] += wx_x*wx_y*local_mesh[l_ID];
				particle_mass[1] += wx_x*wx_y*local_mesh[l_ID + 561];
				particle_mass[2] += wx_x*wx_y*local_mesh[l_ID + 1122];
				particle_mass[3] += wx_x*wx_y*local_mesh[l_ID + 1683];
				particle_mass[4] += wx_x*wx_y*local_mesh[l_ID + 2244];
				particle_mass[5] += wx_x*wx_y*local_mesh[l_ID + 2805];
				particle_mass[6] += wx_x*wx_y*local_mesh[l_ID + 3366];

				l_ID += 1;
			}
			l_ID += 31;
		}

		l_ID   = get_local_id(0) + get_local_id(1)*32;

		particle_mass_data[				rank*ncell_int_padded + wg_ID*512 + l_ID] = particle_mass[0]; 
		particle_mass_data[(rank +   max_np)*ncell_int_padded + wg_ID*512 + l_ID] = particle_mass[1]; 
		particle_mass_data[(rank + 2*max_np)*ncell_int_padded + wg_ID*512 + l_ID] = particle_mass[2]; 
		particle_mass_data[(rank + 3*max_np)*ncell_int_padded + wg_ID*512 + l_ID] = particle_mass[3]; 
		particle_mass_data[(rank + 4*max_np)*ncell_int_padded + wg_ID*512 + l_ID] = particle_mass[4]; 
		particle_mass_data[(rank + 5*max_np)*ncell_int_padded + wg_ID*512 + l_ID] = particle_mass[5]; 
		particle_mass_data[(rank + 6*max_np)*ncell_int_padded + wg_ID*512 + l_ID] = particle_mass[6]; 
	}
}

__kernel void kernel_collect_particles
(
				__global t_real* padded_particle_mass,
				__global t_real* particle_mass,
				__global int*   p_cell,
						 int    nmass,
						 int    max_np,
						 int    ncell_int_padded,
						 int    nparticle
){
	int mass, cell_ID;
	int p_id = get_global_id(0);

	if(p_id < nparticle){
		cell_ID = p_cell[p_id];
		
		particle_mass[7*p_id]     = padded_particle_mass[cell_ID];
		particle_mass[7*p_id + 1] = padded_particle_mass[cell_ID +   max_np*ncell_int_padded];
		particle_mass[7*p_id + 2] = padded_particle_mass[cell_ID + 2*max_np*ncell_int_padded];
		particle_mass[7*p_id + 3] = padded_particle_mass[cell_ID + 3*max_np*ncell_int_padded];
		particle_mass[7*p_id + 4] = padded_particle_mass[cell_ID + 4*max_np*ncell_int_padded];
		particle_mass[7*p_id + 5] = padded_particle_mass[cell_ID + 5*max_np*ncell_int_padded];
		particle_mass[7*p_id + 6] = padded_particle_mass[cell_ID + 6*max_np*ncell_int_padded];
	}
}

