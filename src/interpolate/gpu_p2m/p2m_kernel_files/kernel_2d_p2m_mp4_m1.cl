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
	int cell_pos_x, cell_pos_y, val;
	int wg_x, wg_y;

	p_id = get_global_id(0); 
	
	if(p_id < nparticle) {
        cell_pos_x = (xp[2*p_id] - physmin[0])/h[0]; 
        cell_pos_y = (xp[2*p_id + 1] - physmin[1])/h[1]; 

		wg_x = cell_pos_x/32;
		wg_y = cell_pos_y/16;

		cell_pos_x = cell_pos_x - wg_x*32;
		cell_pos_y = cell_pos_y - wg_y*16;
		cell_ID    = cell_pos_x + 32*cell_pos_y + (wg_x + num_groups[0]*wg_y)*512; 

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
				__global t_real*       sorted_particle_pos
) {
	int p_id, cell_ID, image_group, index, i;
	t_real pos;

	p_id = get_global_id(0); 
	
	if(p_id < nparticle) {
		index = p_cell[p_id];

		pos   = xp[ndim*p_id + rank];
		sorted_particle_pos[index + rank*max_np*ncell_re] = pos;
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
	int   i, rank;
	int2  m_id;
	t_real x, x0, x1, wx_x, wx_y;
	t_real pos_x, pos_y;
	t_real p_mass;
	t_real value[4];

	int l_ID   = get_local_id(0) + get_local_id(1)*32;
	int wg_ID  = get_group_id(0) + get_group_id(1)*get_num_groups(0);

	for(i = l_ID; i < 665; i+=512){
		local_mesh[i] = 0.0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for(rank = 0; rank < max_np; rank++){
		pos_x  = particle_pos_data[             rank*ncell_ext_padded + wg_ID*512 + l_ID];
		pos_y  = particle_pos_data[(rank +   max_np)*ncell_ext_padded + wg_ID*512 + l_ID];

		p_mass = particle_mass_data[rank*ncell_ext_padded + wg_ID*512 + l_ID];

		l_ID = get_local_id(0) + get_local_id(1)*35;

		x1 = (pos_y - physmin[1])/h[1];
		x1 = x1 - floor(x1);
		x0 = (pos_x - physmin[0])/h[0];
		x0 = x0 - floor(x0);
		for(m_id.y = -1; m_id.y < 3; m_id.y++){
			x = fabs(x1 - m_id.y);

			if(x < 1.0)
				wx_y = 1.0 - x*x*(2.5 - 1.5*x);
			else 
				wx_y = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

			for(m_id.x = -1; m_id.x < 3; m_id.x++){
				x = fabs(x0 - m_id.x);

				if(x < 1.0)
					wx_x = 1.0 - x*x*(2.5 - 1.5*x);
				else 
					wx_x = 2.0 + (-4.0 + (2.5 - 0.5*x)*x)*x;

				local_mesh[l_ID] += wx_x*wx_y*p_mass;
				barrier(CLK_LOCAL_MEM_FENCE);
				l_ID += 1;
			}
			l_ID += 31;
		}
		
		l_ID   = get_local_id(0) + get_local_id(1)*32;
	}

	for(i = l_ID; i < 665; i+=512){
		mesh_re[wg_ID*665 + i] = local_mesh[i];
	}
}

void get_mesh_ID_wg_padded(
	__global int* mesh_size_padded,
	__global int* num_groups,
	int  global_id_x,
	int  global_id_y,
	int *mesh_ID_padded,
	int *wg_id_x,
	int *wg_id_y
) {
	int wg_ID, local_id_x, local_id_y, local_ID;

	*wg_id_x = global_id_x/32; 
	*wg_id_y = global_id_y/16;
	wg_ID    = (*wg_id_x) + (*wg_id_y)*num_groups[0];
	
	local_id_x = global_id_x  - (*wg_id_x)*32; 
	local_id_y = global_id_y  - (*wg_id_y)*16; 

	global_id_x = 3 + local_id_x + (*wg_id_x)*35;
	global_id_y = 3 + local_id_y + (*wg_id_y)*19;

	local_id_x = global_id_x  - (*wg_id_x)*35; 
	local_id_y = global_id_y  - (*wg_id_y)*19; 
	local_ID   = local_id_x + local_id_y*35;
	*mesh_ID_padded = wg_ID*665 + local_ID;
}

void get_mesh_ID_to_sew(
	__global int* mesh_size_padded,
	__global int* num_groups,
	int  global_id_x,
	int  global_id_y,
	int  wg_id_x_prev,
	int  wg_id_y_prev,
	int *mesh_ID_padded,
	int *wg_id_x,
	int *wg_id_y
) {
	int wg_ID, local_id_x, local_id_y, local_ID;

	*wg_id_x = global_id_x/32; 
	*wg_id_y = global_id_y/16;
	wg_ID    = (*wg_id_x) + (*wg_id_y)*num_groups[0];

	global_id_x += (wg_id_x_prev - *wg_id_x + 1)*3;
	global_id_y += (wg_id_y_prev - *wg_id_y + 1)*3;
	
	local_id_x = global_id_x  - (*wg_id_x)*32; 
	local_id_y = global_id_y  - (*wg_id_y)*16; 
	local_ID   = local_id_x + local_id_y*35;

	*mesh_ID_padded = wg_ID*665 + local_ID;
}
__kernel void kernel_sew_padded_mesh
(
		const __global t_real* padded_mesh,
			  __global t_real* mesh,
			  __global int*   mesh_size,
			  __global int*   mesh_size_padded,
					   int    nmesh_ext_padded,
			  __global int*   num_groups
) {
	int mesh_id_x, mesh_id_y, mesh_ID, mesh_ID_padded;
	int wg_id_x, wg_id_y, wg_ID;
	int wg_id_x_padded, wg_id_y_padded;
	int local_id_x, local_id_y, local_ID;
	int mesh_ID_padded_temp;

	mesh_id_x = get_global_id(0);
	mesh_id_y = get_global_id(1);
	mesh_ID = mesh_id_x + mesh_id_y*mesh_size[0];

	if((mesh_id_x < mesh_size[0]) && (mesh_id_y < mesh_size[1])) {
		get_mesh_ID_wg_padded(mesh_size_padded, num_groups, mesh_id_x, mesh_id_y, &mesh_ID_padded, &wg_id_x, &wg_id_y);
		mesh[mesh_ID] = padded_mesh[mesh_ID_padded];

		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x + 3, mesh_id_y, wg_id_x, wg_id_y, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded);
		if(wg_id_x != wg_id_x_padded){
			mesh[mesh_ID] += padded_mesh[mesh_ID_padded];
		}
	
		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x, mesh_id_y + 3, wg_id_x, wg_id_y, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded);
		if(wg_id_y != wg_id_y_padded){
			mesh[mesh_ID] += padded_mesh[mesh_ID_padded];
		}

		get_mesh_ID_to_sew(mesh_size_padded, num_groups, mesh_id_x + 3, mesh_id_y + 3, wg_id_x, wg_id_y, &mesh_ID_padded, &wg_id_x_padded, &wg_id_y_padded);
		if((wg_id_x != wg_id_x_padded) && (wg_id_y != wg_id_y_padded)){
			mesh[mesh_ID] += padded_mesh[mesh_ID_padded];
		}
	}
}
