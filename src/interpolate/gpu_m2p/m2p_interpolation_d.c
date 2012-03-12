#include "m2p_interpolation_d.h"
#include <stdio.h>

#ifndef __STAMP
#define __STAMP
static char stamp[] = "***\nmodule '' __FILE__ ''\ncompiled '' __TIMESTAMP__ ''\n***";
#endif

#if __KIND == __SINGLE_PRECISION
typedef float    t_real;
typedef cl_float cl_t_real;
const char* m2p_option_s = "-D__SINGLE_PRECISION -cl-fast-relaxed-math -cl-single-precision-constant -cl-denorms-are-zero -cl-no-signed-zeros";
#elif __KIND == __DOUBLE_PRECISION
typedef double    t_real2;
typedef cl_double cl_t_real2;
const char* m2p_option_d = "-D__DOUBLE_PRECISION -cl-fast-relaxed-math -cl-single-precision-constant -cl-denorms-are-zero -cl-no-signed-zeros";
#endif

#if __KIND == __SINGLE_PRECISION
int ppm_gpu_m2p_initialize_s_
#elif __KIND == __DOUBLE_PRECISION
int ppm_gpu_m2p_initialize_d_
#endif
(
					int*   mesh_size,
          t_real2* minphys_with_ghost,
          t_real2* maxphys_with_ghost,
          int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
){
  cl_int	error 	  = CL_SUCCESS, error_buffer = CL_SUCCESS;
  int    	dim;
	int     different_mesh_size = 0;
	int    	nparticle = *np;
	int		ndim	  = *ndim_in;
	int 	nmass     = *nmass_in;
	char **filename = (char **)malloc(sizeof(char *));

    // Initialize device if not initialized before
	if(!GPU_initialized){
    	clInit(&platform, &device, &context, &queue, &program);
        GPU_initialized = 1;
	}
 
    // Read the specific kernel file if sizes or the kernel has changed
	if((ndim_prev   != ndim)  ||
	   (nmass_prev  != nmass) ||
	   (last_call   != m2p_interpolation) ||
       (last_kernel != *interp_kernel))
	{
		filename[0] = (char *)malloc(sizeof(char)*100);

	   	if(*interp_kernel == MP4_KERNEL){
    	    sprintf(filename[0], "src/interpolate/gpu_m2p/m2p_kernel_files/kernel_%dd_m2p_mp4_m%d.cl", ndim, nmass);
    	}
    	else if(*interp_kernel == BSP2_KERNEL){
        	sprintf(filename[0], "src/interpolate/gpu_m2p/m2p_kernel_files/kernel_%dd_m2p_bsp2_m%d.cl", ndim, nmass);
    	}

        // Build program from the new file
#if __KIND == __SINGLE_PRECISION
    	clBuild(filename, 1, &program, &context, &device, m2p_option_s);	
#elif __KIND == __DOUBLE_PRECISION
    	clBuild(filename, 1, &program, &context, &device, m2p_option_d);	
#endif

		free(filename[0]);
		free(filename);
	}

    // Check if mesh size is different
	for(dim = 0; dim < ndim; dim++){
		if(mesh_size_prev[dim] != mesh_size[dim])	different_mesh_size = 1;
	}

    // if the problem has changed, run a clean-up and reset block sizes
	if((ndim_prev   != ndim)  ||
	   (nmass_prev  != nmass) ||
	   (different_mesh_size) ||
	   (last_call   != m2p_interpolation) ||
       (last_kernel != *interp_kernel))
	   {

		ppm_gpu_clean_up();

		if(ndim == 2){
			for(dim = 0; dim < ndim; dim++){
				BLOCK_SIZE_interior[dim] = BLOCK_SIZE_interior_2d[dim];
			}
		} else if(ndim == 3){
			for(dim = 0; dim < ndim; dim++){
				BLOCK_SIZE_interior[dim] = BLOCK_SIZE_interior_3d[dim];
			}
		}

    	if(*interp_kernel == MP4_KERNEL){
        	for(dim = 0; dim < ndim; dim++){
            	kernel_size[dim] = kernel_size_m2p_mp4[dim];
        	}
    	}
    	else if(*interp_kernel == BSP2_KERNEL){
        	for(dim = 0; dim < ndim; dim++){
            	kernel_size[dim] = kernel_size_m2p_bsp2[dim];
        	}
    	}

        // Calculate number of cells and mesh nodes with and without ghost layers
		ndim_prev		  = ndim;
		ncell_int		  = 1;
		ncell_int_padded  = 1; 
		ncell_ext		  = 1; 
		ncell_ext_padded  = 1; 
    	nmesh_int 	      = 1; 
    	nmesh_int_padded  = 1; 
    	nmesh_ext 	      = 1; 
    	nmesh_ext_padded  = 1; 
    	for (dim = 0; dim < ndim; dim++){
			mesh_size_ext[dim] = mesh_size[dim];
			mesh_size_prev[dim] = mesh_size[dim];
			// compute exterior block sizes
			BLOCK_SIZE_exterior[dim] = BLOCK_SIZE_interior[dim] + 2*kernel_size[dim] + 1;
    		// compute h for domain without ghost layers
			mesh_size_int[dim] = mesh_size_ext[dim] - 2*kernel_size[dim];
        	h_d[dim]  = (maxphys_with_ghost[dim] - minphys_with_ghost[dim])/(mesh_size_ext[dim] - 1);
    		// extend input domain including ghost layers
       		maxphys_d[dim]  = maxphys_with_ghost[dim] - h_d[dim]*kernel_size[dim];
        	minphys_d[dim]  = minphys_with_ghost[dim] + h_d[dim]*kernel_size[dim];
			// compute number of cells in interior domain and exterior domain
			cell_size_int[dim] = mesh_size_int[dim] - 1;
        	cell_size_ext[dim] = cell_size_int[dim] + 2*kernel_size[dim];
			// Calculate group sizes for 2D. In 2D, the group size is 32x16
			num_groups[dim] = divide_roundUp(mesh_size_ext[dim], BLOCK_SIZE_interior[dim]);
			// Calculate mesh sizes when the mesh is padded
			mesh_size_ext_padded[dim] = num_groups[dim]*BLOCK_SIZE_exterior[dim];
			// Calculate cell list sizes for the padded exterior cell list
			cell_size_int_padded[dim] = num_groups[dim]*BLOCK_SIZE_interior[dim];
			// Calculate number of number of cells and mesh points
			ncell_int		  *= cell_size_int[dim];
			ncell_int_padded  *= cell_size_int_padded[dim];
			ncell_ext		  *= cell_size_ext[dim];
			ncell_ext_padded  *= cell_size_ext_padded[dim];
	    	nmesh_int 	      *= mesh_size_int[dim];
	    	nmesh_int_padded  *= mesh_size_int_padded[dim];
	    	nmesh_ext 	      *= mesh_size_ext[dim];
	    	nmesh_ext_padded  *= mesh_size_ext_padded[dim];
    	}

		// CREATE BUFFERS //
    	// Buffer 2: Mass values of particles
    	dev_mass     = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_t_real2)*nparticle*nmass, NULL, &error);
		error_buffer |= error;
    	// Buffer 3: Coordinates of minimum extent of physical domain
		dev_minphys = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real2)*ndim, minphys_d, &error);
		error_buffer |= error;
    	// Buffer 4: Coordinates of maximum extent of physical domain
		dev_maxphys = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real2)*ndim, maxphys_d, &error);
		error_buffer |= error;
    	// Buffer 5: h values
		dev_h = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real2)*ndim, h_d, &error);
		error_buffer |= error;
    	// Buffer 6: Interior cell sizes
		dev_cell_size_int = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, cell_size_int, &error);
		error_buffer |= error;
    	// Buffer 7: Exterior cell sizes
		dev_cell_size_ext = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, cell_size_ext, &error);
		error_buffer |= error;
    	// Buffer 8: Exterior cell sizes when padded
		dev_cell_size_ext_padded = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, cell_size_ext_padded, &error);
		error_buffer |= error;
    	// Buffer 9: Number of groups
		dev_num_groups = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, num_groups, &error);
		error_buffer |= error;
		// Buffer 10: Number of particles per cell
		dev_np_cell  = clCreateBuffer(context, CL_MEM_READ_WRITE,    sizeof(cl_int)*ncell_int_padded, NULL, &error);
		error_buffer |= error;
		// Buffer 11: Indices of cells that particles belong to
		dev_p_cell  = clCreateBuffer(context, CL_MEM_READ_WRITE,    sizeof(cl_int)*nparticle, NULL, &error);
		error_buffer |= error;
    	// Buffer 13: Mesh sizes
    	dev_mesh_size_int = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, mesh_size_int, &error);
		error_buffer |= error;
    	// Buffer 14: Mesh sizes, external
    	dev_mesh_size_ext = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, mesh_size_ext, &error);
		error_buffer |= error;
    	// Buffer 15: Mesh sizes when padded
    	dev_mesh_size_ext_padded = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, mesh_size_ext_padded, &error);
		error_buffer |= error;
    	// Buffer 17: Mesh points, padded
    	dev_mesh_padded = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_t_real2)*nmesh_ext_padded*nmass, NULL, &error);
		error_buffer |= error;

    	for(dim = 0; dim < ndim; dim++){
			// store the mesh size
			mesh_size_prev[dim] = mesh_size[dim];
		}

        // Store the dimension, number of properties and the kernel
		ndim_prev         = ndim;
		nmass_prev        = nmass;
		last_call 		  = m2p_interpolation;
		last_kernel 	  = *interp_kernel;
#if __KIND == __SINGLE_PRECISION
		ptr2clean_up      = &ppm_gpu_m2p_clean_up_s;
#elif __KIND == __DOUBLE_PRECISION
		ptr2clean_up      = &ppm_gpu_m2p_clean_up_d;
#endif
	}

	return error_buffer;
}

#if __KIND == __SINGLE_PRECISION
int ppm_gpu_m2p_interpolation_s_
#elif __KIND == __DOUBLE_PRECISION
int ppm_gpu_m2p_interpolation_d_
#endif
(
                    t_real2* xp,
                    t_real2* mass,
                    t_real2* mesh,
                    t_real2* minphys_with_ghost,
                    t_real2* maxphys_with_ghost,
					int*   mesh_size,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
) {
    cl_int  error;
    int     n, i, j, dim;
	int nparticle 			= *np;
	int	ndim	  			= *ndim_in;
	int nmass    		    = *nmass_in;
  int idim;
    struct timeval t, t_total;
    t_real2  t_elapsed, t_elapsed_total;
  // gpu timers
  cl_ulong start, end;
  float executionTimeInMilliseconds;

#if __KIND == __SINGLE_PRECISION
	error = ppm_gpu_m2p_initialize_s_(mesh_size, minphys_with_ghost, maxphys_with_ghost, np, ndim_in, nmass_in, interp_kernel);
#elif __KIND == __DOUBLE_PRECISION
	error = ppm_gpu_m2p_initialize_d_(mesh_size, minphys_with_ghost, maxphys_with_ghost, np, ndim_in, nmass_in, interp_kernel);
#endif
   	checkError(error, "clCreateBuffer while allocating dev_offset", __FUNCTION__, __LINE__, stamp);

  tic(&t_total);

    // Buffer 1: Positions of particles
    dev_xp_const = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_t_real2)*nparticle*ndim, NULL, &error);
    checkError(error, "clCreateBuffer while allocating devCoorX", __FUNCTION__, __LINE__, stamp);
    // Buffer 16: Mesh points to be copied onto device
    dev_mesh_const = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_t_real2)*nmesh_ext*nmass, NULL, &error);
    checkError(error, "clCreateBuffer while allocating dev_offset", __FUNCTION__, __LINE__, stamp);

    error = clEnqueueWriteBuffer(queue, dev_xp_const, CL_FALSE, 0, sizeof(cl_t_real2)*nparticle*ndim, xp, 0, NULL, &event_write_xp);
    checkError(error, "clCreateBuffer while allocating dev_offset", __FUNCTION__, __LINE__, stamp);
	error = clEnqueueWriteBuffer(queue, dev_mesh_const, CL_FALSE, 0, sizeof(cl_t_real2)*nmesh_ext*nmass, mesh, 0, NULL, &event_write_mesh);
    checkError(error, "clCreateBuffer while allocating dev_offset", __FUNCTION__, __LINE__, stamp);

    clFinish(queue);
  tic(&t);

// INITIALIZING TO ZERO //
	// Initialize np_cell array to 0
#if __KIND == __SINGLE_PRECISION
	error = m2p_initialize_np_cell_s();
#elif __KIND == __DOUBLE_PRECISION
	error = m2p_initialize_np_cell_d();
#endif
    checkError(error, "Error in m2p_initialize_np_cell", __FUNCTION__, __LINE__, stamp);

	// Initialize padded mesh array to 0
#if __KIND == __SINGLE_PRECISION
	error = initialize_padded_mesh_s(nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = initialize_padded_mesh_d(nmass);
#endif
    checkError(error, "Error in initialize_padded_mesh", __FUNCTION__, __LINE__, stamp);

// GET CELL INDICES OF PARTICLES AND COUNT NUMBER OF PARTICLES IN CELLS //
#if __KIND == __SINGLE_PRECISION
    error = m2p_get_particle_positions_s(nparticle, ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
    error = m2p_get_particle_positions_d(nparticle, ndim, nmass);
#endif
    checkError(error, "Error in m2p_get_particle_positions", __FUNCTION__, __LINE__, stamp);

// REDUCE NP_CELL TO GET MAXIMUM NUMBER OF PARTICLES PER CELL //
#if __KIND == __SINGLE_PRECISION
	error = m2p_reduce_np_cell_s();
#elif __KIND == __DOUBLE_PRECISION
	error = m2p_reduce_np_cell_d();
#endif
    checkError(error, "Error in m2p_reduce_np_cell", __FUNCTION__, __LINE__, stamp);

// COMPUTE SIZE OF PARTICLE DATA BUFFER AND ALLOCATE IT ON THE DEVICE
#if __KIND == __SINGLE_PRECISION
	error = m2p_allocate_data_buffers_s(ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = m2p_allocate_data_buffers_d(ndim, nmass);
#endif
    checkError(error, "Error in m2p_allocate_data_buffers", __FUNCTION__, __LINE__, stamp);

// DUPLICATE MESH
#if __KIND == __SINGLE_PRECISION
	error = duplicate_mesh_s(ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = duplicate_mesh_d(ndim, nmass);
#endif
    checkError(error, "Error in duplicate_mesh", __FUNCTION__, __LINE__, stamp);

// PLACE PARTICLES ACCORDINGLY //
#if __KIND == __SINGLE_PRECISION
	error = place_particle_pos_s(nparticle, ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = place_particle_pos_d(nparticle, ndim, nmass);
#endif
    checkError(error, "Error in place_particle_pos", __FUNCTION__, __LINE__, stamp);

// INTERPOLATE MESH POINTS ONTO PARTICLES //
#if __KIND == __SINGLE_PRECISION
	error = m2p_interpolate_s(ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = m2p_interpolate_d(ndim, nmass);
#endif
    checkError(error, "Error in m2p_interpolate", __FUNCTION__, __LINE__, stamp);

// COLLECT PARTICLES FROM PADDED PARTICLE ARRAY //
#if __KIND == __SINGLE_PRECISION
	error = collect_particles_s(nparticle, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = collect_particles_d(nparticle, nmass);
#endif
    checkError(error, "Error in collect_particles", __FUNCTION__, __LINE__, stamp);

   clFinish(queue);

    t_elapsed  = toc(&t);
    t_elapsed *= 1000.0;

    
    // get gpu timing for clist
    clGetEventProfilingInfo(event_get_np_cell, CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_get_np_cell, CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 
    executionTimeInMilliseconds = (end - start) * 1.0e-6f;
    printf("GPU M2P (double) clist kernel time: %lf ms\n", executionTimeInMilliseconds);
    // get gpu timing for part place
    executionTimeInMilliseconds = 0.0;
	  for(idim = 0; idim < ndim; idim++){
    clGetEventProfilingInfo(event_place_particle_pos[idim], CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_place_particle_pos[idim], CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 
    executionTimeInMilliseconds += (end - start) * 1.0e-6f;
    }
    printf("GPU M2P (double) part_place (%d dims) kernel time: %lf ms\n", ndim,executionTimeInMilliseconds);
    // get gpu timing for interpolation
    clGetEventProfilingInfo(event_m2p_interpolate, CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_m2p_interpolate, CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 

    executionTimeInMilliseconds = (end - start) * 1.0e-6f;
    printf("GPU M2P (double) interp kernel time: %lf ms\n", executionTimeInMilliseconds);
    // get gpu timing for part_collect
    clGetEventProfilingInfo(event_collect_particles, CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_collect_particles, CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 

    executionTimeInMilliseconds = (end - start) * 1.0e-6f;
    printf("GPU M2P (double) collect part kernel time: %lf ms\n", executionTimeInMilliseconds);


// READ BACK RESULTS //
	error = clEnqueueWaitForEvents(queue, 1, &event_collect_particles);

    error = clEnqueueReadBuffer(queue, dev_mass, CL_TRUE, 0, sizeof(t_real2)*nparticle*nmass, mass, 0, NULL, &event_read_mesh);
    checkError(error, "clEnqueueReadBuffer while reading back", __FUNCTION__, __LINE__, stamp);

    error |= clFinish(queue);

    t_elapsed_total  = toc(&t_total);
    t_elapsed_total *= 1000.0;

    printf("GPU Time consumtion: %lf ms including data transfer: %lf ms\n", t_elapsed, t_elapsed_total);

// CLEAN-UP INTERMEDIATE BUFFERS
	error |= clReleaseMemObject(dev_xp_const);
	error |= clReleaseMemObject(dev_mesh_const);
	error |= clReleaseMemObject(dev_sorted_particle_position);
	error |= clReleaseMemObject(dev_sorted_particle_mass);

    return 0;
}

#if __KIND == __SINGLE_PRECISION
int m2p_initialize_np_cell_s(){
#elif __KIND == __DOUBLE_PRECISION
int m2p_initialize_np_cell_d(){
#endif
	size_t 		localWorkSize, globalWorkSize, num_wg;
	int    		n, error;
	cl_kernel	kernel;

    kernel = clCreateKernel(program, "kernel_clist_init_to_zero_int", &error);
    checkError(error, "clCreateKernel 'clist_init_kernel' function", __FUNCTION__, __LINE__, stamp);
        
    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_np_cell);
	error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&ncell_int_padded);
	
    // Set work-group sizes
	localWorkSize  = 512;
	globalWorkSize = ncell_int_padded;
	if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the GPU
	error |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_np_cell);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int initialize_padded_mesh_s(int nmass)
#elif __KIND == __DOUBLE_PRECISION
int initialize_padded_mesh_d(int nmass)
#endif
{
	size_t 		localWorkSize, globalWorkSize, num_wg;
	int    		n, error;
	cl_kernel	kernel;
	int size_mesh_padded = ncell_int_padded*nmass;

    kernel = clCreateKernel(program, "kernel_clist_init_to_zero_real", &error);
    
    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_mesh_padded);
	error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&size_mesh_padded);
	
    // Set work-group sizes
	localWorkSize  = 512;
	globalWorkSize = size_mesh_padded;
	if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the GPU
	error |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_padded_mesh);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int m2p_get_particle_positions_s
#elif __KIND == __DOUBLE_PRECISION
int m2p_get_particle_positions_d
#endif
(
						int nparticle,
						int ndim,
						int nmass
){
	size_t 		localWorkSize, globalWorkSize, num_wg;
	int    		n, error;
	cl_kernel	kernel;

	kernel = clCreateKernel(program, "kernel_clist_get_particle_positions", &error);

    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_xp_const);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&nparticle);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&ndim);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&nmass);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&ncell_int_padded);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_h);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_minphys);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_num_groups); 
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_np_cell);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_p_cell);

    // Set work-group sizes
	localWorkSize  = 512;	
	globalWorkSize = nparticle;
	if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

	clEnqueueWaitForEvents(queue, 1, &event_init_np_cell);
	clEnqueueWaitForEvents(queue, 1, &event_write_xp);
    
    // Launch kernel on the GPU
    error |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_np_cell);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int m2p_reduce_np_cell_s()
#elif __KIND == __DOUBLE_PRECISION
int m2p_reduce_np_cell_d()
#endif
{
	int 	    n, error, error_func = CL_SUCCESS;
	size_t 	  localWorkSize, globalWorkSize, num_wg;
	cl_kernel ck_init_max_np, ck_init_max_np2, ck_init_max_np3, ck_get_max_np, ck_get_max_np2, ck_get_max_np3; 
	size_t 	  size_max_np = 0, size_max_np2 = 0, size_max_np3 = 0;

  // Buffer 19: Sorted particle position 
	size_max_np = (ncell_int_padded + 1023)/1024;
	dev_max_np  = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_int)*size_max_np, NULL, &error); 

  if(size_max_np > 1){
    // Buffer 20: Sorted particle position 
    size_max_np2 = (size_max_np + 1023)/1024;
    dev_max_np2  = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_int)*size_max_np2, NULL, &error); 
  }
	if(size_max_np2 > 1){
	  size_max_np3 = (size_max_np2 + 1023)/1024;
	  dev_max_np3  = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_int)*size_max_np3, NULL, &error);
  } 
// Initialize max_np to 0
  ck_init_max_np = clCreateKernel(program, "kernel_clist_init_to_zero_int", &error);
        
  // Set kernel arguments
	n = 0;
  error |= clSetKernelArg(ck_init_max_np, n++, sizeof(cl_mem), (void *)&dev_max_np);
	error |= clSetKernelArg(ck_init_max_np, n++, sizeof(int), (void *)&size_max_np);
	
  // Set work-group sizes
	localWorkSize  = 512;
	globalWorkSize = size_max_np;
	if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
  num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
  globalWorkSize = num_wg*localWorkSize;

  // Launch kernel on the GPU
	error |= clEnqueueNDRangeKernel(queue, ck_init_max_np, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_max_np);
	error_func |= error;
	
// Initialize max_np2 to 0
  if(size_max_np > 1){
    ck_init_max_np2 = clCreateKernel(program, "kernel_clist_init_to_zero_int", &error);
        
    // Set kernel arguments
	  n = 0;
    error |= clSetKernelArg(ck_init_max_np2, n++, sizeof(cl_mem), (void *)&dev_max_np2);
	  error |= clSetKernelArg(ck_init_max_np2, n++, sizeof(int), (void *)&size_max_np2);
	
    // Set work-group sizes
	  localWorkSize  = 512;
	  globalWorkSize = size_max_np2;
	  if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the GPU
	  error |= clEnqueueNDRangeKernel(queue, ck_init_max_np2, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_max_np2);
	  error_func |= error;
  }

// Initialize max_np3 to 0
  if(size_max_np2 > 1){
    ck_init_max_np3 = clCreateKernel(program, "kernel_clist_init_to_zero_int", &error);
        
    // Set kernel arguments
	  n = 0;
    error |= clSetKernelArg(ck_init_max_np3, n++, sizeof(cl_mem), (void *)&dev_max_np3);
	  error |= clSetKernelArg(ck_init_max_np3, n++, sizeof(int),    (void *)&size_max_np3);
	
    // Set work-group sizes
	  localWorkSize  = 512;
	  globalWorkSize = size_max_np3;
	  if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the GPU
	  error |= clEnqueueNDRangeKernel(queue, ck_init_max_np2, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_max_np3);
	  error_func |= error;
  }

// Reduce large blocks
	ck_get_max_np = clCreateKernel(program, "clist_get_max_np_blockwise", &error);

  // Set kernel arguments
	n = 0;
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(cl_mem), (void *)&dev_np_cell);
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(int), (void *)&ncell_int_padded);
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(cl_int)*1536, NULL);
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(cl_mem), (void *)&dev_max_np);

  // Set work-group sizes
	localWorkSize = 512;
	num_wg = size_max_np;
	globalWorkSize = localWorkSize*size_max_np;

	error |= clEnqueueWaitForEvents(queue, 1, &event_get_np_cell);
	error |= clEnqueueWaitForEvents(queue, 1, &event_init_max_np);

  // Launch kernel on the GPU
  error |= clEnqueueNDRangeKernel(queue, ck_get_max_np, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_max_np);
	error_func |= error;

// Reduce smaller blocks
  if(size_max_np > 1){
  	ck_get_max_np2 = clCreateKernel(program, "clist_get_max_np_blockwise", &error);

    // Set kernel arguments
	  n = 0;
    error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(cl_mem), (void *)&dev_max_np);
    error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(int), (void *)&size_max_np);
    error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(cl_int)*1536, NULL);
    error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(cl_mem), (void *)&dev_max_np2);

    // Set work-group sizes
	  localWorkSize = 512;
	  num_wg = size_max_np2;
	  globalWorkSize = localWorkSize*size_max_np2;

	  error |= clEnqueueWaitForEvents(queue, 1, &event_get_max_np);
	  error |= clEnqueueWaitForEvents(queue, 1, &event_init_max_np2);

    // Launch kernel on the GPU
    error |= clEnqueueNDRangeKernel(queue, ck_get_max_np2, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_max_np2);
  }
// Reduce even smaller blocks
  if(size_max_np2 > 1){
  	ck_get_max_np3 = clCreateKernel(program, "clist_get_max_np_blockwise", &error);

    // Set kernel arguments
	  n = 0;
    error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(cl_mem), (void *)&dev_max_np2);
    error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(int), (void *)&size_max_np2);
    error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(cl_int)*1536, NULL);
    error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(cl_mem), (void *)&dev_max_np3);

    // Set work-group sizes
	  localWorkSize = 512;
	  num_wg = size_max_np3;
	  globalWorkSize = localWorkSize*size_max_np3;

	  error |= clEnqueueWaitForEvents(queue, 1, &event_get_max_np2);
	  error |= clEnqueueWaitForEvents(queue, 1, &event_init_max_np3);

    // Launch kernel on the GPU
    error |= clEnqueueNDRangeKernel(queue, ck_get_max_np3, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_max_np3);
  }
// READ BACK MAXIMUM NUMBER OF PARTICLES
  if(size_max_np2 > 1){
	  error |= clEnqueueWaitForEvents(queue, 1, &event_get_max_np3);
    error |= clEnqueueReadBuffer(queue, dev_max_np3, CL_TRUE, 0, sizeof(int)*1, &max_np, 0, NULL, &event_read_max_np);
  } else if(size_max_np > 1){
	  error |= clEnqueueWaitForEvents(queue, 1, &event_get_max_np2);
    error |= clEnqueueReadBuffer(queue, dev_max_np2, CL_TRUE, 0, sizeof(int)*1, &max_np, 0, NULL, &event_read_max_np);
  } else{
	  error |= clEnqueueWaitForEvents(queue, 1, &event_get_max_np);
    error |= clEnqueueReadBuffer(queue, dev_max_np, CL_TRUE, 0, sizeof(int)*1, &max_np, 0, NULL, &event_read_max_np);
  }

	error |= clReleaseMemObject(dev_max_np);
  if(size_max_np2 > 1){
	  error |= clReleaseMemObject(dev_max_np3);
	}
  if(size_max_np > 1){
	  error |= clReleaseMemObject(dev_max_np2);
	}
  error_func |= error;

	return error_func;
}

#if __KIND == __SINGLE_PRECISION
int m2p_allocate_data_buffers_s
#elif __KIND == __DOUBLE_PRECISION
int m2p_allocate_data_buffers_d
#endif
(
					int ndim,
					int nmass
){
	int error = CL_SUCCESS, error_func = CL_SUCCESS;

	error |= clEnqueueWaitForEvents(queue, 1, &event_read_max_np);
	n_sorted_particle_position    = max_np*ndim*ncell_int_padded + 1536; // Extra padding in order to avoid if branch
	size_sorted_particle_position = sizeof(t_real2)*n_sorted_particle_position;
	
	n_sorted_particle_mass    = max_np*nmass*ncell_int_padded + 1536; // Extra padding in order to avoid if branch
	size_sorted_particle_mass = sizeof(t_real2)*n_sorted_particle_mass;

    // Buffer 21: Sorted particle position 
	dev_sorted_particle_position = clCreateBuffer(context, CL_MEM_READ_WRITE, size_sorted_particle_position, NULL, &error);
	error_func |= error;

    // Buffer 22: Sorted particle mass
	dev_sorted_particle_mass = clCreateBuffer(context, CL_MEM_READ_WRITE, size_sorted_particle_mass, NULL, &error);
	error_func |= error;

	return error_func;
}

#if __KIND == __SINGLE_PRECISION
int duplicate_mesh_s
#elif __KIND == __DOUBLE_PRECISION
int duplicate_mesh_d
#endif
(
				int ndim,
				int nmass
){
	size_t 		localWorkSizeND[3], globalWorkSizeND[3], num_wgND[3];
	size_t 		localWorkSize, globalWorkSize, num_wg;
	int    		n, dim, error;
	cl_kernel	kernel;

   	kernel = clCreateKernel(program, "kernel_clist_duplicate_mesh", &error);

    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_mesh_const);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_mesh_padded);
	error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_mesh_size_ext);
	error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_num_groups);
	error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&nmesh_ext_padded);
	error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&nmass);

	for(dim = 0; dim < ndim; dim++){
		localWorkSizeND[dim] = BLOCK_SIZE_exterior[dim];
		globalWorkSizeND[dim] = num_groups[dim]*localWorkSizeND[dim];
	}

	globalWorkSize = 1; 
	for(dim = 0; dim < ndim; dim++){
		globalWorkSize *= globalWorkSizeND[dim];
	}
	
    // Set work-group sizes
	localWorkSize = 512;
	num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
	globalWorkSize = num_wg*localWorkSize;
	
	error |= clEnqueueWaitForEvents(queue, 1, &event_write_mesh);

    // Launch kernel on the GPU
	error |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_duplicate_mesh);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int place_particle_pos_s
#elif __KIND == __DOUBLE_PRECISION
int place_particle_pos_d
#endif
(
					int nparticle,
					int ndim,
					int nmass
){
	size_t 		localWorkSize, globalWorkSize, num_wg;
	int    		n, dim, error = CL_SUCCESS, error_func = CL_SUCCESS;
	cl_kernel	kernel[3];
	// Place particles' positions, one position in each iteration.
	// This is better as we use caches more efficiently and 
	// reduce the drawback of uncoalesced writes to global memory.
	
	error |= clEnqueueWaitForEvents(queue, 1, &event_write_xp);
	error |= clEnqueueWaitForEvents(queue, 1, &event_get_np_cell);

	for(dim = 0; dim < ndim; dim++){
		kernel[dim] = clCreateKernel(program, "kernel_clist_place_particle_pos", &error);
		error_func |= error;

        // Set kernel arguments
		n = 0;
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(cl_mem), (void *)&dev_xp_const);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(cl_mem), (void *)&dev_mass);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(int), (void *)&nparticle);
    	error |= clSetKernelArg(kernel[dim], n++, sizeof(int), (void *)&ndim);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(int), (void *)&nmass);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(int), (void *)&dim);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(int), (void *)&ncell_int_padded);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(int), (void *)&max_np);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(cl_mem), (void *)&dev_p_cell);
	    error |= clSetKernelArg(kernel[dim], n++, sizeof(cl_mem), (void *)&dev_sorted_particle_position);

        // Set work-group sizes
		localWorkSize  = 512;	
		globalWorkSize = nparticle;
		if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
	    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
	    globalWorkSize = num_wg*localWorkSize;

        // Launch kernel on the GPU
	    error |= clEnqueueNDRangeKernel(queue, kernel[dim], 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_place_particle_pos[dim]);
		
		error_func |= error;
	}

	return error_func;
}

#if __KIND == __SINGLE_PRECISION
int m2p_interpolate_s
#elif __KIND == __DOUBLE_PRECISION
int m2p_interpolate_d
#endif
(
				int ndim,
				int nmass
){
	size_t 		localWorkSizeND[3], globalWorkSizeND[3], num_wgND[3];
	int    		n, dim, error;
	cl_kernel	kernel;

   	kernel = clCreateKernel(program, "kernel_interpolate", &error);
        
    // Set kernel arguments
    n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_sorted_particle_position);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_sorted_particle_mass);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_mesh_padded);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_t_real2)*(1715*nmass), NULL);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&ndim);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&nmass);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&max_np);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&ncell_int_padded);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&nmesh_ext_padded);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_minphys);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_h);
    
    for(dim = 0; dim < ndim; dim++) {
		localWorkSizeND[dim]  = BLOCK_SIZE_interior[dim];
       	globalWorkSizeND[dim] = num_groups[dim]*localWorkSizeND[dim];
   	}    

	error |= clEnqueueWaitForEvents(queue, 1, &event_duplicate_mesh);
	error |= clEnqueueWaitForEvents(queue, 1, &event_write_mesh);
    for(dim = 0; dim < ndim; dim++) {
		error |= clEnqueueWaitForEvents(queue, 1, &event_place_particle_pos[dim]);
	}

    // Launch kernel on the GPU
	error |= clEnqueueNDRangeKernel(queue, kernel, ndim, NULL, globalWorkSizeND, localWorkSizeND, 0, NULL, &event_m2p_interpolate);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int collect_particles_s
#elif __KIND == __DOUBLE_PRECISION
int collect_particles_d
#endif
(
					int nparticle,
					int nmass
){
	size_t 		localWorkSize, globalWorkSize, num_wg;
	int    		n, error = CL_SUCCESS;
	cl_kernel	kernel;

   	kernel = clCreateKernel(program, "kernel_collect_particles", &error);
 
    // Set kernel arguments
    n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_sorted_particle_mass);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_mass);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_p_cell);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&nmass);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&max_np);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&ncell_int_padded);
    error |= clSetKernelArg(kernel, n++, sizeof(cl_int), (void *)&nparticle);
    
    localWorkSize  = 512;
    num_wg = (nparticle + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

	error |= clEnqueueWaitForEvents(queue, 1, &event_m2p_interpolate);

    // Launch kernel on the GPU
	error |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_collect_particles);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int ppm_gpu_m2p_clean_up_s(){
#elif __KIND == __DOUBLE_PRECISION
int ppm_gpu_m2p_clean_up_d(){
#endif
	int error = CL_SUCCESS;

    // Release memory objects
	error |= clReleaseMemObject(dev_mass);
	error |= clReleaseMemObject(dev_minphys);
	error |= clReleaseMemObject(dev_maxphys);
	error |= clReleaseMemObject(dev_h);
	error |= clReleaseMemObject(dev_cell_size_int);
	error |= clReleaseMemObject(dev_cell_size_ext);
	error |= clReleaseMemObject(dev_cell_size_ext_padded);
	error |= clReleaseMemObject(dev_num_groups);
	error |= clReleaseMemObject(dev_np_cell);
	error |= clReleaseMemObject(dev_p_cell);
	error |= clReleaseMemObject(dev_mesh_size_int);
	error |= clReleaseMemObject(dev_mesh_size_ext);
	error |= clReleaseMemObject(dev_mesh_size_ext_padded);
	error |= clReleaseMemObject(dev_mesh_padded);

	return error;
}

