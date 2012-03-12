#include "p2m_interpolation_s.h"
#include <stdio.h>

#ifndef __STAMP
#define __STAMP
static char stamp[] = "***\nmodule '' __FILE__ ''\ncompiled '' __TIMESTAMP__ ''\n***";
#endif

#if __KIND == __SINGLE_PRECISION
typedef float    t_real;
typedef cl_float cl_t_real;
const char* p2m_option_s = "-D__SINGLE_PRECISION -cl-fast-relaxed-math -cl-single-precision-constant -cl-denorms-are-zero -cl-no-signed-zeros";
#elif __KIND == __DOUBLE_PRECISION
typedef double    t_real2;
typedef cl_double cl_t_real2;
const char* p2m_option_d = "-D__DOUBLE_PRECISION -cl-fast-relaxed-math -cl-single-precision-constant -cl-denorms-are-zero -cl-no-signed-zeros";
#endif

#if __KIND == __SINGLE_PRECISION
int ppm_gpu_p2m_initialize_s_
#elif __KIND == __DOUBLE_PRECISION
int ppm_gpu_p2m_initialize_d_
#endif
(
                    int*   mesh_size,
                    t_real* minphys_wo_ghost,
                    t_real* maxphys_wo_ghost,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
){
	int dim;
	int error;
	int different_mesh_size = 0;
	int nparticle   = *np;
	int ndim        = *ndim_in;
	int nmass       = *nmass_in;
	char **filename = (char **)malloc(sizeof(char *));

    // Initialize GPU if it has never been initialized before
	if(!GPU_initialized){ 
    	clInit(&platform, &device, &context, &queue, &program);
		GPU_initialized = 1;
	}

    // Read the file if any of #dimension, #properties or kernel
	// is different than the last call.
	// NOTE: ndim_prev and nmass_prev are initialized to 0 from the 
	// beginning, so the first call will satisfy the if-condition.
	if((ndim_prev   != ndim)  ||
	   (nmass_prev  != nmass) ||
	   (last_call   != p2m_interpolation) ||
       (last_kernel != *interp_kernel))
	   {
		filename[0] = (char *)malloc(sizeof(char)*100);

		if(*interp_kernel == MP4_KERNEL){
  			sprintf(filename[0], "src/interpolate/gpu_p2m/p2m_kernel_files/kernel_%dd_p2m_mp4_m%d.cl", ndim, nmass);
		}
		else if(*interp_kernel == BSP2_KERNEL){
  			sprintf(filename[0], "src/interpolate/gpu_p2m/p2m_kernel_files/kernel_%dd_p2m_bsp2_m%d.cl", ndim, nmass);
		}

        // Build program from the file
#if __KIND == __SINGLE_PRECISION
    	clBuild(filename, 1, &program, &context, &device, p2m_option_s);	
#elif __KIND == __DOUBLE_PRECISION
    	clBuild(filename, 1, &program, &context, &device, p2m_option_d);	
#endif

		free(filename[0]);
		free(filename);
	}

    // Store mesh sizes in order to detect necessity for reallocation in next calls
	for(dim = 0; dim < ndim; dim++){
		if(mesh_size_prev[dim] != mesh_size[dim])	different_mesh_size = 1;
	}


    // If the task is not the same, clean-up and re-build
	if((ndim_prev   != ndim)  ||
	   (nmass_prev  != nmass) ||
	   (different_mesh_size) ||
	   (last_call   != p2m_interpolation) ||
       (last_kernel != *interp_kernel))
	   {

		ppm_gpu_clean_up();

        // Set interior block sizes depending on the dimension
		if(ndim == 2){
	   	 	for(dim = 0; dim < ndim; dim++){
				BLOCK_SIZE_interior[dim] = BLOCK_SIZE_interior_2d[dim];
			}
		}
		else if(ndim == 3){
    		for(dim = 0; dim < ndim; dim++){
				BLOCK_SIZE_interior[dim] = BLOCK_SIZE_interior_3d[dim];
			}
		}

        // Set kernel sizes depending on the choice of interpolation scheme
		if(*interp_kernel == MP4_KERNEL){
    		for(dim = 0; dim < ndim; dim++){
				kernel_size[dim] = kernel_size_p2m_mp4[dim];
			}
		}
		else if(*interp_kernel == BSP2_KERNEL){
    		for(dim = 0; dim < ndim; dim++){
				kernel_size[dim] = kernel_size_p2m_bsp2[dim];
			}
		}

        // Compute number of cells in domains with and without ghost layers
		interp_lmem_size  = nmass;
		ncell_ext		  = 1;
		ncell_int		  = 1;
		ncell_ext_padded  = 1;
    	nmesh 			  = 1;
    	nmesh_padded 	  = 1;
    	for(dim = 0; dim < ndim; dim++){
    		// compute exterior block size 
			BLOCK_SIZE_exterior[dim] = BLOCK_SIZE_interior[dim] + 2*kernel_size[dim] - 1;
			// compute local memory size
			interp_lmem_size *= BLOCK_SIZE_exterior[dim];
    		// compute h for domain without ghost layers
        	h_s[dim]  = (maxphys_wo_ghost[dim] - minphys_wo_ghost[dim])/(mesh_size[dim] - 1);
    		// extend input domain including ghost layers
        	maxphys_s[dim]  = maxphys_wo_ghost[dim] + h_s[dim]*kernel_size[dim];
        	minphys_s[dim]  = minphys_wo_ghost[dim] - h_s[dim]*kernel_size[dim];
			// compute number of cells in interior domain and exterior domain
			cell_size_int[dim] = mesh_size[dim] - 1;
        	cell_size_ext[dim] = cell_size_int[dim] + 2*kernel_size[dim];
			// calculate group sizes for 2D. In 2D, the group size is 32x16
			num_groups[dim] = divide_roundUp(cell_size_ext[dim], BLOCK_SIZE_interior[dim]);
			// calculate cell list sizes for the padded exterior cell list
			cell_size_ext_padded[dim] = num_groups[dim]*BLOCK_SIZE_interior[dim];
			// calculate mesh sizes when the mesh is padded
			mesh_size_padded[dim] = num_groups[dim]*BLOCK_SIZE_exterior[dim];
			// calculate number of number of cells and mesh points
        	ncell_ext 	      *= cell_size_ext[dim];
        	ncell_int 		  *= cell_size_int[dim];
        	ncell_ext_padded  *= cell_size_ext_padded[dim];
    	    nmesh             *= mesh_size[dim];
	        nmesh_padded      *= mesh_size_padded[dim];
    	}

		// CREATE BUFFERS //
    	// Buffer 3: Coordinates of minimum extent of physical domain
		dev_minphys = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real)*ndim, minphys_s, &error);
    	// Buffer 4: Coordinates of maximum extent of physical domain
		dev_maxphys = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real)*ndim, maxphys_s, &error);
    	// Buffer 5: h values
		dev_h = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real)*ndim, h_s, &error);
    	// Buffer 6: Exterior cell sizes
		dev_cell_size_ext = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, cell_size_ext, &error);
    	// Buffer 7: Exterior cell sizes when padded
		dev_cell_size_ext_padded = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, cell_size_ext_padded, &error);
    	// Buffer 8: Number of groups
		dev_num_groups = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*ndim, num_groups, &error);
		// Buffer 9: Number of particles per cell
		dev_np_cell  = clCreateBuffer(context, CL_MEM_READ_WRITE,    sizeof(cl_int)*ncell_ext_padded, NULL, &error);
		// Buffer 10: Indices of cells that particles belong to
		dev_p_cell  = clCreateBuffer(context, CL_MEM_READ_WRITE,    sizeof(cl_int)*nparticle, NULL, &error);
    	// Buffer 12: Mesh sizes
    	dev_mesh_size = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real)*ndim, mesh_size, &error);
    	// Buffer 13: Mesh sizes when padded
    	dev_mesh_size_padded = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_t_real)*ndim, mesh_size_padded, &error);
    	// Buffer 14: Mesh points
    	dev_mesh = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_t_real)*nmesh*nmass, NULL, &error);
    	// Buffer 15: Mesh points, padded
    	dev_mesh_padded = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_t_real)*nmesh_padded*nmass, NULL, &error);

    	for(dim = 0; dim < ndim; dim++){
			// store the mesh size
			mesh_size_prev[dim] = mesh_size[dim];
		}

        // Store dimensionality, number of properties, type of interpolation and the pointer to clean-up routine
		ndim_prev    = ndim;
		nmass_prev   = nmass;
		last_call    = p2m_interpolation;
		last_kernel  = *interp_kernel;
#if __KIND == __SINGLE_PRECISION
		ptr2clean_up = &ppm_gpu_p2m_clean_up_s;
#elif __KIND == __DOUBLE_PRECISION
		ptr2clean_up = &ppm_gpu_p2m_clean_up_d;
#endif
	}

	return 0;
}

#if __KIND == __SINGLE_PRECISION
int ppm_gpu_p2m_interpolation_s_
#elif __KIND == __DOUBLE_PRECISION
int ppm_gpu_p2m_interpolation_d_
#endif
(
                    t_real* xp,
                    t_real* mass,
                    t_real* mesh,
                    t_real* minphys_wo_ghost,
                    t_real* maxphys_wo_ghost,
                    int*   mesh_size,
                    int*   np,
					int*   ndim_in,
					int*   nmass_in,
					int*   interp_kernel
) {
    cl_int      error;
    int         n, i, j, dim, mass_val;
    
	int			nparticle = *np;
	int 		ndim 	  = *ndim_in;
	int			nmass	  = *nmass_in;
	cl_kernel	kernel;
  int idim;
	struct timeval t;
	struct timeval t_total;
	t_real t_elapsed, t_elapsed_total;
  // gpu timers
  cl_ulong start, end;
  float executionTimeInMilliseconds;

    // This call will run initialization routines if necessary
#if __KIND == __SINGLE_PRECISION
	ppm_gpu_p2m_initialize_s_(mesh_size, minphys_wo_ghost, maxphys_wo_ghost, np, ndim_in, nmass_in, interp_kernel);
#elif __KIND == __DOUBLE_PRECISION
	ppm_gpu_p2m_initialize_d_(mesh_size, minphys_wo_ghost, maxphys_wo_ghost, np, ndim_in, nmass_in, interp_kernel);
#endif
	
  tic(&t_total);

// CREATE BUFFERS //
    // Buffer 1: Positions of particles
    dev_xp_const 		= clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_t_real)*nparticle*ndim, NULL, &error);
    // Buffer 2: Mass values of particles
    dev_mass_const = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_t_real)*nparticle*nmass, NULL, &error);

	error = clEnqueueWriteBuffer(queue, dev_mass_const, CL_FALSE, 0, sizeof(cl_t_real)*nparticle*nmass, mass, 0, NULL, &event_write_mass);
    checkError(error, "clEnqueueWriteBuffer while writing dev_mass_const", __FUNCTION__, __LINE__, stamp);
	error = clEnqueueWriteBuffer(queue, dev_xp_const, CL_FALSE, 0, sizeof(cl_t_real)*nparticle*ndim, xp, 0, NULL, &event_write_xp);
    checkError(error, "clEnqueueWriteBuffer while writing dev_xp_const", __FUNCTION__, __LINE__, stamp);

    clFinish(queue);
  tic(&t);

// INITIALIZING TO ZERO //
#if __KIND == __SINGLE_PRECISION
	error = p2m_initialize_np_cell_s();
#elif __KIND == __DOUBLE_PRECISION
	error = p2m_initialize_np_cell_d();
#endif
    checkError(error, "Error in p2m_initialize_np_cell", __FUNCTION__, __LINE__, stamp);

// GET CELL INDICES OF PARTICLES AND COUNT NUMBER OF PARTICLES IN CELLS //
#if __KIND == __SINGLE_PRECISION
	error = p2m_get_particle_positions_s(nparticle, ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = p2m_get_particle_positions_d(nparticle, ndim, nmass);
#endif
    checkError(error, "Error in p2m_get_particle_positions", __FUNCTION__, __LINE__, stamp);

// REDUCE NP_CELL TO GET MAXIMUM NUMBER OF PARTICLES PER CELL //
#if __KIND == __SINGLE_PRECISION
    error = p2m_reduce_np_cell_s();
#elif __KIND == __DOUBLE_PRECISION
    error = p2m_reduce_np_cell_d();
#endif
    checkError(error, "Error in p2m_reduce_np_cell", __FUNCTION__, __LINE__, stamp);
    
// COMPUTE SIZE OF PARTICLE DATA BUFFER AND ALLOCATE IT ON THE DEVICE
#if __KIND == __SINGLE_PRECISION
	error = p2m_allocate_data_buffers_s(ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = p2m_allocate_data_buffers_d(ndim, nmass);
#endif
    checkError(error, "Error in p2m_allocate_data_buffer", __FUNCTION__, __LINE__, stamp);

// INITIALIZE PARTICLE DATA BUFFER TO 0
#if __KIND == __SINGLE_PRECISION
	error = init_sorted_particle_mass_s();
#elif __KIND == __DOUBLE_PRECISION
	error = init_sorted_particle_mass_d();
#endif
    checkError(error, "Error in init_sorted_particle_mass", __FUNCTION__, __LINE__, stamp);

// INITIALIZE PADDED MESH TO 0
#if __KIND == __SINGLE_PRECISION
    error = init_padded_mesh_s(nmass);
#elif __KIND == __DOUBLE_PRECISION
    error = init_padded_mesh_d(nmass);
#endif
    checkError(error, "Error in init_padded_mesh", __FUNCTION__, __LINE__, stamp);

// PLACE PARTICLES ACCORDINGLY //
#if __KIND == __SINGLE_PRECISION
    error = place_particle_positions_s(nparticle, ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
    error = place_particle_positions_d(nparticle, ndim, nmass);
#endif
    checkError(error, "Error in place_particle_positions", __FUNCTION__, __LINE__, stamp);
#if __KIND == __SINGLE_PRECISION
	error = place_particle_mass_s(nparticle, ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = place_particle_mass_d(nparticle, ndim, nmass);
#endif
    checkError(error, "Error in place_particle_mass", __FUNCTION__, __LINE__, stamp);

// INTERPOLATE PARTICLES ONTO MESH POINTS. ONE PARTICLES PER CELL IN EACH ITERATION //
#if __KIND == __SINGLE_PRECISION
	error = p2m_interpolate_s(ndim, nmass);
#elif __KIND == __DOUBLE_PRECISION
	error = p2m_interpolate_d(ndim, nmass);
#endif
    checkError(error, "Error in p2m_interpolate", __FUNCTION__, __LINE__, stamp);

// SEW PADDED MESH TO GET THE ACTUAL ONE //
#if __KIND == __SINGLE_PRECISION
	error = sew_padded_mesh_s(mesh_size, ndim);
#elif __KIND == __DOUBLE_PRECISION
	error = sew_padded_mesh_d(mesh_size, ndim);
#endif
    checkError(error, "Error in sew_padded_mesh", __FUNCTION__, __LINE__, stamp);
    clFinish(queue);

    t_elapsed  = toc(&t);
    t_elapsed *= 1000.0;
    
    // get gpu timing for clist
    clGetEventProfilingInfo(event_get_np_cell, CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_get_np_cell, CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 
    executionTimeInMilliseconds = (end - start) * 1.0e-6f;
    printf("GPU P2M (single) clist kernel time: %lf ms\n", executionTimeInMilliseconds);
    // get gpu timing for part place
    executionTimeInMilliseconds = 0.0;
	  for(idim = 0; idim < ndim; idim++){
    clGetEventProfilingInfo(event_place_particle_pos[idim], CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_place_particle_pos[idim], CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 
    executionTimeInMilliseconds += (end - start) * 1.0e-6f;
    }
    printf("GPU P2M (single) part_place (%d dims) kernel time: %lf ms\n", ndim,executionTimeInMilliseconds);
    // get gpu timing for interpolation
    clGetEventProfilingInfo(event_p2m_interpolate, CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_p2m_interpolate, CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 
    executionTimeInMilliseconds = (end - start) * 1.0e-6f;
    printf("GPU P2M (single) interp kernel time: %lf ms\n", executionTimeInMilliseconds);
    // get gpu timing for sewing
    clGetEventProfilingInfo(event_sew_padded_mesh, CL_PROFILING_COMMAND_END, 
                                sizeof(cl_ulong), &end, NULL); 
    clGetEventProfilingInfo(event_sew_padded_mesh, CL_PROFILING_COMMAND_START,  
                                sizeof(cl_ulong), &start, NULL); 
    executionTimeInMilliseconds = (end - start) * 1.0e-6f;
    printf("GPU P2M (single) sewing kernel time: %lf ms\n", executionTimeInMilliseconds);


// READ BACK RESULTS //
	error = clEnqueueWaitForEvents(queue, 1, &event_sew_padded_mesh);
    error = clEnqueueReadBuffer(queue, dev_mesh, CL_FALSE, 0, sizeof(t_real)*nmesh*nmass, mesh, 0, NULL, &event_read_mesh);
    checkError(error, "clEnqueueReadBuffer while reading back", __FUNCTION__, __LINE__, stamp);

	error = clEnqueueWaitForEvents(queue, 1, &event_read_mesh);
    clFinish(queue);

    t_elapsed_total  = toc(&t_total);
    t_elapsed_total *= 1000.0;

    printf("P2M GPU Time consumtion: %lf ms including data transfer: %lf ms \n", t_elapsed, t_elapsed_total);

// RELEASE XP and MASS MEMORY OBJECTS
	error |= clReleaseMemObject(dev_xp_const);
	error |= clReleaseMemObject(dev_mass_const);
	error |= clReleaseMemObject(dev_sorted_particle_position);
	error |= clReleaseMemObject(dev_sorted_particle_mass);
    checkError(error, "Error in releasing memory objects", __FUNCTION__, __LINE__, stamp);

    clFinish(queue);

    return 0;
}

#if __KIND == __SINGLE_PRECISION
int p2m_initialize_np_cell_s()
#elif __KIND == __DOUBLE_PRECISION
int p2m_initialize_np_cell_d()
#endif
{
	size_t localWorkSize, globalWorkSize, num_wg;
	int n, error;
	cl_kernel	kernel;

	// Initialize np_cell array to 0
    kernel = clCreateKernel(program, "kernel_clist_init_to_zero_int", &error);
        
    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_np_cell);
	error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&ncell_ext_padded);
	
    // Set work-group sizes
	localWorkSize  = 512;
	globalWorkSize = ncell_ext_padded;
	if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the device
	error |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_np_cell);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int	p2m_get_particle_positions_s
#elif __KIND == __DOUBLE_PRECISION
int	p2m_get_particle_positions_d
#endif
(
						int nparticle, 
						int ndim, 
						int nmass 
){
	cl_kernel kernel;
	int n, error;
	size_t localWorkSize, globalWorkSize, num_wg;

	kernel = clCreateKernel(program, "kernel_clist_get_particle_positions", &error);

    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(kernel, n++, sizeof(cl_mem), (void *)&dev_xp_const);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&nparticle);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&ndim);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&nmass);
    error |= clSetKernelArg(kernel, n++, sizeof(int), (void *)&ncell_ext_padded);
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

	error |= clEnqueueWaitForEvents(queue, 1, &event_write_xp);
	error |= clEnqueueWaitForEvents(queue, 1, &event_init_np_cell);

    // Launch kernel on the device
    error |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_np_cell);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int p2m_reduce_np_cell_s()
#elif __KIND == __DOUBLE_PRECISION
int p2m_reduce_np_cell_d()
#endif

{
	int 	  n, error;
	size_t 	  localWorkSize, globalWorkSize, num_wg;
	cl_kernel ck_init_max_np, ck_init_max_np2, ck_init_max_np3, ck_get_max_np, ck_get_max_np2, ck_get_max_np3; 
	size_t 	  size_max_np = 0, size_max_np2 = 0, size_max_np3 = 0;

//////// Initialize max_np array to 0
	size_max_np = (ncell_ext_padded + 1023)/1024;
	dev_max_np  = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_int)*size_max_np, NULL, &error); 

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

  // Launch kernel on the device
	error |= clEnqueueNDRangeKernel(queue, ck_init_max_np, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_max_np);

//////// Initialize max_np2 array to 0
	if(size_max_np > 1){
	size_max_np2 = (size_max_np + 1023)/1024;
	dev_max_np2  = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_int)*size_max_np2, NULL, &error); 

  ck_init_max_np2 = clCreateKernel(program, "kernel_clist_init_to_zero_int", &error);
  //checkError(error, "clCreateKernel 'clist_init_kernel' function", __FUNCTION__, __LINE__, stamp);
        
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

  // Launch kernel on the device
	error |= clEnqueueNDRangeKernel(queue, ck_init_max_np2, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_max_np2);
	}
//////// Initialize max_np3 array to 0
	if(size_max_np2 > 1){
	  size_max_np3 = (size_max_np2 + 1023)/1024;
	  dev_max_np3  = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_int)*size_max_np3, NULL, &error); 

    ck_init_max_np3 = clCreateKernel(program, "kernel_clist_init_to_zero_int", &error);
    //checkError(error, "clCreateKernel 'clist_init_kernel' function", __FUNCTION__, __LINE__, stamp);
        
    // Set kernel arguments
	  n = 0;
    error |= clSetKernelArg(ck_init_max_np3, n++, sizeof(cl_mem), (void *)&dev_max_np3);
	  error |= clSetKernelArg(ck_init_max_np3, n++, sizeof(int), (void *)&size_max_np3);
	
    // Set work-group sizes
	  localWorkSize  = 512;
	  globalWorkSize = size_max_np3;
	  if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the device
	  error |= clEnqueueNDRangeKernel(queue, ck_init_max_np3, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_max_np3);
	}
///////////////
	// Reduce large blocks
	ck_get_max_np = clCreateKernel(program, "clist_get_max_np_blockwise", &error);

  // Set kernel arguments
	n = 0;
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(cl_mem), (void *)&dev_np_cell);
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(int), (void *)&ncell_ext_padded);
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(cl_int)*1536, NULL);
  error |= clSetKernelArg(ck_get_max_np, n++, sizeof(cl_mem), (void *)&dev_max_np);

  // Set work-group sizes
	localWorkSize = 512;
	num_wg = size_max_np;
	globalWorkSize = localWorkSize*size_max_np;

	error |= clEnqueueWaitForEvents(queue, 1, &event_get_np_cell);
	error |= clEnqueueWaitForEvents(queue, 1, &event_init_max_np);

  // Launch kernel on the device
  error |= clEnqueueNDRangeKernel(queue, ck_get_max_np, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_max_np);

///////////////
	if(size_max_np > 1){  // Reduce smaller blocks
	    ck_get_max_np2 = clCreateKernel(program, "clist_get_max_np_blockwise", &error);

      // Set kernel arguments
	    n = 0;
      error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(cl_mem), (void *)&dev_max_np);
      error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(int), (void *)&size_max_np);
      error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(cl_int)*1536, NULL);
      error |= clSetKernelArg(ck_get_max_np2, n++, sizeof(cl_mem), (void *)&dev_max_np2);

      // Calculate work-group sizes
	    localWorkSize = 512;
	    num_wg = size_max_np2;
	    globalWorkSize = localWorkSize*size_max_np2;

	    error |= clEnqueueWaitForEvents(queue, 1, &event_init_max_np2);
	    error |= clEnqueueWaitForEvents(queue, 1, &event_get_max_np);

      // Launch kernel on the device
      error |= clEnqueueNDRangeKernel(queue, ck_get_max_np2, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_max_np2);
	}
///////////////
	if(size_max_np2 > 1){  // Reduce even smaller blocks
	    ck_get_max_np3 = clCreateKernel(program, "clist_get_max_np_blockwise", &error);

      // Set kernel arguments
	    n = 0;
      error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(cl_mem), (void *)&dev_max_np2);
      error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(int), (void *)&size_max_np2);
      error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(cl_int)*1536, NULL);
      error |= clSetKernelArg(ck_get_max_np3, n++, sizeof(cl_mem), (void *)&dev_max_np3);

      // Calculate work-group sizes
	    localWorkSize = 512;
	    num_wg = size_max_np3;
	    globalWorkSize = localWorkSize*size_max_np2;

	    error |= clEnqueueWaitForEvents(queue, 1, &event_init_max_np3);
	    error |= clEnqueueWaitForEvents(queue, 1, &event_get_max_np2);

      // Launch kernel on the device
      error |= clEnqueueNDRangeKernel(queue, ck_get_max_np3, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_get_max_np3);
	}
// READ BACK MAXIMUM NUMBER OF PARTICLES
  // Choose buffer which holds the maximum value in the first index
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
	} else if(size_max_np > 1){
		error |= clReleaseMemObject(dev_max_np2);
	}
  
	return error;
}

#if __KIND == __SINGLE_PRECISION
int	p2m_allocate_data_buffers_s
#elif __KIND == __DOUBLE_PRECISION
int	p2m_allocate_data_buffers_d
#endif
(
	int ndim,
	int nmass
){
	int error = CL_SUCCESS;

	error |= clEnqueueWaitForEvents(queue, 1, &event_read_max_np);

    // Calculate sizes of cell arrays which will store particle positions and strengths
	n_sorted_particle_position = max_np*ndim*ncell_ext_padded + 1536; // Extra padding in order to avoid if branch
	size_sorted_particle_position = sizeof(t_real)*n_sorted_particle_position;
	
	n_sorted_particle_mass = max_np*nmass*ncell_ext_padded + 1536; // Extra padding in order to avoid if branch
	size_sorted_particle_mass = sizeof(t_real)*n_sorted_particle_mass;
	
    // Allocate cell arrays which will store particle positions and strengths
	dev_sorted_particle_position = clCreateBuffer(context, CL_MEM_READ_WRITE, size_sorted_particle_position, NULL, &error);
	dev_sorted_particle_mass = clCreateBuffer(context, CL_MEM_READ_WRITE, size_sorted_particle_mass, NULL, &error);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int	init_sorted_particle_mass_s(){
#elif __KIND == __DOUBLE_PRECISION
int	init_sorted_particle_mass_d(){
#endif
	size_t localWorkSize, globalWorkSize, num_wg;
	int n, error;
	cl_kernel ck_init_sorted_particle_mass;

    ck_init_sorted_particle_mass = clCreateKernel(program, "kernel_clist_init_to_zero_real", &error);

    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(ck_init_sorted_particle_mass, n++, sizeof(cl_mem), (void *)&dev_sorted_particle_mass);
	error |= clSetKernelArg(ck_init_sorted_particle_mass, n++, sizeof(int), (void *)&n_sorted_particle_mass);
	
    // Set work-group sizes
	localWorkSize  = 512;
	globalWorkSize = n_sorted_particle_mass;
	if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the device
	error |= clEnqueueNDRangeKernel(queue, ck_init_sorted_particle_mass, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_sorted_particle_mass);
	
	return error;
}

#if __KIND == __SINGLE_PRECISION
int init_padded_mesh_s(int nmass)
#elif __KIND == __DOUBLE_PRECISION
int init_padded_mesh_d(int nmass)
#endif
{
	int n, error;
	int size_padded_mesh;
	cl_kernel ck_init_padded_mesh;	

    ck_init_padded_mesh = clCreateKernel(program, "kernel_clist_init_to_zero_real", &error);

    // Set kernel arguments
	n = 0;
    error |= clSetKernelArg(ck_init_padded_mesh, n++, sizeof(cl_mem), (void *)&dev_mesh_padded);
    
	size_padded_mesh = nmesh_padded*nmass;
	error |= clSetKernelArg(ck_init_padded_mesh, n++, sizeof(int), (void *)&size_padded_mesh);
	
    // Set work-group sizes
	localWorkSize  = 512;
	globalWorkSize = size_padded_mesh;
	if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
    globalWorkSize = num_wg*localWorkSize;

    // Launch kernel on the device
	error |= clEnqueueNDRangeKernel(queue, ck_init_padded_mesh, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_init_padded_mesh);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int place_particle_positions_s
#elif __KIND == __DOUBLE_PRECISION
int place_particle_positions_d
#endif
(
						int nparticle,
						int ndim,
						int nmass
){
	int    dim, n, error;
	size_t localWorkSize, globalWorkSize, num_wg;
	cl_kernel	ck_place_particle_pos[3];

	// Place particles' positions, one position in each iteration. This is better as we use caches more efficiently and 
	// reduce the drawback of uncoalesced writes to global memory.
	for(dim = 0; dim < ndim; dim++){
		ck_place_particle_pos[dim] = clCreateKernel(program, "kernel_clist_place_particle_pos", &error);

        // Set kernel arguments
		n = 0;
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(cl_mem), (void *)&dev_xp_const);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(cl_mem), (void *)&dev_mass_const);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(int), (void *)&nparticle);
    	error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(int), (void *)&ndim);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(int), (void *)&nmass);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(int), (void *)&dim);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(int), (void *)&ncell_ext_padded);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(int), (void *)&max_np);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(cl_mem), (void *)&dev_p_cell);
	    error |= clSetKernelArg(ck_place_particle_pos[dim], n++, sizeof(cl_mem), (void *)&dev_sorted_particle_position);

        // Set work-group sizes
		localWorkSize  = 512;	
		globalWorkSize = nparticle;
		if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
	    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
	    globalWorkSize = num_wg*localWorkSize;

        // Launch kernel on the device
	    error |= clEnqueueNDRangeKernel(queue, ck_place_particle_pos[dim], 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_place_particle_pos[dim]);
	}

	return error;
}

#if __KIND == __SINGLE_PRECISION
int	place_particle_mass_s
#elif __KIND == __DOUBLE_PRECISION
int	place_particle_mass_d
#endif
(
						int nparticle,
						int ndim,
						int nmass
){
	int mass_val, n, error = CL_SUCCESS;
	size_t localWorkSize, globalWorkSize, num_wg;
	cl_kernel ck_place_particle_mass[7];

	error |= clEnqueueWaitForEvents(queue, 1, &event_write_mass);
	error |= clEnqueueWaitForEvents(queue, 1, &event_init_sorted_particle_mass);

	for(mass_val = 0; mass_val < nmass; mass_val++){
		ck_place_particle_mass[mass_val] = clCreateKernel(program, "kernel_clist_place_particle_mass", &error);

        // Set kernel arguments
		n = 0;
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(cl_mem), (void *)&dev_xp_const);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(cl_mem), (void *)&dev_mass_const);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(int), (void *)&nparticle);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(int), (void *)&ndim);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(int), (void *)&nmass);
 	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(int), (void *)&mass_val);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(int), (void *)&ncell_ext_padded);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(int), (void *)&max_np);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(cl_mem), (void *)&dev_p_cell);
	    error |= clSetKernelArg(ck_place_particle_mass[mass_val], n++, sizeof(cl_mem), (void *)&dev_sorted_particle_mass);

        // Set work-group sizes
		localWorkSize  = 512;	
		globalWorkSize = nparticle;
		if(localWorkSize > globalWorkSize) localWorkSize = globalWorkSize;
	    num_wg = (globalWorkSize + localWorkSize - 1)/localWorkSize;
	    globalWorkSize = num_wg*localWorkSize;

        // Launch kernel on the device
	    error |= clEnqueueNDRangeKernel(queue, ck_place_particle_mass[mass_val], 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &event_place_particle_mass[mass_val]);
	}

	return error;
}

#if __KIND == __SINGLE_PRECISION
int	p2m_interpolate_s
#elif __KIND == __DOUBLE_PRECISION
int	p2m_interpolate_d
#endif
(
					int ndim,
					int nmass
){
	int dim, mass_val, n, error;
	size_t localWorkSizeND[3], globalWorkSizeND[3], num_wgND[3];
	cl_kernel	ck_interpolate;

	for(dim = 0; dim < ndim; dim++){
		error |= clEnqueueWaitForEvents(queue, 1, &event_place_particle_pos[dim]);
	}

	for(mass_val = 0; mass_val < nmass; mass_val++){
		error |= clEnqueueWaitForEvents(queue, 1, &event_place_particle_mass[mass_val]);
	}

   	ck_interpolate = clCreateKernel(program, "kernel_interpolate", &error);
       
    // Set kernel arguments
    n = 0;
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_mem), (void *)&dev_sorted_particle_position);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_mem), (void *)&dev_sorted_particle_mass);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_t_real)*interp_lmem_size, NULL);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_int), (void *)&ndim);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_int), (void *)&nmass);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_int), (void *)&max_np);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_int), (void *)&ncell_ext_padded);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_int), (void *)&nmesh_padded);
 	error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_mem), (void *)&dev_minphys);
    error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_mem), (void *)&dev_h);
	error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_mem), (void *)&dev_mesh_size_padded);
	error |= clSetKernelArg(ck_interpolate, n++, sizeof(cl_mem), (void *)&dev_mesh_padded);

    // Calculate work-group sizes
    for(dim = 0; dim < ndim; dim++){
    	localWorkSizeND[dim]  = BLOCK_SIZE_interior[dim];
       	num_wgND[dim]         = num_groups[dim];
       	globalWorkSizeND[dim] = num_wgND[dim]*localWorkSizeND[dim];
   	}    

    // Launch kernel on the device
	error |= clEnqueueNDRangeKernel(queue, ck_interpolate, ndim, NULL, globalWorkSizeND, localWorkSizeND, 0, NULL, &event_p2m_interpolate);
  
  printf("p2m (single) interpolate max_np: %d\n",max_np);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int	sew_padded_mesh_s
#elif __KIND == __DOUBLE_PRECISION
int	sew_padded_mesh_d
#endif
(
					int *mesh_size,
					int  ndim
){
	int  dim, n, error = CL_SUCCESS;
	size_t localWorkSizeND[3], globalWorkSizeND[3], num_wgND[3];
	cl_kernel	ck_sew_padded_mesh;

	error |= clEnqueueWaitForEvents(queue, 1, &event_p2m_interpolate);

    ck_sew_padded_mesh = clCreateKernel(program, "kernel_sew_padded_mesh", &error);
        
    // Set kernel arguments
    n = 0;
    error |= clSetKernelArg(ck_sew_padded_mesh, n++, sizeof(cl_mem), (void *)&dev_mesh_padded);
    error |= clSetKernelArg(ck_sew_padded_mesh, n++, sizeof(cl_mem), (void *)&dev_mesh);
    error |= clSetKernelArg(ck_sew_padded_mesh, n++, sizeof(cl_mem), (void *)&dev_mesh_size);
    error |= clSetKernelArg(ck_sew_padded_mesh, n++, sizeof(cl_mem), (void *)&dev_mesh_size_padded);
    error |= clSetKernelArg(ck_sew_padded_mesh, n++, sizeof(cl_int), (void *)&nmesh_padded);
    error |= clSetKernelArg(ck_sew_padded_mesh, n++, sizeof(cl_mem), (void *)&dev_num_groups);
    
    // Calculate work-group sizes
    for(dim = 0; dim < ndim; dim++) {
    	localWorkSizeND[dim]  = BLOCK_SIZE_interior[dim];
        num_wgND[dim]         = (mesh_size[dim] + localWorkSizeND[dim] - 1)/localWorkSizeND[dim];
        globalWorkSizeND[dim] = num_wgND[dim]*localWorkSizeND[dim];
    }    

    // Launch kernel on the device
	error |= clEnqueueNDRangeKernel(queue, ck_sew_padded_mesh, ndim, NULL, globalWorkSizeND, localWorkSizeND, 0, NULL, &event_sew_padded_mesh);

	return error;
}

#if __KIND == __SINGLE_PRECISION
int ppm_gpu_p2m_clean_up_s(){
#elif __KIND == __DOUBLE_PRECISION
int ppm_gpu_p2m_clean_up_d(){
#endif
	int error = CL_SUCCESS;

	error |= clReleaseMemObject(dev_minphys);
	error |= clReleaseMemObject(dev_maxphys);
	error |= clReleaseMemObject(dev_h);
	error |= clReleaseMemObject(dev_cell_size_ext);
	error |= clReleaseMemObject(dev_cell_size_ext_padded);
	error |= clReleaseMemObject(dev_num_groups);
	error |= clReleaseMemObject(dev_np_cell);
	error |= clReleaseMemObject(dev_p_cell);
	error |= clReleaseMemObject(dev_mesh_size);
	error |= clReleaseMemObject(dev_mesh_size_padded);
	error |= clReleaseMemObject(dev_mesh);
	error |= clReleaseMemObject(dev_mesh_padded);

	return error;
}

