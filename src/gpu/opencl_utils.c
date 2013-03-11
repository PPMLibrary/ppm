
/**
 * Copyright (c) 2011 MOSAIC Group (ETH Zurich)
 * 
 * Author: Omar Awile
 *
 * Note:
 * This software contains source code provided by NVIDIA Corporation.
 **/
#include "opencl_utils.h"

static char stamp_opencl_utils[] = "***\nmodule '' __FILE__ ''\ncompiled '' __TIMESTAMP__ ''\n***";

//const char* option = "-cl-nv-maxrregcount=21";

void tic(struct timeval* t)
{
    gettimeofday(t,NULL);
}

double toc(struct timeval* t)
{
    struct timeval delta,now;
    gettimeofday(&now,NULL);
    timersub(&now,t,&delta);
    return delta.tv_sec + delta.tv_usec/1000000.0; 
}

void clReadCode(char* fileName, char** code, size_t* len)
{
    FILE* fp = NULL;
    size_t rsize, size = 0;

    fp = fopen(fileName,"r");
	if(fp != NULL){
	    fseek(fp, 0, SEEK_END);
	    size = ftell(fp);
	    fseek(fp, 0, SEEK_SET);
	    *code = (char*)malloc(sizeof(char)*size);
	    rsize = fread(*code,1,size,fp);
	    if (rsize != size)
	        printf("Help!! couldn't read as much as I should have\n");
	    fclose(fp);
	    *len = size;
	}
	else{
		printf("******\nCouldn't find the kernel file: %s\n******\n", fileName);
		exit(0);
	}

}

void print_info(const char *errinfo,
                const void *private_info,
                size_t cb,
                void *user_data)
{
    printf("LOG message: %s \n", errinfo);

    return;
}

int clInit(cl_platform_id *platform, cl_device_id *device, cl_context *context, cl_command_queue *queue, cl_program *program)
{
    void (*ptr_print_info)(const char*, const void*, size_t, void*);
    ptr_print_info = &print_info;

    const char *errinfo;
    const void *private_info;
    size_t cb;
    void *user_data;

    // initialize OpenCL
    cl_int error = 0;
    error = clGetPlatformIDs(1, platform, NULL);
    if (error != CL_SUCCESS) {
        printf("Error getting platform id: %s\n", clErrorString(error));
    }

#if __GPU == AMD_CPU
    error = clGetDeviceIDs(*platform, CL_DEVICE_TYPE_CPU, 1, device, NULL);
#else
    error = clGetDeviceIDs(*platform, CL_DEVICE_TYPE_GPU, 1, device, NULL);
#endif
    if (error != CL_SUCCESS) {
        printf("Error getting device id: %s\n", clErrorString(error));
    }

    *context = clCreateContext(NULL,1, device, ptr_print_info, NULL, &error);
    if (error != CL_SUCCESS) {
        printf("Error creating context: %s\n", clErrorString(error));
    }
    
    *queue = clCreateCommandQueue(*context, *device, 
             CL_QUEUE_PROFILING_ENABLE, &error);
    //*queue = clCreateCommandQueue(*context, *device, 0, &error);
    if (error != CL_SUCCESS) {
        printf("Error creating cmd queue: %s\n",clErrorString(error));
    }
    return error;
}

int clGetWorkGroupSize(cl_context *context, cl_device_id *device, 
    int dim, size_t *workgroup_size, cl_ulong *lmemsize)
{
    cl_program program;
    cl_int error = 0;
    cl_kernel	kernel;
    char** cl_code;
    size_t code_lens[1] = {41};
    size_t wgsize;
    size_t wgsizemult;
    cl_code = (char**) malloc(sizeof(char*));
    cl_code[0] = (char*) malloc(sizeof(char)*41);
    cl_code[0] = "__kernel void ppm_wgroup_size_test() { }";
    program = clCreateProgramWithSource(*context, 1, 
            (const char**)cl_code, code_lens, &error);
    if (error != CL_SUCCESS) {
        printf("Error creating program: %s\n",clErrorString(error));
    }
    error = clBuildProgram(program, 1, device, NULL, NULL, NULL);
  

    // Initialize np_cell array to 0
    kernel = clCreateKernel(program, "ppm_wgroup_size_test", &error);

    clGetKernelWorkGroupInfo (kernel,*device,
            CL_KERNEL_WORK_GROUP_SIZE,sizeof(wgsize), &wgsize, NULL);
    
    clGetKernelWorkGroupInfo (kernel,*device,
            CL_KERNEL_LOCAL_MEM_SIZE,sizeof(cl_ulong), lmemsize, NULL);

    workgroup_size[0] = 32;
    if (dim == 3) {
      workgroup_size[1] = 2;
      wgsize = wgsize / workgroup_size[0];
      wgsize /= 2;
      while (wgsize > workgroup_size[1]) {
        wgsize /= 2;
        workgroup_size[1] *= 2;
      }
      workgroup_size[2] = wgsize;
    } else {
      workgroup_size[1] = wgsize / workgroup_size[0];
    }
    return error;
}

int clBuild(char **filename, int fcount, cl_program *program, cl_context *context, cl_device_id *device, const char* option)
{   
    // load and compile CL program
    cl_int error = 0;
    char **cl_code;
    size_t *code_lens;
    int i;

    char* build_log;
    size_t log_len;

    cl_code = (char**)malloc(sizeof(char*)*fcount);
    code_lens = (size_t*)malloc(sizeof(size_t)*fcount);
   
    for (i=0;i<fcount;i++)
        clReadCode(filename[i],&(cl_code[i]),&(code_lens[i]));
    
    *program = clCreateProgramWithSource(*context, fcount, 
            (const char**)cl_code, code_lens, &error);
    if (error != CL_SUCCESS) {
        printf("Error creating program: %s\n",clErrorString(error));
    }
    error = clBuildProgram(*program, 1, device, option, NULL, NULL);
    //checkError(error, "clBuildProgram", __FUNCTION__, __LINE__, stamp_opencl_utils);

    // Show the build log
    error = clGetProgramBuildInfo(*program, *device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_len);
    checkError(error, "clGetBuildProgramInfo", __FUNCTION__, __LINE__, stamp_opencl_utils);

    build_log = (char*)malloc(sizeof(char)*(log_len + 1));
    error = clGetProgramBuildInfo(*program, *device, CL_PROGRAM_BUILD_LOG, log_len,
            build_log, NULL);
    checkError(error, "clGetBuildProgramInfo", __FUNCTION__, __LINE__, stamp_opencl_utils);
    build_log[log_len] = '\0';
    //printf("%s\n",build_log);
    free(build_log);

    return CL_SUCCESS;
}

const char* clErrorString(cl_int error)
{
    static const char* errorString[] = {
        "CL_SUCCESS",
        "CL_DEVICE_NOT_FOUND",
        "CL_DEVICE_NOT_AVAILABLE",
        "CL_COMPILER_NOT_AVAILABLE",
        "CL_MEM_OBJECT_ALLOCATION_FAILURE",
        "CL_OUT_OF_RESOURCES",
        "CL_OUT_OF_HOST_MEMORY",
        "CL_PROFILING_INFO_NOT_AVAILABLE",
        "CL_MEM_COPY_OVERLAP",
        "CL_IMAGE_FORMAT_MISMATCH",
        "CL_IMAGE_FORMAT_NOT_SUPPORTED",
        "CL_BUILD_PROGRAM_FAILURE",
        "CL_MAP_FAILURE",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "CL_INVALID_VALUE",
        "CL_INVALID_DEVICE_TYPE",
        "CL_INVALID_PLATFORM",
        "CL_INVALID_DEVICE",
        "CL_INVALID_CONTEXT",
        "CL_INVALID_QUEUE_PROPERTIES",
        "CL_INVALID_COMMAND_QUEUE",
        "CL_INVALID_HOST_PTR",
        "CL_INVALID_MEM_OBJECT",
        "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
        "CL_INVALID_IMAGE_SIZE",
        "CL_INVALID_SAMPLER",
        "CL_INVALID_BINARY",
        "CL_INVALID_BUILD_OPTIONS",
        "CL_INVALID_PROGRAM",
        "CL_INVALID_PROGRAM_EXECUTABLE",
        "CL_INVALID_KERNEL_NAME",
        "CL_INVALID_KERNEL_DEFINITION",
        "CL_INVALID_KERNEL",
        "CL_INVALID_ARG_INDEX",
        "CL_INVALID_ARG_VALUE",
        "CL_INVALID_ARG_SIZE",
        "CL_INVALID_KERNEL_ARGS",
        "CL_INVALID_WORK_DIMENSION",
        "CL_INVALID_WORK_GROUP_SIZE",
        "CL_INVALID_WORK_ITEM_SIZE",
        "CL_INVALID_GLOBAL_OFFSET",
        "CL_INVALID_EVENT_WAIT_LIST",
        "CL_INVALID_EVENT",
        "CL_INVALID_OPERATION",
        "CL_INVALID_GL_OBJECT",
        "CL_INVALID_BUFFER_SIZE",
        "CL_INVALID_MIP_LEVEL",
        "CL_INVALID_GLOBAL_WORK_SIZE",
    };

    const int errorCount = sizeof(errorString) / sizeof(errorString[0]);

    const int index = -error;

    return (index >= 0 && index < errorCount) ? errorString[index] : "Unspecified Error";
}

void checkError(cl_int error, const char* string, const char* funct_name, int line, char* stamp)
{   
    if(error != CL_SUCCESS){
        printf("***\n%s failed in function '%s' at line %d due to %s\n%s\n", string, funct_name, line,
                clErrorString(error), stamp);
        exit(0);
    }
}

int divide_roundUp(int toBeRounded, int divider){
	return (toBeRounded + divider - 1)/divider;
}
