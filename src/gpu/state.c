#include "state.h"
#include <stdlib.h> 

int GPU_initialized = 0;
int last_call	    = NONE;

int (*ptr2clean_up)() = NULL;

// Generic clean-up function to deallocate
// buffers on the GPU. Add the specific function
// call, as you add a new functionality.
int	ppm_gpu_clean_up(){
	int error = 0;

	if(*ptr2clean_up != NULL)
		error = (*ptr2clean_up)();

	return error;
}
