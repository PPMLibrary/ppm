#ifndef STATE_H
#define STATE_H

#define	NONE 0
#define p2m_interpolation 1
#define m2p_interpolation 2

extern int GPU_initialized;
extern int last_call;

extern int(*ptr2clean_up)();
int	ppm_gpu_clean_up();

#endif
