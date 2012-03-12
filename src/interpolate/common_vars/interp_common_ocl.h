#ifndef INTERP_COMMON_OCL_H
#define INTERP_COMMON_OCL_H

#include "../../gpu/opencl_utils.h"

extern cl_platform_id   platform;
extern cl_device_id     device;
extern cl_context       context;
extern cl_command_queue queue;
extern cl_program       program;

extern cl_mem      dev_xp_const;
extern cl_mem      dev_mass_const;
extern cl_mem      dev_mass;
extern cl_mem      dev_minphys;
extern cl_mem	    dev_maxphys;
extern cl_mem  	dev_h;

extern cl_mem      dev_mesh_size;
extern cl_mem      dev_mesh_size_padded;
extern cl_mem      dev_mesh_size_int;
extern cl_mem      dev_mesh_size_int_padded;
extern cl_mem      dev_mesh_size_ext;
extern cl_mem      dev_mesh_size_ext_padded;

extern cl_mem      dev_num_groups;

extern cl_mem      dev_cell_size_int;
extern cl_mem      dev_cell_size_int_padded;
extern cl_mem      dev_cell_size_ext;
extern cl_mem      dev_cell_size_ext_padded;

extern cl_mem		dev_sorted_particle_position;
extern cl_mem		dev_sorted_particle_mass;

extern cl_mem      dev_max_np;
extern cl_mem      dev_max_np2;
extern cl_mem      dev_max_np3;

extern cl_mem		dev_np_cell;
extern cl_mem		dev_p_cell;

extern cl_mem      dev_mesh;
extern cl_mem      dev_mesh_const;
extern cl_mem      dev_mesh_padded;

extern cl_event    event_write_xp;
extern cl_event    event_write_mass;
extern cl_event    event_write_mesh;
extern cl_event    event_init_np_cell;
extern cl_event    event_init_padded_mesh;
extern cl_event    event_init_sorted_particle_mass;
extern cl_event    event_get_np_cell;
extern cl_event    event_init_max_np;
extern cl_event    event_init_max_np2;
extern cl_event    event_init_max_np3;
extern cl_event    event_get_max_np;
extern cl_event    event_get_max_np2;
extern cl_event    event_get_max_np3;
extern cl_event    event_read_max_np;
extern cl_event    event_read_max_np2;
extern cl_event    event_read_max_np3;
extern cl_event    event_duplicate_mesh;
extern cl_event    event_place_particle_pos[3];
extern cl_event    event_place_particle_mass[7];
extern cl_event    event_m2p_interpolate;
extern cl_event    event_p2m_interpolate;
extern cl_event    event_collect_particles;
extern cl_event    event_sew_padded_mesh;
extern cl_event    event_read_mesh;

#endif
