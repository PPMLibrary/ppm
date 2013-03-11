#include "interp_common_ocl.h"

cl_platform_id   platform;
cl_device_id     device;
cl_context       context;
cl_command_queue queue;
cl_program       program;

cl_mem      test_buffer; //Allocate for testing purposes only!

cl_mem      dev_xp_const;
cl_mem      dev_mass_const;
cl_mem      dev_mass;
cl_mem      dev_minphys;
cl_mem	    dev_maxphys;
cl_mem  	dev_h;

cl_mem      dev_mesh_size;
cl_mem      dev_mesh_size_padded;
cl_mem      dev_mesh_size_int;
cl_mem      dev_mesh_size_int_padded;
cl_mem      dev_mesh_size_ext;
cl_mem      dev_mesh_size_ext_padded;

cl_mem      dev_num_groups;

cl_mem      dev_cell_size_int;
cl_mem      dev_cell_size_int_padded;
cl_mem      dev_cell_size_ext;
cl_mem      dev_cell_size_ext_padded;

cl_mem		dev_sorted_particle_position;
cl_mem		dev_sorted_particle_mass;

cl_mem      dev_max_np;
cl_mem      dev_max_np2;
cl_mem      dev_max_np3;

cl_mem		dev_np_cell;
cl_mem		dev_p_cell;

cl_mem      dev_mesh;
cl_mem      dev_mesh_const;
cl_mem      dev_mesh_padded;

cl_event    event_write_xp;
cl_event    event_write_mass;
cl_event    event_write_mesh;
cl_event    event_init_np_cell;
cl_event    event_init_padded_mesh;
cl_event    event_init_sorted_particle_mass;
cl_event    event_get_np_cell;
cl_event    event_init_max_np;
cl_event    event_init_max_np2;
cl_event    event_init_max_np3;
cl_event    event_get_max_np;
cl_event    event_get_max_np2;
cl_event    event_get_max_np3;
cl_event    event_read_max_np;
cl_event    event_read_max_np2;
cl_event    event_read_max_np3;
cl_event    event_duplicate_mesh;
cl_event    event_place_particle_pos[3];
cl_event    event_place_particle_mass[7];
cl_event    event_m2p_interpolate;
cl_event    event_p2m_interpolate;
cl_event    event_collect_particles;
cl_event    event_sew_padded_mesh;
cl_event    event_read_mesh;

