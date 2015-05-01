      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_core
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_core
      !!! This is the global user module. It contains all
      !!! user-callable routines of the entire ppm core library.
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_data
         USE ppm_module_substart
         USE ppm_module_substop
         USE ppm_module_write
         USE ppm_module_error
         USE ppm_module_map
         USE ppm_module_core_util
         USE ppm_module_alloc
         USE ppm_module_loadbal, ONLY : ppm_loadbal_inquire,          &
         &   ppm_estimate_procspeed,ppm_get_cost,ppm_set_decomp_cost, &
         &   ppm_set_proc_speed
         USE ppm_module_topo
         USE ppm_module_io, ONLY : ppm_io,ppm_doio,ppm_io_open,ppm_io_close, &
         &   ppm_io_inquire,ppm_io_read_ascii,ppm_io_read_binary,ppm_io_write_ascii, &
         &   ppm_io_write_binary,ppm_io_delete,ppm_io_set_unit,ppm_io_unused_unit
         USE ppm_module_neighlist, ONLY : ppm_clist_destroy,ppm_neighlist_MkNeighIdx, &
         &   ppm_neighlist_clist,ppm_neighlist_vlist,ppm_t_clist
         USE ppm_module_kdtree, ONLY : kdtree_s,kdtree_d,kdtree_result_s,kdtree_result_d
         USE ppm_module_inl_k_vlist, ONLY : kdtree_n_nearest,kdtree_n_nearest_around_point, &
         &   kdtree_r_nearest,kdtree_r_nearest_around_point,kdtree_r_count, &
         &   kdtree_r_count_around_point
         USE ppm_module_tree, ONLY : ppm_tree
         USE ppm_module_mesh
         USE ppm_module_io_vtk, ONLY : ppm_vtk_particles,ppm_vtk_mesh_2d, &
         &   ppm_vtk_mesh_3d
         USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo !,ppm_t_clist
         USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch,ppm_t_equi_mesh
         USE ppm_module_particles_typedef, ONLY : ppm_t_particles_s,ppm_t_particles_d, &
         &   ppm_t_neighlist_s,ppm_t_neighlist_d,particles_stats_s,particles_stats_d, &
         &   ppm_t_part_prop_s,ppm_t_part_prop_d
         USE ppm_module_field_typedef, ONLY : ppm_t_discr_info,ppm_t_field
         USE ppm_module_operator_typedef, ONLY : ppm_t_operator, &
         &   ppm_c_operator,ppm_t_dcop_s,ppm_t_dcop_d
         USE ppm_module_interfaces, ONLY : ppm_t_main_abstr,ppm_t_discr_kind, &
         &   ppm_t_discr_data,ppm_t_operator_,ppm_t_operator_discr_,ppm_t_operator_discr, &
         &   ppm_t_discr_info_,ppm_t_field_,ppm_t_subpatch_,ppm_t_mesh_maplist, &
         &   ppm_t_equi_mesh_,ppm_t_neighlist_s_,ppm_t_neighlist_d_,particles_stats_s_, &
         &   particles_stats_d_,ppm_t_particles_s_,ppm_t_particles_d_, &
         &   ppm_param_part_init_cartesian,ppm_param_part_init_random, &
         &   ppm_mesh_ghosts,ppm_mesh_reqput,ppm_mesh_cartesian, &
         &   ppm_param_length_meshflags,ppm_part_ghosts,ppm_part_partial, &
         &   ppm_part_reqput,ppm_part_areinside,ppm_part_cartesian, &
         &   ppm_part_neighlists,ppm_part_global_index,ppm_param_length_partflags, &
         &   ppm_ppt_ghosts,ppm_ppt_partial,ppm_ppt_reqput,ppm_ppt_map_parts, &
         &   ppm_ppt_map_ghosts,ppm_param_length_pptflags,ppm_ops_inc_ghosts, &
         &   ppm_ops_interp,ppm_ops_iscomputed,ppm_ops_vector,ppm_param_length_opsflags
         USE ppm_module_timestats, ONLY : ppm_t_tstats
      !   USE ppm_module_rmsh
         USE ppm_module_ctrl, ONLY : arg, arg_group, parse_args,      &
         &   disable_help, disable_ctrl,set_ctrl_name,                &
#ifdef __F2003
         &    INTEGER_func, LONGINT_func, SINGLE_func, DOUBLE_func,   &
         &    LOGICAL_func, STRING_func, COMPLEX_func, DCOMPLEX_func, &
         &    INTEGER_ARRAY_func, LONGINT_ARRAY_func,                 &
         &    SINGLE_ARRAY_func, DOUBLE_ARRAY_func,                   &
         &    LOGICAL_ARRAY_func, STRING_ARRAY_func,                  &
         &    COMPLEX_ARRAY_func, DCOMPLEX_ARRAY_func,                &
#endif
         &    reset, add_cmd, ctrl_file_name, break_help,             &
         &    find_arg, find_flag, arg_count,                         &
         &    enabling_flag, disabling_flag, exit_gracefully
         USE ppm_module_options, ONLY : ppm_t_options_op
         USE ppm_module_mapping_typedef

      END MODULE ppm_module_core
