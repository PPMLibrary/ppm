      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                     ppm_param
      !
      !  Purpose      : Define global and user-accessible parameters
      !
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

      !-------------------------------------------------------------------------
      !  Define the precision
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_kind_double     = KIND(1.0D0) !SELECTED_REAL_KIND(P=15, R=307)
      INTEGER, PARAMETER :: ppm_kind_single     = KIND(1.0E0) !SELECTED_REAL_KIND(P= 6, R= 37)
      INTEGER, PARAMETER :: ppm_integer         = 13
      INTEGER, PARAMETER :: ppm_logical         = 17
      INTEGER, PARAMETER :: ppm_char            = 256
      INTEGER, PARAMETER :: ppm_kind_int32      = 4 !SELECTED_INT_KIND(6)
      INTEGER, PARAMETER :: ppm_kind_int64      = 8 !SELECTED_INT_KIND(12)
      !-------------------------------------------------------------------------
      !  Define the supported types
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_type_none                 = 0
      INTEGER, PARAMETER :: ppm_type_int                  = 1
      INTEGER, PARAMETER :: ppm_type_longint              = 2
      INTEGER, PARAMETER :: ppm_type_real                 = 3
      INTEGER, PARAMETER :: ppm_type_real_single          = 4
      INTEGER, PARAMETER :: ppm_type_comp                 = 5
      INTEGER, PARAMETER :: ppm_type_comp_single          = 6
      INTEGER, PARAMETER :: ppm_type_logical              = 7
      INTEGER, PARAMETER :: ppm_type_char                 = 8
      !-------------------------------------------------------------------------
      !  Mapping options
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_map_pop             = 1
      INTEGER, PARAMETER :: ppm_param_map_push            = 2
      INTEGER, PARAMETER :: ppm_param_map_global          = 3
      INTEGER, PARAMETER :: ppm_param_map_send            = 4
      INTEGER, PARAMETER :: ppm_param_map_ghost_get       = 5
      INTEGER, PARAMETER :: ppm_param_map_ghost_put       = 6
      INTEGER, PARAMETER :: ppm_param_map_partial         = 7
      INTEGER, PARAMETER :: ppm_param_map_remap           = 8
      INTEGER, PARAMETER :: ppm_param_map_cancel          = 9
      INTEGER, PARAMETER :: ppm_param_map_init            = 10
      !-------------------------------------------------------------------------
      !  Connection options
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_connect_distribute  = 1
      INTEGER, PARAMETER :: ppm_param_connect_send        = 2
      INTEGER, PARAMETER :: ppm_param_connect_prune       = 3
      !-------------------------------------------------------------------------
      !  Define data receive options (map_pop)
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_pop_replace         = 1
      INTEGER, PARAMETER :: ppm_param_pop_add             = 2
      !-------------------------------------------------------------------------
      !  Define ppm-internal PP interaction kernels
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_kernel_laplace2d_p2 =  1
      INTEGER, PARAMETER :: ppm_param_kernel_laplace3d_p2 =  2
      INTEGER, PARAMETER :: ppm_param_kernel_sph2d_p2     =  3
      INTEGER, PARAMETER :: ppm_param_kernel_dx_sph2d_p2  =  4
      INTEGER, PARAMETER :: ppm_param_kernel_dy_sph2d_p2  =  5
      INTEGER, PARAMETER :: ppm_param_kernel_ddx_sph2d_p2 =  6
      INTEGER, PARAMETER :: ppm_param_kernel_ddy_sph2d_p2 =  7
      INTEGER, PARAMETER :: ppm_param_kernel_dxdy_sph2d_p2=  8
      INTEGER, PARAMETER :: ppm_param_kernel_fast3d       =  9
      INTEGER, PARAMETER :: ppm_param_kernel_fast3d_dx    = 10
      INTEGER, PARAMETER :: ppm_param_kernel_fast3d_dy    = 11
      INTEGER, PARAMETER :: ppm_param_kernel_fast3d_dz    = 12
      INTEGER, PARAMETER :: ppm_param_kernel_fast3d_lap   = 13
      !-------------------------------------------------------------------------
      !  Define decision heuristics for dynamic load balancing
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_loadbal_sar         = 1
      !-------------------------------------------------------------------------
      !  Define statistics update methods
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_update_replace      = 1
      INTEGER, PARAMETER :: ppm_param_update_average      = 2
      INTEGER, PARAMETER :: ppm_param_update_expfavg      = 3
      !-------------------------------------------------------------------------
      !  Define domain decomposition techniques
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_decomp_tree         = 1
      INTEGER, PARAMETER :: ppm_param_decomp_pruned_cell  = 2
      INTEGER, PARAMETER :: ppm_param_decomp_bisection    = 3
      INTEGER, PARAMETER :: ppm_param_decomp_xpencil      = 4
      INTEGER, PARAMETER :: ppm_param_decomp_ypencil      = 5
      INTEGER, PARAMETER :: ppm_param_decomp_zpencil      = 6
      INTEGER, PARAMETER :: ppm_param_decomp_cuboid       = 7
      INTEGER, PARAMETER :: ppm_param_decomp_user_defined = 8
      INTEGER, PARAMETER :: ppm_param_decomp_xy_slab      = 10
      INTEGER, PARAMETER :: ppm_param_decomp_xz_slab      = 11
      INTEGER, PARAMETER :: ppm_param_decomp_yz_slab      = 12
      INTEGER, PARAMETER :: ppm_param_decomp_cartesian    = 13
      !-------------------------------------------------------------------------
      !  Define tree types
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_tree_bin            = 1
      INTEGER, PARAMETER :: ppm_param_tree_quad           = 2
      INTEGER, PARAMETER :: ppm_param_tree_oct            = 3
      !-------------------------------------------------------------------------
      !  Define subdomain-to-processor assignment schemes
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_assign_internal     = 1
      INTEGER, PARAMETER :: ppm_param_assign_metis_cut    = 2
      INTEGER, PARAMETER :: ppm_param_assign_metis_comm   = 3
      INTEGER, PARAMETER :: ppm_param_assign_user_defined = 4
      !-------------------------------------------------------------------------
      !  Define particle-mesh schemes
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_assign_ngp          = 0
      INTEGER, PARAMETER :: ppm_param_assign_cic          = 1
      INTEGER, PARAMETER :: ppm_param_assign_tcs          = 2
      INTEGER, PARAMETER :: ppm_param_assign_mp4          = 3
      INTEGER, PARAMETER :: ppm_param_assign_m3p6         = 4
      !-------------------------------------------------------------------------
      !  Define mesh operations
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_mesh_refine         = 0
      INTEGER, PARAMETER :: ppm_param_mesh_coarsen        = 1
      !-------------------------------------------------------------------------
      !  Define boundary conditions
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_bcdef_freespace     = 0
      INTEGER, PARAMETER :: ppm_param_bcdef_periodic      = 1
      INTEGER, PARAMETER :: ppm_param_bcdef_symmetry      = 2
      INTEGER, PARAMETER :: ppm_param_bcdef_antisymmetry  = 3
      INTEGER, PARAMETER :: ppm_param_bcdef_neumann       = 4
      INTEGER, PARAMETER :: ppm_param_bcdef_dirichlet     = 5
      INTEGER, PARAMETER :: ppm_param_bcdef_robin         = 6
      !-------------------------------------------------------------------------
      !  memory allocation options
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_alloc_fit           = 1
      INTEGER, PARAMETER :: ppm_param_alloc_fit_preserve  = 2
      INTEGER, PARAMETER :: ppm_param_alloc_grow          = 3
      INTEGER, PARAMETER :: ppm_param_alloc_grow_preserve = 4
      INTEGER, PARAMETER :: ppm_param_dealloc             = 5
      !-------------------------------------------------------------------------
      !  I/O parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_io_read             = 1
      INTEGER, PARAMETER :: ppm_param_io_write            = 2
      INTEGER, PARAMETER :: ppm_param_io_read_write       = 3
      INTEGER, PARAMETER :: ppm_param_io_replace          = 4
      INTEGER, PARAMETER :: ppm_param_io_append           = 5
      INTEGER, PARAMETER :: ppm_param_io_ascii            = 6
      INTEGER, PARAMETER :: ppm_param_io_binary           = 7
      INTEGER, PARAMETER :: ppm_param_io_distributed      = 8
      INTEGER, PARAMETER :: ppm_param_io_centralized      = 9
      INTEGER, PARAMETER :: ppm_param_io_same             = 10
      INTEGER, PARAMETER :: ppm_param_io_root             = 11
      INTEGER, PARAMETER :: ppm_param_io_sum              = 12
      INTEGER, PARAMETER :: ppm_param_io_split            = 13
      INTEGER, PARAMETER :: ppm_param_io_concat           = 14
      INTEGER, PARAMETER :: ppm_param_io_single           = 15
      INTEGER, PARAMETER :: ppm_param_io_double           = 16
      INTEGER, PARAMETER :: ppm_param_io_hdf5             = 17
      !-------------------------------------------------------------------------
      !  error severity levels
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_error_notice              = -1
      INTEGER, PARAMETER :: ppm_error_warning             = -2
      INTEGER, PARAMETER :: ppm_error_error               = -3
      INTEGER, PARAMETER :: ppm_error_fatal               = -4
      !-------------------------------------------------------------------------
      !  Toplogy parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_topo_part           = 1
      INTEGER, PARAMETER :: ppm_param_topo_field          = 2
      INTEGER, PARAMETER :: ppm_param_id_internal         = 3
      INTEGER, PARAMETER :: ppm_param_id_user             = 4
      !-------------------------------------------------------------------------
      !  Misc parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_undefined           = -1
      INTEGER, PARAMETER :: ppm_param_success             = 0
      INTEGER, PARAMETER :: ppm_param_topo_undefined      = -1
      !-------------------------------------------------------------------------
      !  RMSH parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_rmsh_kernel_bsp2    = 1
      INTEGER, PARAMETER :: ppm_param_rmsh_kernel_mp4     = 2
      INTEGER, PARAMETER :: ppm_param_rmsh_kernel_m3p6    = 3
      !-------------------------------------------------------------------------
      !  Operators
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_op_fd               = 1
      INTEGER, PARAMETER :: ppm_param_op_pse              = 2
      INTEGER, PARAMETER :: ppm_param_op_dcpse            = 3

