      !-------------------------------------------------------------------------
      !  Module       :                     ppm_param
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Define global and user-accessible parameters
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_param.h,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.45  2006/09/04 18:34:54  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.44  2006/05/15 12:40:47  pchatela
      !  Added Hammer 3-4-7-12 quadrature rules
      !
      !  Revision 1.43  2005/08/05 09:28:08  ivos
      !  Added STS ode scheme.
      !
      !  Revision 1.42  2005/03/10 01:49:19  ivos
      !  Added order parameters for GMM.
      !
      !  Revision 1.41  2005/01/13 13:23:26  michaebe
      !  added m3p6
      !
      !  Revision 1.40  2005/01/13 13:19:08  ivos
      !  Added ppm_param_assign_user_defined.
      !
      !  Revision 1.39  2004/10/29 16:01:36  kotsalie
      !  RED BLACK SOR PARAMETER
      !
      !  Revision 1.38  2004/10/13 16:43:06  davidch
      !  added constants for 3d sph kernels
      !
      !  Revision 1.37  2004/10/01 15:14:37  hiebers
      !  added ppm_param_ode_rk4
      !
      !  Revision 1.36  2004/09/28 14:17:33  kotsalie
      !  Added parameter for 4th order
      !
      !  Revision 1.35  2004/09/24 15:03:54  ivos
      !  Added params for slab decomposition types.
      !
      !  Revision 1.34  2004/09/22 18:21:53  kotsalie
      !  Added parameters for the mg
      !
      !  Revision 1.33  2004/09/22 10:32:54  ivos
      !  Added ppm_param_tree_bin, ppm_param_tree_quad, ppm_param_tree_oct.
      !
      !  Revision 1.32  2004/09/17 12:05:09  ivos
      !  cosmetics.
      !
      !  Revision 1.31  2004/08/31 12:15:18  ivos
      !  Added ppm_param_id_internal and ppm_param_id_user.
      !
      !  Revision 1.30  2004/08/13 15:29:19  michaebe
      !  added another ode scheme
      !
      !  Revision 1.29  2004/08/12 10:45:26  michaebe
      !  inserted kernel paramters
      !
      !  Revision 1.28  2004/07/29 12:56:59  hiebers
      !  Added ppm_param_kernel_*sph2d*
      !
      !  Revision 1.27  2004/07/27 09:05:32  michaebe
      !  added user-accessbile ode parameters
      !
      !  Revision 1.26  2004/07/27 07:06:04  ivos
      !  Added BEM parameters.
      !
      !  Revision 1.25  2004/07/23 12:54:04  ivos
      !  Added parameters for laplacian kernels in 2d and 3d.
      !
      !  Revision 1.24  2004/07/20 16:43:51  ivos
      !  Added data update and loadbalance parameters.
      !
      !  Revision 1.23  2004/07/16 11:43:04  ivos
      !  Added ppm_param_topo_field and ppm_param_topo_part.
      !
      !  Revision 1.22  2004/05/11 14:53:46  ivos
      !  Added ppm_param_io_single and ppm_param_io_double.
      !
      !  Revision 1.21  2004/05/11 13:31:39  oingo
      !  Added connection options
      !
      !  Revision 1.20  2004/05/06 07:44:08  ivos
      !  Added I/O parameters.
      !
      !  Revision 1.19  2004/04/13 15:18:56  oingo
      !  Added ppm_param_decomp_null
      !  Removed ppm_param_map_ring
      !
      !  Revision 1.18  2004/04/07 11:48:26  oingo
      !  Added ppm_param_map_ring
      !
      !  Revision 1.17  2004/04/05 11:01:55  ivos
      !  Added ppm_param_pop_add and ppm_param_pop_replace.
      !
      !  Revision 1.16  2004/04/02 15:24:47  ivos
      !  Added ppm_param_map_init.
      !
      !  Revision 1.15  2004/04/01 14:11:06  ivos
      !  added ppm_param_map_ghost_init.
      !
      !  Revision 1.14  2004/02/25 14:45:36  ivos
      !  Added new METIS assignment types ppm_param_assign_nodal_cut, 
      !  _nodal_comm,
      !  _dual_cut and _dual_comm.
      !
      !  Revision 1.13  2004/02/18 14:53:58  walther
      !  Now receiving the ghost particles is named ppm_param_map_ghost_get and
      !  sending ghost values back to their host processor is 
      !  ppm_param_map_ghost_put.
      !
      !  Revision 1.12  2004/02/18 14:30:54  walther
      !  Renamed parameters: ppm_param_map_update is now ppm_param_map_partial,
      !  ppm_param_map_ghost is used in both the field and particle ghost 
      !  mapping
      !  and refers to the ghost VALUE; therefore I added the 
      !  ppm_param_map_ghostpart
      !  as the parameter to map the ghost POSITION (particles only).
      !
      !  Revision 1.11  2004/02/12 17:07:34  ivos
      !  Added parameters to define sub2proc assignment types.
      !
      !  Revision 1.10  2004/02/12 08:33:27  michaebe
      !  inserted a _io_ for the append,replace,overwrite params to make them
      !  more style conforming
      !
      !  Revision 1.9  2004/02/10 12:32:08  michaebe
      !  added some miscellaneous params (overwrite,replace,append,insert)
      !
      !  Revision 1.8  2004/02/06 15:43:29  walther
      !  Added ppm_param_map_use_symmetry.
      !
      !  Revision 1.7  2004/02/05 09:02:59  walther
      !  Changed header.
      !
      !  Revision 1.6  2004/02/04 17:21:41  ivos
      !  Added mesh operation params ppm_param_mesh_refine and
      !  ppm_param_mesh_coarsen.
      !
      !  Revision 1.5  2004/01/29 13:56:12  ivos
      !  Added decomposition methods for meshes.
      !
      !  Revision 1.4  2004/01/09 15:56:03  ivos
      !  Added error level constants.
      !
      !  Revision 1.3  2003/12/12 17:13:43  ivos
      !  Added map type ppm_param_map_cancel = 10.
      !
      !  Revision 1.2  2003/11/24 15:02:54  ivos
      !  Added map type ppm_param_map_partial = 10
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define the precision
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_kind_double     = KIND(1.0D0)
      INTEGER, PARAMETER :: ppm_kind_single     = KIND(1.0E0)
      INTEGER, PARAMETER :: ppm_integer         = 13
      INTEGER, PARAMETER :: ppm_logical         = 17
      INTEGER, PARAMETER :: ppm_char            = 256

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
      INTEGER, PARAMETER :: ppm_param_update_overwrite    = 1 
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
      INTEGER, PARAMETER :: ppm_param_decomp_null         = 9
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
      INTEGER, PARAMETER :: ppm_param_assign_nodal_cut    = 2
      INTEGER, PARAMETER :: ppm_param_assign_nodal_comm   = 3
      INTEGER, PARAMETER :: ppm_param_assign_dual_cut     = 4
      INTEGER, PARAMETER :: ppm_param_assign_dual_comm    = 5
      INTEGER, PARAMETER :: ppm_param_assign_user_defined = 6
      
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

      !------------------------------------------------------------------------
      !  Define equation ppm solves 
      !------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_eq_poisson          = 1 

      !------------------------------------------------------------------------
      !  Define order for finite difference
      !------------------------------------------------------------------------ 
      INTEGER, PARAMETER :: ppm_param_order_1             = 1
      INTEGER, PARAMETER :: ppm_param_order_2             = 2
      INTEGER, PARAMETER :: ppm_param_order_3             = 3
      INTEGER, PARAMETER :: ppm_param_order_4             = 4

      !------------------------------------------------------------------------
      !  Define smoother
      !------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_smooth_rbsor        = 1 

      !-------------------------------------------------------------------------
      !  Define boundary conditions
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_bcdef_periodic      = 1
      INTEGER, PARAMETER :: ppm_param_bcdef_freespace     = 2
      INTEGER, PARAMETER :: ppm_param_bcdef_symmetry      = 3
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

      !-------------------------------------------------------------------------
      !  error severity levels
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_error_notice              = -4
      INTEGER, PARAMETER :: ppm_error_warning             = -3
      INTEGER, PARAMETER :: ppm_error_error               = -2
      INTEGER, PARAMETER :: ppm_error_fatal               = -1

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

      !-------------------------------------------------------------------------
      !  RMSH parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_rmsh_kernel_bsp2    = 1
      INTEGER, PARAMETER :: ppm_param_rmsh_kernel_mp4     = 2
      INTEGER, PARAMETER :: ppm_param_rmsh_kernel_m3p6    = 3
      
