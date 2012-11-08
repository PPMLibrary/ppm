      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_loadbal
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


#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4

      MODULE ppm_module_loadbal
      !!! This module contains all routines needed for
      !!! dynamic load balancing.
      !!!
      !!! Currently only the SOR (stop at rise) method is provided as the
      !!! heuristic to estimate the time when the domain needs to be
      !!! redecomposed.

         USE ppm_module_topo_typedef
         USE ppm_module_particles_typedef
         USE ppm_module_operator_typedef
         USE ppm_module_interfaces
         USE ppm_module_data_loadbal
         USE ppm_module_data
         USE ppm_module_alloc
         !----------------------------------------------------------------------
         !  Define interface to load balance inquiry routine
         !----------------------------------------------------------------------
         INTERFACE ppm_loadbal_sendsub
            MODULE PROCEDURE loadbal_sendsub
         END INTERFACE

         INTERFACE ppm_loadbal_recvsub
            MODULE PROCEDURE loadbal_recvsub
         END INTERFACE
          INTERFACE ppm_loadbal_map_subpart
             MODULE PROCEDURE loadbal_map_subpart_s
             MODULE PROCEDURE loadbal_map_subpart_d
          END INTERFACE
         !----------------------------------------------------------------------
         !  Define interface to load balance inquiry routine
         !----------------------------------------------------------------------
!         INTERFACE ppm_loadbal_inquire
!            MODULE PROCEDURE loadbal_inq_s
!            MODULE PROCEDURE loadbal_inq_d
!         END INTERFACE
!
!         INTERFACE ppm_loadbal_inquire_sar
!            MODULE PROCEDURE loadbal_inquire_sar_s
!            MODULE PROCEDURE loadbal_inquire_sar_d
!         END INTERFACE

!         INTERFACE ppm_loadbal_inquire_dlb
!            MODULE PROCEDURE loadbal_inquire_dlb_s
!            MODULE PROCEDURE loadbal_inquire_dlb_d
!         END INTERFACE

!         INTERFACE ppm_loadbal_do_dlb
!            MODULE PROCEDURE loadbal_do_dlb
!         END INTERFACE
         !----------------------------------------------------------------------
         !  Define interface to processor speed estimator
         !----------------------------------------------------------------------
         INTERFACE ppm_estimate_procspeed
            MODULE PROCEDURE est_procspeed_s
            MODULE PROCEDURE est_procspeed_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to cost inquiry routine
         !----------------------------------------------------------------------
         INTERFACE ppm_get_cost
            MODULE PROCEDURE ppm_get_cost_s
            MODULE PROCEDURE ppm_get_cost_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to decomposition cost update routine
         !----------------------------------------------------------------------
         INTERFACE ppm_set_decomp_cost
            MODULE PROCEDURE set_dcost_s
            MODULE PROCEDURE set_dcost_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to processor speed setter
         !----------------------------------------------------------------------
         INTERFACE ppm_set_proc_speed
            MODULE PROCEDURE ppm_set_proc_speed_s
            MODULE PROCEDURE ppm_set_proc_speed_d
         END INTERFACE
         
         !----------------------------------------------------------------------
         !  Balancing circuit DLB model
         !----------------------------------------------------------------------
         INTERFACE ppm_loadbal_bc
            MODULE PROCEDURE loadbal_bc_s
            MODULE PROCEDURE loadbal_bc_d
         END INTERFACE
         !----------------------------------------------------------------------
         !  Choose a sequence of subdomains to send to an underloaded processor
         !----------------------------------------------------------------------
         INTERFACE ppm_loadbal_choose_sub
            MODULE PROCEDURE loadbal_choose_sub_s
            MODULE PROCEDURE loadbal_choose_sub_d
         END INTERFACE
         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "loadbal/ppm_loadbal_sendsub.f"
#include "loadbal/ppm_loadbal_recvsub.f"

#define __KIND __SINGLE_PRECISION
#include "loadbal/ppm_estimate_proc_speed.f"
#include "loadbal/ppm_loadbal_bc.f"
#include "loadbal/ppm_loadbal_choose_sub.f"
#include "loadbal/ppm_loadbal_map_subpart.f"
#include "loadbal/ppm_set_decomp_cost.f"
#include "loadbal/ppm_set_proc_speed.f"
#include "loadbal/ppm_get_cost.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "loadbal/ppm_estimate_proc_speed.f"
#include "loadbal/ppm_loadbal_bc.f"
#include "loadbal/ppm_loadbal_choose_sub.f"
#include "loadbal/ppm_loadbal_map_subpart.f"
#include "loadbal/ppm_set_decomp_cost.f"
#include "loadbal/ppm_set_proc_speed.f"
#include "loadbal/ppm_get_cost.f"
#undef __KIND

      END MODULE ppm_module_loadbal
