      !--*- f90 -*--------------------------------------------------------------
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

      MODULE ppm_module_loadbal
      !!! This module contains all routines needed for
      !!! dynamic load balancing.
      !!!
      !!! Currently only the SOR (stop at rise) method is provided as the
      !!! heuristic to estimate the time when the domain needs to be
      !!! redecomposed.

         USE ppm_module_topo_typedef
         !----------------------------------------------------------------------
         !  Define interface to load balance inquiry routine
         !----------------------------------------------------------------------
         INTERFACE ppm_loadbal_inquire
            MODULE PROCEDURE loadbal_inq_s
            MODULE PROCEDURE loadbal_inq_d
         END INTERFACE

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
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "loadbal/ppm_loadbal_inquire.f"
#include "loadbal/ppm_estimate_proc_speed.f"
#include "loadbal/ppm_get_cost.f"
#include "loadbal/ppm_set_decomp_cost.f"
#include "loadbal/ppm_set_proc_speed.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "loadbal/ppm_loadbal_inquire.f"
#include "loadbal/ppm_estimate_proc_speed.f"
#include "loadbal/ppm_get_cost.f"
#include "loadbal/ppm_set_decomp_cost.f"
#include "loadbal/ppm_set_proc_speed.f"
#undef __KIND


      END MODULE ppm_module_loadbal
