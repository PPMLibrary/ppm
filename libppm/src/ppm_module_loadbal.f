      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_loadbal
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
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
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "loadbal/ppm_loadbal_inquire.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "loadbal/ppm_estimate_proc_speed.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "loadbal/ppm_estimate_proc_speed.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "loadbal/ppm_get_cost.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "loadbal/ppm_get_cost.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "loadbal/ppm_set_decomp_cost.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "loadbal/ppm_set_decomp_cost.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "loadbal/ppm_set_proc_speed.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "loadbal/ppm_set_proc_speed.f"
#undef __KIND

      END MODULE ppm_module_loadbal
