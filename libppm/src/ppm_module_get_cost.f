      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_get_cost
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 concerned with load balancing. They are callable by
      !                 the user.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_get_cost.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:37  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_get_cost

         !----------------------------------------------------------------------
         !  Define interface to cost inquiry routine
         !----------------------------------------------------------------------
         INTERFACE ppm_get_cost
            MODULE PROCEDURE ppm_get_cost_s
            MODULE PROCEDURE ppm_get_cost_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_get_cost.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_get_cost.f"
#undef __KIND

      END MODULE ppm_module_get_cost
