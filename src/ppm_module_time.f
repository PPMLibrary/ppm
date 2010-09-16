      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                  ppm_module_time
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_time
      !!! This module provides the routines
      !!! concerned with load balancing. They are callable by the user.
         !----------------------------------------------------------------------
         !  Define interfaces to the timing routine
         !----------------------------------------------------------------------
         INTERFACE ppm_time
            MODULE PROCEDURE ppm_time_s
            MODULE PROCEDURE ppm_time_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "util/ppm_time.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_time.f"
#undef __KIND

      END MODULE ppm_module_time
