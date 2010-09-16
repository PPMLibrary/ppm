      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_substart
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

      MODULE ppm_module_substart
      !!! This module provides utility routines.
         !----------------------------------------------------------------------
         !  Define interfaces to substart
         !----------------------------------------------------------------------
         INTERFACE substart
            MODULE PROCEDURE substart_s
            MODULE PROCEDURE substart_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __KIND __SINGLE_PRECISION
#include "util/substart.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/substart.f"
#undef __KIND

      END MODULE ppm_module_substart
