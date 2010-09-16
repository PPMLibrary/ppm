      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_substop
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

      MODULE ppm_module_substop
      !!! This module provides utility routines.
         !----------------------------------------------------------------------
         !  Define interfaces to substop
         !----------------------------------------------------------------------
         INTERFACE substop
            MODULE PROCEDURE substop_s
            MODULE PROCEDURE substop_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __KIND __SINGLE_PRECISION
#include "util/substop.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/substop.f"
#undef __KIND

      END MODULE ppm_module_substop
