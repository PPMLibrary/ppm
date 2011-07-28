      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_copy_image_subs
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

      MODULE ppm_module_copy_image_subs
      !!! This module provides the decomposition routines.
         !----------------------------------------------------------------------
         !  Define interface to periodic image copy routine
         !----------------------------------------------------------------------
         INTERFACE ppm_copy_image_subs
            MODULE PROCEDURE copy_imgsubs_s
            MODULE PROCEDURE copy_imgsubs_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_copy_image_subs.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_copy_image_subs.f"
#undef __KIND

      END MODULE ppm_module_copy_image_subs
