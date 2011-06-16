      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_scale_domain
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

      MODULE ppm_module_scale_domain
      !!! This module provides routines for scaling the domain
         !----------------------------------------------------------------------
         !  Define interfaces to the domain scaling routine
         !----------------------------------------------------------------------
         INTERFACE ppm_scale_domain
            MODULE PROCEDURE ppm_scale_domain_s
            MODULE PROCEDURE ppm_scale_domain_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "util/ppm_scale_domain.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_scale_domain.f"
#undef __KIND

      END MODULE ppm_module_scale_domain
