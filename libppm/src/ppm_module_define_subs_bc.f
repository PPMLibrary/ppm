      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_define_subs_bc
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

      MODULE ppm_module_define_subs_bc
      !!! This module provides the decomposition routines.
         !----------------------------------------------------------------------
         !  Define interface to boundary condition definition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_define_subs_bc
            MODULE PROCEDURE define_subsbc_s
            MODULE PROCEDURE define_subsbc_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_define_subs_bc.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_define_subs_bc.f"
#undef __KIND

      END MODULE ppm_module_define_subs_bc
