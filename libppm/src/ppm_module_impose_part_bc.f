      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_impose_part_bc
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2

      MODULE ppm_module_impose_part_bc
      !!! This module includes the source code user callable
      !!! ppm_impose_part_bc routine to map the particles back
      !!! inside the computational domain for periodic BC.
         !----------------------------------------------------------------------
         !  Define interface to ppm_map_field_ghost
         !----------------------------------------------------------------------
         INTERFACE ppm_impose_part_bc
             MODULE PROCEDURE ppm_impose_part_bc_s
             MODULE PROCEDURE ppm_impose_part_bc_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "map/ppm_impose_part_bc.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_impose_part_bc.f"
#undef __KIND

      END MODULE ppm_module_impose_part_bc
