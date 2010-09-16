      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_topo_check
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

      MODULE ppm_module_topo_check
      !!! This module provides the routines to check topologies -
      !!! callable from outside.
         !----------------------------------------------------------------------
         !  Define interfaces to topology checking routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_check
            MODULE PROCEDURE ppm_topo_check_s
            MODULE PROCEDURE ppm_topo_check_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_check.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_check.f"
#undef __KIND

      END MODULE ppm_module_topo_check
