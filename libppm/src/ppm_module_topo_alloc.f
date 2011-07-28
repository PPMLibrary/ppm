      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                  ppm_module_topo_alloc
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

      MODULE ppm_module_topo_alloc
      !!! This module contains the topology allocation and deallocation
      !!! subroutines.

         !----------------------------------------------------------------------
         !  Define interfaces to topology allocation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_alloc
            MODULE PROCEDURE ppm_topo_alloc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to topology deallocation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_dealloc
            MODULE PROCEDURE ppm_topo_dealloc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "topo/ppm_topo_alloc.f"

#include "topo/ppm_topo_dealloc.f"

      END MODULE ppm_module_topo_alloc
