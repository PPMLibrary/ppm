      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_topo_copy
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

      MODULE ppm_module_topo_copy
      !!! This module provides the routines to copy topology structures
         !----------------------------------------------------------------------
         !  Define interfaces to topology copy routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_copy
            MODULE PROCEDURE ppm_topo_copy
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "topo/ppm_topo_copy.f"

      END MODULE ppm_module_topo_copy
