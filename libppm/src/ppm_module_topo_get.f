      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_topo_get
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

      MODULE ppm_module_topo_get
      !!! This module provides the routines to read an internally stored topology.
         !----------------------------------------------------------------------
         !  Define interfaces to topology store routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_get
            MODULE PROCEDURE ppm_topo_get
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "topo/ppm_topo_get.f"

      END MODULE ppm_module_topo_get
