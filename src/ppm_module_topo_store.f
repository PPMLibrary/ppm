      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_topo_store
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

      MODULE ppm_module_topo_store
      !!! This module provides the routine that stores newly created topologies.
         !----------------------------------------------------------------------
         !  Define interfaces to topology store routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_store
            MODULE PROCEDURE ppm_topo_store_s
            MODULE PROCEDURE ppm_topo_store_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_store.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_store.f"
#undef __KIND

      END MODULE ppm_module_topo_store
