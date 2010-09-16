      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_mktopo
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

      MODULE ppm_module_mktopo
      !!! This module provides the interface to the routines that create
      !!! topologies.
         !----------------------------------------------------------------------
         !  Define interfaces to the main topology routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_mktopo
            MODULE PROCEDURE ppm_topo_mkpart_s
            MODULE PROCEDURE ppm_topo_mkpart_d
            MODULE PROCEDURE ppm_topo_mkfield_s
            MODULE PROCEDURE ppm_topo_mkfield_d
            MODULE PROCEDURE ppm_topo_mkgeom_s
            MODULE PROCEDURE ppm_topo_mkgeom_d
            MODULE PROCEDURE ppm_topo_mktree_s
            MODULE PROCEDURE ppm_topo_mktree_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_mkpart.f"
#include "topo/ppm_topo_mkfield.f"
#include "topo/ppm_topo_mkgeom.f"
#include "topo/ppm_topo_mktree.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_mkpart.f"
#include "topo/ppm_topo_mkfield.f"
#include "topo/ppm_topo_mkgeom.f"
#include "topo/ppm_topo_mktree.f"
#undef __KIND

      END MODULE ppm_module_mktopo
