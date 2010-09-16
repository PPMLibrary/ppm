      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_topo_box2subs
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

      MODULE ppm_module_topo_box2subs
      !!! This module provides the routine that converts tree boxes to subs.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_box2subs
            MODULE PROCEDURE ppm_topo_box2subs_s
            MODULE PROCEDURE ppm_topo_box2subs_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_box2subs.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_box2subs.f"
#undef __KIND

      END MODULE ppm_module_topo_box2subs
