      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_find_neigh
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

      MODULE ppm_module_find_neigh
      !!! This module provides the decomposition routines.
         !----------------------------------------------------------------------
         !  Define interface to neighbor finding routine
         !----------------------------------------------------------------------
         INTERFACE ppm_find_neigh
            MODULE PROCEDURE ppm_find_neigh_s
            MODULE PROCEDURE ppm_find_neigh_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_find_neigh.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_find_neigh.f"
#undef __KIND

      END MODULE ppm_module_find_neigh
