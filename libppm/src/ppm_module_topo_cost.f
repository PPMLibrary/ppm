      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_topo_cost
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

      MODULE ppm_module_topo_cost
      !!! This module provides the subdomain cost calculation routine.
         !----------------------------------------------------------------------
         !  Define interface to topology cost computation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_cost
            MODULE PROCEDURE ppm_topo_cost_s
            MODULE PROCEDURE ppm_topo_cost_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_cost.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_cost.f"
#undef __KIND

      END MODULE ppm_module_topo_cost
