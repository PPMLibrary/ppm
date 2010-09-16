      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_topo_subs2proc
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

      MODULE ppm_module_topo_subs2proc
      !!! This module provides the topology subdomain to processor assignment
      !!! routines.
         !----------------------------------------------------------------------
         !  Define interface to internal sub-to-proc assignment routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_subs2proc
            MODULE PROCEDURE ppm_topo_subs2proc_s
            MODULE PROCEDURE ppm_topo_subs2proc_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_subs2proc.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_subs2proc.f"
#undef __KIND

      END MODULE ppm_module_topo_subs2proc
