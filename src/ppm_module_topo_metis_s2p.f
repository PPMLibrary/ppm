      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_topo_metis_s2p
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

      MODULE ppm_module_topo_metis_s2p
      !!! This module provides the decomposition routines.
         !----------------------------------------------------------------------
         !  Define interface to METIS-based sub-to-proc assignment
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_metis_s2p
            MODULE PROCEDURE ppm_topo_metis_s2p_s
            MODULE PROCEDURE ppm_topo_metis_s2p_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "topo/ppm_topo_metis_s2p.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "topo/ppm_topo_metis_s2p.f"
#undef __KIND

      END MODULE ppm_module_topo_metis_s2p
