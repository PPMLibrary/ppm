      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_check_topoid
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_check_id
      !!! This module provides the utility routines,
      !!! which check topology and mesh IDs for their validity.
         !----------------------------------------------------------------------
         !  Define interface to the topoid check routine
         !----------------------------------------------------------------------
         INTERFACE ppm_check_topoid
            MODULE PROCEDURE ppm_check_topoid
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to the meshid check routine
         !----------------------------------------------------------------------
         INTERFACE ppm_check_meshid
            MODULE PROCEDURE ppm_check_meshid
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "topo/ppm_check_topoid.f"

#include "topo/ppm_check_meshid.f"

      END MODULE ppm_module_check_id
