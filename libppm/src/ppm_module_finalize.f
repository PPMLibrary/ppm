      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_finalize
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_finalize
      !!! This module provides the finalization routine callable from outside
         !----------------------------------------------------------------------
         !  Define interface to ppm_finalize
         !----------------------------------------------------------------------
         INTERFACE ppm_finalize
            MODULE PROCEDURE ppm_finalize
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_finalize.f"

      END MODULE ppm_module_finalize
