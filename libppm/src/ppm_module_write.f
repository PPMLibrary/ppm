      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                  ppm_module_write
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_write
      !!! This module provides the utility write routine.
         !----------------------------------------------------------------------
         !  Define interface to the write routine
         !----------------------------------------------------------------------
         INTERFACE ppm_write
            MODULE PROCEDURE ppm_write
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "io/ppm_write.f"

      END MODULE ppm_module_write
