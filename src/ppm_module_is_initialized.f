      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_is_initialized
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_is_initialized
      !!! This module provides the routines
      !!! callable from the outside.
         !----------------------------------------------------------------------
         !  Define interface to global status inquiry routine
         !----------------------------------------------------------------------
         INTERFACE ppm_is_initialized
            MODULE PROCEDURE ppm_is_initialized
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "util/ppm_is_initialized.f"

      END MODULE ppm_module_is_initialized
