      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                  ppm_module_init
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_init
      !!! This module provides the init routine - callable from the outside.
         !----------------------------------------------------------------------
         !  Define interface to ppm_init 
         !----------------------------------------------------------------------
         INTERFACE ppm_init
            MODULE PROCEDURE ppm_init
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_init.f"

      END MODULE ppm_module_init
