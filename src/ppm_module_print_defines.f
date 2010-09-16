      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_print_defines
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_print_defines
      !!! This module provides utility routines.

         !----------------------------------------------------------------------
         !  Define interface to the print defines routine
         !----------------------------------------------------------------------
         INTERFACE ppm_print_defines
            MODULE PROCEDURE ppm_print_defines
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "util/ppm_print_defines.f"

      END MODULE ppm_module_print_defines
