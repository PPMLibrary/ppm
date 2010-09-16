      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                   ppm_module_log
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_log
      !!! This module provides the utility routines.
         !----------------------------------------------------------------------
         !  Define interface to the log write rountine
         !----------------------------------------------------------------------
         INTERFACE ppm_log
            MODULE PROCEDURE ppm_log
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "util/ppm_log.f"

      END MODULE ppm_module_log
