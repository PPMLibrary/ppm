      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_error
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_error
      !!! This module provides the error routines.
         !----------------------------------------------------------------------
         !  Header file for error codes
         !----------------------------------------------------------------------
         INCLUDE 'ppm_error.h'

         !----------------------------------------------------------------------
         !  Define interface to the error handler routine
         !----------------------------------------------------------------------
         INTERFACE ppm_error
            MODULE PROCEDURE ppm_error
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "util/ppm_error.f"

      END MODULE ppm_module_error
