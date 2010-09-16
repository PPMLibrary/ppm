      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_collapse_list
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_util_collapse_list
      !!! This module provides the utility routines.
         !----------------------------------------------------------------------
         !  Define interface to the list inversion routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_collapse_list
            MODULE PROCEDURE ppm_util_collapse_list
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "util/ppm_util_collapse_list.f"

      END MODULE ppm_module_util_collapse_list
