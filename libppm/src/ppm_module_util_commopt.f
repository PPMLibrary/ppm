      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_commopt
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_util_commopt
      !!! This module provides the utility routines.
         !----------------------------------------------------------------------
         !  Define interface to the communication optimization routine
         !----------------------------------------------------------------------
        INTERFACE ppm_util_commopt
            MODULE PROCEDURE ppm_util_commopt
        END INTERFACE

        INTERFACE ppm_util_commopt_cart
            MODULE PROCEDURE ppm_util_commopt_cart
        END INTERFACE
         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "util/ppm_util_commopt.f"

#include "util/ppm_util_commopt_cart.f"
      END MODULE ppm_module_util_commopt
