      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_tree_cutpos
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_tree_cutpos
      !!! This module provides the
      !!! routine that determines the optimal position where to cut a box.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_cutpos
            MODULE PROCEDURE ppm_tree_cutpos_s
            MODULE PROCEDURE ppm_tree_cutpos_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree_cutpos.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree_cutpos.f"
#undef __KIND

      END MODULE ppm_module_tree_cutpos
