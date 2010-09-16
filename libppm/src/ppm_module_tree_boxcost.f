      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_tree_boxcost
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

      MODULE ppm_module_tree_boxcost
      !!! This module provides the
      !!! routine that computes the costs of tree boxes.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_boxcost
            MODULE PROCEDURE ppm_tree_boxcost_s
            MODULE PROCEDURE ppm_tree_boxcost_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree_boxcost.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree_boxcost.f"
#undef __KIND

      END MODULE ppm_module_tree_boxcost
