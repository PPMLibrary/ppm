      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                   ppm_module_tree
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
#define __TREE             3
#define __DECOMP           4

      MODULE ppm_module_tree
      !!! This module provides the tree routine.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree
            MODULE PROCEDURE ppm_tree_ts
            MODULE PROCEDURE ppm_tree_td
            MODULE PROCEDURE ppm_tree_ds
            MODULE PROCEDURE ppm_tree_dd
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __TYPE __TREE
#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#undef __TYPE

#define __TYPE __DECOMP
#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree.f"
#undef __KIND
#undef __TYPE

      END MODULE ppm_module_tree
