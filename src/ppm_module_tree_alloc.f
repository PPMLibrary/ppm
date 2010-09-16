      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_tree_alloc
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

      MODULE ppm_module_tree_alloc
      !!! This module provides the routines
      !!! that allocate tree data structures.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_alloc
            MODULE PROCEDURE ppm_tree_alloc_ts
            MODULE PROCEDURE ppm_tree_alloc_td
            MODULE PROCEDURE ppm_tree_alloc_ds
            MODULE PROCEDURE ppm_tree_alloc_dd
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __TYPE __TREE
#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree_alloc.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree_alloc.f"
#undef __KIND
#undef __TYPE

#define __TYPE __DECOMP
#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree_alloc.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree_alloc.f"
#undef __KIND
#undef __TYPE

      END MODULE ppm_module_tree_alloc
