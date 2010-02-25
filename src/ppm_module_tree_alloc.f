      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_tree_alloc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 that allocate tree data structures.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_tree_alloc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/22 10:32:06  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
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
#include "ppm_tree_alloc.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_tree_alloc.f"
#undef __KIND
#undef __TYPE

#define __TYPE __DECOMP
#define __KIND __SINGLE_PRECISION
#include "ppm_tree_alloc.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_tree_alloc.f"
#undef __KIND
#undef __TYPE

      END MODULE ppm_module_tree_alloc
