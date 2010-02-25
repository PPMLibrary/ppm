      !-------------------------------------------------------------------------
      !  Module       :                   ppm_module_tree
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the tree
      !                 routine.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_tree.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/08/31 11:24:29  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.2  2004/12/03 17:14:03  ivos
      !  Added lhbx_cut and lpdx_cut.
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

      MODULE ppm_module_tree

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
#include "ppm_tree.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_tree.f"
#undef __KIND
#undef __TYPE

#define __TYPE __DECOMP
#define __KIND __SINGLE_PRECISION
#include "ppm_tree.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_tree.f"
#undef __KIND
#undef __TYPE

      END MODULE ppm_module_tree
