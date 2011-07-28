      !-------------------------------------------------------------------------
      !  Module       :           ppm_module_decomp_pruned_cell
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 decomposition routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_decomp_pruned_cell.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/04 13:35:03  walther
      !  Moved the large local arrays in ppm_decomp_prune_cell to the module.
      !
      !  Revision 1.1  2004/07/26 07:29:33  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
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

      MODULE ppm_module_decomp_pruned_cell
         !----------------------------------------------------------------------
         !  Local work arrays used in ppm_decomp_pruned_cell
         !----------------------------------------------------------------------
         INTEGER , DIMENSION(:), POINTER, PRIVATE :: npbx,npbxg

         !----------------------------------------------------------------------
         !  Define interface to pruned cell decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_pruned_cell
            MODULE PROCEDURE ppm_decomp_pruned_cell_s
            MODULE PROCEDURE ppm_decomp_pruned_cell_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_decomp_pruned_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_decomp_pruned_cell.f"
#undef __KIND

      END MODULE ppm_module_decomp_pruned_cell
