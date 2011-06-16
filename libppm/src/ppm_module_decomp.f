      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_decomp
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

      MODULE ppm_module_decomp
      !!! This module provides the PPM decomposition routines.
      !!!
      !!! PPM core provides 4 different decomposition routines that are called
      !!! within the topology routines.

         !----------------------------------------------------------------------
         !  Local work arrays used in ppm_decomp_pruned_cell
         !----------------------------------------------------------------------
         INTEGER , DIMENSION(:), POINTER, PRIVATE :: npbx,npbxg

         !----------------------------------------------------------------------
         !  Define interface to box split routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_boxsplit
            MODULE PROCEDURE decomp_bsplit_s
            MODULE PROCEDURE decomp_bsplit_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to cartesian decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_cartesian
            MODULE PROCEDURE decomp_cart_s
            MODULE PROCEDURE decomp_cart_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to pruned cell decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_pruned_cell
            MODULE PROCEDURE decomp_pcell_s
            MODULE PROCEDURE decomp_pcell_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to tree decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_tree
            MODULE PROCEDURE decomp_tree_s
            MODULE PROCEDURE decomp_tree_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_boxsplit.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_boxsplit.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_cartesian.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_cartesian.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_pruned_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_pruned_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "decomp/ppm_decomp_tree.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "decomp/ppm_decomp_tree.f"
#undef __KIND

END MODULE ppm_module_decomp
