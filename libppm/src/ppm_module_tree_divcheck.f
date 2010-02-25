      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_tree_divcheck
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 routine that checks if a box is subdivisible.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_tree_divcheck.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/22 10:32:10  ivos
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

      MODULE ppm_module_tree_divcheck

         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_divcheck
            MODULE PROCEDURE ppm_tree_divcheck_s
            MODULE PROCEDURE ppm_tree_divcheck_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_tree_divcheck.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_tree_divcheck.f"
#undef __KIND

      END MODULE ppm_module_tree_divcheck
