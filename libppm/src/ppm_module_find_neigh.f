      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_find_neigh
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
      !  $Log: ppm_module_find_neigh.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 08:56:11  ivos
      !  Renamed.
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

      MODULE ppm_module_find_neigh

         !----------------------------------------------------------------------
         !  Define interface to neighbor finding routine
         !----------------------------------------------------------------------
         INTERFACE ppm_find_neigh
            MODULE PROCEDURE ppm_find_neigh_s
            MODULE PROCEDURE ppm_find_neigh_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_find_neigh.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_find_neigh.f"
#undef __KIND

      END MODULE ppm_module_find_neigh
