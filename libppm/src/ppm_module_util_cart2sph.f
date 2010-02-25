      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_util_cart2sph
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 that convert cartesian to spherical coordinates.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_cart2sph.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/02/09 15:27:26  polasekb
      !  initial implementation
      !
      !  Revision 1.1  2004/12/02 10:07:57  polasekb
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

      MODULE ppm_module_util_cart2sph

         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_util_cart2sph
            MODULE PROCEDURE ppm_util_cart2sph_s
            MODULE PROCEDURE ppm_util_cart2sph_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_util_cart2sph.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_cart2sph.f"
#undef __KIND

      END MODULE ppm_module_util_cart2sph
