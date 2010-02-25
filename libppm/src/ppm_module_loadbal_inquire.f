      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_loadbal_inquire
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 concerned with load balancing. They are callable by
      !                 the user.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_loadbal_inquire.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:42  ivos
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

      MODULE ppm_module_loadbal_inquire

         !----------------------------------------------------------------------
         !  Define interface to load balance inquiry routine
         !----------------------------------------------------------------------
         INTERFACE ppm_loadbal_inquire
            MODULE PROCEDURE ppm_loadbal_inquire_s
            MODULE PROCEDURE ppm_loadbal_inquire_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_loadbal_inquire.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_loadbal_inquire.f"
#undef __KIND

      END MODULE ppm_module_loadbal_inquire
