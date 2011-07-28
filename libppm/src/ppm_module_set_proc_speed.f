      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_set_proc_speed
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
      !  $Log: ppm_module_set_proc_speed.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:30:07  ivos
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

      MODULE ppm_module_set_proc_speed

         !----------------------------------------------------------------------
         !  Define interface to processor speed setter
         !----------------------------------------------------------------------
         INTERFACE ppm_set_proc_speed
            MODULE PROCEDURE ppm_set_proc_speed_s
            MODULE PROCEDURE ppm_set_proc_speed_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_set_proc_speed.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_set_proc_speed.f"
#undef __KIND

      END MODULE ppm_module_set_proc_speed
