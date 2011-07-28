      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_error
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_error.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:35  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_error

         !----------------------------------------------------------------------
         !  Header file for error codes
         !----------------------------------------------------------------------
         INCLUDE 'ppm_error.h'

         !----------------------------------------------------------------------
         !  Define interface to the error handler routine
         !----------------------------------------------------------------------
         INTERFACE ppm_error
            MODULE PROCEDURE ppm_error
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_error.f"

      END MODULE ppm_module_error
