      !-------------------------------------------------------------------------
      !  Module       :                   ppm_module_log
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
      !  $Log: ppm_module_log.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:43  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_log

         !----------------------------------------------------------------------
         !  Define interface to the log write rountine
         !----------------------------------------------------------------------
         INTERFACE ppm_log
            MODULE PROCEDURE ppm_log
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_log.f"

      END MODULE ppm_module_log
