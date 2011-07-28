      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_check_topoid
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routine which checks topology IDs for their
      !                 validity.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_check_topoid.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/08/31 12:13:46  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_check_topoid

         !----------------------------------------------------------------------
         !  Define interface to the topoid check routine
         !----------------------------------------------------------------------
         INTERFACE ppm_check_topoid
            MODULE PROCEDURE ppm_check_topoid
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_check_topoid.f"

      END MODULE ppm_module_check_topoid
