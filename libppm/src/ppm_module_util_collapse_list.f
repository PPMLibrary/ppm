      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_util_collapse_list
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
      !  $Log: ppm_module_util_collapse_list.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/02/09 16:29:53  pchatela
      !  Initial insertion
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_util_collapse_list

         !----------------------------------------------------------------------
         !  Define interface to the list inversion routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_collapse_list
            MODULE PROCEDURE ppm_util_collapse_list
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_util_collapse_list.f"

      END MODULE ppm_module_util_collapse_list
