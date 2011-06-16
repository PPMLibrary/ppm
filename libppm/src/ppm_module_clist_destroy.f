      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_clist_destroy
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for neighbor
      !                 search routines (cell lists, Verlet lists).
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_clist_destroy.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 12:02:07  ivos
      !  REnamed due to compiler name size limit.
      !
      !  Revision 1.1  2004/07/26 07:30:01  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschebngraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_clist_destroy

         !----------------------------------------------------------------------
         !  Define interface to ppm_clist_destroy
         !----------------------------------------------------------------------
         INTERFACE ppm_clist_destroy
            MODULE PROCEDURE ppm_clist_destroy
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_clist_destroy.f"

      END MODULE ppm_module_clist_destroy
