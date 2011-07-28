      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_map_part_load
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the loading of
      !                 the stored mapping. This routine is user callable.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map_part_load.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2006/10/10 20:57:23  walther
      !  *** empty log message ***
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define the module
      !-------------------------------------------------------------------------
      MODULE ppm_module_map_part_load

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_load
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_load
            MODULE PROCEDURE ppm_map_part_load
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_part_load.f"

      END MODULE ppm_module_map_part_load
