      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_map_part_store
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the storing of
      !                 mapping.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map_part_store.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2006/10/10 20:58:10  walther
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
      MODULE ppm_module_map_part_store

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_store
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_store
            MODULE PROCEDURE ppm_map_part_store
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_part_store.f"

      END MODULE ppm_module_map_part_store
