      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_map_connect_prune
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the mapping
      !                 routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map_connect_prune.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/11 15:24:33  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 08:24:45  ivos
      !  Check-in after module cleanup. All ppm_map_part_connect* were
      !  renamed to ppm_map_connect*. Also the modules.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_map_connect_prune

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:)  , POINTER :: id_temp,id_inv
         INTEGER, DIMENSION(:,:), POINTER :: cd_local

         PRIVATE :: id_temp,id_inv,cd_local

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_connect_prune
         !----------------------------------------------------------------------
         INTERFACE ppm_map_connect_prune
            MODULE PROCEDURE ppm_map_connect_prune
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_connect_prune.f"

      END MODULE ppm_module_map_connect_prune
