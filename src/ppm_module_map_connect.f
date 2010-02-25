      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_map_connect
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
      !  $Log: ppm_module_map_connect.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
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
     
      MODULE ppm_module_map_connect

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_connect
         !----------------------------------------------------------------------
         INTERFACE ppm_map_connect
            MODULE PROCEDURE ppm_map_connect
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_connect.f"

      END MODULE ppm_module_map_connect
