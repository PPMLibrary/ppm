      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_map_part_ghost_put
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
      !  $Log: ppm_module_map_part_ghost_put.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2006/04/06 14:57:07  walther
      !  the overloading of ppm_map_part_ghost_put is no longer necessary;
      !  now only one ppm_map_part_ghost_put (as the only arguments is info).
      !
      !  Revision 1.1  2004/07/26 07:29:52  ivos
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
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2

      MODULE ppm_module_map_part_ghost_put

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ghost_put
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ghost_put
            MODULE PROCEDURE ppm_map_part_ghost_put
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_part_ghost_put.f"
      END MODULE ppm_module_map_part_ghost_put
