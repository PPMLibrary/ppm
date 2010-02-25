      !-------------------------------------------------------------------------
      !  Module       :                   ppm_module_map
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all user-callable routines for
      !                 mapping particles and fields.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.24  2006/10/10 21:27:33  walther
      !  Now using ppm_module_map_part_load and ppm_module_map_part_store.
      !
      !  Revision 1.23  2004/08/03 09:13:57  ivos
      !  Added ppm_module_impose_part_bc.
      !
      !  Revision 1.22  2004/07/26 15:38:49  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.21  2004/07/26 13:40:29  ivos
      !  Initial implementation. These are meta-modules for the user-
      !  callable functions. Only these mod files will be given away
      !  to the user.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_map

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_map_part
         USE ppm_module_map_part_ghost
         USE ppm_module_map_part_ring_shift
         USE ppm_module_map_field
         USE ppm_module_map_field_ghost
         USE ppm_module_map_connect
         USE ppm_module_impose_part_bc
         USE ppm_module_map_part_load
         USE ppm_module_map_part_store   
         
      END MODULE ppm_module_map
