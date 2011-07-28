      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                   ppm_module_map
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_map
      !!! This module contains all user-callable routines for mapping
      !!! particles and fields.
      !!!
      !!! You may either `USE` this module or directly the modules for the
      !!! subroutines you will be using.
         ! TODO: map_field_init added from Petros codebase (this might
         ! need revising, maybe transparent handling as done for ghost_init)
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------

         USE ppm_module_map_connect
         USE ppm_module_map_field_ghost
         USE ppm_module_map_field_global
         USE ppm_module_map_field
         USE ppm_module_map_part_util
         USE ppm_module_map_part_ghost
         USE ppm_module_map_part_global
         USE ppm_module_map_part_partial
         USE ppm_module_map_part
         USE ppm_module_impose_part_bc
         
      END MODULE ppm_module_map
