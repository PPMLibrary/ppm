      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_mesh
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_mesh
      !!! This module contains all user-callable routines for meshes

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_mesh_alloc
         USE ppm_module_mesh_block_intersect
         USE ppm_module_mesh_define
         USE ppm_module_mesh_derive
         USE ppm_module_mesh_finalize
         USE ppm_module_mesh_on_subs
         USE ppm_module_mesh_store


      END MODULE ppm_module_mesh
