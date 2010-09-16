      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_user
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_user
      !!! This is the global user module. It contains all
      !!! user-callable routines of the entire ppm core library.
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_map
         USE ppm_module_user_util
         USE ppm_module_loadbal
         USE ppm_module_topo
         USE ppm_module_io
         USE ppm_module_neighlist
         USE ppm_module_tree
         USE ppm_module_mesh
         USE ppm_module_rmsh

      END MODULE ppm_module_user
