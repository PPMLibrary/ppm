      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_topo
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_topo
      !!! This module contains all user-callable routines
      !!! needed to create and manage PPM topologies.
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_mktopo
         USE ppm_module_topo_check
         USE ppm_module_mesh_define
         USE ppm_module_scale_domain
         USE ppm_module_topo_get

      END MODULE ppm_module_topo
