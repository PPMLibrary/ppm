      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_user_util
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_user_util
      !!! This module contains all user-callable utility routines.
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_init
         USE ppm_module_finalize
         USE ppm_module_find_duplicates
         USE ppm_module_is_initialized
         USE ppm_module_time

      END MODULE ppm_module_user_util
