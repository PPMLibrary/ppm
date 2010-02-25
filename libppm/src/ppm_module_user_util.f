      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_user_util
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all user-callable utility
      !                 routines.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_user_util.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 13:40:31  ivos
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

      MODULE ppm_module_user_util

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_init
         USE ppm_module_finalize
         USE ppm_module_find_duplicates
         USE ppm_module_get_revision
         USE ppm_module_is_initialized
         USE ppm_module_time
         
      END MODULE ppm_module_user_util
