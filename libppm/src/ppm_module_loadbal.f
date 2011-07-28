      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_loadbal
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all user-callable routines
      !                 needed for dynamic load balancing.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_loadbal.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2004/07/26 13:40:31  ivos
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

      MODULE ppm_module_loadbal

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_get_cost
         USE ppm_module_set_proc_speed
         USE ppm_module_estimate_proc_speed
         USE ppm_module_loadbal_inquire
         USE ppm_module_set_decomp_cost
         
      END MODULE ppm_module_loadbal
