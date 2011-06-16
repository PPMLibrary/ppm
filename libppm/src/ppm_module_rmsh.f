      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_rmsh
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all user-callable routines
      !                 needed for interpolation and remeshing
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_rmsh.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/12/02 10:03:13  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_rmsh

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_rmsh_comp_weights
         USE ppm_module_rmsh_create_part
         USE ppm_module_rmsh_remesh
         USE ppm_module_interp_m2p
         USE ppm_module_interp_p2m
         
      END MODULE ppm_module_rmsh
