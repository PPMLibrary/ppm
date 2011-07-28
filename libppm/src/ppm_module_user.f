      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_user
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This is the global user module. It contains all
      !                 user-callable routines of the entire ppm library.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_user.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.25  2006/09/04 18:34:53  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.24  2006/05/12 15:49:41  pchatela
      !  gmm and hamjac back in place...
      !
      !  Revision 1.23  2006/05/12 15:45:11  pchatela
      !  Removed modules hamjac and gmm for release 1.0
      !
      !  Revision 1.22  2006/02/03 09:36:38  ivos
      !  Added module_fmm
      !
      !  Revision 1.21  2005/07/25 00:35:28  ivos
      !  Added ppm_module_hamjac
      !
      !  Revision 1.20  2005/03/10 01:49:48  ivos
      !  Added ppm_module_gmm.
      !
      !  Revision 1.19  2004/12/02 10:02:50  ivos
      !  Added ppm_module_tree and ppm_module_rmsh.
      !
      !  Revision 1.18  2004/09/22 18:32:31  kotsalie
      !  changed multigrid to mg
      !
      !  Revision 1.17  2004/07/26 13:40:28  ivos
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

      MODULE ppm_module_user

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_comp_part
         USE ppm_module_map
         USE ppm_module_user_util
         USE ppm_module_loadbal
         USE ppm_module_topo
         USE ppm_module_user_io
         USE ppm_module_neighlist
         USE ppm_module_tree
         USE ppm_module_rmsh
         
      END MODULE ppm_module_user
