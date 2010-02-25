      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_comp_part
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all user-callable routines to
      !                 compute kernel-PP interactions on a set of particles.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_comp_part.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 13:40:30  ivos
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

      MODULE ppm_module_comp_part

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_comp_pp_verlet
         USE ppm_module_comp_pp_cell
         USE ppm_module_comp_pp_ring
         USE ppm_module_comp_pp_correct
         USE ppm_module_comp_pp_mk_table
         
      END MODULE ppm_module_comp_part
