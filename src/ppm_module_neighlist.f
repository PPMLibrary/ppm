      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_neighlist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the user-callable functions 
      !                 for the neighbor list routines.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_neighlist.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/07/26 15:38:50  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.5  2004/07/26 13:40:33  ivos
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

      MODULE ppm_module_neighlist

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_data_neighlist
         USE ppm_module_neighlist_clist
         USE ppm_module_clist_destroy
         USE ppm_module_neighlist_vlist
         USE ppm_module_neighlist_MkNeighIdx
         
      END MODULE ppm_module_neighlist
