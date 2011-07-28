      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_data_part
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the mapping
      !                 routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_part.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:28:19  ivos
      !  Initial implementation. Originated from splitting the old ppm
      !  modules.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_data_part

         !----------------------------------------------------------------------
         !  Declare the target topoid for the mapping
         !----------------------------------------------------------------------
         INTEGER                       :: ppm_target_topoid = -1

      END MODULE ppm_module_data_part
