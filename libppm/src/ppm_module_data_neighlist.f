      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_data_neighlist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for neighbor
      !                 search routines (cell lists, Verlet lists).
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_neighlist.f,v $
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
      !  ETH Zentrum, Hirschebngraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_data_neighlist

         !----------------------------------------------------------------------
         !  Define data TYPEs
         !----------------------------------------------------------------------
         ! Pointer to cell list (needed to make lists of cell lists)
         TYPE ppm_type_ptr_to_clist
             ! particle index list
             INTEGER, DIMENSION(:), POINTER    :: lpdx
             ! first particle in each cell
             INTEGER, DIMENSION(:), POINTER    :: lhbx
         END TYPE

      END MODULE ppm_module_data_neighlist
