      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_neighlist_vlist
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
      !  $Log: ppm_module_neighlist_vlist.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/11 15:26:20  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 07:30:03  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschebngraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_neighlist_vlist

         USE ppm_module_data_neighlist, ONLY: ppm_type_ptr_to_clist
         PRIVATE :: ppm_type_ptr_to_clist,clist

         !----------------------------------------------------------------------
         !  Temporary cell list memory
         !----------------------------------------------------------------------
         TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER :: clist

         !----------------------------------------------------------------------
         !  Define interface to ppm_neighlist_vlist
         !----------------------------------------------------------------------
         INTERFACE ppm_neighlist_vlist
            MODULE PROCEDURE ppm_neighlist_vlist_d
            MODULE PROCEDURE ppm_neighlist_vlist_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_neighlist_vlist.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_neighlist_vlist.f"
#undef  __KIND

      END MODULE ppm_module_neighlist_vlist
