      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_map_part_get_sub
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the mapping
      !                 routine ppm_map_part_get_sub.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map_part_get_sub.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2005/02/07 15:44:37  walther
      !  First release - not tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6

      MODULE ppm_module_map_part_get_sub
         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER                         :: nlist1,nlist2,nlist3,nlist4
         INTEGER , DIMENSION(:), POINTER :: ilist1,ilist2,ilist3,ilist4

         !----------------------------------------------------------------------
         !  Define interfaces to the ppm_map_part
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_get_sub
            MODULE PROCEDURE ppm_map_part_get_sub_s
            MODULE PROCEDURE ppm_map_part_get_sub_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_get_sub.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_get_sub.f"
#undef __KIND

      END MODULE ppm_module_map_part_get_sub
