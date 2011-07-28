      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_map_part_util
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
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

      MODULE ppm_module_map_part_util
      !!! This module provides various particle mapping utility and helper
      !!! routines and their work-arrays.

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: plist_des,plist_act,plist_exc
         INTEGER, DIMENSION(:), POINTER :: srlist1,srlist2
         INTEGER                        :: nlist1,nlist2,nlist3,nlist4
         INTEGER ,DIMENSION(:), POINTER :: ilist1,ilist2,ilist3,ilist4

         PRIVATE :: plist_des,plist_act,plist_exc,srlist1,srlist2
         PRIVATE :: nlist1,nlist2,nlist3,nlist4
         PRIvATE :: ilist1,ilist2,ilist3,ilist4
         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_eqdistrib
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_eqdistrib
            MODULE PROCEDURE ppm_map_part_eqdistrib_d
            MODULE PROCEDURE ppm_map_part_eqdistrib_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_send
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_cancel
            MODULE PROCEDURE ppm_map_part_cancel
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to the ppm_map_part_get_sub
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_get_sub
            MODULE PROCEDURE ppm_map_part_get_sub_s
            MODULE PROCEDURE ppm_map_part_get_sub_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_load
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_load
            MODULE PROCEDURE ppm_map_part_load
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ring_shift
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ring_shift
            MODULE PROCEDURE ppm_map_part_ring_shift_s
            MODULE PROCEDURE ppm_map_part_ring_shift_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_store
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_store
            MODULE PROCEDURE ppm_map_part_store
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_eqdistrib.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_eqdistrib.f"
#undef __KIND

#include "map/ppm_map_part_cancel.f"

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_get_sub.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_get_sub.f"
#undef __KIND

#include "map/ppm_map_part_load.f"

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_ring_shift.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_ring_shift.f"
#undef __KIND

#include "map/ppm_map_part_store.f"

      END MODULE ppm_module_map_part_util
