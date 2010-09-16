      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_map_part_pop
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

      MODULE ppm_module_map_part
      !!! This module provides the basic mapping routines for particles, namely
      !!! the push, sned and pop routine.
      !!!
      !!! The module also holds the needed work arrays for these routines.
         
         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: sends,recvs
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: sendd,recvd
         INTEGER, DIMENSION(:), POINTER   :: nsend,nrecv,psend,precv
         INTEGER, DIMENSION(:,:), POINTER :: pp,qq
         INTEGER                          :: old_nsendlist

         PRIVATE :: sends,recvs,sendd,recvd,nsend,nrecv,psend,precv,pp,qq
         PRIVATE :: old_nsendlist

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_pop
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_pop
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_pop_1dd
            MODULE PROCEDURE ppm_map_part_pop_1ds
            MODULE PROCEDURE ppm_map_part_pop_1di
            MODULE PROCEDURE ppm_map_part_pop_1dl
            MODULE PROCEDURE ppm_map_part_pop_1ddc
            MODULE PROCEDURE ppm_map_part_pop_1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_pop_2dd
            MODULE PROCEDURE ppm_map_part_pop_2ds
            MODULE PROCEDURE ppm_map_part_pop_2di
            MODULE PROCEDURE ppm_map_part_pop_2dl
            MODULE PROCEDURE ppm_map_part_pop_2ddc
            MODULE PROCEDURE ppm_map_part_pop_2dsc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_push
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_push
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_push_1dd
            MODULE PROCEDURE ppm_map_part_push_1ds
            MODULE PROCEDURE ppm_map_part_push_1di
            MODULE PROCEDURE ppm_map_part_push_1dl
            MODULE PROCEDURE ppm_map_part_push_1ddc
            MODULE PROCEDURE ppm_map_part_push_1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_push_2dd
            MODULE PROCEDURE ppm_map_part_push_2ds
            MODULE PROCEDURE ppm_map_part_push_2di
            MODULE PROCEDURE ppm_map_part_push_2dl
            MODULE PROCEDURE ppm_map_part_push_2ddc
            MODULE PROCEDURE ppm_map_part_push_2dsc
         END INTERFACE
         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_send
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_send
            MODULE PROCEDURE ppm_map_part_send
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "map/ppm_map_part_push.f"
#undef __KIND
#undef __DIM

#include "map/ppm_map_part_send.f"

      END MODULE ppm_module_map_part
