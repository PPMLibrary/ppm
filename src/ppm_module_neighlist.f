      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_neighlist
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_neighlist
      !!! This module provides neighbor
      !!! search routines (cell lists, Verlet lists).
         
         USE ppm_module_data_neighlist, ONLY: ppm_type_ptr_to_clist
         
         !----------------------------------------------------------------------
         !  Temporary cell list memory
         !----------------------------------------------------------------------
         TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER :: clist
         PRIVATE :: clist
         
         !----------------------------------------------------------------------
         !  Define interface to ppm_clist_destroy
         !----------------------------------------------------------------------
         INTERFACE ppm_clist_destroy
            MODULE PROCEDURE ppm_clist_destroy
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_neighlist_MkNeighIdx
         !----------------------------------------------------------------------
         INTERFACE ppm_neighlist_MkNeighIdx
            MODULE PROCEDURE ppm_neighlist_MkNeighIdx
         END INTERFACE

         INTERFACE ppm_neighlist_clist
            MODULE PROCEDURE ppm_neighlist_clist_d
            MODULE PROCEDURE ppm_neighlist_clist_s
         END INTERFACE

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

#include "neighlist/ppm_clist_destroy.f"

#include "neighlist/ppm_neighlist_MkNeighIdx.f"

#define __KIND __SINGLE_PRECISION
#include "neighlist/ppm_neighlist_clist.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "neighlist/ppm_neighlist_clist.f"
#undef  __KIND

#define __KIND __SINGLE_PRECISION
#include "neighlist/ppm_neighlist_vlist.f"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "neighlist/ppm_neighlist_vlist.f"
#undef  __KIND

      END MODULE ppm_module_neighlist
