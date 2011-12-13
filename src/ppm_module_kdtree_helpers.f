!!
!!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!!
!! Licensed under the Academic Free License version 1.1 found in file LICENSE
!! with additional provisions found in that same file.
!!
!module kdtree2_precision_module
  
  !integer, parameter :: sp = kind(0.0)
  !integer, parameter :: dp = kind(0.0d0)

  !private :: sp, dp

  !!
  !! You must comment out exactly one
  !! of the two lines.  If you comment
  !! out kdkind = sp then you get single precision
  !! and if you comment out kdkind = dp 
  !! you get double precision.
  !!

  !!integer, parameter :: kdkind = sp  
  !integer, parameter :: kdkind = dp  
  !public :: kdkind

!end module kdtree2_precision_module

!module kdtree2_priority_queue_module
module ppm_module_kdtree_helpers
  !
  ! maintain a priority queue (PQ) of data, pairs of 'priority/payload', 
  ! implemented with a binary heap.  This is the type, and the 'dis' field
  ! is the priority.
  !
  USE ppm_module_typedef

#define __KIND __SINGLE_PRECISION
#define kdkind ppm_kind_single
#define  DTYPE(a) a/**/_s
#include "kdtree/ppm_module_kdtree_helpers.inc"
#undef  __KIND
#undef  DTYPE
#undef  kdkind

#define __KIND __DOUBLE_PRECISION
#define kdkind ppm_kind_double
#define  DTYPE(a) a/**/_d
#include "kdtree/ppm_module_kdtree_helpers.inc"
#undef  __KIND
#undef  DTYPE
#undef  kdkind

  INTERFACE pq_create
      MODULE PROCEDURE pq_create_s
      MODULE PROCEDURE pq_create_d
  END INTERFACE
  INTERFACE pq_insert
      MODULE PROCEDURE pq_insert_s
      MODULE PROCEDURE pq_insert_d
  END INTERFACE
  INTERFACE pq_extract_max
      MODULE PROCEDURE pq_extract_max_s
      MODULE PROCEDURE pq_extract_max_d
  END INTERFACE
  INTERFACE pq_replace_max
      MODULE PROCEDURE pq_replace_max_s
      MODULE PROCEDURE pq_replace_max_d
  END INTERFACE
  INTERFACE pq_max
      MODULE PROCEDURE pq_max_s
      MODULE PROCEDURE pq_max_d
  END INTERFACE
  INTERFACE pq_maxpri
      MODULE PROCEDURE pq_maxpri_s
      MODULE PROCEDURE pq_maxpri_d
  END INTERFACE

  public :: pq_create
  public :: pq_delete, pq_insert
  public :: pq_extract_max, pq_max, pq_replace_max, pq_maxpri
  private

contains


#define __KIND __SINGLE_PRECISION
#define kdkind ppm_kind_single
#define  DTYPE(a) a/**/_s
#include "kdtree/ppm_kdtree_helpers.f"
#undef  __KIND
#undef  DTYPE
#undef  kdkind

#define __KIND __DOUBLE_PRECISION
#define kdkind ppm_kind_double
#define  DTYPE(a) a/**/_d
#include "kdtree/ppm_kdtree_helpers.f"
#undef  __KIND
#undef  DTYPE
#undef  kdkind

end module ppm_module_kdtree_helpers



