MODULE ppm_module_kdtree
  USE ppm_module_typedef
  USE ppm_module_kdtree_helpers
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported. 
  !
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module. 

  INTERFACE kdtree2_create
      MODULE PROCEDURE kdtree2_create_s
      MODULE PROCEDURE kdtree2_create_d
  END INTERFACE
  INTERFACE kdtree2_destroy
      MODULE PROCEDURE kdtree2_destroy_s
      MODULE PROCEDURE kdtree2_destroy_d
  END INTERFACE
  INTERFACE kdtree2_n_nearest
      MODULE PROCEDURE kdtree2_n_nearest_s
      MODULE PROCEDURE kdtree2_n_nearest_d
  END INTERFACE
  INTERFACE kdtree2_n_nearest_around_point
      MODULE PROCEDURE kdtree2_n_nearest_around_point_s
      MODULE PROCEDURE kdtree2_n_nearest_around_point_d
  END INTERFACE
  INTERFACE kdtree2_r_nearest
      MODULE PROCEDURE kdtree2_r_nearest_s
      MODULE PROCEDURE kdtree2_r_nearest_d
  END INTERFACE
  INTERFACE kdtree2_r_nearest_around_point
      MODULE PROCEDURE kdtree2_r_nearest_around_point_s
      MODULE PROCEDURE kdtree2_r_nearest_around_point_d
  END INTERFACE
  INTERFACE kdtree2_sort_results
      MODULE PROCEDURE kdtree2_sort_results_s
      MODULE PROCEDURE kdtree2_sort_results_d
  END INTERFACE
  INTERFACE kdtree2_r_count
      MODULE PROCEDURE kdtree2_r_count_s
      MODULE PROCEDURE kdtree2_r_count_d
  END INTERFACE
  INTERFACE kdtree2_r_count_around_point
      MODULE PROCEDURE kdtree2_r_count_around_point_s
      MODULE PROCEDURE kdtree2_r_count_around_point_d
  END INTERFACE
  INTERFACE kdtree2_n_nearest_brute_force
      MODULE PROCEDURE kdtree2_n_nearest_brute_force_s
      MODULE PROCEDURE kdtree2_n_nearest_brute_force_d
  END INTERFACE
  INTERFACE kdtree2_r_nearest_brute_force
      MODULE PROCEDURE kdtree2_r_nearest_brute_force_s
      MODULE PROCEDURE kdtree2_r_nearest_brute_force_d
  END INTERFACE

  !----------------------------------------------------------------
  !
  !-------------DATA TYPE, CREATION, DELETION---------------------
  PUBLIC :: kdtree2_s, kdtree2_result_s, tree_node_s
  PUBLIC :: kdtree2_d, kdtree2_result_d, tree_node_d
  PUBLIC :: kdtree2_create, kdtree2_destroy
  !---------------------------------------------------------------
  !-------------------SEARCH ROUTINES-----------------------------
  PUBLIC :: kdtree2_n_nearest,kdtree2_n_nearest_around_point
  ! Return fixed number of nearest neighbors around arbitrary vector,
  ! or extant point in dataset, with decorrelation window. 
  !
  PUBLIC :: kdtree2_r_nearest, kdtree2_r_nearest_around_point
  ! Return points within a fixed ball of arb vector/extant point 
  !
  PUBLIC :: kdtree2_sort_results
  ! Sort, in order of increasing distance, rseults from above.
  !
  PUBLIC :: kdtree2_r_count, kdtree2_r_count_around_point 
  ! Count points within a fixed ball of arb vector/extant point 
  !
  PUBLIC :: kdtree2_n_nearest_brute_force, kdtree2_r_nearest_brute_force
  ! brute force of kdtree2_[n|r]_nearest
  !----------------------------------------------------------------


  INTEGER, PARAMETER :: bucket_size = 12
  ! The maximum number of points to keep in a terminal node.

#define __KIND __SINGLE_PRECISION
#define kdkind ppm_kind_single
#define  DTYPE(a) a/**/_s
#include "kdtree/ppm_module_kdtree.inc"
#undef  __KIND
#undef  DTYPE
#undef  kdkind

#define __KIND __DOUBLE_PRECISION
#define kdkind ppm_kind_double
#define  DTYPE(a) a/**/_d
#include "kdtree/ppm_module_kdtree.inc"
#undef  __KIND
#undef  DTYPE
#undef  kdkind

  PRIVATE
  ! everything else is private.

  TYPE(tree_search_record_s), SAVE, TARGET :: sr_s   ! A GLOBAL VARIABLE for search
  TYPE(tree_search_record_d), SAVE, TARGET :: sr_d   ! A GLOBAL VARIABLE for search

CONTAINS

#define __KIND __SINGLE_PRECISION
#define kdkind ppm_kind_single
#define  DTYPE(a) a/**/_s
#define  sr sr_s
#include "kdtree/ppm_kdtree.f"
#undef  sr
#undef  __KIND
#undef  DTYPE
#undef  kdkind

#define __KIND __DOUBLE_PRECISION
#define kdkind ppm_kind_double
#define  DTYPE(a) a/**/_d
#define  sr sr_d
#include "kdtree/ppm_kdtree.f"
#undef  sr
#undef  __KIND
#undef  DTYPE
#undef  kdkind

END MODULE ppm_module_kdtree

