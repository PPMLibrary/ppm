      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_data_tree
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      MODULE ppm_module_data_tree
      !!! This module contains all data structures and
      !!! definitions that are used in the tree routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this modules should not be directly accessed
      !!! by the user. They are managed interally by the library.
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Data TYPEs
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !  Global data
         !----------------------------------------------------------------------
         ! switches for tree type: do we have particles and/or a mesh?
         LOGICAL                          :: have_particles,have_mesh
         ! Ranked particle lists in all tree boxes
         INTEGER, DIMENSION(:,:), POINTER :: tree_lhbx
         INTEGER, DIMENSION(:  ), POINTER :: tree_lpdx,lhbx_cut,lpdx_cut

         !----------------------------------------------------------------------
         !  Work arrays
         !----------------------------------------------------------------------
         ! list of all current tree boxes and number of subdivisions per box
         INTEGER , DIMENSION(:  ), POINTER       :: boxlist,ndiv
         ! global number of mesh cells and local number of cells per box
         INTEGER , DIMENSION(:,:), POINTER       :: Nmc,Nm_box
         ! boxID for each particle and number of particles per box
         INTEGER , DIMENSION(:), POINTER         :: cbox,npbx
         ! particle-based costs of all boxes
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: pcst_d
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: pcst_s
#ifdef __MPI
         ! accumulated costs from all processors
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: pcsum_d
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: pcsum_s
#endif

      END MODULE ppm_module_data_tree
