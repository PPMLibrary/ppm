      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_data_tree
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are GLOBAL to all tree routines.
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_tree.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2006/07/03 12:56:56  ivos
      !  Added comments to explain all the variables.
      !
      !  Revision 1.5  2005/08/31 13:35:20  ivos
      !  Moved daclaration of icut from here to ppm_tree since it is used
      !  as argument to subroutines and is always small (ppm_dim).
      !
      !  Revision 1.4  2005/08/31 11:24:29  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.3  2005/02/01 13:22:57  ivos
      !  Moved declarations of lhbx_cut and lpdx_cut to module_data_tree.
      !
      !  Revision 1.2  2004/12/03 17:13:44  ivos
      !  Added tree_lhbx and tree_lpdx.
      !
      !  Revision 1.1  2004/09/22 10:32:10  ivos
      !  Initial implementation.
      !
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      MODULE ppm_module_data_tree

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
