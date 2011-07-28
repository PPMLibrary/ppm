      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_util_eigen_3sym
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 that compute Eigenvalues and Eigenvectors of
      !                 symmetric 3x3 matrices.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_eigen_3sym.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/17 14:07:59  ivos
      !  Initial implementation.
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_util_eigen_3sym

         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_util_eigen_3sym
            MODULE PROCEDURE ppm_util_eigen_3sym_s
            MODULE PROCEDURE ppm_util_eigen_3sym_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_util_eigen_3sym.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_eigen_3sym.f"
#undef __KIND

      END MODULE ppm_module_util_eigen_3sym
