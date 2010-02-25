      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_util_eigen_2sym
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 that compute Eigenvalues and Eigenvectors of
      !                 symmetric 2x2 matrices.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_eigen_2sym.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/09/17 14:07:58  ivos
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

      MODULE ppm_module_util_eigen_2sym

         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_util_eigen_2sym
            MODULE PROCEDURE ppm_util_eigen_2sym_s
            MODULE PROCEDURE ppm_util_eigen_2sym_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_util_eigen_2sym.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_eigen_2sym.f"
#undef __KIND

      END MODULE ppm_module_util_eigen_2sym
