      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_io_read_binary
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the IO routines.
      !                 It also included those routines and provides
      !                 INTERFACEs.
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_io_read_binary.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:40  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
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

      MODULE ppm_module_io_read_binary

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_read_binary
         !----------------------------------------------------------------------
         INTERFACE ppm_io_read_binary
             MODULE PROCEDURE ppm_io_read_binarys
             MODULE PROCEDURE ppm_io_read_binaryd
             MODULE PROCEDURE ppm_io_read_binaryi
             MODULE PROCEDURE ppm_io_read_binaryl
             MODULE PROCEDURE ppm_io_read_binarysc
             MODULE PROCEDURE ppm_io_read_binarydc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_io_read_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io_read_binary.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io_read_binary.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io_read_binary.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io_read_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io_read_binary.f"
#undef __KIND

      END MODULE ppm_module_io_read_binary
