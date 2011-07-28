      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_io_write_ascii
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
      !  $Log: ppm_module_io_write_ascii.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:41  ivos
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

      MODULE ppm_module_io_write_ascii

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_write_ascii
         !----------------------------------------------------------------------
         INTERFACE ppm_io_write_ascii
             MODULE PROCEDURE ppm_io_write_asciis
             MODULE PROCEDURE ppm_io_write_asciid
             MODULE PROCEDURE ppm_io_write_asciii
             MODULE PROCEDURE ppm_io_write_asciil
             MODULE PROCEDURE ppm_io_write_asciisc
             MODULE PROCEDURE ppm_io_write_asciidc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_io_write_ascii.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io_write_ascii.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io_write_ascii.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io_write_ascii.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io_write_ascii.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io_write_ascii.f"
#undef __KIND

      END MODULE ppm_module_io_write_ascii
