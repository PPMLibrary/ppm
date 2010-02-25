      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_comp_pp_correct
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 concerned with computing kernel interactions
      !                 within a set of particles.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_comp_pp_correct.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:27  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
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
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __SINGLE_PRECISION_COMPLEX 3
#define __DOUBLE_PRECISION_COMPLEX 4

#define __INTERNAL                 5
#define __USER_FUNCTION            6
#define __LOOKUP_TABLE             7

      MODULE ppm_module_comp_pp_correct

         !----------------------------------------------------------------------
         !  Define interfaces to the discretisation correction routine
         !----------------------------------------------------------------------
         INTERFACE ppm_comp_pp_correct
            MODULE PROCEDURE ppm_comp_pp_correct_si
            MODULE PROCEDURE ppm_comp_pp_correct_di
            MODULE PROCEDURE ppm_comp_pp_correct_sci
            MODULE PROCEDURE ppm_comp_pp_correct_dci
            MODULE PROCEDURE ppm_comp_pp_correct_su
            MODULE PROCEDURE ppm_comp_pp_correct_du
            MODULE PROCEDURE ppm_comp_pp_correct_scu
            MODULE PROCEDURE ppm_comp_pp_correct_dcu
            MODULE PROCEDURE ppm_comp_pp_correct_st
            MODULE PROCEDURE ppm_comp_pp_correct_dt
            MODULE PROCEDURE ppm_comp_pp_correct_sct
            MODULE PROCEDURE ppm_comp_pp_correct_dct
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KERNEL __INTERNAL
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_correct.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __USER_FUNCTION
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_correct.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __LOOKUP_TABLE
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_correct.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_correct.f"
#undef __KIND
#undef __KERNEL

      END MODULE ppm_module_comp_pp_correct
