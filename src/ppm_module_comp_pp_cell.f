      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_comp_pp_cell
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
      !  $Log: ppm_module_comp_pp_cell.f,v $
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

      MODULE ppm_module_comp_pp_cell

         !----------------------------------------------------------------------
         !  Define interfaces to the cell list versions
         !----------------------------------------------------------------------
         INTERFACE ppm_comp_pp_cell
            MODULE PROCEDURE ppm_comp_pp_cell_si
            MODULE PROCEDURE ppm_comp_pp_cell_di
            MODULE PROCEDURE ppm_comp_pp_cell_sci
            MODULE PROCEDURE ppm_comp_pp_cell_dci
            MODULE PROCEDURE ppm_comp_pp_cell_su
            MODULE PROCEDURE ppm_comp_pp_cell_du
            MODULE PROCEDURE ppm_comp_pp_cell_scu
            MODULE PROCEDURE ppm_comp_pp_cell_dcu
            MODULE PROCEDURE ppm_comp_pp_cell_st
            MODULE PROCEDURE ppm_comp_pp_cell_dt
            MODULE PROCEDURE ppm_comp_pp_cell_sct
            MODULE PROCEDURE ppm_comp_pp_cell_dct
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KERNEL __INTERNAL
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __USER_FUNCTION
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __LOOKUP_TABLE
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND
#undef __KERNEL

      END MODULE ppm_module_comp_pp_cell
