      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_doio
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
      !  $Log: ppm_module_doio.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2004/12/02 10:02:31  ivos
      !  bugfix: ppm_kind_* are now declared PRIVATE to avoid name clashes
      !  with user programs that include ppm_param.h.
      !
      !  Revision 1.2  2004/11/11 15:24:32  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 07:29:34  ivos
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

      MODULE ppm_module_doio

         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------
#include "ppm_define.h"

         !----------------------------------------------------------------------
         !  Work arrays
         !----------------------------------------------------------------------
         REAL(ppm_kind_single)   , DIMENSION(:)    , POINTER :: rbuf_s
         REAL(ppm_kind_double)   , DIMENSION(:)    , POINTER :: rbuf_d
         COMPLEX(ppm_kind_single), DIMENSION(:)    , POINTER :: rbuf_sc
         COMPLEX(ppm_kind_double), DIMENSION(:)    , POINTER :: rbuf_dc
         INTEGER                 , DIMENSION(:)    , POINTER :: rbuf_i
         LOGICAL                 , DIMENSION(:)    , POINTER :: rbuf_l

         PRIVATE :: rbuf_s,rbuf_d,rbuf_sc,rbuf_dc,rbuf_i,rbuf_l

#ifdef __MPI
         REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: abuffer_s
         REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: abuffer_d
         COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: abuffer_sc
         COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: abuffer_dc
         INTEGER                 , DIMENSION(:), POINTER :: abuffer_i
         LOGICAL                 , DIMENSION(:), POINTER :: abuffer_l
         REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: bbuffer_s
         REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: bbuffer_d
         COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: bbuffer_sc
         COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: bbuffer_dc
         INTEGER                 , DIMENSION(:), POINTER :: bbuffer_i
         LOGICAL                 , DIMENSION(:), POINTER :: bbuffer_l

         PRIVATE :: abuffer_s,abuffer_d,abuffer_sc,abuffer_dc
         PRIVATE :: abuffer_i,abuffer_l
         PRIVATE :: bbuffer_s,bbuffer_d,bbuffer_sc,bbuffer_dc
         PRIVATE :: bbuffer_i,bbuffer_l
#endif

         !----------------------------------------------------------------------
         !  Define interface to ppm_doio
         !----------------------------------------------------------------------
         INTERFACE ppm_doio
             MODULE PROCEDURE ppm_doio_s
             MODULE PROCEDURE ppm_doio_d
             MODULE PROCEDURE ppm_doio_i
             MODULE PROCEDURE ppm_doio_l
             MODULE PROCEDURE ppm_doio_sc
             MODULE PROCEDURE ppm_doio_dc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_doio.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_doio.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_doio.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_doio.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_doio.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_doio.f"
#undef __KIND

      END MODULE ppm_module_doio
