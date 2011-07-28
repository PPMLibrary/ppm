      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_io
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
      !  $Log: ppm_module_io.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.10  2004/12/02 10:02:31  ivos
      !  bugfix: ppm_kind_* are now declared PRIVATE to avoid name clashes
      !  with user programs that include ppm_param.h.
      !
      !  Revision 1.9  2004/11/11 15:24:31  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.8  2004/07/26 07:29:39  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !  Revision 1.7  2004/06/02 13:41:59  ivos
      !  Declared ppm_kind_single and ppm_kind_double PRIVATE to avoid
      !  naming conflicts with user programs that include ppm_param.h
      !
      !  Revision 1.6  2004/05/27 10:44:54  ivos
      !  Moved the buffers here and added INTERFACEs and includes for ppm_doio,
      !  ppm_io_write_ascii, ppm_io_write_binary, ppm_io_read_ascii and
      !  ppm_io_read_binary.
      !
      !  Revision 1.5  2004/05/13 11:40:02  ivos
      !  Added ppm_io_unused_unit.
      !
      !  Revision 1.4  2004/05/11 14:54:39  ivos
      !  Added ppm_io_delete.
      !
      !  Revision 1.3  2004/05/06 10:42:17  ivos
      !  Added interfaces to ppm_io_open, ppm_io_close, ppm_io_inquire and
      !  ppm_io.
      !
      !  Revision 1.2  2004/05/06 07:47:05  ivos
      !  Renamed ppm_set_unit to ppm_io_set_unit.
      !
      !  Revision 1.1  2004/02/18 17:44:18  walther
      !  Initial implementation of the io module for the ppm library.
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
#define __CHARACTER                7

      MODULE ppm_module_io

         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Work arrays
         !----------------------------------------------------------------------
         REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: abuf_s
         REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: abuf_d
         COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: abuf_sc
         COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: abuf_dc
         INTEGER                 , DIMENSION(:), POINTER :: abuf_i
         LOGICAL                 , DIMENSION(:), POINTER :: abuf_l

         !----------------------------------------------------------------------
         !  Define interface to ppm_io
         !----------------------------------------------------------------------
         INTERFACE ppm_io
             MODULE PROCEDURE ppm_io_0ds
             MODULE PROCEDURE ppm_io_0dd
             MODULE PROCEDURE ppm_io_0di
             MODULE PROCEDURE ppm_io_0dl
             MODULE PROCEDURE ppm_io_0dsc
             MODULE PROCEDURE ppm_io_0ddc
             MODULE PROCEDURE ppm_io_0dc

             MODULE PROCEDURE ppm_io_1ds
             MODULE PROCEDURE ppm_io_1dd
             MODULE PROCEDURE ppm_io_1di
             MODULE PROCEDURE ppm_io_1dl
             MODULE PROCEDURE ppm_io_1dsc
             MODULE PROCEDURE ppm_io_1ddc

             MODULE PROCEDURE ppm_io_2ds
             MODULE PROCEDURE ppm_io_2dd
             MODULE PROCEDURE ppm_io_2di
             MODULE PROCEDURE ppm_io_2dl
             MODULE PROCEDURE ppm_io_2dsc
             MODULE PROCEDURE ppm_io_2ddc

             MODULE PROCEDURE ppm_io_3ds
             MODULE PROCEDURE ppm_io_3dd
             MODULE PROCEDURE ppm_io_3di
             MODULE PROCEDURE ppm_io_3dl
             MODULE PROCEDURE ppm_io_3dsc
             MODULE PROCEDURE ppm_io_3ddc

             MODULE PROCEDURE ppm_io_4ds
             MODULE PROCEDURE ppm_io_4dd
             MODULE PROCEDURE ppm_io_4di
             MODULE PROCEDURE ppm_io_4dl
             MODULE PROCEDURE ppm_io_4dsc
             MODULE PROCEDURE ppm_io_4ddc

             MODULE PROCEDURE ppm_io_5ds
             MODULE PROCEDURE ppm_io_5dd
             MODULE PROCEDURE ppm_io_5di
             MODULE PROCEDURE ppm_io_5dl
             MODULE PROCEDURE ppm_io_5dsc
             MODULE PROCEDURE ppm_io_5ddc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM 0
#define __KIND __SINGLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#define __KIND __CHARACTER
#include "ppm_io.f"
#undef __KIND
#undef __DIM

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#undef __DIM
         
#define __DIM 3
#define __KIND __SINGLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#undef __DIM
         
#define __DIM 4
#define __KIND __SINGLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#undef __DIM
         
#define __DIM 5
#define __KIND __SINGLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io.f"
#undef __KIND
#undef __DIM 

      END MODULE ppm_module_io
