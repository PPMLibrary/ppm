      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_io
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
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
      !!! This module provides the interfaces to all PPM I/O routines.
      !!!
      !!! This module contains all data structures and definitions that
      !!! are private to the I/O routines.
         
         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------
#include "ppm_define.h"

         !----------------------------------------------------------------------
         !  Work arrays
         !----------------------------------------------------------------------

         REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: abuf_s
         REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: abuf_d
         COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: abuf_sc
         COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: abuf_dc
         INTEGER                 , DIMENSION(:), POINTER :: abuf_i
         LOGICAL                 , DIMENSION(:), POINTER :: abuf_l
         
         PRIVATE :: abuf_s, abuf_d, abuf_sc, abuf_dc, abuf_i, abuf_l

         REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: rbuf_s
         REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: rbuf_d
         COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: rbuf_sc
         COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: rbuf_dc
         INTEGER                 , DIMENSION(:), POINTER :: rbuf_i
         LOGICAL                 , DIMENSION(:), POINTER :: rbuf_l

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
         !  Define interface to ppm_io
         !----------------------------------------------------------------------
         INTERFACE ppm_io
             !!! Performs parallel I/O of scalars and vectors
             !!! of Fortran intrinsic types
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
         !  Define interface to ppm_doio
         !----------------------------------------------------------------------
         INTERFACE ppm_doio
         !!! Internal routine that performs actual I/O
             MODULE PROCEDURE ppm_doio_s
             MODULE PROCEDURE ppm_doio_d
             MODULE PROCEDURE ppm_doio_i
             MODULE PROCEDURE ppm_doio_l
             MODULE PROCEDURE ppm_doio_sc
             MODULE PROCEDURE ppm_doio_dc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_open
         !----------------------------------------------------------------------
         INTERFACE ppm_io_open
             !!! Opens a file
             MODULE PROCEDURE ppm_io_open
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_close
         !----------------------------------------------------------------------
         INTERFACE ppm_io_close
             !!! Closes a file
             MODULE PROCEDURE ppm_io_close
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_inquire
         !----------------------------------------------------------------------
         INTERFACE ppm_io_inquire
             !!! Inquires the state of a I/O Unit
             MODULE PROCEDURE ppm_io_inquire
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_read_ascii
         !----------------------------------------------------------------------
         INTERFACE ppm_io_read_ascii
             !!! Reads ASCII arrays from an I/O Unit
             MODULE PROCEDURE io_readascii_s
             MODULE PROCEDURE io_readascii_d
             MODULE PROCEDURE io_readascii_i
             MODULE PROCEDURE io_readascii_l
             MODULE PROCEDURE io_readascii_sc
             MODULE PROCEDURE io_readascii_dc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_read_binary
         !----------------------------------------------------------------------
         INTERFACE ppm_io_read_binary
             !!! Reads binary arrays from an I/O Unit
             MODULE PROCEDURE io_readbin_s
             MODULE PROCEDURE io_readbin_d
             MODULE PROCEDURE io_readbin_i
             MODULE PROCEDURE io_readbin_l
             MODULE PROCEDURE io_readbin_sc
             MODULE PROCEDURE io_readbin_dc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_write_ascii
         !----------------------------------------------------------------------
         INTERFACE ppm_io_write_ascii
             !!! Writes ASCII arrays to an I/O Unit
             MODULE PROCEDURE io_writeascii_s
             MODULE PROCEDURE io_writeascii_d
             MODULE PROCEDURE io_writeascii_i
             MODULE PROCEDURE io_writeascii_l
             MODULE PROCEDURE io_writeascii_sc
             MODULE PROCEDURE io_writeascii_dc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_write_binary
         !----------------------------------------------------------------------
         INTERFACE ppm_io_write_binary
             !!! Writes binary arrays to an I/O Unit
             MODULE PROCEDURE io_writebin_s
             MODULE PROCEDURE io_writebin_d
             MODULE PROCEDURE io_writebin_i
             MODULE PROCEDURE io_writebin_l
             MODULE PROCEDURE io_writebin_sc
             MODULE PROCEDURE io_writebin_dc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_delete
         !----------------------------------------------------------------------
         INTERFACE ppm_io_delete
         !!! Deletes a file
             MODULE PROCEDURE ppm_io_delete
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_set_unit
         !----------------------------------------------------------------------
         INTERFACE ppm_io_set_unit
         !!! Sets I/O units for stout, stderr and log file
             MODULE PROCEDURE ppm_io_set_unit
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_unused_unit
         !----------------------------------------------------------------------
         INTERFACE ppm_io_unused_unit
         !!! Routine that finds the first unused I/O unit
             MODULE PROCEDURE ppm_io_unused_unit
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM 0
#define __KIND __SINGLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __CHARACTER
#include "io/ppm_io.f"
#undef __KIND
#undef __DIM

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#undef __DIM
         
#define __DIM 3
#define __KIND __SINGLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#undef __DIM
         
#define __DIM 4
#define __KIND __SINGLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#undef __DIM
         
#define __DIM 5
#define __KIND __SINGLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io.f"
#undef __KIND
#undef __DIM 

#define __KIND __SINGLE_PRECISION
#include "io/ppm_doio.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_doio.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_doio.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_doio.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_doio.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_doio.f"
#undef __KIND

#include "io/ppm_io_open.f"

#include "io/ppm_io_close.f"

#include "io/ppm_io_inquire.f"

#define __KIND __SINGLE_PRECISION
#include "io/ppm_io_read_ascii.f"
#include "io/ppm_io_read_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io_read_ascii.f"
#include "io/ppm_io_read_binary.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io_read_ascii.f"
#include "io/ppm_io_read_binary.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io_read_ascii.f"
#include "io/ppm_io_read_binary.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io_read_ascii.f"
#include "io/ppm_io_read_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io_read_ascii.f"
#include "io/ppm_io_read_binary.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "io/ppm_io_write_ascii.f"
#include "io/ppm_io_write_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "io/ppm_io_write_ascii.f"
#include "io/ppm_io_write_binary.f"
#undef __KIND
#define __KIND __INTEGER
#include "io/ppm_io_write_ascii.f"
#include "io/ppm_io_write_binary.f"
#undef __KIND
#define __KIND __LOGICAL
#include "io/ppm_io_write_ascii.f"
#include "io/ppm_io_write_binary.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "io/ppm_io_write_ascii.f"
#include "io/ppm_io_write_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "io/ppm_io_write_ascii.f"
#include "io/ppm_io_write_binary.f"
#undef __KIND

#include "io/ppm_io_delete.f"

#include "io/ppm_io_set_unit.f"

#include "io/ppm_io_unused_unit.f"

      END MODULE ppm_module_io
