      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_data_io
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
      !  $Log: ppm_module_data_io.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/07/26 11:48:10  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.1  2004/07/26 07:28:17  ivos
      !  Initial implementation. Originated from splitting the old ppm
      !  modules.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_io

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_single, ppm_kind_double 
         PRIVATE :: ppm_kind_single, ppm_kind_double

         !----------------------------------------------------------------------
         !  Mode (parallel or serial) for each unit (internal numbering)
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER          :: ppm_io_mode 

         !----------------------------------------------------------------------
         !  IO format (ASCII or binary) for each unit (internal numbering)
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER          :: ppm_io_format

         !----------------------------------------------------------------------
         !  Inverse list of currently used unit numbers by ppm_module_io.
         !  Index is Fortran Unit number, value is index to above lists or
         !  0 if that particular Fortran unit is not used by ppm_io.
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER          :: ppm_io_unit

         !----------------------------------------------------------------------
         !  I/O buffers for single and double precision types
         !----------------------------------------------------------------------
         REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: sbuffer
         REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: dbuffer
         COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: scbuffer
         COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: dcbuffer

      END MODULE ppm_module_data_io
