      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_data_io
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_io
      !!! This module contains data structures and definitions that
      !!! are used in the PPM I/O routines.
      !!!
      !!! [NOTE]
      !!! The members of this modules should not be accessed directly by the PPM
      !!! client developer. They are managed internally by the library.
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
