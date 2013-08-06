      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_data_io
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
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

         IMPLICIT NONE

         PRIVATE :: ppm_kind_single, ppm_kind_double

         !----------------------------------------------------------------------
         !  Mode (parallel or serial) for each unit (internal numbering)
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: ppm_io_mode  => NULL()

         !----------------------------------------------------------------------
         !  IO format (ASCII or binary) for each unit (internal numbering)
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: ppm_io_format  => NULL()

         !----------------------------------------------------------------------
         !  Inverse list of currently used unit numbers by ppm_module_io.
         !  Index is Fortran Unit number, value is index to above lists or
         !  0 if that particular Fortran unit is not used by ppm_io.
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: ppm_io_unit  => NULL()

         !----------------------------------------------------------------------
         !  I/O buffers for single and double precision types
         !----------------------------------------------------------------------
         REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: sbuffer  => NULL()
         REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: dbuffer  => NULL()
         COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: scbuffer  => NULL()
         COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: dcbuffer  => NULL()

      END MODULE ppm_module_data_io
