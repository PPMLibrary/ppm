      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                  ppm_module_mpi
      !-------------------------------------------------------------------------
      ! Copyright (c) 2014 CSE Lab (ETH Zurich), MOSAIC Group (MPI-CBG Dresden),
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------

      MODULE ppm_module_mpi
      !----------------------------------------------------------------------
      ! Module for MPI calls
      !----------------------------------------------------------------------
#ifdef __MPI
!           USE mpi
!           USE mpi_f08
#endif
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
#ifdef __MPI
          INCLUDE "mpif.h"
#else
          ! TYPES
          INTEGER, PARAMETER :: MPI_LOGICAL          = 1
          INTEGER, PARAMETER :: MPI_CHARACTER        = 1
          INTEGER, PARAMETER :: MPI_INTEGER          = 1
          INTEGER, PARAMETER :: MPI_REAL             = 1
          INTEGER, PARAMETER :: MPI_REAL8            = 1
          INTEGER, PARAMETER :: MPI_REAL16           = 1
          INTEGER, PARAMETER :: MPI_DOUBLE_PRECISION = 1
          ! OPERANDS
          INTEGER, PARAMETER :: MPI_LAND = 1
          INTEGER, PARAMETER :: MPI_MAX  = 1
          INTEGER, PARAMETER :: MPI_MIN  = 1
          INTEGER, PARAMETER :: MPI_SUM  = 1
          ! GLOBAL
          INTEGER, PARAMETER :: MPI_UNDEFINED   = -1
          INTEGER, PARAMETER :: MPI_COMM_NULL   =  0
          INTEGER, PARAMETER :: MPI_COMM_SELF   =  1
          INTEGER, PARAMETER :: MPI_COMM_WORLD  =  1
          INTEGER, PARAMETER :: MPI_STATUS_SIZE =  1
#endif

      END MODULE ppm_module_mpi