      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_data
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

      MODULE ppm_module_data
      !!! Declares global data types and variables.
      !!!
      !!! [NOTE]
      !!! Most of the declared variables in this module should not be accessed
      !!! directly by the PPM client developer, they are used internally in the
      !!! library.

         IMPLICIT NONE
         !----------------------------------------------------------------------
         !
         !----------------------------------------------------------------------
         INCLUDE 'ppm_param.h'

         !----------------------------------------------------------------------
         !  Global TYPEs
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !  internal particle lists
         !----------------------------------------------------------------------
         INTEGER,               DIMENSION(:,:), POINTER :: list_sub => NULL()
         INTEGER,               DIMENSION(:  ), POINTER :: store_info => NULL()
         INTEGER,               DIMENSION(4)            :: ppm_rmsh_kernelsize

         !----------------------------------------------------------------------
         !  Precision
         !----------------------------------------------------------------------
         INTEGER                                        :: ppm_kind
         !!! Precision
         INTEGER                                        :: ppm_mpi_kind
         !!! MPI Precision

         !----------------------------------------------------------------------
         !  Dimensionality
         !----------------------------------------------------------------------
         INTEGER                                        :: ppm_dim
         !!! Dimensionality

         !----------------------------------------------------------------------
         !  Debugging
         !----------------------------------------------------------------------
         INTEGER                                        :: ppm_debug = 0
         !!! are we in Debugging mode?

         !----------------------------------------------------------------------
         !  Has ppm_init been called?
         !----------------------------------------------------------------------
         LOGICAL                                        :: ppm_initialized = .FALSE.
         !!! Has ppm_init been called?

         !----------------------------------------------------------------------
         !  parallel variables
         !----------------------------------------------------------------------
         INTEGER                                        :: ppm_nproc
         INTEGER                                        :: ppm_rank
         INTEGER                                        :: ppm_comm
         ! relative speeds of the processors (for load balancing)
         REAL(ppm_kind_double), DIMENSION(:),   POINTER :: ppm_proc_speed => NULL()

         !----------------------------------------------------------------------
         !  Numerical tolerance. Differences smaller than this are considered
         !  zero.
         !----------------------------------------------------------------------
         REAL(ppm_kind_double)                          :: ppm_myepsd
         REAL(ppm_kind_single)                          :: ppm_myepss

         !----------------------------------------------------------------------
         !  Constants (computed in ppm_init)
         !----------------------------------------------------------------------
         REAL(ppm_kind_double)                          :: ppm_pi_d
         REAL(ppm_kind_single)                          :: ppm_pi_s

         !----------------------------------------------------------------------
         !  I/O Units
         !----------------------------------------------------------------------
         INTEGER                                        :: ppm_stdout = 6
         INTEGER                                        :: ppm_stderr = 0
         INTEGER                                        :: ppm_logfile = -1

         !----------------------------------------------------------------------
         !  Global parameter
         !----------------------------------------------------------------------
         INTEGER, PARAMETER                             :: ppm_big_i=HUGE(1)

      END MODULE ppm_module_data
