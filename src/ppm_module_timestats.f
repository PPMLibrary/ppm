      !-------------------------------------------------------------------------
      !  Module       :                     ppm_module_timestats
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

      MODULE ppm_module_timestats
         USE ppm_module_error,    ONLY : ppm_error, ppm_err_argument
         USE ppm_module_data,     ONLY : ppm_rank, ppm_nproc, ppm_comm,&
                                       & ppm_char, ppm_error_fatal, &
                                       & ppm_kind_single, ppm_kind_double,&
                                       & ppm_error_error
         USE ppm_module_substart, ONLY : substart
         USE ppm_module_substop,  ONLY : substop
         IMPLICIT NONE

         TYPE ppm_t_tstats
             REAL(ppm_kind_double), DIMENSION(:), POINTER :: times => NULL()
             CHARACTER(LEN=ppm_char)                      :: label
             ! only allocated and filled when ppm_tstats_collect is called
             REAL(ppm_kind_double), DIMENSION(:), POINTER :: tmin => NULL()
             REAL(ppm_kind_double), DIMENSION(:), POINTER :: tmax => NULL()
             REAL(ppm_kind_double), DIMENSION(:), POINTER :: tmean => NULL()
         END TYPE

         INTEGER                                      :: ppm_tstats_nsamples = 100
         INTEGER                                      :: ppm_ntstats = 0
         TYPE(ppm_t_tstats), DIMENSION(:), POINTER    :: ppm_tstats => NULL()


         CONTAINS

#include "util/ppm_util_timestats_setup.f"
#include "util/ppm_util_tic.f"
#include "util/ppm_util_toc.f"

      END MODULE ppm_module_timestats
