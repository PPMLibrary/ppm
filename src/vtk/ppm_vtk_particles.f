      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                 ppm_vtk_particles
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_vtk_particles_s(topoid,particles,filename,info,&
      &                              step,with_ghosts)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_vtk_particles_d(topoid,particles,filename,info,&
      &                              step,with_ghosts)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_topo
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    
      !--------------------------------------------------------------------
      !  Arguments
      !--------------------------------------------------------------------
      INTEGER,                         INTENT(IN   ) :: topoid
#if   __KIND == __SINGLE_PRECISION
      TYPE(ppm_t_particles_s), POINTER,INTENT(IN   ) :: particles
#elif __KIND == __DOUBLE_PRECISION
      TYPE(ppm_t_particles_d), POINTER,INTENT(IN   ) :: particles
#endif
      CHARACTER(LEN=*),                INTENT(IN   ) :: filename
      INTEGER,                         INTENT(  OUT) :: info
      INTEGER,               OPTIONAL, INTENT(IN   ) :: step
      LOGICAL,               OPTIONAL, INTENT(IN   ) :: with_ghosts
      !--------------------------------------------------------------------
      !  Variables
      !--------------------------------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: caller='ppm_vtk_particles'
      CHARACTER(LEN=ppm_char)              :: scratch
      CHARACTER(LEN=ppm_char)              :: fname
      REAL(ppm_kind_double)                :: t0
      INTEGER                              :: i, j, k, l, nd, N
      INTEGER                              :: nb_wps, nb_wpv, nb_wpv_field
      INTEGER                              :: nb_wpi
      LOGICAL                              :: ghosts
      TYPE(ppm_t_topo), POINTER            :: topo => NULL()
      REAL(MK), DIMENSION(:,:), POINTER     :: xp  => NULL()
      REAL(MK), DIMENSION(:),   POINTER     :: wp  => NULL()
      INTEGER, DIMENSION(:),   POINTER      :: wpi  => NULL()
      !--------------------------------------------------------------------
      !  Code
      !--------------------------------------------------------------------
      CALL substart(caller,t0,info)
      
      topo => ppm_topo(topoid)%t
      IF (PRESENT(step)) THEN
         WRITE(fname,'(A,A,I0)') &
              filename(1:LEN_TRIM(filename)), '.', step
      ELSE
         fname = filename
      END IF

      ghosts = .FALSE.
      ! print ghosts?
      IF (PRESENT(with_ghosts)) THEN
         ghosts = with_ghosts
      ELSE
         ghosts = .FALSE.
      END IF



      nb_wpi = particles%niprop
      nb_wps = particles%nprop




#ifdef __MPI
      ! write parallel file
      IF (ppm_rank .EQ. 0) THEN
         WRITE(scratch,'(A,A)') fname(1:LEN_TRIM(fname)), '.pvtp'
         OPEN(iUnit, FILE=scratch(1:LEN_TRIM(scratch)), &
              IOSTAT=info, ACTION='WRITE')
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            WRITE(errtxt,'(2A)') 'Failed to open file: ', &
                 scratch(1:LEN_TRIM(scratch))
            CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
            GOTO 9999
         END IF
#define VTK_FILE_TYPE "PPolyData"
#define VTK_PARALLEL
#include "vtk/print_header.f"
         WRITE(iUnit,'(A)') "    <PPointData>"
         DO i=1,nb_wpi
         WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
              particles%iprop(i)%name &
              (1:LEN_TRIM(particles%iprop(i)%name)), &
              "' type='Float64' />"
         END DO
         DO i=1,nb_wps
         WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
              particles%prop(i)%name &
              (1:LEN_TRIM(particles%prop(i)%name)), &
              "' type='Float64' />"
         END DO
         WRITE(iUnit,'(A)') "    </PPointData>"              
         WRITE(iUnit,'(A)') "    <PPoints>"
         WRITE(iUnit,'(A)') "      <PDataArray NumberOfComponents='3' type='Float64' />"
         WRITE(iUnit,'(A)') "    </PPoints>"
         ! find the basename of the file
         DO i=0,ppm_nproc-1
            WRITE(iUnit,'(A,A,A,I0,A)') "    <Piece Source='",     &
                 fname(INDEX(fname, '/', .true.)+1:LEN_TRIM(fname)), &
                 ".", i, ".vtp' />"
         END DO
         ! close
#include "vtk/print_end_header.f"
         CLOSE(iUnit)
      END IF
      ! append rank to name
      WRITE(scratch,'(A,A,I0,A)') fname(1:LEN_TRIM(fname)), &
                                  '.', ppm_rank, '.vtp'
#else
      WRITE(scratch,'(A,A)') fname(1:LEN_TRIM(fname)), '.vtp'
#endif

      ! open output file
      OPEN(iUnit, FILE=scratch(1:LEN_TRIM(scratch)), &
           IOSTAT=info, ACTION='WRITE')
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         WRITE(errtxt,'(2A)') 'Failed to open file: ', &
               scratch(1:LEN_TRIM(scratch))
         CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
         GOTO 9999
      END IF

      ! write data
      IF (ghosts) THEN
         N = particles%mp
      ELSE
         N = particles%np
      END IF
      ! write header
#define VTK_FILE_TYPE "PolyData"
#define VTK_NPOINTS N
#define VTK_NVERTS  N
#include "vtk/print_header.f"

      ! print properties
      IF (nb_wpi .GT. 0 .OR. nb_wps .GT. 0) THEN

         ! print names
         WRITE(iUnit,'(A)',advance='no') "      <PointData" 
         IF (nb_wpi .GT. 0) THEN
            WRITE(iUnit,'(A)',advance='no') " Integers='"
            DO i=1,nb_wpi
               WRITE(iUnit,'(A)',advance='no') &
                    particles%iprop(i)%name &
                    & (1:LEN_TRIM(particles%iprop(i)%name))
               IF (i .LT. nb_wpi) WRITE(iUnit,'(A)',advance='no') " "
            END DO
         END IF
         IF (nb_wps .GT. 0) THEN
            IF (nb_wpi .GT. 0) &
                 WRITE(iUnit,'(A)',advance='no') "'"
            WRITE(iUnit,'(A)',advance='no') " Scalars='"
            DO i=1,nb_wps
               WRITE(iUnit,'(A)',advance='no') &
                    particles%prop(i)%name &
                    (1:LEN_TRIM(particles%prop(i)%name))
               IF (i .LT. nb_wps) WRITE(iUnit,'(A)',advance='no') " "
            END DO
         END IF
         WRITE(iUnit,'(A)') "'>"

         DO k=1,nb_wpi
            wpi => particles%iprop(k)%wp
#define VTK_NAME particles%iprop(k)%name
#define VTK_TYPE "Float64"
#define VTK_INTEGER wpi
#include "vtk/print_data_array.f"
            wpi => particles%iprop(k)%wp
            !wpi => set_wpi(Particles,wpi_l(k),read_only=.TRUE.)
         END DO
         DO k=1,nb_wps
            wp => particles%prop(k)%wp
#define VTK_NAME particles%prop(k)%name
#define VTK_TYPE "Float64"
#define VTK_SCALAR wp
#include "vtk/print_data_array.f"
            wp => particles%prop(k)%wp
            !wp => set_wps(Particles,wps_l(k),read_only=.TRUE.)
         END DO
         WRITE(iUnit,'(A)') "      </PointData>"
      END IF

      ! print point coordinates
      WRITE(iUnit,'(A)') "      <Points>"
      xp => particles%xp
      nd =  ppm_dim
#define VTK_TYPE "Float64"
#define VTK_NDIM "3"
#define VTK_VECTOR xp
#define APPEND_ZEROS
#include "vtk/print_data_array.f"
      xp => particles%xp
      WRITE(iUnit,'(A)') "      </Points>"

      ! create a vertex for every point
      WRITE(iUnit,'(A)') "      <Verts>"
      ! connectivity
      N = N - 1
#define VTK_RANGE N
#define VTK_RANGE_START 0
#define VTK_NAME "connectivity"
#define VTK_TYPE "Int32"
#include "vtk/print_data_array.f"
      ! offsets
      N = N + 1
#define VTK_RANGE N
#define VTK_NAME "offsets"
#define VTK_TYPE "Int32"
#include "vtk/print_data_array.f"
      WRITE(iUnit,'(A)') "      </Verts>"

      ! close
#include "vtk/print_end_header.f"
      ! close file
9999       CONTINUE
      CLOSE(iUnit)
      CALL substop(caller,t0,info)
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_vtk_particles_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_vtk_particles_d
#endif
