      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_io_vtk
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
      MODULE ppm_module_io_vtk
         USE ppm_module_typedef,   ONLY: ppm_char, ppm_error_fatal, &
                                         ppm_kind_single, ppm_kind_double
         USE ppm_module_error,     ONLY: ppm_error, ppm_err_argument
         USE ppm_module_data,      ONLY: ppm_rank, ppm_nproc, ppm_comm
         USE ppm_module_typedef,   ONLY: ppm_error_error
         USE ppm_module_substart,  ONLY: substart
         USE ppm_module_substop,   ONLY: substop
         USE ppm_module_particles, ONLY: ppm_t_particles,  &
                                         get_xp, set_xp,   &
                                         get_wps, set_wps, &
                                         get_wpv, set_wpv

         IMPLICIT NONE

         PUBLIC :: ppm_vtk_particle_cloud, ppm_vtk_grid_like, &
              ppm_vtk_init, ppm_vtk_section, ppm_vtk_data, ppm_vtk_close
         PRIVATE
         !----------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------
#include "ppm_define.h"

#ifdef __MPI
  INCLUDE 'mpif.h'
#endif
         !----------------------------------------------------------------------
         !  New Interface
         !----------------------------------------------------------------------
!          INTERFACE ppm_vtk_particle_cloud
!             MODULE PROCEDURE ppm_vtk_particle_cloud
!          END INTERFACE

!          INTERFACE ppm_vtk_grid_like
!             MODULE PROCEDURE ppm_vtk_grid_like
!          END INTERFACE

         !----------------------------------------------------------------------
         !  Ugly Interface
         !----------------------------------------------------------------------
!          INTERFACE ppm_vtk_init
!             MODULE PROCEDURE ppm_vtk_init_s
!             MODULE PROCEDURE ppm_vtk_init_d
!          END INTERFACE

!          INTERFACE ppm_vtk_data
!             MODULE PROCEDURE ppm_vtk_data_INTEGER
!             MODULE PROCEDURE ppm_vtk_data_LONGINT
!             MODULE PROCEDURE ppm_vtk_data_SINGLE
!             MODULE PROCEDURE ppm_vtk_data_DOUBLE
!          END INTERFACE

         !----------------------------------------------------------------------
         !  Variables
         !----------------------------------------------------------------------
         INTEGER                 :: iUnit = 123
         CHARACTER(LEN=ppm_char) :: errtxt
         CHARACTER(LEN=ppm_char) :: scratch

         !  Ugly interface
         CHARACTER(LEN=ppm_char) :: vtk_type
         CHARACTER(LEN=ppm_char) :: vtk_section
         INTEGER                 :: current_section

       CONTAINS

         !----------------------------------------------------------------------
         !  New Interface
         !----------------------------------------------------------------------

         SUBROUTINE ppm_vtk_particle_cloud(filename, Particles, info, &
              with_ghosts, wps_list, wpv_list)
           !--------------------------------------------------------------------
           !  Arguments
           !--------------------------------------------------------------------
           CHARACTER(LEN=*),                INTENT(IN   ) :: filename
           TYPE(ppm_t_particles), POINTER,  INTENT(IN   ) :: Particles
           INTEGER,                         INTENT(  OUT) :: info
           LOGICAL,               OPTIONAL, INTENT(IN   ) :: with_ghosts
           INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wps_list
           INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: wpv_list
           !--------------------------------------------------------------------
           !  Variables
           !--------------------------------------------------------------------
           CHARACTER(LEN=*), PARAMETER :: caller='ppm_vtk_particle_cloud'
           CHARACTER(LEN=ppm_char)              :: scratch
           REAL(ppm_kind_double)                :: t0
           INTEGER                              :: i, j, k, nd, N
           INTEGER                              :: nb_wps, nb_wpv
           INTEGER, DIMENSION(:),   ALLOCATABLE :: wps_l, wpv_l
           LOGICAL                              :: ghosts
           REAL(8), DIMENSION(:,:), POINTER     :: xp  => NULL()
           REAL(8), DIMENSION(:),   POINTER     :: wp  => NULL()
           !--------------------------------------------------------------------
           !  Code
           !--------------------------------------------------------------------
           CALL substart(caller,t0,info)

           ! print ghosts?
           IF (PRESENT(with_ghosts)) THEN
              ghosts = with_ghosts
           ELSE
              ghosts = .FALSE.
           END IF

           ! create the list of properties to print
           IF(PRESENT(wps_list)) THEN
              nb_wps=SIZE(wps_list)
              ALLOCATE(wps_l(nb_wps),STAT=info)
              wps_l=wps_list
              DO i=1,nb_wps
                 IF (wps_l(i).GT.Particles%max_wpsid) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'property index exceeds size of property array',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
                 IF (Particles%wps_m(wps_l(i)).NE.1) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'trying to printout a property that is not mapped &
                         & to the particles',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
              ENDDO
           ELSE
              !printout all properties i for which wps_m(i) = 1
              nb_wps = 0
              DO i=1,Particles%max_wpsid
                 IF (Particles%wps_m(i).EQ.1) &
                      nb_wps = nb_wps + 1
              ENDDO
              ALLOCATE(wps_l(nb_wps),STAT=info)
              nb_wps = 0
              DO i=1,Particles%max_wpsid
                 IF (Particles%wps_m(i).EQ.1) THEN
                    nb_wps = nb_wps + 1
                    wps_l(nb_wps) = i
                 ENDIF
              ENDDO
           ENDIF
           IF(PRESENT(wpv_list)) THEN
              nb_wpv=SIZE(wpv_list)
              ALLOCATE(wpv_l(nb_wpv),STAT=info)
              wpv_l=wpv_list
              DO i=1,nb_wpv
                 IF (wpv_l(i).GT.Particles%max_wpvid) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'property index exceeds size of property array',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
                 IF (Particles%wpv_m(wpv_l(i)).NE.1) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,   &
                         &  'trying to printout a property that is not mapped &
                         & to the particles',&
                         &  __LINE__,info)
                    GOTO 9999
                 ENDIF
              ENDDO
           ELSE
              !printout all properties i for which wpv_m(i) = 1
              nb_wpv = 0
              DO i=1,Particles%max_wpvid
                 IF (Particles%wpv_m(i).EQ.1) &
                      nb_wpv = nb_wpv + 1
              ENDDO
              ALLOCATE(wpv_l(nb_wpv),STAT=info)
              nb_wpv = 0
              DO i=1,Particles%max_wpvid
                 IF (Particles%wpv_m(i).EQ.1) THEN
                    nb_wpv = nb_wpv + 1
                    wpv_l(nb_wpv) = i
                 ENDIF
              ENDDO
           ENDIF

#ifdef __MPI
           ! write parallel file
           IF (ppm_rank .EQ. 0) THEN
              WRITE(scratch,'(A,A)') filename(1:LEN_TRIM(filename)), '.pvtp'
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
              DO i=1,nb_wps
              WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                   Particles%wps(wps_l(i))%name &
                   (1:LEN_TRIM(Particles%wps(wps_l(i))%name)), &
                   "' type='Float64' />"
              END DO
              DO i=1,nb_wpv
              WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                   Particles%wpv(wpv_l(i))%name &
                   (1:LEN_TRIM(Particles%wpv(wpv_l(i))%name)), &
                   "' type='Float64' />"
              END DO
              WRITE(iUnit,'(A)') "    </PPointData>"              
              WRITE(iUnit,'(A)') "    <PPoints>"
              WRITE(iUnit,'(A)') "      <PDataArray NumberOfComponents='3' type='Float64' />"
              WRITE(iUnit,'(A)') "    </PPoints>"
              DO i=0,ppm_nproc-1
                 WRITE(iUnit,'(A,A,A,I0,A)') "    <Piece Source='", &
                      filename(1:LEN_TRIM(filename)), ".", i, ".vtp' />"
              END DO
              ! close
#include "vtk/print_end_header.f"
              CLOSE(iUnit)
           END IF
           ! append rank to name
           WRITE(scratch,'(A,A,I0,A)') filename(1:LEN_TRIM(filename)), &
                                       '.', ppm_rank, '.vtp'
#else
           WRITE(scratch,'(A,A)') filename(1:LEN_TRIM(filename)), '.vtp'
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
              N = Particles%Mpart
           ELSE
              N = Particles%Npart
           END IF

           ! write header
#define VTK_FILE_TYPE "PolyData"
#define VTK_NPOINTS N
#define VTK_NVERTS  N
#include "vtk/print_header.f"

           ! print properties
           IF (nb_wps .GT. 0 .OR. nb_wpv .GT. 0) THEN

              ! print names
              WRITE(iUnit,'(A)',advance='no') "      <PointData" 
              IF (nb_wps .GT. 0) THEN
                 WRITE(iUnit,'(A)',advance='no') " Scalars='"
                 DO i=1,nb_wps
                    WRITE(iUnit,'(A)',advance='no') &
                         Particles%wps(wps_l(i))%name &
                         (1:LEN_TRIM(Particles%wps(wps_l(i))%name))
                    IF (i .LT. nb_wps) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              IF (nb_wpv .GT. 0) THEN
                 WRITE(iUnit,'(A)',advance='no') "' Vectors='"
                 DO i=1,nb_wpv
                    WRITE(iUnit,'(A)',advance='no') &
                         Particles%wpv(wpv_l(i))%name &
                         (1:LEN_TRIM(Particles%wpv(wpv_l(i))%name))
                    IF (i .LT. nb_wpv) WRITE(iUnit,'(A)',advance='no') " "
                 END DO
              END IF
              WRITE(iUnit,'(A)') "'>"

              ! property values
              DO k=1,nb_wps
                 wp => get_wps(Particles,wps_l(k),with_ghosts=ghosts)
#define VTK_NAME Particles%wps(wps_l(k))%name
#define VTK_TYPE "Float64"
#define VTK_SCALAR wp
#include "vtk/print_data_array.f"
                 wp => set_wps(Particles,wps_l(k),read_only=.TRUE.)
              END DO
              DO k=1,nb_wpv
                 xp => get_wpv(Particles,wpv_l(k),with_ghosts=ghosts)
                 nd = SIZE(xp,1)
#define VTK_NAME Particles%wpv(wpv_l(k))%name
#define VTK_TYPE "Float64"
#define VTK_NDIM "3"
#define VTK_VECTOR xp
#define APPEND_ZEROS
#include "vtk/print_data_array.f"
                 xp => set_wpv(Particles,wpv_l(i),read_only=.TRUE.)
              END DO
              WRITE(iUnit,'(A)') "      </PointData>"
           END IF

           ! print point coordinates
           WRITE(iUnit,'(A)') "      <Points>"
           xp => get_xp(Particles,with_ghosts=ghosts)
           nd = SIZE(xp,1)
#define VTK_TYPE "Float64"
#define VTK_NDIM "3"
#define VTK_VECTOR xp
#define APPEND_ZEROS
#include "vtk/print_data_array.f"
           xp => set_xp(Particles,read_only=.TRUE.)
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
         END SUBROUTINE ppm_vtk_particle_cloud

         SUBROUTINE ppm_vtk_grid_like()
         END SUBROUTINE ppm_vtk_grid_like

         !----------------------------------------------------------------------
         !  Parallel VTK output
         !----------------------------------------------------------------------
         ! todo

         !----------------------------------------------------------------------
         !  Ugly interface
         !----------------------------------------------------------------------

         SUBROUTINE ppm_vtk_init(filename, type, info,              &
                                 version, byte_order, whole_extent, &
                                 origin, spacing, extent,           &
                                 npoints, nverts, nlines, nstrips, npolys)
           ! args
           CHARACTER(LEN=*),                              INTENT(IN   ) :: filename
           CHARACTER(LEN=*),                              INTENT(IN   ) :: type
           INTEGER,                                       INTENT(  OUT) :: info
           CHARACTER(LEN=*),                    OPTIONAL, INTENT(IN   ) :: version
           CHARACTER(LEN=*),                    OPTIONAL, INTENT(IN   ) :: byte_order
           REAL(ppm_kind_single), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: whole_extent
           REAL(ppm_kind_single), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: origin
           REAL(ppm_kind_single), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: spacing
           REAL(ppm_kind_single), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: extent
           INTEGER,                             OPTIONAL, INTENT(IN   ) :: npoints
           INTEGER,                             OPTIONAL, INTENT(IN   ) :: nverts
           INTEGER,                             OPTIONAL, INTENT(IN   ) :: nlines
           INTEGER,                             OPTIONAL, INTENT(IN   ) :: nstrips
           INTEGER,                             OPTIONAL, INTENT(IN   ) :: npolys
           ! vars
           CHARACTER(LEN=*), PARAMETER :: caller = 'ppm_vtk_init'
           INTEGER                     :: i
           CHARACTER(LEN=ppm_char)     :: iversion
           CHARACTER(LEN=ppm_char)     :: ibyte_order
           ! defaults
           IF (.NOT. PRESENT(version)) THEN
              iversion = "0.1"
           ELSE
              iversion = version
           END IF
           IF (.NOT. PRESENT(byte_order)) THEN
              ibyte_order = "LittleEndian"
           ELSE
              ibyte_order = byte_order
           END IF
           ! save type for ppm_vtk_close and reset current_section
           vtk_type = type
           current_section = 0
           ! open output file
           OPEN(iUnit, FILE=filename, IOSTAT=info, ACTION='WRITE')
           IF (info .NE. 0) THEN
              info = ppm_error_fatal
              WRITE(errtxt,'(2A)') 'Failed to open file: ', &
                    filename(1:LEN_TRIM(filename))
              CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
              GOTO 9999
           END IF
           WRITE(iUnit,'(A)')  "<?xml version='1.0' ?>"
           WRITE(iUnit,'(7A)') "<VTKFile type='", type, &
                               "' version='",     iversion(1:LEN_TRIM(iversion)),  &
                               "' byte_order='",  ibyte_order(1:LEN_TRIM(ibyte_order)), "'>"
           WRITE(iUnit,'(2A)',advance='no') "  <", type
           IF (PRESENT(whole_extent)) THEN
              WRITE(iUnit, '(A)', advance='no') " WholeExtent='"
              DO i=LBOUND(whole_extent,1),UBOUND(whole_extent,1)
                 WRITE(scratch, *) whole_extent(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', advance='no') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(whole_extent,1)) &
                      WRITE(iUnit, '(A)', advance='no') " "
              END DO
              WRITE(iUnit, '(A)', advance='no') "'"
           END IF
           IF (PRESENT(origin)) THEN
              WRITE(iUnit, '(A)', advance='no') " Origin='"
              DO i=LBOUND(origin,1),UBOUND(origin,1)
                 WRITE(scratch, *) origin(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', advance='no') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(origin,1)) &
                      WRITE(iUnit, '(A)', advance='no') " "
              END DO
              WRITE(iUnit, '(A)', advance='no') "'"
           END IF
           IF (PRESENT(spacing)) THEN
              WRITE(iUnit, '(A)', advance='no') " Spacing='"
              DO i=LBOUND(spacing,1),UBOUND(spacing,1)
                 WRITE(scratch, *) spacing(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', advance='no') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(spacing,1)) &
                      WRITE(iUnit, '(A)', advance='no') " "
              END DO
              WRITE(iUnit, '(A)', advance='no') "'"
           END IF
           WRITE(iUnit,'(A)') ">"
           WRITE(iUnit,'(A)', advance='no')     "    <Piece"
           IF (PRESENT(extent)) THEN
              WRITE(iUnit, '(A)', advance='no') " Extent='"
              DO i=LBOUND(extent,1),UBOUND(extent,1)
                 WRITE(scratch, *) extent(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', advance='no') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(extent,1)) &
                      WRITE(iUnit, '(A)', advance='no') " "
              END DO
              WRITE(iUnit, '(A)', advance='no') "'"
           ELSE
              IF (PRESENT(npoints)) THEN
                 WRITE(scratch, *) npoints
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', advance='no') " NumberOfPoints='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              END IF
              IF (PRESENT(nverts)) THEN
                 WRITE(scratch, *) nverts
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', advance='no') " NumberOfVerts='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              ENd IF
              IF (PRESENT(nlines)) THEN
                 WRITE(scratch, *) nlines
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', advance='no') " NumberOfLines='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              END IF
              IF (PRESENT(nstrips)) THEN
                 WRITE(scratch, *) nstrips
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', advance='no') " NumberOfStrips='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              END IF
              IF (PRESENT(npolys)) THEN
                 WRITE(scratch, *) npolys
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', advance='no') " NumberOfPolys='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              END IF
           END IF
           WRITE(iUnit, '(A)') ">"
9999       CONTINUE
         END SUBROUTINE ppm_vtk_init

         SUBROUTINE ppm_vtk_section(name, attr)
           CHARACTER(LEN=*),           INTENT(IN   ) :: name
           CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: attr
           IF (current_section .GT. 0) THEN
              ! close previous section
              WRITE(iUnit,'(3A)') "      </", vtk_section(1:LEN_TRIM(vtk_section)), ">"
           END IF
           vtk_section = ADJUSTL(name)
           current_section = current_section + 1
           WRITE(iUnit,'(2A)', advance='no') "      <", &
                vtk_section(1:LEN_TRIM(vtk_section))
           IF (PRESENT(attr)) THEN
              WRITE(iUnit,'(2A)', advance='no') " ", attr(1:LEN_TRIM(attr))
           END IF
           WRITE (iUnit, '(A)') '>'
         END SUBROUTINE ppm_vtk_section

         SUBROUTINE ppm_vtk_data(name, data, info, &
                                 format, type, ncomponents, offset)
           ! args
           CHARACTER(LEN=*),                              INTENT(IN   ) :: name
           REAL(ppm_kind_single), DIMENSION(:),           INTENT(IN   ) :: data
           INTEGER,                                       INTENT(  OUT) :: info
           CHARACTER(LEN=*),                    OPTIONAL, INTENT(IN   ) :: format
           CHARACTER(LEN=*),                    OPTIONAL, INTENT(IN   ) :: type
           INTEGER,                             OPTIONAL, INTENT(IN   ) :: ncomponents
           INTEGER,                             OPTIONAL, INTENT(IN   ) :: offset
           ! vars
           CHARACTER(LEN=*),       PARAMETER :: caller = 'ppm_vtk_data'
           INTEGER                           :: i
           CHARACTER(LEN=ppm_char)           :: iformat
           ! defaults
           IF (.NOT. PRESENT(format)) THEN
              iformat = ADJUSTL('ascii')
           ELSE
              iformat = ADJUSTL(format)
           END IF
           IF (current_section .EQ. 0) THEN
              info = ppm_error_fatal
              WRITE(errtxt,'(A)') 'All vtk data must be inside a section!'
              CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
              GOTO 9999
           END IF
           WRITE(iUnit,'(A)',advance='no') "        <DataArray"
           IF (PRESENT(type)) THEN
              WRITE(iUnit,'(3A)',advance='no') " type='", type(1:LEN_TRIM(type)), "'"
           END IF
           WRITE(iUnit,'(3A)',advance='no') " Name='", name(1:LEN_TRIM(name)), "'"
           IF (PRESENT(ncomponents)) THEN
              WRITE(scratch, *) ncomponents
              scratch = ADJUSTL(scratch)
              WRITE(iUnit,'(3A)',advance='no') " NumberOfComponents='", &
                   scratch(1:LEN_TRIM(scratch)), "'"
           END IF
           WRITE(iUnit,'(3A)',advance='no') " format='", iformat(1:LEN_TRIM(iformat)), "'"
           IF (PRESENT(offset)) THEN
              WRITE(scratch, *) offset
              scratch = ADJUSTL(scratch)
              WRITE(iUnit,'(3A)',advance='no') " offset='", &
                   scratch(1:LEN_TRIM(scratch)), "'"
           END IF
           WRITE(iUnit,'(A)') ">"
           DO i=LBOUND(data,1),UBOUND(data,1)
              WRITE(scratch, *) data(i)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', advance='no') &
                   scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(data,1)) &
                   WRITE(iUnit, '(A)', advance='no') " "
           END DO
           WRITE(iUnit,'(/A)') "        </DataArray>"
9999       CONTINUE
         END SUBROUTINE ppm_vtk_data

         SUBROUTINE ppm_vtk_close()
           IF (current_section .GT. 0) THEN
              ! close open section
              WRITE(iUnit,'(3A)') "      </", vtk_section(1:LEN_TRIM(vtk_section)), ">"
           END IF
#define VTK_FILE_TYPE vtk_type
#include "vtk/print_end_header.f"
#undef VTK_FILE_TYPE
           ! close file
           CLOSE(iUnit)
         END SUBROUTINE ppm_vtk_close

      END MODULE ppm_module_io_vtk
