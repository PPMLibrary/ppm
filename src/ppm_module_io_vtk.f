      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_io_vtk
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

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 4
#define __2D               2
#define __3D               3
#define __VEC              5
#define __SCA              6

      MODULE ppm_module_io_vtk
         USE ppm_module_data,      ONLY: ppm_char, ppm_error_fatal, &
         &   ppm_kind_single, ppm_kind_double,ppm_error_error,      &
         &   ppm_rank,ppm_nproc, ppm_comm
         USE ppm_module_error,     ONLY: ppm_error, ppm_err_argument
         USE ppm_module_substart,  ONLY: substart
         USE ppm_module_substop,   ONLY: substop
         USE ppm_module_interfaces
         USE ppm_module_topo_typedef
         IMPLICIT NONE

         PUBLIC :: ppm_vtk_particles,ppm_vtk_mesh_2d,ppm_vtk_mesh_3d
         PRIVATE
         !----------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------

#ifdef __MPI
         INCLUDE 'mpif.h'
#endif
         !----------------------------------------------------------------------
         !  New Interface
         !----------------------------------------------------------------------
          INTERFACE ppm_vtk_particles
             MODULE PROCEDURE ppm_vtk_particles_s
             MODULE PROCEDURE ppm_vtk_particles_d
          END INTERFACE

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

#define __KIND __SINGLE_PRECISION
#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "vtk/ppm_vtk_particles.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "vtk/ppm_vtk_particles.f"
#undef __KIND

#define __DIM __2D
#include "vtk/ppm_vtk_mesh.f"
#undef __DIM
#define __DIM __3D
#include "vtk/ppm_vtk_mesh.f"
#undef __DIM
         !----------------------------------------------------------------------
         !  Parallel VTK output
         !----------------------------------------------------------------------
         ! todo

         !----------------------------------------------------------------------
         !  Ugly interface
         !----------------------------------------------------------------------

         SUBROUTINE ppm_vtk_init(filename, vtype, info, &
         &          version, byte_order, whole_extent,  &
         &          origin, spacing, extent, npoints,   &
         &          nverts, nlines, nstrips, npolys)
           ! args
           CHARACTER(LEN=*),                              INTENT(IN   ) :: filename
           CHARACTER(LEN=*),                              INTENT(IN   ) :: vtype
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
           vtk_type = vtype
           current_section = 0
           ! open output file
           OPEN(iUnit, FILE=filename, IOSTAT=info, ACTION='WRITE')
           IF (info .NE. 0) THEN
              info = ppm_error_fatal
              WRITE(errtxt,'(2A)') 'Failed to open file: ', &
              & filename(1:LEN_TRIM(filename))
              CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
              GOTO 9999
           END IF
           WRITE(iUnit,'(A)')  "<?xml version='1.0' ?>"
           WRITE(iUnit,'(7A)') "<VTKFile type='", vtype, &
           &                   "' version='",     iversion(1:LEN_TRIM(iversion)),  &
           &                   "' byte_order='",  ibyte_order(1:LEN_TRIM(ibyte_order)), "'>"
           WRITE(iUnit,'(2A)',ADVANCE='NO') "  <", vtype
           IF (PRESENT(whole_extent)) THEN
              WRITE(iUnit, '(A)', ADVANCE='NO') " WholeExtent='"
              DO i=LBOUND(whole_extent,1),UBOUND(whole_extent,1)
                 WRITE(scratch, *) whole_extent(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', ADVANCE='NO') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(whole_extent,1)) WRITE(iUnit, '(A)', ADVANCE='NO') " "
              END DO
              WRITE(iUnit, '(A)', ADVANCE='NO') "'"
           END IF
           IF (PRESENT(origin)) THEN
              WRITE(iUnit, '(A)', ADVANCE='NO') " Origin='"
              DO i=LBOUND(origin,1),UBOUND(origin,1)
                 WRITE(scratch, *) origin(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', ADVANCE='NO') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(origin,1)) &
                      WRITE(iUnit, '(A)', ADVANCE='NO') " "
              END DO
              WRITE(iUnit, '(A)', ADVANCE='NO') "'"
           END IF
           IF (PRESENT(spacing)) THEN
              WRITE(iUnit, '(A)', ADVANCE='NO') " Spacing='"
              DO i=LBOUND(spacing,1),UBOUND(spacing,1)
                 WRITE(scratch, *) spacing(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', ADVANCE='NO') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(spacing,1)) &
                      WRITE(iUnit, '(A)', ADVANCE='NO') " "
              END DO
              WRITE(iUnit, '(A)', ADVANCE='NO') "'"
           END IF
           WRITE(iUnit,'(A)') ">"
           WRITE(iUnit,'(A)', ADVANCE='NO')     "    <Piece"
           IF (PRESENT(extent)) THEN
              WRITE(iUnit, '(A)', ADVANCE='NO') " Extent='"
              DO i=LBOUND(extent,1),UBOUND(extent,1)
                 WRITE(scratch, *) extent(i)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', ADVANCE='NO') &
                      scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(extent,1)) &
                      WRITE(iUnit, '(A)', ADVANCE='NO') " "
              END DO
              WRITE(iUnit, '(A)', ADVANCE='NO') "'"
           ELSE
              IF (PRESENT(npoints)) THEN
                 WRITE(scratch, *) npoints
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfPoints='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              END IF
              IF (PRESENT(nverts)) THEN
                 WRITE(scratch, *) nverts
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfVerts='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              ENd IF
              IF (PRESENT(nlines)) THEN
                 WRITE(scratch, *) nlines
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfLines='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              END IF
              IF (PRESENT(nstrips)) THEN
                 WRITE(scratch, *) nstrips
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfStrips='", &
                      scratch(1:LEN_TRIM(scratch)), "'"
              END IF
              IF (PRESENT(npolys)) THEN
                 WRITE(scratch, *) npolys
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfPolys='", &
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
           WRITE(iUnit,'(2A)', ADVANCE='NO') "      <", &
                vtk_section(1:LEN_TRIM(vtk_section))
           IF (PRESENT(attr)) THEN
              WRITE(iUnit,'(2A)', ADVANCE='NO') " ", attr(1:LEN_TRIM(attr))
           END IF
           WRITE (iUnit, '(A)') '>'
         END SUBROUTINE ppm_vtk_section

         SUBROUTINE ppm_vtk_data(name, data, info, format, &
         &          vtype, ncomponents, offset)
           ! args
           CHARACTER(LEN=*),                              INTENT(IN   ) :: name
           REAL(ppm_kind_single), DIMENSION(:),           INTENT(IN   ) :: data
           INTEGER,                                       INTENT(  OUT) :: info
           CHARACTER(LEN=*),                    OPTIONAL, INTENT(IN   ) :: format
           CHARACTER(LEN=*),                    OPTIONAL, INTENT(IN   ) :: vtype
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
           WRITE(iUnit,'(A)',ADVANCE='NO') "        <DataArray"
           IF (PRESENT(vtype)) THEN
              WRITE(iUnit,'(3A)',ADVANCE='NO') " type='", vtype(1:LEN_TRIM(vtype)), "'"
           END IF
           WRITE(iUnit,'(3A)',ADVANCE='NO') " Name='", name(1:LEN_TRIM(name)), "'"
           IF (PRESENT(ncomponents)) THEN
              WRITE(scratch, *) ncomponents
              scratch = ADJUSTL(scratch)
              WRITE(iUnit,'(3A)',ADVANCE='NO') " NumberOfComponents='", &
                   scratch(1:LEN_TRIM(scratch)), "'"
           END IF
           WRITE(iUnit,'(3A)',ADVANCE='NO') " format='", iformat(1:LEN_TRIM(iformat)), "'"
           IF (PRESENT(offset)) THEN
              WRITE(scratch, *) offset
              scratch = ADJUSTL(scratch)
              WRITE(iUnit,'(3A)',ADVANCE='NO') " offset='", &
                   scratch(1:LEN_TRIM(scratch)), "'"
           END IF
           WRITE(iUnit,'(A)') ">"
           DO i=LBOUND(data,1),UBOUND(data,1)
              WRITE(scratch, *) data(i)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', ADVANCE='NO') &
                   scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(data,1)) &
                   WRITE(iUnit, '(A)', ADVANCE='NO') " "
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
