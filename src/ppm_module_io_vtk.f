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

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4
#define __VEC              5
#define __SCA              6

      MODULE ppm_module_io_vtk
         USE ppm_module_data,      ONLY: ppm_char, ppm_error_fatal, &
                                         ppm_kind_single, ppm_kind_double,&
                                         ppm_error_error
         USE ppm_module_error,     ONLY: ppm_error, ppm_err_argument
         USE ppm_module_data,      ONLY: ppm_rank, ppm_nproc, ppm_comm
         USE ppm_module_substart,  ONLY: substart
         USE ppm_module_substop,   ONLY: substop
         USE ppm_module_interfaces
         USE ppm_module_topo_typedef
         USE ppm_module_mesh_typedef

         IMPLICIT NONE

         PUBLIC :: ppm_vtk_particles, ppm_vtk_fields, &
              ppm_t_field_2ds,ppm_t_field_2dd,ppm_t_field_3ds,ppm_t_field_3dd, &
              ppm_t_particles_s,ppm_t_particles_d, &
              ppm_t_prop_s,ppm_t_prop_d,ppm_t_prop_i
         PRIVATE
         !----------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------

#ifdef __MPI
  INCLUDE 'mpif.h'
#endif
         !----------------------------------------------------------------------
         !  Types
         !----------------------------------------------------------------------

         TYPE ppm_t_field_2dd
             REAL(ppm_kind_double), DIMENSION(:,:,:), POINTER :: fdata => NULL()
             CHARACTER(LEN=ppm_char)                          :: fname
         END TYPE
         TYPE ppm_t_field_3dd
             REAL(ppm_kind_double), DIMENSION(:,:,:,:), POINTER :: fdata => NULL()
             CHARACTER(LEN=ppm_char)                            :: fname
         END TYPE
         
         TYPE ppm_t_field_2ds
             REAL(ppm_kind_single), DIMENSION(:,:,:), POINTER :: fdata => NULL()
             CHARACTER(LEN=ppm_char)                          :: fname
         END TYPE
         TYPE ppm_t_field_3ds
             REAL(ppm_kind_single), DIMENSION(:,:,:,:), POINTER :: fdata => NULL()
             CHARACTER(LEN=ppm_char)                            :: fname
         END TYPE

         TYPE ppm_t_particles_s
             REAL(ppm_kind_single), DIMENSION(:,:), POINTER     :: xp => NULL()
             INTEGER                                            :: np = 0
             INTEGER                                            :: mp = 0
             TYPE(ppm_t_prop_s),    DIMENSiON(:)  , POINTER     :: prop => NULL()
             TYPE(ppm_t_prop_i),    DIMENSION(:)  , POINTER     :: iprop => NULL()
             INTEGER                                            :: nprop = 0
             INTEGER                                            :: niprop = 0
         END TYPE
         TYPE ppm_t_particles_d
             REAL(ppm_kind_double), DIMENSION(:,:), POINTER     :: xp => NULL()
             INTEGER                                            :: np = 0
             INTEGER                                            :: mp = 0
             TYPE(ppm_t_prop_d),    DIMENSiON(:)  , POINTER     :: prop => NULL()
             TYPE(ppm_t_prop_i),    DIMENSION(:)  , POINTER     :: iprop => NULL()
             INTEGER                                            :: nprop = 0
             INTEGER                                            :: niprop = 0
         END TYPE

         TYPE ppm_t_prop_s
             REAL(ppm_kind_single), DIMENSION(:)  , POINTER     :: wp => NULL()
             CHARACTER(LEN=ppm_char)                            :: name
         END TYPE
         TYPE ppm_t_prop_d
             REAL(ppm_kind_double), DIMENSION(:)  , POINTER     :: wp => NULL()
             CHARACTER(LEN=ppm_char)                            :: name
         END TYPE
         TYPE ppm_t_prop_i
             INTEGER              , DIMENSION(:)  , POINTER     :: wp => NULL()
             CHARACTER(LEN=ppm_char)                            :: name
         END TYPE

         !----------------------------------------------------------------------
         !  New Interface
         !----------------------------------------------------------------------
          INTERFACE ppm_vtk_particles
             MODULE PROCEDURE ppm_vtk_particles_s
             MODULE PROCEDURE ppm_vtk_particles_d
          END INTERFACE

          INTERFACE ppm_vtk_fields
             MODULE PROCEDURE ppm_vtk_fields_2ds
             MODULE PROCEDURE ppm_vtk_fields_2dd
             MODULE PROCEDURE ppm_vtk_fields_3ds
             MODULE PROCEDURE ppm_vtk_fields_3dd
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

#define __KIND __DOUBLE_PRECISION
#include "vtk/ppm_vtk_particles.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION
#include "vtk/ppm_vtk_particles.f"
#undef __KIND

#define __DIM __2D
#define __KIND __DOUBLE_PRECISION
#include "vtk/ppm_vtk_fields.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION
#include "vtk/ppm_vtk_fields.f"
#undef __KIND
#undef __DIM
#define __DIM __3D
#define __KIND __DOUBLE_PRECISION
#include "vtk/ppm_vtk_fields.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION
#include "vtk/ppm_vtk_fields.f"
#undef __KIND
#undef __DIM
         !----------------------------------------------------------------------
         !  Parallel VTK output
         !----------------------------------------------------------------------
         ! todo

         !----------------------------------------------------------------------
         !  Ugly interface
         !----------------------------------------------------------------------

         SUBROUTINE ppm_vtk_init(filename, vtype, info,             &
                                 version, byte_order, whole_extent, &
                                 origin, spacing, extent,           &
                                 npoints, nverts, nlines, nstrips, npolys)
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
                    filename(1:LEN_TRIM(filename))
              CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
              GOTO 9999
           END IF
           WRITE(iUnit,'(A)')  "<?xml version='1.0' ?>"
           WRITE(iUnit,'(7A)') "<VTKFile type='", vtype, &
                               "' version='",     iversion(1:LEN_TRIM(iversion)),  &
                               "' byte_order='",  ibyte_order(1:LEN_TRIM(ibyte_order)), "'>"
           WRITE(iUnit,'(2A)',advance='no') "  <", vtype
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
                                 format, vtype, ncomponents, offset)
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
           WRITE(iUnit,'(A)',advance='no') "        <DataArray"
           IF (PRESENT(vtype)) THEN
              WRITE(iUnit,'(3A)',advance='no') " type='", vtype(1:LEN_TRIM(vtype)), "'"
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
