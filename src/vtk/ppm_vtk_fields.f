      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                 ppm_vtk_fields
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

#if   __DIM  == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_vtk_fields_2ds(topoid,meshid,nfields,fields,ghostsize,&
      &                         filename,info,step)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_vtk_fields_2dd(topoid,meshid,nfields,fields,ghostsize,&
      &                         filename,info,step)
#endif
#elif   __DIM  == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_vtk_fields_3ds(topoid,meshid,nfields,fields,ghostsize,&
      &                         filename,info,step)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_vtk_fields_3dd(topoid,meshid,nfields,fields,ghostsize,&
      &                         filename,info,step)
#endif
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
    
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,            INTENT(IN)            :: topoid
      INTEGER,            INTENT(IN)            :: meshid
      INTEGER,            INTENT(IN)            :: nfields
#if   __DIM  == __2D
#if   __KIND == __SINGLE_PRECISION
      TYPE(ppm_t_field_2ds),  DIMENSION(:), POINTER :: fields
#elif __KIND == __DOUBLE_PRECISION
      TYPE(ppm_t_field_2dd),  DIMENSION(:), POINTER :: fields
#endif
#elif   __DIM  == __3D
#if   __KIND == __SINGLE_PRECISION
      TYPE(ppm_t_field_3ds),  DIMENSION(:), POINTER :: fields
#elif __KIND == __DOUBLE_PRECISION
      TYPE(ppm_t_field_3dd),  DIMENSION(:), POINTER :: fields
#endif
#endif
      INTEGER, DIMENSION(:)                     :: ghostsize
      CHARACTER(LEN=*),   INTENT(IN)            :: filename
      INTEGER,            INTENT(OUT)           :: info
      INTEGER,               OPTIONAL, INTENT(IN   ) :: step
      !-------------------------------------------------------------------------
      !  Variables
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), PARAMETER     :: caller = 'ppm_vtk_fields'
      INTEGER                         :: i,j,k,l,isub,ifield
      REAL(MK)                        :: t0
      TYPE(ppm_t_topo), POINTER       :: topo => NULL()
      TYPE(ppm_t_equi_mesh), POINTER  :: mesh => NULL()
      REAL(MK), DIMENSION(3)          :: min_phys,max_phys
      REAL(MK), DIMENSION(3)          :: h
      INTEGER, DIMENSION(6)           :: whole_ext,extent
      CHARACTER(LEN=ppm_char)         :: scratch
      CHARACTER(LEN=ppm_char)         :: imgdatstr
      CHARACTER(LEN=ppm_char)         :: fname

      !--------------------------------------------------------------------
      !  Code
      !--------------------------------------------------------------------
      CALL substart(caller,t0,info)
     
      ! TODO check parameters

      topo => ppm_topo(topoid)%t
      SELECT TYPE (t => ppm_mesh%vec(meshid))
      TYPE IS (ppm_t_equi_mesh)
          mesh => t
      END SELECT

      IF (PRESENT(step)) THEN
         WRITE(fname,'(A,A,I0)') &
              filename(1:LEN_TRIM(filename)), '.', step
      ELSE
         fname = filename
      END IF
      h(:) = 0
      min_phys(:) = 0
      max_phys(:) = 0
      whole_ext(:) = 0
      extent(:) = 0
#if   __KIND == __SINGLE_PRECISION
      min_phys(1:ppm_dim) = topo%min_physs(1:ppm_dim)
      max_phys(1:ppm_dim) = topo%max_physs(1:ppm_dim)
#elif __KIND == __DOUBLE_PRECISION
      min_phys(1:ppm_dim) = topo%min_physd(1:ppm_dim)
      max_phys(1:ppm_dim) = topo%max_physd(1:ppm_dim)
#endif
      h(1:ppm_dim) = (max_phys(1:ppm_dim) - min_phys(1:ppm_dim)) / mesh%Nm(1:ppm_dim)
      DO i=1,2*ppm_dim,2 
        whole_ext(i) = 0
        whole_ext(i+1) = mesh%Nm((i+1)/2)-1
      ENDDO
      ! write parallel file
      IF (ppm_rank .EQ. 0) THEN
          WRITE(scratch,'(A,A)') fname(1:LEN_TRIM(fname)), '.pvti'
          OPEN(iUnit, FILE=scratch(1:LEN_TRIM(scratch)), &
          IOSTAT=info, ACTION='WRITE')
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              WRITE(errtxt,'(2A)') 'Failed to open file: ', &
              scratch(1:LEN_TRIM(scratch))
              CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
              GOTO 9999
          END IF
#define VTK_FILE_TYPE "PImageData"
#define VTK_WHOLE_EXTENT whole_ext
#define VTK_GHOSTLEVEL 0
#define VTK_ORIGIN min_phys
#define VTK_SPACING h
#define VTK_PARALLEL
#include "vtk/print_header.f"
          WRITE(iUnit,'(A)') "    <PPointData Scalars='scalars'>"              
          DO l=1,nfields
              WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
              fields(l)%fname &
              (1:LEN_TRIM(fields(l)%fname)), &
              "' type='Float64' />"
          END DO
          WRITE(iUnit,'(A)') "    </PPointData>"              
          ! find the basename of the file
          DO l=1,topo%nsubs
              WRITE(iUnit,'(A)', advance='no')     "    <Piece"
              WRITE(iUnit, '(A)', advance='no') " Extent='"
              DO i=1,ppm_dim
                  WRITE(scratch, '(I0,A,I0)') mesh%istart(i,l)-1,' ',&
                  &                           mesh%istart(i,l)+mesh%nnodes(i,l)-2
                  scratch = ADJUSTL(scratch)
                  WRITE(iUnit, '(A)', advance='no') scratch(1:LEN_TRIM(scratch))
                  IF (i .LT. 3) WRITE(iUnit, '(A)', advance='no') " "
              END DO
              IF (ppm_dim.NE.3) THEN
                  WRITE(iUnit, '(A)', advance='no') '0 0'
              ENDIF
              WRITE(iUnit, '(A)', advance='no') "'"
              WRITE(iUnit, '(A,A,A,I0,A)')" Source='",     &
              fname(INDEX(fname, '/', .true.)+1:LEN_TRIM(fname)), &
              ".", l, ".vti' />"
          END DO
          ! close
#include "vtk/print_end_header.f"
          CLOSE(iUnit)
      END IF

      DO l=1,topo%nsublist
          isub = topo%isublist(l)
          ! append rank to name
          WRITE(scratch,'(A,A,I0,A)') fname(1:LEN_TRIM(fname)), &
          '.', isub, '.vti'

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
          ! write header
          DO i=2,2*ppm_dim,2
              extent(i-1) = mesh%istart(i/2,isub)-1
              extent(i) = mesh%istart(i/2,isub)+mesh%nnodes(i/2,isub)-2
          ENDDO
#define VTK_FILE_TYPE "ImageData"
#define VTK_WHOLE_EXTENT whole_ext
#define VTK_ORIGIN min_phys
#define VTK_SPACING h
#define VTK_EXTENT extent
#include "vtk/print_header.f"

          WRITE(iUnit,'(A)',advance='no') "      <PointData" 
          WRITE(iUnit,'(A)',advance='no') " Scalars='"
          DO ifield=1,nfields
              WRITE(iUnit,'(A)',advance='no') &
              fields(ifield)%fname(1:LEN_TRIM(fields(ifield)%fname))
              IF (ifield .LT. nfields) WRITE(iUnit,'(A)',advance='no') " "
          ENDDO
          WRITE(iUnit,'(A)') "'>"
          DO ifield=1,nfields
          ! write data
#define VTK_NAME fields(ifield)%fname
#define VTK_TYPE "Float64"
#define VTK_MESH fields(ifield)%fdata
#define VTK_MESH_SUB l
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND mesh%nnodes(1,isub)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND mesh%nnodes(2,isub)
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND mesh%nnodes(3,isub)
#include "vtk/print_data_array.f"
          ENDDO
          WRITE(iUnit,'(A)') "      </PointData>"
         ! close
#include "vtk/print_end_header.f"
         ! close file
         CLOSE(iUnit)
      ENDDO ! subs

9999  CONTINUE

      CALL substop(caller,t0,info)
      RETURN
#if   __DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_vtk_fields_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_vtk_fields_2dd
#endif
#elif   __DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_vtk_fields_3ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_vtk_fields_3dd
#endif
#endif

