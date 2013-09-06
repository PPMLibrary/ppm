      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                 ppm_vtk_fields
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

#if   __DIM  == __2D
      SUBROUTINE ppm_vtk_fields_2d(filename,Mesh,info,step)
#elif   __DIM  == __3D
      SUBROUTINE ppm_vtk_fields_3d(filename,Mesh,info,step)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      INTEGER, PARAMETER :: MK = ppm_kind_double

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh_),    INTENT(INOUT) :: Mesh
      CHARACTER(LEN=*),           INTENT(IN   ) :: filename
      INTEGER,                    INTENT(  OUT) :: info
      INTEGER,          OPTIONAL, INTENT(IN   ) :: step
      !-------------------------------------------------------------------------
      !  Variables
      !-------------------------------------------------------------------------
      INTEGER                         :: i,j,k,l,isub,ifield,icomp
      TYPE(ppm_t_topo), POINTER       :: topo => NULL()
      REAL(MK), DIMENSION(3)          :: min_phys,max_phys
      REAL(MK), DIMENSION(3)          :: h
      INTEGER, DIMENSION(6)           :: whole_ext,extent
      CHARACTER(LEN=ppm_char)         :: scratch
      !CHARACTER(LEN=ppm_char)         :: imgdatstr
      CHARACTER(LEN=ppm_char)         :: fname,vname
      CLASS(ppm_t_subpatch_),POINTER  :: p => NULL()
      CLASS(ppm_t_subpatch_data_),POINTER  :: pdat => NULL()
      CLASS(ppm_t_discr_data), POINTER  :: field => NULL()
#if   __DIM  == __2D
      REAL(MK),DIMENSION(:,:),POINTER   :: fdata2 => NULL()
      REAL(MK),DIMENSION(:,:,:),POINTER :: fdata3 => NULL()
#elif   __DIM  == __3D
      REAL(MK),DIMENSION(:,:,:),POINTER  :: fdata3 => NULL()
      REAL(MK),DIMENSION(:,:,:,:),POINTER:: fdata4 => NULL()
#endif

      start_subroutine("ppm_vtk_fields")
      !-------------------------------------------------------------------------
      !  Code
      !-------------------------------------------------------------------------

      topo => ppm_topo(Mesh%topoid)%t
      check_associated(topo)

      IF (PRESENT(step)) THEN
         WRITE(fname,'(A,A,I0)') filename(1:LEN_TRIM(filename)), '.', step
      ELSE
         fname = filename
      END IF
      h(:) = 0.0_MK
      min_phys(:) = 0.0_MK
      max_phys(:) = 0.0_MK
      whole_ext(:) = 0
      extent(:) = 0

      If (ASSOCIATED(topo%min_physd)) THEN
         min_phys(1:ppm_dim) = topo%min_physd(1:ppm_dim)
         max_phys(1:ppm_dim) = topo%max_physd(1:ppm_dim)
      ELSE
         min_phys(1:ppm_dim) = topo%min_physs(1:ppm_dim)
         max_phys(1:ppm_dim) = topo%max_physs(1:ppm_dim)
      ENDIF

      h(1:ppm_dim) = Mesh%h(1:ppm_dim)

      DO i=1,2*ppm_dim,2
        whole_ext(i) = 0
        whole_ext(i+1) = Mesh%Nm((i+1)/2)-1
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

          !print only the fields of type real
          field => Mesh%mdata%begin()
          DO WHILE (ASSOCIATED(field))
              IF (field%data_type .EQ. ppm_type_real) THEN
                  IF (field%lda.EQ.1) THEN
                      WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                          TRIM(ADJUSTL(field%name)),"' type='Float64' />"
                  ELSE
                      DO icomp=1,field%lda
                          WRITE(iUnit,'(3A,I0,A)') "      <PDataArray Name='", &
                              TRIM(ADJUSTL(field%name)),"_",icomp,"' type='Float64' />"
                      ENDDO
                  ENDIF
              ENDIF
              field => Mesh%mdata%next()
          ENDDO
          WRITE(iUnit,'(A)') "    </PPointData>"
          ! find the basename of the file
          !DO l=1,topo%nsubs

          l = 0
          p => Mesh%subpatch%begin()
          DO WHILE (ASSOCIATED(p))
              l = l+1
              WRITE(iUnit,'(A)', advance='no')     "    <Piece"
              WRITE(iUnit, '(A)', advance='no') " Extent='"
              DO i=1,ppm_dim
                  WRITE(scratch, '(I0,A,I0)') p%istart(i)-1,' ',&
                  &                           p%istart(i)+p%nnodes(i)-2
                  WRITE(iUnit, '(A)', advance='no') TRIM(ADJUSTL(scratch))
                  IF (i .LT. 3) WRITE(iUnit, '(A)', advance='no') " "
              END DO
              IF (ppm_dim.NE.3) THEN
                  WRITE(iUnit, '(A)', advance='no') '0 0'
              ENDIF
              WRITE(iUnit, '(A)', advance='no') "'"
              WRITE(iUnit, '(A,A,A,I0,A)')" Source='",     &
              fname(INDEX(fname, '/', .true.)+1:LEN_TRIM(fname)), &
              ".", l, ".vti' />"
              p => Mesh%subpatch%next()
          END DO
          ! close
#include "vtk/print_end_header.f"
          CLOSE(iUnit)
      END IF

      !DO l=1,topo%nsublist
      l = 0
      p => Mesh%subpatch%begin()
      DO WHILE (ASSOCIATED(p))
          !isub = topo%isublist(l)
          l = l + 1
          ! append rank to name
          WRITE(scratch,'(A,A,I0,A)') fname(1:LEN_TRIM(fname)), &
          '.', l, '.vti'

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
              extent(i-1) = p%istart(i/2)-1
              extent(i) = p%istart(i/2)+p%nnodes(i/2)-2
          ENDDO
#define VTK_FILE_TYPE "ImageData"
#define VTK_WHOLE_EXTENT whole_ext
#define VTK_ORIGIN min_phys
#define VTK_SPACING h
#define VTK_EXTENT extent
#include "vtk/print_header.f"

          WRITE(iUnit,'(A)',advance='no') "      <PointData"
          WRITE(iUnit,'(A)',advance='no') " Scalars='"


          pdat => p%subpatch_data%begin()
          DO WHILE (ASSOCIATED(pdat))
              IF (pdat%discr_data%data_type .EQ. ppm_type_real) THEN
                  IF (pdat%discr_data%lda.EQ.1) THEN
                          WRITE(iUnit,'(A)',advance='no') &
                              TRIM(ADJUSTL(pdat%discr_data%name))
                          WRITE(iUnit,'(A)',advance='no') " "
                  ELSE
                      DO icomp=1,pdat%discr_data%lda
                          WRITE(iUnit,'(2A,I0)',advance='no') &
                              TRIM(ADJUSTL(pdat%discr_data%name)),"_",icomp
                          WRITE(iUnit,'(A)',advance='no') " "
                      ENDDO
                  ENDIF

              ENDIF
              pdat => p%subpatch_data%next()
          ENDDO
          !FIXME
          ! add this?
          ! BACKSPACE(iUnit)

          WRITE(iUnit,'(A)') "'>"

          !write data for all fields discretized here
          pdat => p%subpatch_data%begin()
          DO WHILE (ASSOCIATED(pdat))
              IF (pdat%discr_data%data_type .EQ. ppm_type_real) THEN
                  IF (pdat%discr_data%lda.EQ.1) THEN !scalar
#if   __DIM == __2D
                      fdata2 =>  pdat%data_2d_rd
#elif   __DIM == __3D
                      fdata3 =>  pdat%data_3d_rd
#endif
#define VTK_NAME pdat%discr_data%name
#define VTK_TYPE "Float64"
#if   __DIM == __2D
#define VTK_MESH fdata2
#elif   __DIM == __3D
#define VTK_MESH fdata3
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)
#endif
#undef VTK_MESH_COMPONENT
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)
#include "vtk/print_data_array.f"
                  ELSE !vector field
#if   __DIM == __2D
                      fdata3 =>  pdat%data_3d_rd
#elif   __DIM == __3D
                      fdata4 =>  pdat%data_4d_rd
#endif
                      DO icomp=1,pdat%discr_data%lda
                          WRITE(vname,'(2A,I0)') &
                              TRIM(ADJUSTL(pdat%discr_data%name)),"_",icomp
#define VTK_NAME vname
#define VTK_TYPE "Float64"
#if   __DIM == __2D
#define VTK_MESH fdata3
#elif   __DIM == __3D
#define VTK_MESH fdata4
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)
#endif
#define VTK_MESH_COMPONENT icomp
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)
#include "vtk/print_data_array.f"
                      ENDDO
                  ENDIF
              ENDIF
              pdat => p%subpatch_data%next()
          ENDDO ! pdat
          WRITE(iUnit,'(A)') "      </PointData>"
         ! close
#include "vtk/print_end_header.f"
         ! close file
         CLOSE(iUnit)
         p => Mesh%subpatch%next()
      ENDDO ! subpatch

      end_subroutine()

#if   __DIM == __2D
      END SUBROUTINE ppm_vtk_fields_2d
#elif   __DIM == __3D
      END SUBROUTINE ppm_vtk_fields_3d
#endif

