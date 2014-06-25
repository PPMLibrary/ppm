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
      INTEGER, PARAMETER :: MK  = ppm_kind_double
      INTEGER, PARAMETER :: MKS = ppm_kind_single

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh_), INTENT(INOUT) :: Mesh
      CHARACTER(LEN=*),        INTENT(IN   ) :: filename
      INTEGER,                 INTENT(  OUT) :: info
      INTEGER,       OPTIONAL, INTENT(IN   ) :: step
      !-------------------------------------------------------------------------
      !  Variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

!       CLASS(ppm_t_subpatch_),      POINTER :: p     => NULL()
      CLASS(ppm_t_subpatch_data_), POINTER :: pdat
      CLASS(ppm_t_discr_data),     POINTER :: field

#if   __DIM  == __2D
      REAL(MKS), DIMENSION(:,:),     POINTER :: fdata2_rs
      REAL(MKS), DIMENSION(:,:,:),   POINTER :: fdata3_rs

      REAL(MK),  DIMENSION(:,:),     POINTER :: fdata2_rd
      REAL(MK),  DIMENSION(:,:,:),   POINTER :: fdata3_rd

      INTEGER,   DIMENSION(:,:),     POINTER :: fdata2_i
      INTEGER,   DIMENSION(:,:,:),   POINTER :: fdata3_i
#elif   __DIM  == __3D
      REAL(MKS), DIMENSION(:,:,:),   POINTER :: fdata3_rs
      REAL(MKS), DIMENSION(:,:,:,:), POINTER :: fdata4_rs

      REAL(MK),  DIMENSION(:,:,:),   POINTER :: fdata3_rd
      REAL(MK),  DIMENSION(:,:,:,:), POINTER :: fdata4_rd

      INTEGER,   DIMENSION(:,:,:),   POINTER :: fdata3_i
      INTEGER,   DIMENSION(:,:,:,:), POINTER :: fdata4_i
#endif
      REAL(MK), DIMENSION(3) :: min_phys,max_phys
      REAL(MK), DIMENSION(3) :: h

      INTEGER               :: i,j,k,l,ipatch
      INTEGER               :: isub,sub2proc
      INTEGER               :: ifield,icomp
      INTEGER, DIMENSION(6) :: whole_ext,extent
      INTEGER, DIMENSION(3) :: istart,iend,nc

      CHARACTER(LEN=ppm_char) :: scratch
      CHARACTER(LEN=ppm_char) :: fname,vname

      LOGICAL :: lopen,lexists

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

      h        = 0.0_MK
      min_phys = 0.0_MK
      max_phys = 0.0_MK

      whole_ext = 0
      extent    = 0

      If (ASSOCIATED(topo%min_physd)) THEN
         min_phys(1:ppm_dim) = topo%min_physd(1:ppm_dim)
         max_phys(1:ppm_dim) = topo%max_physd(1:ppm_dim)
      ELSE
         min_phys(1:ppm_dim) = REAL(topo%min_physs(1:ppm_dim),MK)
         max_phys(1:ppm_dim) = REAL(topo%max_physs(1:ppm_dim),MK)
      ENDIF

      h(1:ppm_dim) = Mesh%h(1:ppm_dim)

      DO i=1,2*ppm_dim,2
         whole_ext(i) = 0
         whole_ext(i+1) = Mesh%Nm((i+1)/2)-1
      ENDDO

      ! write parallel file
      IF (ppm_rank.EQ.0) THEN
         WRITE(scratch,'(A,A)') fname(1:LEN_TRIM(fname)), '.pvti'
         OPEN(iUnit, FILE=scratch(1:LEN_TRIM(scratch)), &
         & IOSTAT=info, ACTION='WRITE')
         IF (info.NE.0) THEN
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
            SELECT CASE (field%data_type)
            CASE (ppm_type_real_single)
               IF (field%lda.EQ.1) THEN
                  WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                  & TRIM(ADJUSTL(field%name)),"' type='Float32' />"
               ELSE
                  DO icomp=1,field%lda
                     WRITE(iUnit,'(3A,I0,A)') "      <PDataArray Name='", &
                     & TRIM(ADJUSTL(field%name)),"_",icomp,"' type='Float32' />"
                  ENDDO
               ENDIF

            CASE (ppm_type_real)
               IF (field%lda.EQ.1) THEN
                  WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                  & TRIM(ADJUSTL(field%name)),"' type='Float64' />"
               ELSE
                  DO icomp=1,field%lda
                     WRITE(iUnit,'(3A,I0,A)') "      <PDataArray Name='", &
                     & TRIM(ADJUSTL(field%name)),"_",icomp,"' type='Float64' />"
                  ENDDO
               ENDIF

            CASE (ppm_type_int)
               IF (field%lda.EQ.1) THEN
                  WRITE(iUnit,'(3A)') "      <PDataArray Name='", &
                  & TRIM(ADJUSTL(field%name)),"' type='Int32' />"
               ELSE
                  DO icomp=1,field%lda
                     WRITE(iUnit,'(3A,I0,A)') "      <PDataArray Name='", &
                     & TRIM(ADJUSTL(field%name)),"_",icomp,"' type='Int32' />"
                  ENDDO
               ENDIF

            END SELECT
            field => Mesh%mdata%next()
         ENDDO
         WRITE(iUnit,'(A)') "    </PPointData>"
         ! find the basename of the file
         DO isub=1,topo%nsubs
            DO ipatch=1,Mesh%subpatch_by_sub(isub)%nsubpatch
               SELECT TYPE(p => Mesh%subpatch_by_sub(isub)%vec(ipatch)%t)
               CLASS IS (ppm_t_subpatch_)
                  WRITE(iUnit,'(A)', ADVANCE='NO')     "    <Piece"
                  WRITE(iUnit, '(A)', ADVANCE='NO') " Extent='"
                  DO i=1,ppm_dim
                     !TODO
                     !check the boundary condition
                     !it has been implemented in a wrong way
                     IF (p%bc(2*i).EQ.ppm_param_bcdef_periodic.OR. &
                     &   p%bc(2*i).EQ.-1) THEN
                        WRITE(scratch, '(I0,A,I0)') p%istart(i)-1,' ',&
                        & p%iend(i)
                     ELSE
                        WRITE(scratch, '(I0,A,I0)') p%istart(i)-1,' ',&
                        & p%iend(i)-1
                     ENDIF
                     WRITE(iUnit, '(A)', ADVANCE='NO') TRIM(ADJUSTL(scratch))
                     IF (i .LT. 3) WRITE(iUnit, '(A)', ADVANCE='NO') " "
                  END DO
                  IF (ppm_dim.NE.3) THEN
                     WRITE(iUnit, '(A)', ADVANCE='NO') '0 0'
                  ENDIF
                  WRITE(iUnit, '(A)', ADVANCE='NO') "'"

                  sub2proc=topo%sub2proc(isub)

                  WRITE(iUnit, '(A,A,A,I0,A)')" Source='",              &
                  & fname(INDEX(fname, '/', .TRUE.)+1:LEN_TRIM(fname)), &
                  & ".", sub2proc, ".vti' />"
               END SELECT
            ENDDO ! ipatch
         ENDDO ! isub
         ! close
#include "vtk/print_end_header.f"
         CLOSE(iUnit)
      END IF

      DO l=1,topo%nsublist
         isub=topo%isublist(l)
         IF (Mesh%subpatch_by_sub(isub)%nsubpatch.GT.0) THEN
            sub2proc=topo%sub2proc(isub)
            WRITE(scratch,'(A,A,I0,A)') fname(1:LEN_TRIM(fname)), &
            & '.', sub2proc, '.vti'
            OPEN(iUnit, FILE=scratch(1:LEN_TRIM(scratch)), &
            & IOSTAT=info, ACTION='WRITE')
            IF (info.NE.0) THEN
               info = ppm_error_fatal
               WRITE(errtxt,'(2A)') 'Failed to open file: ', &
               & scratch(1:LEN_TRIM(scratch))
               CALL ppm_error(ppm_err_argument, caller, errtxt, __LINE__, info)
               GOTO 9999
            END IF
#define VTK_FILE_TYPE "ImageData"
#define VTK_WHOLE_EXTENT whole_ext
#define VTK_ORIGIN min_phys
#define VTK_SPACING h
#define VTK_PARALLEL
#include "vtk/print_header.f"
#undef  VTK_PARALLEL
            EXIT
         ELSE
            CYCLE
         ENDIF
      ENDDO

      DO l=1,topo%nsublist
         isub=topo%isublist(l)
         DO ipatch=1,Mesh%subpatch_by_sub(isub)%nsubpatch
            SELECT TYPE(p => Mesh%subpatch_by_sub(isub)%vec(ipatch)%t)
            CLASS IS (ppm_t_subpatch_)
               ! write header
               DO i=2,2*ppm_dim,2
                  extent(i-1) = p%istart(i/2)-1
                  IF (p%bc(i).EQ.ppm_param_bcdef_periodic.OR. &
                  &   p%bc(i).EQ.-1) THEN
                     extent(i) = p%iend(i/2)
                  ELSE
                     extent(i) = p%iend(i/2)-1
                  ENDIF
               ENDDO
#define VTK_EXTENT extent
#include "vtk/print_fields_header1.f"
               WRITE(iUnit,'(A)',ADVANCE='NO') "      <PointData"
               WRITE(iUnit,'(A)',ADVANCE='NO') " Scalars='"

               pdat => p%subpatch_data%begin()
               DO WHILE (ASSOCIATED(pdat))
                  IF (pdat%discr_data%lda.EQ.1) THEN
                     WRITE(iUnit,'(A)',ADVANCE='NO') &
                     & TRIM(ADJUSTL(pdat%discr_data%name))
                     WRITE(iUnit,'(A)',ADVANCE='NO') " "
                  ELSE
                     DO icomp=1,pdat%discr_data%lda
                        WRITE(iUnit,'(2A,I0)',ADVANCE='NO') &
                        & TRIM(ADJUSTL(pdat%discr_data%name)),"_",icomp
                        WRITE(iUnit,'(A)',ADVANCE='NO') " "
                     ENDDO
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
                  SELECT CASE (pdat%discr_data%data_type)
                  CASE (ppm_type_real_single)
                     IF (pdat%discr_data%lda.EQ.1) THEN !scalar
                        nc=0
                        SELECT CASE (p%bc(2))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(1)=1
                        END SELECT
                        SELECT CASE (p%bc(4))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(2)=1
                        END SELECT

#if   __DIM == __2D
                        fdata2_rs =>  pdat%data_2d_rs
#elif   __DIM == __3D
                        SELECT CASE (p%bc(6))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(3)=1
                        END SELECT

                        fdata3_rs =>  pdat%data_3d_rs
#endif

#define VTK_NAME pdat%discr_data%name
#define VTK_TYPE "Float32"
#if   __DIM == __2D
#define VTK_MESH fdata2_rs
#elif   __DIM == __3D
#define VTK_MESH fdata3_rs
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)+nc(3)
#endif
#undef VTK_MESH_COMPONENT
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)+nc(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)+nc(2)
#include "vtk/print_data_array.f"
                     ELSE !vector field
                        nc=0
                        SELECT CASE (p%bc(2))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(1)=1
                        END SELECT
                        SELECT CASE (p%bc(4))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(2)=1
                        END SELECT

#if   __DIM == __2D
                        fdata3_rs =>  pdat%data_3d_rs
#elif   __DIM == __3D
                        SELECT CASE (p%bc(6))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(3)=1
                        END SELECT

                        fdata4_rs =>  pdat%data_4d_rs
#endif
                        DO icomp=1,pdat%discr_data%lda
                           WRITE(vname,'(2A,I0)') &
                           & TRIM(ADJUSTL(pdat%discr_data%name)),"_",icomp
#define VTK_NAME vname
#define VTK_TYPE "Float32"
#if   __DIM == __2D
#define VTK_MESH fdata3_rs
#elif   __DIM == __3D
#define VTK_MESH fdata4_rs
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)+nc(3)
#endif
#define VTK_MESH_COMPONENT icomp
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)+nc(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)+nc(2)
#include "vtk/print_data_array.f"
                        ENDDO ! icomp=1,pdat%discr_data%lda
                     ENDIF ! (pdat%discr_data%lda.EQ.?)

                  CASE (ppm_type_real)
                     IF (pdat%discr_data%lda.EQ.1) THEN !scalar
                        nc=0
                        SELECT CASE (p%bc(2))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(1)=1
                        END SELECT
                        SELECT CASE (p%bc(4))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(2)=1
                        END SELECT

#if   __DIM == __2D
                        fdata2_rd =>  pdat%data_2d_rd
#elif   __DIM == __3D
                        SELECT CASE (p%bc(6))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(3)=1
                        END SELECT

                        fdata3_rd =>  pdat%data_3d_rd
#endif

#define VTK_NAME pdat%discr_data%name
#define VTK_TYPE "Float64"
#if   __DIM == __2D
#define VTK_MESH fdata2_rd
#elif   __DIM == __3D
#define VTK_MESH fdata3_rd
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)+nc(3)
#endif
#undef VTK_MESH_COMPONENT
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)+nc(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)+nc(2)
#include "vtk/print_data_array.f"
                     ELSE !vector field
                        nc=0
                        SELECT CASE (p%bc(2))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(1)=1
                        END SELECT
                        SELECT CASE (p%bc(4))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(2)=1
                        END SELECT

#if   __DIM == __2D
                        fdata3_rd =>  pdat%data_3d_rd
#elif   __DIM == __3D
                        SELECT CASE (p%bc(6))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(3)=1
                        END SELECT

                        fdata4_rd =>  pdat%data_4d_rd
#endif
                        DO icomp=1,pdat%discr_data%lda
                           WRITE(vname,'(2A,I0)') &
                           & TRIM(ADJUSTL(pdat%discr_data%name)),"_",icomp
#define VTK_NAME vname
#define VTK_TYPE "Float64"
#if   __DIM == __2D
#define VTK_MESH fdata3_rd
#elif   __DIM == __3D
#define VTK_MESH fdata4_rd
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)+nc(3)
#endif
#define VTK_MESH_COMPONENT icomp
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)+nc(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)+nc(2)
#include "vtk/print_data_array.f"
                        ENDDO ! icomp=1,pdat%discr_data%lda
                     ENDIF ! (pdat%discr_data%lda.EQ.?)

                  CASE (ppm_type_int)
                     IF (pdat%discr_data%lda.EQ.1) THEN !scalar
                        nc=0
                        SELECT CASE (p%bc(2))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(1)=1
                        END SELECT
                        SELECT CASE (p%bc(4))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(2)=1
                        END SELECT

#if   __DIM == __2D
                        fdata2_i =>  pdat%data_2d_i
#elif   __DIM == __3D
                        SELECT CASE (p%bc(6))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(3)=1
                        END SELECT

                        fdata3_i =>  pdat%data_3d_i
#endif

#define VTK_NAME pdat%discr_data%name
#define VTK_TYPE "Int32"
#if   __DIM == __2D
#define VTK_MESH fdata2_i
#elif   __DIM == __3D
#define VTK_MESH fdata3_i
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)+nc(3)
#endif
#undef VTK_MESH_COMPONENT
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)+nc(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)+nc(2)
#include "vtk/print_data_array.f"
                     ELSE !vector field
                        nc=0
                        SELECT CASE (p%bc(2))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(1)=1
                        END SELECT
                        SELECT CASE (p%bc(4))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(2)=1
                        END SELECT

#if   __DIM == __2D
                        fdata3_i =>  pdat%data_3d_i
#elif   __DIM == __3D
                        SELECT CASE (p%bc(6))
                        CASE (ppm_param_bcdef_periodic,-1)
                           nc(3)=1
                        END SELECT

                        fdata4_i =>  pdat%data_4d_i
#endif
                        DO icomp=1,pdat%discr_data%lda
                           WRITE(vname,'(2A,I0)') &
                           & TRIM(ADJUSTL(pdat%discr_data%name)),"_",icomp
#define VTK_NAME vname
#define VTK_TYPE "Int32"
#if   __DIM == __2D
#define VTK_MESH fdata3_i
#elif   __DIM == __3D
#define VTK_MESH fdata4_i
#define VTK_MESH_KLBOUND 1
#define VTK_MESH_KUBOUND p%nnodes(3)+nc(3)
#endif
#define VTK_MESH_COMPONENT icomp
#define VTK_MESH_ILBOUND 1
#define VTK_MESH_IUBOUND p%nnodes(1)+nc(1)
#define VTK_MESH_JLBOUND 1
#define VTK_MESH_JUBOUND p%nnodes(2)+nc(2)
#include "vtk/print_data_array.f"
                        ENDDO ! icomp=1,pdat%discr_data%lda
                     ENDIF ! (pdat%discr_data%lda.EQ.?)

                  END SELECT ! (pdat%discr_data%data_type)
                  pdat => p%subpatch_data%next()
               ENDDO ! pdat
               WRITE(iUnit,'(A)') "      </PointData>"
#ifndef VTK_PARALLEL
               WRITE(iUnit,'(A)')  "    </Piece>"
#endif
            END SELECT ! (p => Mesh%subpatch_by_sub(isub)%vec(ipatch)%t)
         ENDDO ! ipatch
      ENDDO ! l=1,topo%nsublist

      DO l=1,topo%nsublist
         isub=topo%isublist(l)
         IF (Mesh%subpatch_by_sub(isub)%nsubpatch.GT.0) THEN
            ! close
#include "vtk/print_end_header1.f"
            ! close file
            CLOSE(iUnit)
            EXIT
         ENDIF
      ENDDO

      end_subroutine()

#if   __DIM == __2D
      END SUBROUTINE ppm_vtk_fields_2d
#elif   __DIM == __3D
      END SUBROUTINE ppm_vtk_fields_3d
#endif

