      !-------------------------------------------------------------------------
      !     Subroutine   :                   p2m_interp_bc
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

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_bc_ss_2d(Mesh,Field,field_up,p2m_bcdef,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bc_ds_2d(Mesh,Field,field_up,p2m_bcdef,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_bc_sv_2d(Mesh,Field,field_up,lda,p2m_bcdef,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bc_dv_2d(Mesh,Field,field_up,lda,p2m_bcdef,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_bc_ss_3d(Mesh,Field,field_up,p2m_bcdef,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bc_ds_3d(Mesh,Field,field_up,p2m_bcdef,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_bc_sv_3d(Mesh,Field,field_up,lda,p2m_bcdef,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bc_dv_3d(Mesh,Field,field_up,lda,p2m_bcdef,info)
#endif
#endif
#endif
      !!! Apply boundary conditions in Particle to mesh interpolation
      !!!
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_map
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !------------------------------------------------------------------------!
      !  Arguments
      !------------------------------------------------------------------------!
      CLASS(ppm_t_equi_mesh_)                      :: Mesh
      !!! Mesh
      CLASS(ppm_t_field_)                          :: Field
      !!! Field
#if   __MODE == __SCA
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:    ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:  ) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
#elif __MODE == __VEC
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:  ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
      INTEGER                       , INTENT(IN)   :: lda
      !!! leading dimension of up
#endif
      INTEGER, DIMENSION(:)                        :: p2m_bcdef
      !!! boundary conditions
      !!!    - ppm_param_bcdef_symmetry
      !!!    - ppm_param_bcdef_antisymmetry
      INTEGER                       , INTENT( OUT) :: info
      !!! Returns status, 0 upon success
      !------------------------------------------------------------------------!
      !  Local variables
      !------------------------------------------------------------------------!
      INTEGER :: i,j,k,l
      INTEGER :: xlo,ylo,zlo,xhi,yhi,zhi

      ! aliases
      CLASS(ppm_t_subpatch_), POINTER :: p

      start_subroutine("p2m_interp_bc")

      !  loop over subpatches
      p => Mesh%subpatch%begin()
      subpatch: DO WHILE (ASSOCIATED(p))
         CALL p%get_field(Field,field_up,info)
         or_fail("get_field failed for this subpatch")

#if   __DIME == __2D
#if   __MODE == __SCA
         IF (p%bc(1).GE.0) THEN
            xlo = 1
            ylo = 1
            xhi = 1 + p%ghostsize(1)
            yhi = p%nnodes(2)
            SELECT CASE(p2m_bcdef(1))
            CASE (ppm_param_bcdef_symmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)+field_up(xlo-i+1,j)
                  ENDDO
               ENDDO
            CASE (ppm_param_bcdef_antisymmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)-field_up(xlo-i+1,j)
                  ENDDO
               ENDDO
            END SELECT
         ENDIF
         IF (p%bc(2).GE.0) THEN
            xlo = p%nnodes(1) - p%ghostsize(2)
            ylo = 1
            xhi = p%nnodes(1)
            yhi = p%nnodes(2)
            SELECT CASE(p2m_bcdef(2))
            CASE (ppm_param_bcdef_symmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)+field_up(2*xhi-i,j)
                  ENDDO
               ENDDO
            CASE (ppm_param_bcdef_antisymmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)-field_up(2*xhi-i,j)
                  ENDDO
               ENDDO
            END SELECT
         ENDIF
         IF (p%bc(3).GE.0) THEN
            xlo = 1
            ylo = 1
            xhi = p%nnodes(1)
            yhi = 1 + p%ghostsize(3)
            SELECT CASE(p2m_bcdef(3))
            CASE (ppm_param_bcdef_symmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)+field_up(i,ylo-j+1)
                  ENDDO
               ENDDO
            CASE (ppm_param_bcdef_antisymmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)-field_up(i,ylo-j+1)
                  ENDDO
               ENDDO
            END SELECT
         ENDIF
         IF (p%bc(4).GE.0) THEN
            xlo = 1
            ylo = p%nnodes(2) - p%ghostsize(4)
            xhi = p%nnodes(1)
            yhi = p%nnodes(2)
            SELECT CASE(p2m_bcdef(4))
            CASE (ppm_param_bcdef_symmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)+field_up(i,2*yhi-j)
                  ENDDO
               ENDDO
            CASE (ppm_param_bcdef_antisymmetry)
               DO j=ylo,yhi
                  DO i=xlo,xhi
                     field_up(i,j)=field_up(i,j)-field_up(i,2*yhi-j)
                  ENDDO
               ENDDO
            END SELECT
         ENDIF
#elif __MODE == __VEC
         IF (p%bc(1).GE.0) THEN
            xlo = 1
            ylo = 1
            xhi = 1 + p%ghostsize(1)
            yhi = p%nnodes(2)
            SELECT CASE(p2m_bcdef(1))
            CASE (ppm_param_bcdef_symmetry)
               DO l=1,lda
                  DO j=ylo,yhi
                     DO i=xlo,xhi
                        field_up(l,i,j)=field_up(l,i,j)+field_up(l,xlo-i+1,j)
                     ENDDO
                  ENDDO
               ENDDO
            CASE (ppm_param_bcdef_antisymmetry)
               DO l=1,lda
                  DO j=ylo,yhi
                     DO i=xlo,xhi
                        field_up(l,i,j)=field_up(l,i,j)-field_up(l,xlo-i+1,j)
                     ENDDO
                  ENDDO
               ENDDO
            END SELECT
         ENDIF
         IF (p%bc(2).GE.0) THEN
                  xlo = p%nnodes(1) - p%ghostsize(2)
                  ylo = 1
                  xhi = p%nnodes(1)
                  yhi = p%nnodes(2)
                  SELECT CASE(p2m_bcdef(2))
                  CASE(ppm_param_bcdef_symmetry)
                     DO l=1,lda
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j) = field_up(l,i,j) + &
     &                              field_up(l,2*xhi-i,j)
                        ENDDO
                     ENDDO
                     ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                     DO l=1,lda
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j) = field_up(l,i,j) - &
     &                              field_up(l,2*xhi-i,j)
                        ENDDO
                     ENDDO
                     ENDDO
                  END SELECT
            ENDIF
            IF (p%bc(3).GE.0) THEN
                  xlo = 1
                  ylo = 1
                  xhi = p%nnodes(1)
                  yhi = 1 + p%ghostsize(3)
                  SELECT CASE(p2m_bcdef(3))
                  CASE(ppm_param_bcdef_symmetry)
                     DO l=1,lda
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j) = field_up(l,i,j) + &
     &                              field_up(l,i,ylo-j+1)
                        ENDDO
                     ENDDO
                     ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                     DO l=1,lda
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j) = field_up(l,i,j) - &
     &                              field_up(l,i,ylo-j+1)
                        ENDDO
                     ENDDO
                     ENDDO
                  END SELECT
            ENDIF
            IF (p%bc(4).GE.0) THEN
                  xlo = 1
                  ylo = p%nnodes(2) - p%ghostsize(4)
                  xhi = p%nnodes(1)
                  yhi = p%nnodes(2)
                  SELECT CASE(p2m_bcdef(4))
                  CASE(ppm_param_bcdef_symmetry)
                     DO l=1,lda
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j) = field_up(l,i,j) + &
     &                              field_up(l,i,2*yhi-j)
                        ENDDO
                     ENDDO
                     ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                     DO l=1,lda
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j) = field_up(l,i,j) - &
     &                              field_up(l,i,2*yhi-j)
                        ENDDO
                     ENDDO
                     ENDDO
                 END SELECT
            ENDIF
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
            IF (p%bc(1).GE.0) THEN
               xlo = 1
               ylo = 1
               zlo = 1
               xhi = 1 + p%ghostsize(1)
               yhi = p%nnodes(2)
               zhi = p%nnodes(3)
               SELECT CASE(p2m_bcdef(1))
               CASE(ppm_param_bcdef_symmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) + &
     &                              field_up(xlo-i+1,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
               CASE(ppm_param_bcdef_antisymmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) - &
     &                              field_up(xlo-i+1,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
               END SELECT
            ENDIF
            IF (p%bc(2).GE.0) THEN
               xlo = p%nnodes(1) - p%ghostsize(2)
               ylo = 1
               zlo = 1
               xhi = p%nnodes(1)
               yhi = p%nnodes(2)
               zhi = p%nnodes(3)
               SELECT CASE(p2m_bcdef(2))
               CASE(ppm_param_bcdef_symmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) + &
     &                              field_up(2*xhi-i,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
               CASE(ppm_param_bcdef_antisymmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) - &
     &                              field_up(2*xhi-i,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
               END SELECT
            ENDIF
            IF (p%bc(3).GE.0) THEN
               xlo = 1
               ylo = 1
               zlo = 1
               xhi = p%nnodes(1)
               yhi = 1 + p%ghostsize(3)
               zhi = p%nnodes(3)
               SELECT CASE(p2m_bcdef(3))
               CASE(ppm_param_bcdef_symmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) + &
     &                              field_up(i,ylo-j+1,k)
                        ENDDO
                     ENDDO
                  ENDDO
               CASE(ppm_param_bcdef_antisymmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) - &
     &                              field_up(i,ylo-j+1,k)
                        ENDDO
                     ENDDO
                  ENDDO
               END SELECT
            ENDIF
            IF (p%bc(4).GE.0) THEN
               xlo = 1
               ylo = p%nnodes(2) - p%ghostsize(4)
               zlo = 1
               xhi = p%nnodes(1)
               yhi = p%nnodes(2)
               zhi = p%nnodes(3)
               SELECT CASE(p2m_bcdef(4))
               CASE(ppm_param_bcdef_symmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) + &
     &                              field_up(i,2*yhi-j,k)
                        ENDDO
                     ENDDO
                  ENDDO
               CASE(ppm_param_bcdef_antisymmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) - &
     &                              field_up(i,2*yhi-j,k)
                        ENDDO
                     ENDDO
                  ENDDO
               END SELECT
            ENDIF
            IF (p%bc(5).GE.0) THEN
               xlo = 1
               ylo = 1
               zlo = 1
               xhi = p%nnodes(1)
               yhi = p%nnodes(2)
               zhi = 1 + p%ghostsize(5)
               SELECT CASE(p2m_bcdef(5))
               CASE(ppm_param_bcdef_symmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) + &
     &                              field_up(i,j,zlo-k+1)
                        ENDDO
                     ENDDO
                  ENDDO
               CASE(ppm_param_bcdef_antisymmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) - &
     &                              field_up(i,j,zlo-k+1)
                        ENDDO
                     ENDDO
                  ENDDO
               END SELECT
            ENDIF
            IF (p%bc(6).GE.0) THEN
               xlo = 1
               ylo = 1
               zlo = p%nnodes(3) - p%ghostsize(6)
               xhi = p%nnodes(1)
               yhi = p%nnodes(2)
               zhi = p%nnodes(3)
               SELECT CASE(p2m_bcdef(6))
               CASE(ppm_param_bcdef_symmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) + &
     &                              field_up(i,j,2*zhi-k)
                        ENDDO
                     ENDDO
                  ENDDO
               CASE(ppm_param_bcdef_antisymmetry)
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(i,j,k) = field_up(i,j,k) - &
     &                              field_up(i,j,2*zhi-k)
                        ENDDO
                     ENDDO
                  ENDDO
               END SELECT
            ENDIF
#elif __MODE == __VEC
            IF (p%bc(1).GE.0) THEN
                  xlo = 1
                  ylo = 1
                  zlo = 1
                  xhi = 1 + p%ghostsize(1)
                  yhi = p%nnodes(2)
                  zhi = p%nnodes(3)
                  SELECT CASE(p2m_bcdef(1))
                  CASE(ppm_param_bcdef_symmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) + &
     &                              field_up(l,xlo-i+1,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) - &
     &                              field_up(l,xlo-i+1,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  END SELECT
            ENDIF
            IF (p%bc(2).GE.0) THEN
                  xlo = p%nnodes(1) - p%ghostsize(2)
                  ylo = 1
                  zlo = 1
                  xhi = p%nnodes(1)
                  yhi = p%nnodes(2)
                  zhi = p%nnodes(3)
                  SELECT CASE(p2m_bcdef(2))
                  CASE(ppm_param_bcdef_symmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) + &
     &                              field_up(l,2*xhi-i,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) - &
     &                              field_up(l,2*xhi-i,j,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  END SELECT
            ENDIF
            IF (p%bc(3).GE.0) THEN
                  xlo = 1
                  ylo = 1
                  zlo = 1
                  xhi = p%nnodes(1)
                  yhi = 1 + p%ghostsize(3)
                  zhi = p%nnodes(3)
                  SELECT CASE(p2m_bcdef(3))
                  CASE(ppm_param_bcdef_symmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) + &
     &                              field_up(l,i,ylo-j+1,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) - &
     &                              field_up(l,i,ylo-j+1,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  END SELECT
            ENDIF
            IF (p%bc(4).GE.0) THEN
                  xlo = 1
                  ylo = p%nnodes(2) - p%ghostsize(4)
                  zlo = 1
                  xhi = p%nnodes(1)
                  yhi = p%nnodes(2)
                  zhi = p%nnodes(3)
                  SELECT CASE(p2m_bcdef(4))
                  CASE(ppm_param_bcdef_symmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) + &
     &                              field_up(l,i,2*yhi-j,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) - &
     &                              field_up(l,i,2*yhi-j,k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                 END SELECT
            ENDIF
            IF (p%bc(5).GE.0) THEN
                  xlo = 1
                  ylo = 1
                  zlo = 1
                  xhi = p%nnodes(1)
                  yhi = p%nnodes(2)
                  zhi = 1 + p%ghostsize(5)
                  SELECT CASE(p2m_bcdef(5))
                  CASE(ppm_param_bcdef_symmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) + &
     &                              field_up(l,i,j,zlo-k+1)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) - &
     &                              field_up(l,i,j,zlo-k+1)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  END SELECT
            ENDIF
            IF (p%bc(6).GE.0) THEN
                  xlo = 1
                  ylo = 1
                  zlo = p%nnodes(3) - p%ghostsize(6)
                  xhi = p%nnodes(1)
                  yhi = p%nnodes(2)
                  zhi = p%nnodes(3)
                  SELECT CASE(p2m_bcdef(6))
                  CASE(ppm_param_bcdef_symmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) + &
     &                              field_up(l,i,j,2*zhi-k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  CASE(ppm_param_bcdef_antisymmetry)
                  DO l=1,lda
                  DO k= zlo,zhi
                     DO j=ylo,yhi
                        DO i=xlo,xhi
                           field_up(l,i,j,k) = field_up(l,i,j,k) - &
     &                      field_up(l,i,j,2*zhi-k)
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  END SELECT
            ENDIF
#endif
#endif

          p => Mesh%subpatch%next()
        ENDDO subpatch       ! loop over subpatches

      end_subroutine()
      RETURN

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bc_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bc_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bc_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bc_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bc_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bc_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bc_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bc_dv_3d
#endif
#endif
#endif
