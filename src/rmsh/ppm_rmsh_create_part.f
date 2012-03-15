      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_rmsh_create_part
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
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_sss_2d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_ssv_2d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_svs_2d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_svv_2d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,vp,lda2)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_sss_3d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_ssv_3d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_svs_3d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_svv_3d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,vp,lda2)
#endif
#endif
#endif
#elif __KIND == __DOUBLE_PRECISION
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dss_2d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dsv_2d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dvs_2d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dvv_2d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,vp,lda2)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dss_3d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dsv_3d(topoid,meshid,xp,Np,up,&
           & field_up,cutoff,info,resetpos,field_wp,wp,vp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dvs_3d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,field_wp,wp,vp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dvv_3d(topoid,meshid,xp,Np,up,lda,&
           & field_up,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,vp,lda2)
#endif
#endif
#endif
#endif
      !!! Creates particles up from field_up
      !!!
      !!! [NOTE]
      !!! Should take square cutoff for comparison
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_data_rmsh
      USE ppm_module_data_mesh
      USE ppm_module_data
      USE ppm_module_check_id

      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(mk),  DIMENSION(:,:)       , POINTER        :: xp
      !!! Particle positions
      INTEGER, INTENT(out)                             :: np
      !!! New number of particles
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER        :: up
      !!! New master particle values
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#endif
      !!! Particle values on a field e.g. as generated by rmsh_remesh
#elif __MODE == __VEC
      INTEGER                         , INTENT(in)     :: lda
      !!! Leading dimension
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: up
      !!! New master particle values
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_up
#endif
      !!! Particle values on a field e.g. as generated by rmsh_remesh
#endif
      INTEGER                        ,  INTENT(IN   )  :: topoid
      !!! Topology ID
      INTEGER                        ,  INTENT(IN   )  :: meshid
      !!! Mesh ID
      REAL(mk), DIMENSION(2)         ,  INTENT(in   )  :: cutoff
      !!! Lower (element 1) and upper (element 2) bound of particle
      !!! strengths. Only particles with strengths in this band
      !!! will be created.
      INTEGER                        ,  INTENT(  OUT)  :: info
      !!! Returns 0 upon success
      LOGICAL , OPTIONAL             ,  INTENT(in   )  :: resetpos
      !!! Reset the particle positions? (default .false.)
#if   __MODE == __VEC
      REAL(mk), OPTIONAL, DIMENSION(:), POINTER        :: cutoff_weights
      !!! Gives the weights of the linear
      !!! combination of the vector elements. The cutoff is then applied
      !!! to the weighted sum.
#endif
#if   __MODE2 == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER, OPTIONAL :: wp
      !!! New slave particle values. Has to be present if field_wp is.
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER, OPTIONAL :: field_wp
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER, OPTIONAL :: field_wp
#endif
      !!! A slave field for which to create particles as well where field_up 
      !!! is within cutoff. 
#elif __MODE2 == __VEC
      INTEGER                         , INTENT(in)     :: lda2
      !!! Has to be present if field_wp is and if it is a vector field.
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: wp
      !!! New slave particle values. Has to be present if field_wp is.
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_wp
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_wp
#endif
      !!! A slave field for which to create particles as well
      !!! where field_up is within cutoff.
#endif
      REAL(mk),  DIMENSION(:)         , POINTER, OPTIONAL :: vp
      !!! Particle volumes
      !-------------------------------------------------------------------------
      !  Locals
      !-------------------------------------------------------------------------
      INTEGER                                          :: isub, isubl
      INTEGER                                          :: dim
      INTEGER ,  DIMENSION(:,:),        POINTER        :: istart => NULL()
      INTEGER ,  DIMENSION(:,:),        POINTER        :: ndata  => NULL()
      LOGICAL                                          :: reset
      INTEGER ,  DIMENSION(2  )                        :: ldu, ldl
      INTEGER                                          :: iopt
      REAL(mk) , DIMENSION(:,:),        POINTER        :: min_sub => NULL()
      REAL(mk) , DIMENSION(:),          POINTER        :: min_phys => NULL()
      REAL(mk) , DIMENSION(:),          POINTER        :: max_phys => NULL()
      REAL(mk)                                         :: dx,dy,dz
      INTEGER,   DIMENSION(:  ),        POINTER        :: nm => NULL()
      INTEGER                                          :: inp, nnx, nny, nnz
      INTEGER                                          :: startx, starty, startz
      INTEGER                                          :: i,j,k
      INTEGER                                          :: ilowbc, ihighbc
      INTEGER                                          :: jlowbc, jhighbc
      INTEGER                                          :: klowbc, khighbc
      REAL(mk)                                         :: xbc_factor, ybc_factor
      REAL(mk)                                         :: zbc_factor
      REAL(mk)                                         :: strength
      LOGICAL                                          :: lok,lslave
      LOGICAL                                          :: with_vol
      TYPE(ppm_t_equi_mesh), POINTER                   :: p_mesh => NULL()
      TYPE(ppm_t_topo)     , POINTER                   :: topo   => NULL()
#if  __MODE == __VEC
      LOGICAL                                          :: with_weighting
#endif
      CALL substart('ppm_rmsh_create_part',t0,info)

      dim = ppm_dim
      !-------------------------------------------------------------------------
      !  Check Arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check topoid
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_check_topoid(topoid,lok,info)
         IF (info .NE. 0 .OR. .NOT.lok) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_create_part',  &
                 & 'topo_id is invalid!',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF

      !-------------------------------------------------------------------------
      !  Check meshid
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_check_meshid(topoid,meshid,lok,info)
         IF (info .NE. 0 .OR. .NOT.lok) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_create_part',  &
                 & 'mesh_id is invalid!',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF

      !-------------------------------------------------------------------------
      !  Get meshid
      !-------------------------------------------------------------------------
      !meshid = ppm_topo(topo_id)%t%mesh(mesh_id)

      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------
      topo   => ppm_topo(topoid)%t
      SELECT TYPE (t => ppm_mesh%vec(meshid))
      TYPE IS (ppm_t_equi_mesh)
          p_mesh => t
      END SELECT
      istart => p_mesh%istart
      nm     => p_mesh%nm
      ndata  => p_mesh%nnodes
#if   __KIND == __SINGLE_PRECISION
      min_sub  => topo%min_subs
      min_phys => topo%min_physs
      max_phys => topo%max_physs
#elif __KIND == __DOUBLE_PRECISION
      min_sub  => topo%min_subd
      min_phys => topo%min_physd
      max_phys => topo%max_physd
#endif

      dx = (max_phys(1)-min_phys(1))/REAL(nm(1)-1,mk)
      dy = (max_phys(2)-min_phys(2))/REAL(nm(2)-1,mk)
      IF(dim.EQ.3) dz = (max_phys(3)-min_phys(3))/REAL(nm(3)-1,mk)

#if  __MODE == __VEC
      !-------------------------------------------------------------------------
      !  check the cutoff weights if present
      !-------------------------------------------------------------------------
      with_weighting = PRESENT(cutoff_weights)
      IF(ppm_debug.GT.0) THEN
         IF(with_weighting) THEN
            IF(.NOT.ASSOCIATED(cutoff_weights)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_create_part', &
                    & 'cutoff_weights is not allocated',__LINE__,info)
               GOTO 9999
            ELSE
               DO i=1,lda
                  IF(cutoff_weights(i).LT.0.0_mk&
                       &.OR.cutoff_weights(i).GT.1.0_mk) THEN
                     info = ppm_error_error
                     CALL ppm_error(ppm_err_argument,'ppm_rmsh_create_part',&
                          & '0 <= cutoff_weights <= 1 is violated.'&
                          &,__LINE__,info)
                     GOTO 9999
                  END IF
               END DO
            END IF
         END IF
      END IF
#endif

      !-------------------------------------------------------------------------
      !  check whether to reset the positions
      !-------------------------------------------------------------------------
      IF(PRESENT(resetpos)) THEN
         reset = resetpos
      ELSE
         reset = .FALSE.
      END IF
      
      !-------------------------------------------------------------------------
      !  check whether to reset the positions
      !-------------------------------------------------------------------------
      with_vol = PRESENT(vp)
      
      !-------------------------------------------------------------------------
      !  Check whether we have slave fields
      !-------------------------------------------------------------------------
#if   __MODE2 == __SCA
      lslave = .FALSE.
      IF (PRESENT(field_wp)) lslave = .TRUE.
#elif __MODE2 == __VEC
      lslave = .TRUE.
#endif

      !-------------------------------------------------------------------------
      !  Loop over the subs and compute the number of particles that will be
      !  created
      !-------------------------------------------------------------------------
      inp = 0
      DO isub=1,topo%nsublist
         isubl = topo%isublist(isub)
         !----------------------------------------------------------------------
         !  check the boundary conditions
         !----------------------------------------------------------------------
         startx = 1
         IF(topo%subs_bc(1,isubl).EQ.1&
              &.AND.topo%bcdef(1).EQ.ppm_param_bcdef_freespace) THEN
            startx = 2
         END IF
         IF(topo%subs_bc(2,isub).EQ.1&
              &.AND.topo%bcdef(2).NE.ppm_param_bcdef_periodic &
              &.AND.topo%bcdef(2).NE.ppm_param_bcdef_freespace) THEN
            nnx = ndata(1,isubl)
         ELSE
            nnx = ndata(1,isubl)-1
         END IF

         starty = 1
         IF(topo%subs_bc(3,isubl).EQ.1&
              &.AND.topo%bcdef(3).EQ.ppm_param_bcdef_freespace) THEN
            starty = 2
         END IF
         IF(topo%subs_bc(4,isubl).EQ.1&
              &.AND.topo%bcdef(4).NE.ppm_param_bcdef_periodic &
              &.AND.topo%bcdef(4).NE.ppm_param_bcdef_freespace) THEN
            nny = ndata(2,isubl)
         ELSE
            nny = ndata(2,isubl)-1
         END IF

#if   __DIME == __3D
         startz = 1
         IF(topo%subs_bc(5,isubl).EQ.1&
              &.AND.topo%bcdef(5).EQ.ppm_param_bcdef_freespace) THEN
            startz = 2
         END IF
         IF(topo%subs_bc(6,isubl).EQ.1&
              &.AND.topo%bcdef(6).NE.ppm_param_bcdef_periodic &
              &.AND.topo%bcdef(6).NE.ppm_param_bcdef_freespace) THEN
            nnz = ndata(3,isubl)
         ELSE
            nnz = ndata(3,isubl)-1
         END IF
#elif __DIME == __2D
         startz = 1
         nnz = 1
#endif

#if   __DIME == __3D
         DO k=startz,nnz
            DO j=starty,nny
               DO i=startx,nnx
                  !-------------------------------------------------------------
                  !  compute strength
                  !-------------------------------------------------------------
#if   __MODE == __SCA
                  strength = ABS(field_up(i,j,k,isub))
#elif __MODE == __VEC
                  IF(with_weighting) THEN
                     strength = SQRT(SUM(&
                     & (field_up(1:lda,i,j,k,isub)*cutoff_weights(1:lda))**2 &
                     &  ))
                  ELSE
                     strength = SQRT(SUM(field_up(1:lda,i,j,k,isub)**2))
                  END IF
#endif
                  IF((strength.GT.cutoff(1)).AND.(strength.LT.cutoff(2))) THEN
                     inp = inp + 1
                  END IF
               END DO
            END DO
         END DO
#elif __DIME == __2D

      DO j=starty,nny
            DO i=startx,nnx
               !----------------------------------------------------------------
               !  compute strength
               !----------------------------------------------------------------
#if   __MODE == __SCA
               strength = ABS(field_up(i,j,isub))
#elif __MODE == __VEC
               IF(with_weighting) THEN
                  strength = SQRT(SUM(&
                       & (field_up(1:lda,i,j,isub)*cutoff_weights(1:lda))**2 &
                       &  ))
               ELSE
                  strength = SQRT(SUM(field_up(1:lda,i,j,isub)**2))
               END IF
#endif
               IF((strength.GT.cutoff(1)).AND.(strength.LT.cutoff(2))) THEN
                  inp = inp + 1
               END IF
            END DO
         END DO
#endif
      END DO
      np = inp

      !-------------------------------------------------------------------------
      !  grow the up array such as to accomodate inp elements
      !-------------------------------------------------------------------------
#if   __MODE == __SCA
      ldl(1) = 1
      ldu(1) = np
#elif __MODE == __VEC
      ldl(1) = 1
      ldu(1) = lda
      ldl(2) = 1
      ldu(2) = np
#endif
      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(up,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_rmsh_create_part', &
              & 'could not allocate new particles (up)',__LINE__,info)
         GOTO 9999
      END IF

      IF (lslave) THEN
#if   __MODE2 == __SCA
          ldl(1) = 1
          ldu(1) = np
#elif __MODE2 == __VEC
          ldl(1) = 1
          ldu(1) = lda2
          ldl(2) = 1
          ldu(2) = np
#endif
          iopt = ppm_param_alloc_fit
          CALL ppm_alloc(wp,ldu,iopt,info)
          IF(info.NE.0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_rmsh_create_part', &
                  & 'could not allocate new particles (wp)',__LINE__,info)
             GOTO 9999
          END IF
      END IF

      !-------------------------------------------------------------------------
      !  grow the xp array if necessary
      !-------------------------------------------------------------------------
      IF(reset) THEN
         ldu(1) = dim
         ldu(2) = np

         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(xp,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_rmsh_create_part', &
                 & 'could not allocate new particles (xp)',__LINE__,info)
            GOTO 9999
         END IF
      END IF
      
      !-------------------------------------------------------------------------
      !  grow the vp array if necessary
      !-------------------------------------------------------------------------
      IF(with_vol) THEN
         ldu(1) = np

         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(vp,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_rmsh_create_part', &
                 & 'could not allocate new volumes (vp)',__LINE__,info)
            GOTO 9999
         END IF
      END IF
      
      inp = 0
      !-------------------------------------------------------------------------
      !  loop over subs again and create particles
      !-------------------------------------------------------------------------
      DO isub=1,topo%nsublist
         isubl = topo%isublist(isub)
         !----------------------------------------------------------------------
         !  check the boundary conditions
         !----------------------------------------------------------------------
         startx = 1
         IF(topo%subs_bc(1,isubl).EQ.1&
              &.AND.topo%bcdef(1).EQ.ppm_param_bcdef_freespace) THEN
            startx = 2
         END IF
         IF(topo%subs_bc(2,isubl).EQ.1&
              &.AND.topo%bcdef(2).NE.ppm_param_bcdef_periodic &
              &.AND.topo%bcdef(2).NE.ppm_param_bcdef_freespace) THEN
            nnx = ndata(1,isubl)
         ELSE
            nnx = ndata(1,isubl)-1
         END IF

         starty = 1
         IF(topo%subs_bc(3,isubl).EQ.1&
              &.AND.topo%bcdef(3).EQ.ppm_param_bcdef_freespace) THEN
            starty = 2
         END IF
         IF(topo%subs_bc(4,isubl).EQ.1&
              &.AND.topo%bcdef(4).NE.ppm_param_bcdef_periodic &
              &.AND.topo%bcdef(4).NE.ppm_param_bcdef_freespace) THEN
            nny = ndata(2,isubl)
         ELSE
            nny = ndata(2,isubl)-1
         END IF
#if   __DIME == __3D

         startz = 1
         IF(topo%subs_bc(5,isubl).EQ.1&
              &.AND.topo%bcdef(5).EQ.ppm_param_bcdef_freespace) THEN
            startz = 2
         END IF
         IF(topo%subs_bc(6,isubl).EQ.1&
              &.AND.topo%bcdef(6).NE.ppm_param_bcdef_periodic &
              &.AND.topo%bcdef(6).NE.ppm_param_bcdef_freespace) THEN
            nnz = ndata(3,isubl)
         ELSE
            nnz = ndata(3,isubl)-1
         END IF
#elif __DIME == __2D
         startz = 1
         nnz = 1
#endif
         
#if   __DIME == __3D
         ilowbc = 0
         IF(topo%subs_bc(1,isubl).EQ.1&
              &.AND.((topo%bcdef(1).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(1).NE.ppm_param_bcdef_freespace))) THEN
            ilowbc = 1 
         END IF
         ihighbc = ndata(1,isubl)+1
         IF(topo%subs_bc(2,isubl).EQ.1&
              &.AND.((topo%bcdef(2).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(2).NE.ppm_param_bcdef_freespace))) THEN
            ihighbc = ndata(1,isubl)
         END IF
         jlowbc = 0
         IF(topo%subs_bc(3,isubl).EQ.1&
              &.AND.((topo%bcdef(3).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(3).NE.ppm_param_bcdef_freespace))) THEN
            jlowbc = 1 
         END IF
         jhighbc = ndata(2,isubl)+1
         IF(topo%subs_bc(4,isubl).EQ.1&
              &.AND.((topo%bcdef(4).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(4).NE.ppm_param_bcdef_freespace))) THEN
            jhighbc = ndata(2,isubl)
         END IF
         klowbc = 0
         IF(topo%subs_bc(5,isubl).EQ.1&
              &.AND.((topo%bcdef(5).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(5).NE.ppm_param_bcdef_freespace))) THEN
            klowbc = 1 
         END IF
         khighbc = ndata(3,isubl)+1
         IF(topo%subs_bc(6,isubl).EQ.1&
              &.AND.((topo%bcdef(6).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(6).NE.ppm_param_bcdef_freespace))) THEN
            khighbc = ndata(3,isubl)
         END IF

         DO k=startz,nnz
            IF (k.EQ.klowbc .OR. k.EQ.khighbc) THEN
               zbc_factor = 0.5_mk
            ELSE
               zbc_factor = 1.0_mk
            END IF
            DO j=starty,nny
               IF (j.EQ.jlowbc .OR. j.EQ.jhighbc) THEN
                  ybc_factor = 0.5_mk
               ELSE
                  ybc_factor = 1.0_mk
               END IF
               DO i=startx,nnx
                  IF (i.EQ.ilowbc .OR. i.EQ.ihighbc) THEN
                     xbc_factor = 0.5_mk
                  ELSE
                     xbc_factor = 1.0_mk
                  END IF
                  !-------------------------------------------------------------
                  !  compute strength
                  !-------------------------------------------------------------
#if   __MODE == __SCA
                  strength = ABS(field_up(i,j,k,isub))
#elif __MODE == __VEC
                  IF(with_weighting) THEN
                     strength = SQRT(SUM(&
                     & (field_up(1:lda,i,j,k,isub)*cutoff_Weights(1:lda))**2 &
                     &  ))
                  ELSE
                     strength = SQRT(SUM(field_up(1:lda,i,j,k,isub)**2))
                  END IF
#endif
                  IF((strength.GT.cutoff(1)).AND.(strength.LT.cutoff(2))) THEN
                     inp = inp + 1
#if   __MODE == __SCA
                     up(inp) = xbc_factor*ybc_factor*zbc_factor* &
     &                         field_up(i,j,k,isub)
#elif __MODE == __VEC
                     up(1:lda,inp) = xbc_factor*ybc_factor*zbc_factor* &
     &                               field_up(1:lda,i,j,k,isub)
#endif
                     IF (lslave) THEN
#if   __MODE2 == __SCA
                         wp(inp) = xbc_factor*ybc_factor*zbc_factor* &
     &                             field_wp(i,j,k,isub)
#elif __MODE2 == __VEC
                         wp(1:lda2,inp) = xbc_factor*ybc_factor*zbc_factor* &
     &                                    field_wp(1:lda2,i,j,k,isub)
#endif
                     ENDIF
                     IF(reset) THEN
                        xp(1,inp) = REAL(i-1,mk)*dx + min_sub(1,isubl)
                        xp(2,inp) = REAL(j-1,mk)*dy + min_sub(2,isubl)
                        xp(3,inp) = REAL(k-1,mk)*dz + min_sub(3,isubl)
                     END IF
                  END IF
               END DO
            END DO
         END DO
#elif __DIME == __2D
         ilowbc = 0
         IF(topo%subs_bc(1,isubl).EQ.1&
              &.AND.((topo%bcdef(1).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(1).NE.ppm_param_bcdef_freespace))) THEN
            ilowbc = 1 
         END IF
         ihighbc = ndata(1,isubl)+1
         IF(topo%subs_bc(2,isubl).EQ.1&
              &.AND.((topo%bcdef(2).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(2).NE.ppm_param_bcdef_freespace))) THEN
            ihighbc = ndata(1,isubl)
         END IF
         jlowbc = 0
         IF(topo%subs_bc(3,isubl).EQ.1&
              &.AND.((topo%bcdef(3).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(3).NE.ppm_param_bcdef_freespace))) THEN
            jlowbc = 1 
         END IF
         jhighbc = ndata(2,isubl)+1
         IF(topo%subs_bc(4,isubl).EQ.1&
              &.AND.((topo%bcdef(4).NE.ppm_param_bcdef_periodic)&
              &  .AND.(topo%bcdef(4).NE.ppm_param_bcdef_freespace))) THEN
            jhighbc = ndata(2,isubl)
         END IF
         DO j=starty,nny
            IF (j.EQ.jlowbc .OR. j.EQ.jhighbc) THEN
               ybc_factor = 0.5_mk
            ELSE
               ybc_factor = 1.0_mk
            END IF
            DO i=startx,nnx
               IF (i.EQ.ilowbc .OR. i.EQ.ihighbc) THEN
                  xbc_factor = 0.5_mk
               ELSE
                  xbc_factor = 1.0_mk
               END IF
               !----------------------------------------------------------------
               !  compute strength
               !----------------------------------------------------------------
#if   __MODE == __SCA
               strength = ABS(field_up(i,j,isub))
#elif __MODE == __VEC
               IF(with_weighting) THEN
                  strength = SQRT(SUM((field_up(1:lda,i,j,isub)* &
     &                       cutoff_Weights(1:lda))**2))
               ELSE
                  strength = SQRT(SUM(field_up(1:lda,i,j,isub)**2))
               END IF
#endif
               IF((strength.GT.cutoff(1)).AND.(strength.LT.cutoff(2))) THEN
                  inp = inp + 1
#if   __MODE == __SCA
                  up(inp) = xbc_factor*ybc_factor*field_up(i,j,isub)
#elif __MODE == __VEC
                  up(1:lda,inp) = xbc_factor*ybc_factor* &
     &                            field_up(1:lda,i,j,isub)
#endif
                  IF (lslave) THEN
#if   __MODE2 == __SCA
                      wp(inp) = xbc_factor*ybc_factor*field_wp(i,j,isub)
#elif __MODE2 == __VEC
                      wp(1:lda2,inp) = xbc_factor*ybc_factor* &
     &                                 field_wp(1:lda2,i,j,isub)
#endif
                  ENDIF
                  IF(reset) THEN
                     xp(1,inp) = REAL(i-1,mk)*dx + min_sub(1,isubl)
                     xp(2,inp) = REAL(j-1,mk)*dy + min_sub(2,isubl)
                  END IF
                  IF(with_vol) THEN
                     vp(inp) = xbc_factor*ybc_factor
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
#endif
      END DO

      !-------------------------------------------------------------------------
      !  done
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_rmsh_create_part',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        ! check xp
         IF(ASSOCIATED(xp)) THEN
            IF(SIZE(xp,1).LT.dim) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument, &
                    & 'ppm_rmsh_create_part', &
                    & 'leading dimension of xp insufficient',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDIF
         ! check field_up
         IF(.NOT.ASSOCIATED(field_up)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'field_up empty. Run ppm_rmsh_remesh first.',__LINE__,info)
            GOTO 8888
         ENDIF
#if   __MODE2 == __VEC
         IF(lda2 .LT. 1) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'lda2 must be > 0.',__LINE__,info)
            GOTO 8888
         ENDIF
#elif __MODE2 == __SCA
         IF(PRESENT(wp).AND..NOT.PRESENT(field_wp)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'field_wp must be present if wp is.',__LINE__,info)
            GOTO 8888
         ENDIF
         IF(PRESENT(field_wp).AND..NOT.PRESENT(wp)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'wp must be present if field_wp is.',__LINE__,info)
            GOTO 8888
         ENDIF
#endif
 8888    CONTINUE
      END SUBROUTINE check

#if   __KIND == __SINGLE_PRECISION
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_sss_2d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_ssv_2d
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_svs_2d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_svv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_sss_3d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_ssv_3d
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_svs_3d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_svv_3d
#endif
#endif
#endif
#elif __KIND == __DOUBLE_PRECISION
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_dss_2d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_dsv_2d
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_dvs_2d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_dvv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_dss_3d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_dsv_3d
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      END SUBROUTINE ppm_rmsh_create_part_dvs_3d
#elif __MODE2 == __VEC
      END SUBROUTINE ppm_rmsh_create_part_dvv_3d
#endif
#endif
#endif
#endif
