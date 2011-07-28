      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_rmsh_create_part
      !-------------------------------------------------------------------------
      !
      !  Purpose      : creates particles up from field_up
      !
      !  Input        : xp(:,:)        (F) particle positions
      !                 field_up([[:],:],:,:,:)
      !                                (F) particle values on a field
      !                                    e.g. as generated by rmsh_remesh
      !                 topo_id        (I) user topology id
      !                 mesh_id        (I) user mesh id
      !                 cutoff(2)      (F) Lower (element 1) and upper
      !                                    (element 2) bound of particle
      !                                    strengths. Only particles
      !                                    with strengths in this band
      !                                    will be created.
      !                 resetpos       (L) OPTIONAL: whether to reset the
      !                                    particle positions (default .false.)
      !                 cutoff_weights(:) (F) OPTIONAL, only in the
      !                                    vector case. Gives the
      !                                    weights of the linear
      !                                    combination of the vector
      !                                    elements. The cutoff is then
      !                                    applied to the weighted sum.
      !                 field_wp([[:],:],:,:,:)
      !                                (F) OPTIONAL. A slave field for
      !                                    which to create particles as
      !                                    well where field_up is within
      !                                    cutoff.  
      !                 lda2           (I) lda of the slave field. OPTIONAL
      !                                    has to be present if field_wp is 
      !                                    and if it is a vector field.
      !
      !  Output       : np             (I) new number of particles
      !                 up([:,]:)      (F) new master particle values
      !                 wp([:,]:)      (F) new slave particle values.
      !                                    OPTIONAL. Has to be present 
      !                                    if field_wp is. 
      !                 info           (I) RETURN status
      !
      !  Remarks      : should take square cutoff for comparison
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_rmsh_create_part.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2006/02/03 09:34:04  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.7  2005/06/21 17:36:09  ivos
      !  bugfix in overloading.
      !
      !  Revision 1.6  2005/06/21 01:07:57  ivos
      !  Added the OPTIONAL slave arrays to create multiple particle
      !  properties (maybe by calling this several times). Overloaded for
      !  scalar and vector versions.
      !
      !  Revision 1.5  2004/12/08 19:36:14  ivos
      !  Added upper cutoff. Now particles within a band can be created
      !  (useful for level sets). cutoff is now a 2-array.
      !
      !  Revision 1.4  2004/08/27 08:09:57  michaebe
      !  added optional weighting vector for particle killing criterion
      !
      !  Revision 1.3  2004/08/18 15:03:33  michaebe
      !  bugfix for #0000021 (size vs. ubound)
      !
      !  Revision 1.2  2004/08/12 08:10:56  michaebe
      !  fixed grid spacing bug
      !
      !  Revision 1.1  2004/08/11 08:58:57  michaebe
      !  Initial implementaiton
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_sss_2d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_ssv_2d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_svs_2d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_svv_2d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,lda2)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_sss_3d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_ssv_3d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_svs_3d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_svv_3d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,lda2)
#endif
#endif
#endif

#elif __KIND == __DOUBLE_PRECISION          
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dss_2d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dsv_2d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dvs_2d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dvv_2d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,lda2)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dss_3d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dsv_3d(xp,np,up,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,field_wp,wp,lda2)
#endif
#elif __MODE == __VEC
#if   __MODE2 == __SCA
      SUBROUTINE ppm_rmsh_create_part_dvs_3d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,field_wp,wp)
#elif __MODE2 == __VEC
      SUBROUTINE ppm_rmsh_create_part_dvv_3d(xp,np,up,lda,field_up,&
           & topo_id,mesh_id,cutoff,info,resetpos,cutoff_weights,  &
           & field_wp,wp,lda2)
#endif
#endif
#endif
#endif          

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
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      
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
      INTEGER, INTENT(out)                             :: np
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER        :: up
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#endif
#elif __MODE == __VEC
      INTEGER                         , INTENT(in)     :: lda
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: up
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_up
#endif     
#endif     
      INTEGER                        ,  INTENT(IN   )  :: topo_id, mesh_id
      REAL(mk), DIMENSION(2)         ,  INTENT(in   )  :: cutoff
      INTEGER                        ,  INTENT(  OUT)  :: info
      LOGICAL , OPTIONAL             ,  INTENT(in   )  :: resetpos
#if   __MODE == __VEC
      REAL(mk), OPTIONAL, DIMENSION(:), POINTER        :: cutoff_weights
#endif
#if   __MODE2 == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER, OPTIONAL :: wp
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER, OPTIONAL :: field_wp
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER, OPTIONAL :: field_wp
#endif
#elif __MODE2 == __VEC
      INTEGER                         , INTENT(in)     :: lda2
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: wp
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_wp
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_wp
#endif     
#endif     
      
      !-------------------------------------------------------------------------
      !  Locals
      !-------------------------------------------------------------------------
      INTEGER                                          :: isub, isubl
      INTEGER                                          :: dim
      INTEGER                                          :: meshid, topoid
      INTEGER ,  DIMENSION(:,:),        POINTER        :: istart
      INTEGER ,  DIMENSION(:,:),        POINTER        :: ndata
      LOGICAL                                          :: reset
      INTEGER ,  DIMENSION(2  )                        :: ldu, ldl
      INTEGER                                          :: iopt
      REAL(mk) , DIMENSION(:,:,:),      POINTER        :: min_sub
      REAL(mk) , DIMENSION(:,:  ),      POINTER        :: min_phys, max_phys
      REAL(mk)                                         :: dx,dy,dz
      INTEGER,   DIMENSION(:  ),        POINTER        :: nm
      INTEGER                                          :: inp, nnx, nny, nnz
      INTEGER                                          :: i,j,k
      REAL(mk)                                         :: strength
      LOGICAL                                          :: lok,lslave
#if  __MODE == __VEC
      LOGICAL                                          :: with_weighting
#endif      
      CALL substart('ppm_rmsh_create_part',t0,info)

      dim = ppm_dim
      !-------------------------------------------------------------------------
      !  Check Arguments
      !-------------------------------------------------------------------------
      IF(ppm_debug .GT. 0) THEN
         ! check xp
         IF(ASSOCIATED(xp)) THEN
            IF(SIZE(xp,1).LT.dim) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument, &
                    & 'ppm_rmsh_create_part', &
                    & 'leading dimension of xp insufficient',__LINE__,info)
               GOTO 9999
            END IF
         END IF
         ! check field_up
         IF(.NOT.ASSOCIATED(field_up)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'field_up empty. Run ppm_rmsh_remesh first.',__LINE__,info)
            GOTO 9999
         END IF
#if   __MODE2 == __VEC
         IF(lda2 .LT. 1) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'lda2 must be > 0.',__LINE__,info)
            GOTO 9999
         END IF
#elif __MODE2 == __SCA
         IF(PRESENT(wp).AND..NOT.PRESENT(field_wp)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'field_wp must be present if wp is.',__LINE__,info)
            GOTO 9999
         END IF
         IF(PRESENT(field_wp).AND..NOT.PRESENT(wp)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument, &
                 & 'ppm_rmsh_create_part', &
                 & 'wp must be present if field_wp is.',__LINE__,info)
            GOTO 9999
         END IF
#endif

      END IF

      !-------------------------------------------------------------------------
      !  Check topoid
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_check_topoid(ppm_param_id_user,topo_id,lok,info)
         IF (info .NE. 0 .OR. .NOT.lok) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
                 & 'topo_id is invalid!',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      
      !-------------------------------------------------------------------------
      !  Get the internal topoid
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)

      !-------------------------------------------------------------------------
      !  Check meshid
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_check_meshid(ppm_param_id_user,mesh_id,topoid,lok,info)
         IF (info .NE. 0 .OR. .NOT.lok) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
                 & 'mesh_id is invalid!',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF

      !-------------------------------------------------------------------------
      !  Get the internal meshid
      !-------------------------------------------------------------------------
      meshid = ppm_meshid(topoid)%internal(mesh_id)
      
      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------
      istart => ppm_cart_mesh(meshid,topoid)%istart
      nm     => ppm_cart_mesh(meshid,topoid)%nm
      ndata  => ppm_cart_mesh(meshid,topoid)%nnodes
#if   __KIND == __SINGLE_PRECISION
      min_sub  => ppm_min_subs
      min_phys => ppm_min_physs
      max_phys => ppm_max_physs
#elif __KIND == __DOUBLE_PRECISION
      min_sub => ppm_min_subd
      min_phys => ppm_min_physd
      max_phys => ppm_max_physd
#endif

      dx = (max_phys(1,topoid)-min_phys(1,topoid))/REAL(nm(1)-1,mk)
      dy = (max_phys(2,topoid)-min_phys(2,topoid))/REAL(nm(2)-1,mk)
      IF(dim.EQ.3) dz = (max_phys(3,topoid)-min_phys(3,topoid))/REAL(nm(3)-1,mk)

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
      DO isub=1,ppm_nsublist(topoid)
         isubl = ppm_isublist(isub,topoid)
         !----------------------------------------------------------------------
         !  check the boundary conditions
         !----------------------------------------------------------------------
         IF(ppm_subs_bc(2,isubl,topoid).EQ.1&
              &.AND.ppm_bcdef(2,topoid).NE.ppm_param_bcdef_periodic) THEN
            nnx = ndata(1,isubl)
         ELSE
            nnx = ndata(1,isubl)-1
         END IF

         IF(ppm_subs_bc(4,isubl,topoid).EQ.1&
              &.AND.ppm_bcdef(4,topoid).NE.ppm_param_bcdef_periodic) THEN
            nny = ndata(2,isubl)
         ELSE
            nny = ndata(2,isubl)-1
         END IF

#if   __DIME == __3D
         IF(ppm_subs_bc(6,isubl,topoid).EQ.1&
              &.AND.ppm_bcdef(6,topoid).NE.ppm_param_bcdef_periodic) THEN
            nnz = ndata(3,isubl)
         ELSE
            nnz = ndata(3,isubl)-1
         END IF
#elif __DIME == __2D
         nnz = 1
#endif
         
#if   __DIME == __3D
         DO k=1,nnz
            DO j=1,nny
               DO i=1,nnx
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

      DO j=1,nny
            DO i=1,nnx
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
      
      inp = 0
      !-------------------------------------------------------------------------
      !  loop over subs again and create particles
      !-------------------------------------------------------------------------
      DO isub=1,ppm_nsublist(topoid)
         isubl = ppm_isublist(isub,topoid)
         !----------------------------------------------------------------------
         !  check the boundary conditions
         !----------------------------------------------------------------------
         IF(ppm_subs_bc(2,isubl,topoid).EQ.1&
              &.AND.ppm_bcdef(2,topoid).NE.ppm_param_bcdef_periodic) THEN
            nnx = ndata(1,isubl)
         ELSE
            nnx = ndata(1,isubl)-1
         END IF

         IF(ppm_subs_bc(4,isubl,topoid).EQ.1&
              &.AND.ppm_bcdef(4,topoid).NE.ppm_param_bcdef_periodic) THEN
            nny = ndata(2,isubl)
         ELSE
            nny = ndata(2,isubl)-1
         END IF
#if   __DIME == __3D
         IF(ppm_subs_bc(6,isubl,topoid).EQ.1&
              &.AND.ppm_bcdef(6,topoid).NE.ppm_param_bcdef_periodic) THEN
            nnz = ndata(3,isubl)
         ELSE
            nnz = ndata(3,isubl)-1
         END IF
#elif __DIME == __2D
         nnz = 1
#endif
         
#if   __DIME == __3D
         DO k=1,nnz
            DO j=1,nny
               DO i=1,nnx
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
                     up(inp)       = field_up(i,j,k,isub)
#elif __MODE == __VEC
                     up(1:lda,inp) = field_up(1:lda,i,j,k,isub)
#endif
                     IF (lslave) THEN
#if   __MODE2 == __SCA
                         wp(inp)        = field_wp(i,j,k,isub)
#elif __MODE2 == __VEC
                         wp(1:lda2,inp) = field_wp(1:lda2,i,j,k,isub)
#endif
                     ENDIF
                     IF(reset) THEN
                        xp(1,inp) = REAL(i-1,mk)*dx + min_sub(1,isubl,topoid)
                        xp(2,inp) = REAL(j-1,mk)*dy + min_sub(2,isubl,topoid)
                        xp(3,inp) = REAL(k-1,mk)*dz + min_sub(3,isubl,topoid)
                     END IF
                  END IF
               END DO
            END DO
         END DO
#elif __DIME == __2D
         DO j=1,nny
            DO i=1,nnx
               !----------------------------------------------------------------
               !  compute strength
               !----------------------------------------------------------------
#if   __MODE == __SCA
               strength = ABS(field_up(i,j,isub))
#elif __MODE == __VEC
               IF(with_weighting) THEN
                  strength = SQRT(SUM(&
                    & (field_up(1:lda,i,j,isub)*cutoff_Weights(1:lda))**2 &
                    &  ))
               ELSE
                  strength = SQRT(SUM(field_up(1:lda,i,j,isub)**2))
               END IF
#endif
               IF((strength.GT.cutoff(1)).AND.(strength.LT.cutoff(2))) THEN
                  inp = inp + 1
#if   __MODE == __SCA
                  up(inp)       = field_up(i,j,isub)
#elif __MODE == __VEC
                  up(1:lda,inp) = field_up(1:lda,i,j,isub)
#endif
                  IF (lslave) THEN
#if   __MODE2 == __SCA
                      wp(inp)        = field_wp(i,j,isub)
#elif __MODE2 == __VEC
                      wp(1:lda2,inp) = field_wp(1:lda2,i,j,isub)
#endif
                  ENDIF
                  IF(reset) THEN
                     xp(1,inp) = REAL(i-1,mk)*dx + min_sub(1,isubl,topoid)
                     xp(2,inp) = REAL(j-1,mk)*dy + min_sub(2,isubl,topoid)
                  END IF
               END IF
            END DO
         END DO
#endif
      END DO

      !-------------------------------------------------------------------------
      !  done
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_rmsh_create_part',t0,info)
      RETURN
      
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
