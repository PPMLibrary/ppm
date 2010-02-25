      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_rmsh_remesh
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine recompute the new particles values, thx to
      !                 the precomputed weights...
      !                 
      !
      !  Input        : xp(:,:)      (F) the position of the particles
      !                 Np           (I) the number of particles (on processor)
      !                 topo_id      (I) topology identifier (user)
      !                                  of target
      !                 mesh_id      (I) Id of the mesh (user)
      !                 up(:[,:])    (F) The array to be remeshed
      !                 kernel       (I) Choice of the kernel used to
      !                                  compute the weights
      !                 ghostsize(:) (I) Number of ghost particles
      !                 wx1,2,3_user(I) OPTIONAL :
      !                                 If present, these weights are used,
      !                                 instead of internal weights
      !
      !
      !  Output       : info         (I) return status
      !                 field_up(:,:,:,:) (I) mesh on which particles have been 
      !                                  distributed (if not present field)
      !
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    : 
      !-------------------------------------------------------------------------
      !  $Log: ppm_rmsh_remesh.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.12  2005/03/10 01:51:13  ivos
      !  Resolved CVS conflict.
      !
      !  Revision 1.11  2004/09/04 15:38:44  michaebe
      !  unrolled remeshing for lda=4
      !
      !  Revision 1.10  2004/08/27 09:23:08  michaebe
      !  unrolled loop for lda=3 dim=3 krnl=mp4; makes it more than 2x 
      !  as fast @G5
      !
      !  Revision 1.9  2004/08/24 07:34:10  michaebe
      !  removed quotes from comments
      !
      !  Revision 1.8  2004/08/18 15:03:33  michaebe
      !  bugfix for #0000021 (size vs. ubound)
      !
      !  Revision 1.7  2004/08/16 12:09:42  michaebe
      !  included tup as tmp variable
      !
      !  Revision 1.6  2004/08/16 11:46:36  michaebe
      !  removed a write statemetn
      !
      !  Revision 1.5  2004/08/16 11:04:19  michaebe
      !  included handling of np=0 cases
      !
      !  Revision 1.4  2004/08/13 12:03:03  michaebe
      !  removed epsilon from left neighbor calc
      !
      !  Revision 1.3  2004/08/12 08:11:14  michaebe
      !  fixed grid spacing bug
      !
      !  Revision 1.2  2004/08/10 16:11:20  michaebe
      !  Countless and major changes..
      !
      !  Revision 1.1  2004/08/09 12:00:00  michaebe
      !  revised initinal implementation (P. Gautier).
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
 
      !-------------------------------------------------------------------------
      !  Have to modify field_up to be optional. if its not present then an
      !  internal field is allocated; Nfortunately a little bit of a pain cause
      !  the overloading is not going to get a grip on whats going on then.
      !-------------------------------------------------------------------------

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ss_2d(xp,Np,up,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ds_2d(xp,Np,up,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_sv_2d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_dv_2d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ss_3d(xp,Np,up,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ds_3d(xp,Np,up,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_sv_3d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_dv_3d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#endif
#endif       

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_map_field_ghost
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_data_rmsh

      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: xp
      INTEGER  , DIMENSION(:)         , INTENT(IN   )  :: ghostsize
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER        :: up
      REAL(mk)                                         :: tup
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#endif
#elif __MODE == __VEC
      INTEGER                         , INTENT(in)     :: lda
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: up
      REAL(mk) , DIMENSION(lda)                        :: tup
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_up
#endif     
#endif     
      INTEGER                         , INTENT(IN   )  :: Np
      INTEGER                         , INTENT(IN   )  :: topo_id, mesh_id
      INTEGER                         , INTENT(IN   )  :: kernel
      INTEGER                         , INTENT(  OUT)  :: info
      REAL(MK), OPTIONAL, DIMENSION(:,:,:),  POINTER   :: wx1_user,&
     &                                                    wx2_user,&
     &                                                    wx3_user

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)             :: len_phys
      REAL(MK), DIMENSION(:,:,:)   , POINTER   :: wx1,wx2,wx3
      REAL(MK), DIMENSION(:,:)     , POINTER   :: min_phys,max_phys
      REAL(MK), DIMENSION(ppm_dim)             :: dxi,dx
      REAL(MK)                                 :: dv1,dv2,dv3,epsilon
      INTEGER,  DIMENSION(ppm_dim)             :: Nc
      INTEGER,  DIMENSION(ppm_dim+2)           :: ldu,ldl
      REAL(MK)                                 :: kernel_support
      INTEGER                                  :: i,j,k,ii,jj,kk,iidec,l
      INTEGER                                  :: kkdec,ip1,ip2,ip3
      INTEGER                                  :: nb_sub,isub,dim,npart,isubl
      INTEGER                                  :: ifrom,ito,ip,max_partnumber
      INTEGER                                  :: ipart,nlist2,nlist1
      INTEGER                                  :: jjdec,iopt,idom,kkk
      INTEGER                                  :: ndata1_max,ndata2_max
      INTEGER                                  :: ndata3_max,nsubs,idim
      INTEGER,  DIMENSION(ppm_dim)             :: Nm
      INTEGER,  DIMENSION(:,:),      POINTER   :: ndata
      INTEGER,  DIMENSION(:,:),      POINTER   :: istart
      INTEGER,  DIMENSION(6)                   :: bcdef
      LOGICAL                                  :: consistent
      INTEGER                                  :: maptype, topoid, meshid
#if __MODE == __SCA
      INTEGER                                  :: lda = 1
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart('ppm_rmsh_remesh',t0,info)

      !-------------------------------------------------------------------------
      !  This is a hack! Somehow having trouble with constructors in the
      !  ppm_module_data_rmsh module
      !-------------------------------------------------------------------------
      ppm_rmsh_kernelsize = (/1,2,2,4/)

      dim = ppm_dim
      consistent = .TRUE.

      !-------------------------------------------------------------------------
      !  Check if arguments consistent
      !-------------------------------------------------------------------------
      IF(ppm_debug.GT.0) THEN
         IF(ppm_dim.EQ.3) THEN
            IF(   (PRESENT(wx1_user)).AND.&
     &            (PRESENT(wx2_user)).AND.&
     &            (PRESENT(wx3_user))) THEN
               consistent = .TRUE.
            ELSEIF(PRESENT(wx1_user).AND.(.NOT.(PRESENT(wx2_user))).AND. &
     &          (.NOT.(PRESENT(wx3_user)))) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             ' wx2_user or wx3_user weights lacking',__LINE__,info)
               GOTO 9999
            ELSEIF(PRESENT(wx2_user).AND.(.NOT.(PRESENT(wx3_user))).AND. &
     &         (.NOT.(PRESENT(wx1_user)))) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             ' wx1_user or wx3_user weights lacking',__LINE__,info)
               GOTO 9999
            ELSEIF(PRESENT(wx3_user).AND.(.NOT.(PRESENT(wx1_user))).AND. &
     &         (.NOT.(PRESENT(wx2_user)))) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             ' wx1_user or wx2_user weights lacking',__LINE__,info)
               GOTO 9999
            ENDIF
         ELSE
           IF(PRESENT(wx3_user)) THEN
               !----------------------------------------------------------------
               !  We dont care, maybe issue a warning..
               !----------------------------------------------------------------
            ELSEIF (PRESENT(wx1_user).AND..NOT.PRESENT(wx2_user)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'wx2_user not supplied',__LINE__,info)
               GOTO 9999
            ELSEIF (PRESENT(wx2_user).AND..NOT.PRESENT(wx1_user)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'wx1_user not supplied',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF

         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
         IF (Np .GT. 0) THEN
            IF (SIZE(xp,2) .LT. Np) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'not enough particles contained in xp',__LINE__,info)
               GOTO 9999
            ENDIF
            IF (SIZE(xp,1) .LT. ppm_dim) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'leading dimension of xp insufficient',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
         IF (Np .LE. 0) THEN
            IF (Np .LT. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'particles must be specified',__LINE__,info)
               GOTO 9999
            END IF
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  Check if kernel is implemented
         !----------------------------------------------------------------------
         IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &          'wrong kernel definition',__LINE__,info)
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  This can never happen, if it does then the library is sick
         !----------------------------------------------------------------------
         kernel_support = ppm_rmsh_kernelsize(kernel)*2
        !kernel_support = ppm_rmsh_kernelsize(kernel)*2
        !IF (.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) & 
        !     & .OR.(kernel_support.EQ.6))) THEN
        !   info = ppm_error_error
        !   CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
        !        &     'wrong kernel support',__LINE__,info)
        !   GOTO 9999
        !END IF
        
         IF (SIZE(up).LT.Np) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &         'wrong up array size',__LINE__,info)
            GOTO 9999
         ENDIF
         !----------------------------------------------------------------------
         !  This is also wrong: Ghostsize can be bigger...
         !----------------------------------------------------------------------
         !DO idim=1,dim
         !   IF (ghostsize(idim).NE.ppm_rmsh_kernelsize(kernel)) THEN
         !      info = ppm_error_error
         !      CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
         !           &     'wrong size of ghostsize',__LINE__,info)
         !      GOTO 9999
         !   END IF
         !END DO
         DO idim=1,dim
            IF (ghostsize(idim).LT.(ppm_rmsh_kernelsize(kernel)-1)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'ghostsize too small for this kernel',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF
     
      !-------------------------------------------------------------------------
      !  Check meshid and topoid validity
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         !----------------------------------------------------------------------
         !  Check that we have at least one topology - in that case the arrays
         !  must be associated, and we do not have to check this explicitly
         !  later
         !----------------------------------------------------------------------
         IF (topo_id .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &          'topo_id cannot be negative!',__LINE__,info)
            GOTO 9999
         ENDIF
        
         !----------------------------------------------------------------------
         !  check that the user topoid is less that the maximum length of the
         !  ppm_internal_topoid array
         !----------------------------------------------------------------------
         IF (topo_id .GT. UBOUND(ppm_internal_topoid,1)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &          'topo_id too large',__LINE__,info)
            GOTO 9999
         ENDIF
        
         !----------------------------------------------------------------------
         !  Check that the topology is defined
         !----------------------------------------------------------------------
         IF (topo_id .GE. 0) THEN
            IF (ppm_internal_topoid(topo_id) .EQ. -HUGE(topo_id)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'Desired topology is not defined',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF
     
     
      !-------------------------------------------------------------------------
      !  Get the internal topoid
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)
      IF(ppm_debug.GT.0) THEN
         IF (mesh_id.LT.LBOUND(ppm_meshid(topoid)%internal,1).OR.&
     &      mesh_id.GT.UBOUND(ppm_meshid(topoid)%internal,1)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &          'mesh_id too large',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (mesh_id .GT. 0) THEN
            IF (ppm_meshid(topoid)%internal(mesh_id) .EQ.    &
     &               -HUGE(mesh_id)) THEN
               info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &              'Desired mesh is not defined',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Get the internal meshid
      !-------------------------------------------------------------------------
      meshid = ppm_meshid(topoid)%internal(mesh_id)
     
      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------
      istart => ppm_cart_mesh(meshid,topoid)%istart
 
      IF(ppm_debug.GT.0) THEN   
         !----------------------------------------------------------------------
         !  Check weights again
         !----------------------------------------------------------------------
         SELECT CASE(ppm_dim)
           
         CASE(2)
            IF(PRESENT(wx1_user)) THEN
               IF(   .NOT.ASSOCIATED(wx1_user).OR.&
     &               .NOT.ASSOCIATED(wx2_user) ) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &                'cannot remesh with empty weights',__LINE__,info)
                  GOTO 9999
               END IF
            END IF
           
         CASE(3)
            IF(PRESENT(wx1_user)) THEN
               IF(   .NOT.ASSOCIATED(wx1_user).OR.&
     &               .NOT.ASSOCIATED(wx2_user).OR.&
     &               .NOT.ASSOCIATED(wx3_user)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &                'cannot remesh with empty weights',__LINE__,info)
                  GOTO 9999
               END IF
            END IF
           
         END SELECT

         IF(.NOT.PRESENT(wx1_user)) THEN
#if   __KIND == __SINGLE_PRECISION
            IF(.NOT.ASSOCIATED(wx1_s)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'compute weights before remeshing',__LINE__,info)
               GOTO 9999
            END IF
#elif __KIND == __DOUBLE_PRECISION
            IF(.NOT.ASSOCIATED(wx1_d)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'compute weights before remeshing',__LINE__,info)
               GOTO 9999
            END IF
#endif
         END IF
      END IF

      IF(np.EQ.0) GOTO 5555
      !-------------------------------------------------------------------------
      !  Done with Argument checks, do aliasing
      !-------------------------------------------------------------------------
      IF(PRESENT(wx1_user)) THEN
         wx1 => wx1_user
         wx2 => wx2_user
         IF(ppm_dim.EQ.3) wx3 => wx3_user
      ELSE
#if   __KIND == __SINGLE_PRECISION
         wx1 => wx1_s
         wx2 => wx2_s
         IF(ppm_dim.EQ.3) wx3 => wx3_s
#elif __KIND == __DOUBLE_PRECISION
         wx1 => wx1_d  
         wx2 => wx2_d
         IF(ppm_dim.EQ.3) wx3 => wx3_d
#endif
      END IF

 5555 CONTINUE
      !-------------------------------------------------------------------------
      !  Some more aliasing.
      !  WARNING! Did he get the internal/external stuff right?
      !-------------------------------------------------------------------------
      nsubs                = ppm_nsublist(topoid)
      ! ndata(1:dim,1:nsubs) = ppm_cart_mesh(topoid,meshid)%nnodes
      ndata                => ppm_cart_mesh(meshid,topoid)%nnodes
      Nm(1:dim)            = ppm_cart_mesh(meshid,topoid)%Nm
      bcdef(1:2*ppm_dim)   = ppm_bcdef(1:2*ppm_dim,topoid)

      !-------------------------------------------------------------------------
      !  Allocate memory for field if necessary
      !-------------------------------------------------------------------------
      ndata1_max = MAXVAL(ndata(1,ppm_isublist(1:nsubs,topoid)))
      ndata2_max = MAXVAL(ndata(2,ppm_isublist(1:nsubs,topoid)))
      IF(ppm_dim.EQ.3) ndata3_max=MAXVAL(ndata(3,ppm_isublist(1:nsubs,topoid)))
      ! WARNING: distinguish between sca and vector allocation
      ! WARNING: distinguish between 2d and 3d allocation!!
      IF(ASSOCIATED(field_up)) THEN
         iopt    = ppm_param_alloc_grow
      ELSE
         iopt    = ppm_param_alloc_fit
      END IF
#if   __DIME == __3D
#if   __MODE == __VEC
      !-------------------------------------------------------------------------
      !  Allocation for 3D vector
      !-------------------------------------------------------------------------
      ldl(1)     = 1
      ldl(2)     = 1 - ghostsize(1)
      ldl(3)     = 1 - ghostsize(2)
      ldl(4)     = 1 - ghostsize(3)
      ldl(5)     = 1
      ldu(1)     = lda
      ldu(2)     = ndata1_max + ghostsize(1)
      ldu(3)     = ndata2_max + ghostsize(2)
      ldu(4)     = ndata3_max + ghostsize(3)
      ldu(5)     = nsubs
#elif __MODE == __SCA
      !-------------------------------------------------------------------------
      !  Allocation for 3D scalar
      !-------------------------------------------------------------------------
      ldl(1)     = 1 - ghostsize(1)
      ldl(2)     = 1 - ghostsize(2)
      ldl(3)     = 1 - ghostsize(3)
      ldl(4)     = 1
      ldu(1)     = ndata1_max + ghostsize(1)
      ldu(2)     = ndata2_max + ghostsize(2)
      ldu(3)     = ndata3_max + ghostsize(3)
      ldu(4)     = nsubs
#endif
#elif __DIME == __2D
#if   __MODE == __VEC
      !-------------------------------------------------------------------------
      !  Allocation for 2D vector
      !-------------------------------------------------------------------------
      ldl(1)     = 1
      ldl(2)     = 1 - ghostsize(1)
      ldl(3)     = 1 - ghostsize(2)
      ldl(4)     = 1
      ldu(1)     = lda
      ldu(2)     = ndata1_max + ghostsize(1)
      ldu(3)     = ndata2_max + ghostsize(2)
      ldu(4)     = nsubs
#elif __MODE == __SCA
      !-------------------------------------------------------------------------
      !  Allocation for 2D vector
      !-------------------------------------------------------------------------
      ldl(1)     = 1 - ghostsize(1)
      ldl(2)     = 1 - ghostsize(2)
      ldl(3)     = 1
      ldu(1)     = ndata1_max + ghostsize(1)
      ldu(2)     = ndata2_max + ghostsize(2)
      ldu(3)     = nsubs
#endif     
#endif        
      CALL ppm_alloc(field_up,ldl,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_rmsh_remesh', &
     &       'field_up',__LINE__,info)
         GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Initialize the field
      !-------------------------------------------------------------------------
      field_up = 0.0_mk
      IF(np.EQ.0) GOTO 9999
      !-------------------------------------------------------------------------
      !  Recover particle lists
      !-------------------------------------------------------------------------
      !  dont need to do this as they are in ppm_module_data_rmsh
      
      !-------------------------------------------------------------------------
      !  Pre computing operations
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_phys => ppm_min_physs
      max_phys => ppm_max_physs
#elif __KIND == __DOUBLE_PRECISION
      min_phys => ppm_min_physd
      max_phys => ppm_max_physd
#endif  
      DO i = 1,ppm_dim
         Nc(i)       = Nm(i)-1
         len_phys(i) = max_phys(i,topoid) - min_phys(i,topoid)
         dx(i)       = len_phys(i)/REAL(Nc(i),MK)
      END DO
     
      dxi = 1.0_mk/dx
      dv1 = dx(1)
      dv2 = dx(2)
      if(ppm_dim.eq.3) dv3 = dx(3)

      !-------------------------------------------------------------------------
      !  Start computation, taking out the select case for speed
      !-------------------------------------------------------------------------
      SELECT CASE(kernel)
         !----------------------------------------------------------------------
         ! M prime 4
         !----------------------------------------------------------------------
      CASE(ppm_param_rmsh_kernel_mp4)
         DO isub = 1,nsubs
            isubl = ppm_isublist(isub,topoid)
            DO ip=1,store_info(isub)
               ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1,topoid) &
     &                    )*dxi(1))+2-istart(1,isubl)
               ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2,topoid) &
     &                    )*dxi(2))+2-istart(2,isubl)
#if   __DIME == __3D
               ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3,topoid) &
     &                    )*dxi(3))+2-istart(3,isubl)
#endif
#if   __MODE == __VEC
               tup(1:lda) = up(1:lda,list_sub(isub,ip))
#elif __MODE == __SCA
               tup        = up(list_sub(isub,ip))
#endif              
               !----------------------------------------------------------------
               !  flipped upside down for ppm conformance
               !----------------------------------------------------------------
               if(lda.eq.3.and.ppm_dim.eq.3) then
#if     __DIME == __3D
               kkdec = 0
               DO kk=ip3-1,ip3+2
                  kkdec = kkdec + 1                 
#endif
                  jjdec = 0
                  DO jj=ip2-1,ip2+2
                     jjdec = jjdec + 1
 
                     iidec = 0
                     DO ii=ip1-1,ip1+2
                        iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
!                       field_up(1:lda,ii,jj,kk,isub)&
!                            &=field_up(1:lda,ii,jj,kk,isub)&
!                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
!                            &*wx3(kkdec,isub,ip)*tup(1:lda)
                        field_up(1,ii,jj,kk,isub) &
     &                      =field_up(1,ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup(1)
                        field_up(2,ii,jj,kk,isub) &
     &                      =field_up(2,ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup(2)
                        field_up(3,ii,jj,kk,isub) &
     &                      =field_up(3,ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup(3)
                       
                       
#elif __MODE == __SCA
                        field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup
#endif                          
#elif __DIME == __2D
#if   __MODE == __VEC                       
                        field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *tup(1:lda)
#elif __MODE == __SCA
                        field_up(ii,jj,isub)=field_up(ii,jj,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *tup
#endif
#endif                       
                     END DO
                  END DO
#if  __DIME == __3D
               END DO
#endif                 
            ELSEIF (lda.eq.4.and.ppm_dim.eq.3) then
#if     __DIME == __3D
               kkdec = 0
               DO kk=ip3-1,ip3+2
                  kkdec = kkdec + 1
#endif
                  jjdec = 0
                  DO jj=ip2-1,ip2+2
                     jjdec = jjdec + 1

                     iidec = 0
                     DO ii=ip1-1,ip1+2
                        iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
!                       field_up(1:lda,ii,jj,kk,isub)&
!                            &=field_up(1:lda,ii,jj,kk,isub)&
!                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
!                            &*wx3(kkdec,isub,ip)*tup(1:lda)
                        field_up(1,ii,jj,kk,isub)&
                             &=field_up(1,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(1)
                        field_up(2,ii,jj,kk,isub)&
                             &=field_up(2,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(2)
                        field_up(3,ii,jj,kk,isub)&
                             &=field_up(3,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(3)
                        field_up(4,ii,jj,kk,isub)&
                             &=field_up(4,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(4)
                       
                       
#elif __MODE == __SCA
                        field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup
#endif                          
#elif __DIME == __2D
#if   __MODE == __VEC                       
                        field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub)&
                          &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                          &*tup(1:lda)
#elif __MODE == __SCA
                        field_up(ii,jj,isub)=field_up(ii,jj,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*tup
#endif
#endif                       
                     END DO
                  END DO
#if  __DIME == __3D
               END DO
#endif                 
            ELSE
#if     __DIME == __3D
              kkdec = 0
              DO kk=ip3-1,ip3+2
                 kkdec = kkdec + 1
                 
#endif
                 jjdec = 0
                 DO jj=ip2-1,ip2+2
                    jjdec = jjdec + 1

                    iidec = 0
                    DO ii=ip1-1,ip1+2
                       iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
                       DO l=1,lda
                          field_up(l,ii,jj,kk,isub)&
                               &=field_up(l,ii,jj,kk,isub)&
                               &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                               &*wx3(kkdec,isub,ip)*tup(l)
                       END DO
                       
#elif __MODE == __SCA
                       field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*wx3(kkdec,isub,ip)*tup
#endif                          
#elif __DIME == __2D
#if   __MODE == __VEC                       
                       field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub)&
                         &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                         &*tup(1:lda)
#elif __MODE == __SCA
                       field_up(ii,jj,isub)=field_up(ii,jj,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*tup
#endif
#endif                       
                    END DO
                 END DO
#if  __DIME == __3D
              END DO
#endif                 
           END IF
           END DO           
        END DO
        !-----------------------------------------------------------------------
        ! BS2
        !-----------------------------------------------------------------------
     CASE(ppm_param_rmsh_kernel_bsp2)
        DO isub = 1,nsubs
           isubl = ppm_isublist(isub,topoid)
           DO ip=1,store_info(isub)
              ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1,topoid)&
                   &)*dxi(1))+2-istart(1,isubl)
              ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2,topoid)&
                   &)*dxi(2))+2-istart(2,isubl)
#if   __DIME == __3D
              ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3,topoid)&
                   &)*dxi(3))+2-istart(3,isubl)
#endif              
              !-----------------------------------------------------------------
              !  flipped upside down for ppm conformance
              !-----------------------------------------------------------------
#if     __DIME == __3D
              kkdec = 0
              DO kk=ip3,ip3+1
                 kkdec = kkdec + 1
#endif
                 jjdec = 0
                 DO jj=ip2,ip2+1
                    jjdec = jjdec + 1

                    iidec = 0
                    DO ii=ip1,ip1+1
                       iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
                       field_up(1:lda,ii,jj,kk,isub)&
                            &=field_up(1:lda,ii,jj,kk,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*wx3(kkdec,isub,ip)*up(1:lda,list_sub(isub,ip))
#elif __MODE == __SCA
                       field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*wx3(kkdec,isub,ip)*up(list_sub(isub,ip))
#endif                          
#elif __DIME == __2D
#if   __MODE == __VEC                       
                       field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*up(1:lda,list_sub(isub,ip))
#elif __MODE == __SCA
                       field_up(ii,jj,isub)=field_up(ii,jj,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*up(list_sub(isub,ip))
#endif
#endif                       
                    END DO
                 END DO
                 
#if  __DIME == __3D
              END DO
#endif                 
           END DO
        END DO
     END SELECT

     !--------------------------------------------------------------------------
     !  Now map the ghosts in order to get consistent values at the border of
     !  the subdomains.
     !--------------------------------------------------------------------------
     maptype = ppm_param_map_init
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#endif     
     IF (info .NE. 0) GOTO 9999
         maptype = ppm_param_map_ghost_put
#if   __MODE == __SCA     
         CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
     &      info)
#elif __MODE == __VEC     
         CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
     &      info)
#endif     

     IF (info .NE. 0) GOTO 9999
     maptype = ppm_param_map_push
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#endif     

     IF (info .NE. 0) GOTO 9999
     maptype = ppm_param_map_send
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#endif     
     IF (info .NE. 0) GOTO 9999
     maptype = ppm_param_map_pop
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
          & info)
#endif     

     !--------------------------------------------------------------------------
     !  Return 
     !--------------------------------------------------------------------------
     9999 CONTINUE
     CALL substop('ppm_rmsh_remesh',t0,info)
     RETURN
     
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ss_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_sv_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ss_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_sv_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_dv_3d
#endif
#endif
#endif
