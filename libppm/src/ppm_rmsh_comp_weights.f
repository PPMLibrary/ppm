      !-------------------------------------------------------------------------
      !     Subroutine   :                   ppm_rmsh_comp_weights
      !-------------------------------------------------------------------------
      !     
      !     Purpose      : This routine computes the weights of the particles 
      !                    which need to be remeshed in 2D.
      !      
      !     Input        : xp(:,:)            (F) particle positions
      !                    Np                 (I) number of particles.
      !                    topo_id            (I) topology identifier (user
      !                                           numbering)
      !                                           of target
      !                    mesh_id            (I) id of the mesh (user)
      !                    kernel             (I) Choice of the kernel used to
      !                                           compute the weights.
      !     
      !     Input/Output : wx1,2,3_user(:,:,:)(F) Optional :
      !                                           if present, this routine does
      !                                           nothing, because these 
      !                                           weights are used after for 
      !                                           remeshing
      !     Output       : info               (I) return status
      !      
      !     Remarks      : 
      !     
      !     References   :
      !     
      !     Revisions    :
      !-------------------------------------------------------------------------
      !     $Log: ppm_rmsh_comp_weights.f,v $
      !     Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !     CBL version of the PPM library
      !
      !     Revision 1.18  2006/02/03 09:34:04  ivos
      !     Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !     local subs in topo_store. Several mapping routines however need the
      !     info about all (global) subs.
      !     Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !     occurrences.
      !
      !     Revision 1.17  2005/08/22 17:12:35  michaebe
      !     fixed wtime bug
      !
      !     Revision 1.16  2005/08/22 11:44:02  michaebe
      !     removed include of mpif.h
      !
      !     Revision 1.15  2005/06/04 00:23:12  michaebe
      !     put the deleted revision history back in.
      !
      !     Revision 1.14  2005/05/31 11:58:18  pchatela
      !     Fixed bsp2 (min corner missing in interpolation formula)
      !
      !     Revision 1.13  2005/03/10 01:51:13  ivos
      !     Resolved CVS conflict.
      !
      !     Revision 1.12  2005/03/07 22:46:55  michaebe
      !     forgot a bracket, sorry.
      !
      !     Revision 1.11  2005/03/07 22:03:47  michaebe
      !     fixed dble and min phys bug
      !
      !     Revision 1.10  2004/11/09 10:06:14  ivos
      !     CVS fix.
      !
      !     Revision 1.8  2004/09/29 07:27:50  michaebe
      !     put a if ppm_debug in front of a ppm_write statement
      !
      !     Revision 1.7  2004/08/18 15:03:33  michaebe
      !     bugfix for #0000021 (size vs. ubound)
      !
      !     Revision 1.6  2004/08/18 12:53:31  michaebe
      !     took the select case for the kernel out of the loop
      !
      !     Revision 1.5  2004/08/16 11:04:00  michaebe
      !     included handling of np=0 cases
      !
      !     Revision 1.4  2004/08/13 12:03:27  michaebe
      !     removed epsilon from left neighbor calc
      !
      !     Revision 1.3  2004/08/12 08:10:27  michaebe
      !     fixed some pointer nullification issues.
      !
      !     Revision 1.2  2004/08/10 16:11:36  michaebe
      !     Countless and major changes
      !
      !     Revision 1.1  2004/08/09 11:06:28  michaebe
      !     revised initial implementation (P. Gautier)
      !
      !-------------------------------------------------------------------------
      !     Parallel Particle Mesh Library (PPM)
      !     Institute of Computational Science
      !     ETH Zentrum, Hirschengraben 84
      !     CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_comp_weights_s(xp,Np, &
      & topo_id,mesh_id,kernel,info, &
      & wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_comp_weights_d(xp,Np, &
      & topo_id,mesh_id,kernel,info, &
      & wx1_user,wx2_user,wx3_user)
#endif
      !-------------------------------------------------------------------------
      !  INCLUDES
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_data_rmsh
      USE ppm_module_write
      USE ppm_module_util_time
      IMPLICIT NONE
            
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
     !--------------------------------------------------------------------------
     ! Arguments     
     !--------------------------------------------------------------------------
     REAL(MK), DIMENSION(:,:)       , POINTER       :: xp
     INTEGER                        , INTENT(IN   ) :: Np
     INTEGER                        , INTENT(IN   ) :: topo_id,mesh_id
     INTEGER                        , INTENT(IN   ) :: kernel
     INTEGER                        , INTENT(  OUT) :: info
     REAL(MK), OPTIONAL, DIMENSION(:,:,:), POINTER  :: wx1_user,&
          &                                            wx2_user,&
          &                                            wx3_user

     !--------------------------------------------------------------------------
     ! Local variables 
     !--------------------------------------------------------------------------
     INTEGER,  DIMENSION(:,:)     , POINTER :: istart
     INTEGER,  DIMENSION(:)       , POINTER :: ilist1,ilist2
     REAL(MK), DIMENSION(:,:)     , POINTER :: min_phys,max_phys
     REAL   ,  DIMENSION(ppm_dim)           :: dxi,dx
     REAL   ,  DIMENSION(ppm_dim)           :: len_phys
     REAL(MK)                               :: x1,x2,x3,epsilon
     INTEGER                                :: kernel_support
     INTEGER,  DIMENSION(ppm_dim+2)         :: ldu,ldl
     INTEGER,  DIMENSION(ppm_dim)           :: Nc
     INTEGER                                :: i,j,k,ii,jj,kk,iidec
     INTEGER                                :: jjdec,nb_sub,npart,ipart
     INTEGER                                :: kkdec,ip1,nbpt_z,nlist1
     INTEGER                                :: ip2,ip3,nbpt_x,nbpt_y,iface
     INTEGER                                :: isub,ifrom,ito,ip,iopt,isubl
     INTEGER                                :: max_partnumber,idom,nlist2,idoml
     INTEGER, DIMENSION(1)                  :: lda
     INTEGER, DIMENSION(ppm_dim)            :: Nm
     INTEGER                                :: nsubs
     INTEGER, DIMENSION(6)                  :: bcdef
     INTEGER                                :: topoid, meshid
     LOGICAL                                :: internal_weights
     ! aliases
     REAL(mk), DIMENSION(:,:,:),    POINTER :: min_sub, max_sub
     REAL(MK), DIMENSION(:,:,:)   , POINTER :: wx1,wx2,wx3
     REAL(mk)                               :: myeps
     REAL(mk)                               :: tim1s, tim1e
     CHARACTER(len=256)                     :: msg

     
     !--------------------------------------------------------------------------
     !  Initialise 
     !--------------------------------------------------------------------------

     !--------------------------------------------------------------------------
     !  This is a hack! Somehow having trouble with constructors in the
     !  ppm_module_data_rmsh module
     !--------------------------------------------------------------------------
     ppm_rmsh_kernelsize = (/1,2,2,4/)
     
     NULLIFY(wx1,wx2,wx3)
     
     CALL substart('ppm_rmsh_comp_weights',t0,info)

     internal_weights = .FALSE.

     IF(ppm_dim.EQ.3) THEN
        IF(PRESENT(wx1_user).AND.(.NOT.(PRESENT(wx2_user))).AND. &
             & (.NOT.(PRESENT(wx3_user)))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &  ' wx2_user or wx3_user weights lacking',__LINE__,info)
           GOTO 9999
        ELSEIF(PRESENT(wx2_user).AND.(.NOT.(PRESENT(wx3_user))).AND. &
             & (.NOT.(PRESENT(wx1_user)))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &  ' wx1_user or wx3_user weights lacking',__LINE__,info)
           GOTO 9999
        ELSEIF(PRESENT(wx3_user).AND.(.NOT.(PRESENT(wx1_user))).AND. &
             & (.NOT.(PRESENT(wx2_user)))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &  ' wx1_user or wx2_user weights lacking',__LINE__,info)
           GOTO 9999
        END IF
     ELSEIF (ppm_dim.EQ.2) THEN
        IF(PRESENT(wx1_user).AND.(.NOT.PRESENT(wx2_user)).OR.&
             &PRESENT(wx2_user).AND.(.NOT.PRESENT(wx1_user))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights', &
                &  'a wx_user array is missing',__LINE__,info)
           GOTO 9999
        END IF
     END IF
        
     !--------------------------------------------------------------------------
     !  Check arguments
     !--------------------------------------------------------------------------
     IF (ppm_debug .GT. 0) THEN 
        IF (Np .GT. 0) THEN
           IF (SIZE(xp,2) .LT. Np) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &  'not enough particles contained in xp',__LINE__,info)
              GOTO 9999
           ENDIF
           IF (SIZE(xp,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &  'leading dimension of xp insufficient',__LINE__,info)
              GOTO 9999
           ENDIF
        ENDIF
        IF (Np .LE. 0) THEN
           IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &  'particles must be specified',__LINE__,info)
              GOTO 9999
           END IF
           GOTO 9999
        END IF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &     'wrong kernel definition',__LINE__,info)
           GOTO 9999
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) & 
             &  .OR.(kernel_support.EQ.6))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &     'wrong kernel support',__LINE__,info)
           GOTO 9999
        ENDIF
     ENDIF
     
     !--------------------------------------------------------------------------
     !  Check meshid and topoid validity
     !--------------------------------------------------------------------------
     IF (ppm_debug .GT. 0) THEN
        !-----------------------------------------------------------------------
        !  Check that we have at least one topology - in that case the arrays
        !  must be associated, and we do not have to check this explicitly
        !  later
        !-----------------------------------------------------------------------
        IF (topo_id .LT. 0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                & 'topo_id cannot be negative!',__LINE__,info)
           GOTO 9999
        ENDIF
        
        !-----------------------------------------------------------------------
        !  check that the user topoid is less that the maximum length of the
        !  ppm_internal_topoid array
        !-----------------------------------------------------------------------
        IF (topo_id .GT. UBOUND(ppm_internal_topoid,1)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                & 'topo_id too large',__LINE__,info)
           GOTO 9999
        ENDIF
        
        !-----------------------------------------------------------------------
        !  Check that the topology is defined
        !-----------------------------------------------------------------------
        IF (topo_id .GE. 0) THEN
           IF (ppm_internal_topoid(topo_id) .EQ. -HUGE(topo_id)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   & 'Desired topology is not defined',__LINE__,info)
              GOTO 9999
           ENDIF
        ENDIF
     ENDIF

     IF(np.EQ.0) GOTO 9999
     
     !--------------------------------------------------------------------------
     !  Get the internal topoid
     !--------------------------------------------------------------------------
     topoid = ppm_internal_topoid(topo_id)
     IF(ppm_debug.GT.0) THEN
        IF (mesh_id.LT.LBOUND(ppm_meshid(topoid)%internal,1).OR.&
             & mesh_id.GT.UBOUND(ppm_meshid(topoid)%internal,1)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                & 'mesh_id too large',__LINE__,info)
           GOTO 9999
        ENDIF
        IF (mesh_id .GT. 0) THEN
           IF (ppm_meshid(topoid)%internal(mesh_id) .EQ.    &
                &               -HUGE(mesh_id)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   & 'Desired mesh is not defined',__LINE__,info)
              GOTO 9999
           ENDIF
        ENDIF
     END IF
     !--------------------------------------------------------------------------
     !  Get the internal meshid
     !--------------------------------------------------------------------------
     meshid = ppm_meshid(topoid)%internal(mesh_id)

     !--------------------------------------------------------------------------
     !  Get istart
     !--------------------------------------------------------------------------
     istart => ppm_cart_mesh(meshid,topoid)%istart
     
     !--------------------------------------------------------------------------
     !  Assignment of the useful arrays/scalar
     !--------------------------------------------------------------------------
     Nm(1:ppm_dim) = ppm_cart_mesh(meshid,topoid)%Nm
     bcdef(1:(2*ppm_dim)) = ppm_bcdef(1:(2*ppm_dim),topoid)
     nsubs = ppm_nsublist(topoid)
    
     !--------------------------------------------------------------------------
     !  Alloc memory for particle lists
     !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
     !  or a finalize routine.
     !--------------------------------------------------------------------------
     iopt   = ppm_param_alloc_fit
     ldu(1) = Np
     CALL ppm_alloc(ilist1,ldu,iopt,info)
     IF (info .NE. 0) THEN
       info = ppm_error_fatal
       CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights_3d',     &
     &               'particle list 1 ILIST1',__LINE__,info)
       GOTO 9999
     ENDIF
     CALL ppm_alloc(ilist2,ldu,iopt,info)
     IF (info .NE. 0) THEN
       info = ppm_error_fatal
       CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights_3d',     &
     &               'particle list 2 ILIST2',__LINE__,info)
       GOTO 9999
     ENDIF

     iopt   = ppm_param_alloc_fit
     ldu(1) = nsubs
     CALL ppm_alloc(store_info,ldu,iopt,info)
     IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights_3d',     &
     &                'store_info allocation : problem',__LINE__,info)
        GOTO 9999
     ENDIF

     !--------------------------------------------------------------------------
     !  Mesh spacing
     !--------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
     min_phys => ppm_min_physs
     max_phys => ppm_max_physs
#elif __KIND == __DOUBLE_PRECISION
     min_phys => ppm_min_physd
     max_phys => ppm_max_physd
#endif
     
     DO i = 1,ppm_dim
        Nc(i)       = Nm(i) - 1
        len_phys(i) = max_phys(i,topoid) - min_phys(i,topoid)
        dx(i)       = len_phys(i)/REAL(Nc(i),MK)
     ENDDO
     dxi     = 1.0_mk/dx
     epsilon = 0.000001_mk
     
     !--------------------------------------------------------------------------
     !  Initialize the particle list
     !--------------------------------------------------------------------------
     nlist1     = 0
     store_info = 0

     DO ipart=1,Np
        nlist1         = nlist1 + 1
        ilist1(nlist1) = ipart
     ENDDO
#if   __KIND == __SINGLE_PRECISION
     myeps = ppm_myepss
     min_sub => ppm_min_subs
     max_sub => ppm_max_subs
#elif __KIND == __DOUBLE_PRECISION
     myeps = ppm_myepsd
     min_sub => ppm_min_subd
     max_sub => ppm_max_subd
#endif

     !--------------------------------------------------------------------------
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !--------------------------------------------------------------------------
     DO idom = ppm_nsublist(topoid),1,-1   
        idoml = ppm_isublist(idom,topoid)
        !-----------------------------------------------------------------------
        !  Loop over the remaining particles 
        !-----------------------------------------------------------------------
        nlist2 = 0
        npart = 0
        DO i=1,nlist1
           ipart = ilist1(i)
           !--------------------------------------------------------------------
           !  If the particle is inside the current subdomain, assign it
           !--------------------------------------------------------------------
           IF(ppm_dim.EQ.3) THEN
              !-----------------------------------------------------------------
              ! the particle is in the closure of the subdomain
              !-----------------------------------------------------------------
              IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
                   &xp(3,ipart).GE.min_sub(3,idoml,topoid) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml,topoid) .AND. &
                   &xp(3,ipart).LE.max_sub(3,idoml,topoid) ) ) THEN

                 IF( (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(3,ipart).LT.max_sub(3,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(6,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

                    npart = npart + 1
                    store_info(idom) = npart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           ELSEIF (ppm_dim.EQ.2) THEN
              
              IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml,topoid) ) ) THEN
                 
                 IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                    
                    npart = npart + 1
                    store_info(idom) = npart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           END IF
        END DO
        !-----------------------------------------------------------------------
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !-----------------------------------------------------------------------
        IF (nlist2.NE.nlist1) THEN
           nlist1 = nlist2
           DO i=1,nlist1
              ilist1(i) = ilist2(i)
           ENDDO
        ENDIF
        
        !-----------------------------------------------------------------------
        !  Exit if the list is empty
        !-----------------------------------------------------------------------
        IF (nlist1.EQ.0) EXIT
     END DO
     !--------------------------------------------------------------------------
     !  Check that we sold all the particles 
     !--------------------------------------------------------------------------
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_rmsh_comp_weights',  &
             &       'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF
     
     !--------------------------------------------------------------------------
     ! select the particle for with several domain computation (french)
     !--------------------------------------------------------------------------
     
     max_partnumber = 0
     DO idom=1,ppm_nsublist(topoid)
        IF(store_info(idom).GE.max_partnumber) THEN
           max_partnumber = store_info(idom)
        END IF
     END DO
     iopt   = ppm_param_alloc_fit
     ldu(1) = ppm_nsublist(topoid)
     ldu(2) = max_partnumber
     
     CALL ppm_alloc(list_sub,ldu,iopt,info)
     IF(info.NE.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights_3d',     &
     &        'problem in internal allocation',__LINE__,info)
        GOTO 9999
     END IF
     
     list_sub=0
     
     !--------------------------------------------------------------------------
     !  Initialize the particle list
     !--------------------------------------------------------------------------
     nlist1     = 0
     DO ipart=1,Np
        nlist1         = nlist1 + 1
        ilist1(nlist1) = ipart
     ENDDO
     
     !--------------------------------------------------------------------------
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !--------------------------------------------------------------------------
     DO idom = ppm_nsublist(topoid),1,-1
        idoml = ppm_isublist(idom,topoid)
        !-----------------------------------------------------------------------
        !  loop over the remaining particles 
        !-----------------------------------------------------------------------
        nlist2 = 0
        npart = 0

        DO i=1,nlist1
           ipart = ilist1(i)
           
           !--------------------------------------------------------------------
           !  If the particle is inside the current subdomain, assign it
           !--------------------------------------------------------------------
           IF(ppm_dim.EQ.3) THEN
              IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
                   &xp(3,ipart).GE.min_sub(3,idoml,topoid) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml,topoid) .AND. &
                   &xp(3,ipart).LE.max_sub(3,idoml,topoid) ) ) THEN
                 
                 IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(3,ipart).LT.max_sub(3,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(6,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
                    
                    
                    npart = npart + 1
                    list_sub(idom,npart) = ipart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           ELSEIF (ppm_dim.EQ.2) THEN

              IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml,topoid) ) ) THEN
                                    
                 IF(   (xp(1,ipart).LT.max_sub(1,idom,topoid) .OR.  &
                      & (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
                      & (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                    
                    npart = npart + 1
                    list_sub(idom,npart) = ipart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           END IF
           
        END DO
        !-----------------------------------------------------------------------
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !-----------------------------------------------------------------------
        IF (nlist2.NE.nlist1) THEN
           nlist1 = nlist2
           DO i=1,nlist1
              ilist1(i) = ilist2(i)
           ENDDO
        ENDIF
        
        !-----------------------------------------------------------------------
        !  Exit if the list is empty
        !-----------------------------------------------------------------------
        IF (nlist1.EQ.0) EXIT
     END DO

     !--------------------------------------------------------------------------
     !  Check that we sold all the particles 
     !--------------------------------------------------------------------------
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_rmsh_comp_weights_3d',  &
             &       'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF
     
     !--------------------------------------------------------------------------
     !  Allocate and alias the weights if we need them.
     !--------------------------------------------------------------------------
     max_partnumber = 0
     DO idom = 1,ppm_nsublist(topoid)
        IF(store_info(idom).GE.max_partnumber) THEN
           max_partnumber = store_info(idom)  
        END IF
     END DO

     iopt = ppm_param_alloc_fit
     ldl(1) = 1
     ldl(2) = 1
     ldl(3) = 1 - ppm_rmsh_kernelsize(kernel)
     ldu(1) = ppm_rmsh_kernelsize(kernel)*2
     ldu(2) = ppm_nsublist(topoid)
     ldu(3) = max_partnumber + ppm_rmsh_kernelsize(kernel)
     
     IF(((PRESENT(wx1_user)).AND.PRESENT(wx2_user).AND. &
          & PRESENT(wx3_user).AND.ppm_dim.EQ.3).OR.&
          &(PRESENT(wx1_user).AND.PRESENT(wx2_user).AND.&
          ppm_dim.eq.2)) THEN
        internal_weights = .FALSE.

        IF (.NOT.(ASSOCIATED(wx1_user)).OR. &
             &   .NOT.(ASSOCIATED(wx2_user))) THEN
#if   __KIND == __SINGLE_PRECISION
           CALL ppm_alloc(wx1_s,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx1_user => wx1_s
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_alloc(wx1_d,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx1_user => wx1_d
#endif
#if   __KIND == __SINGLE_PRECISION
           CALL ppm_alloc(wx2_s,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx2_user => wx2_s
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_alloc(wx2_d,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx2_user => wx2_d
#endif
        END IF
        IF(PRESENT(wx3_user)) THEN
           IF(ppm_dim.EQ.3.AND.&
                &.NOT.ASSOCIATED(wx3_user)) THEN
#if   __KIND == __SINGLE_PRECISION
              CALL ppm_alloc(wx3_s,ldl,ldu,iopt,info)
              IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                      &            'pb in weights allocation',__LINE__,info)
                 GOTO 9999
              END IF
              
              wx3_user => wx3_s
#elif __KIND == __DOUBLE_PRECISION
              CALL ppm_alloc(wx3_d,ldl,ldu,iopt,info)
              IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                  &            'pb in weights allocation',__LINE__,info)
                 GOTO 9999
              END IF
              wx3_user => wx3_d
#endif       
           END IF
        END IF
        wx1 => wx1_user
        wx2 => wx2_user
        IF(ppm_dim.EQ.3) wx3 => wx3_user
     ELSE
#if   __KIND == __SINGLE_PRECISION
        CALL ppm_alloc(wx1_s,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
            &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx1 => wx1_s
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_alloc(wx1_d,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx1 => wx1_d
#endif       
#if   __KIND == __SINGLE_PRECISION
        CALL ppm_alloc(wx2_s,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx2 => wx2_s
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_alloc(wx2_d,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx2 => wx2_d
#endif
        IF(ppm_dim.EQ.3) THEN
#if   __KIND == __SINGLE_PRECISION
           CALL ppm_alloc(wx3_s,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx3 => wx3_s
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_alloc(wx3_d,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx3 => wx3_d
#endif       
        END IF
     END IF
     !--------------------------------------------------------------------------
     !  Initialize Weights
     !--------------------------------------------------------------------------
     CALL ppm_util_time(tim1s)
     wx1 = 0.0_mk
     wx2 = 0.0_mk
     IF(ppm_dim.EQ.3) wx3 = 0.0_mk
     !--------------------------------------------------------------------------
     !  Beginning of the computation
     !--------------------------------------------------------------------------
     !  loop over subs
     SELECT CASE(kernel)
        
     CASE(ppm_param_rmsh_kernel_bsp2)
        
        DO isub = 1,ppm_nsublist(topoid)
           !  loop over particles therein
           DO ip = 1,store_info(isub)
              isubl = ppm_isublist(isub,topoid)
              
              ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1,topoid) &
     &             )*dxi(1))+2-istart(1,isubl)
              
              ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2,topoid) &
     &             )*dxi(2))+2-istart(2,isubl)
              IF(ppm_dim.EQ.3) THEN
                 ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3,topoid) &
     &             )*dxi(3))+2-istart(3,isubl)
              END IF
              
              iidec = 0
              jjdec = 0
              kkdec = 0
              
              !-----------------------------------------------------------------
              !  Bsplines 2
              !-----------------------------------------------------------------
              
              DO ii    = ip1,ip1+1
                 iidec = iidec + 1
                 x1 = ABS(REAL(ii-2+istart(1,isubl),mk)-( &
     &               xp(1,list_sub(isub,ip))-min_phys(1,topoid))*dxi(1))
                 
                 IF(x1.LE.1.0_mk) THEN            
                    wx1(iidec,isub,ip) =  1.0_mk - x1
                 END IF
              END DO
              
              ! Loop over y direction
              DO jj    = ip2,ip2+1
                 jjdec = jjdec + 1
                 x2 = ABS(REAL(jj-2+istart(2,isubl),mk)-( &
     &               xp(2,list_sub(isub,ip))-min_phys(2,topoid))*dxi(2))
                 
                 IF(x2.LE.1.0_mk) THEN     
                    wx2(jjdec,isub,ip) =  1.0_mk - x2
                 END IF
              END DO
              IF(ppm_dim.EQ.3) THEN              
                 ! Loop over z direction
                 DO kk    = ip3,ip3+1
                    kkdec = kkdec + 1
                    x3 = ABS(REAL(kk-2+istart(3,isubl),mk)- (&
     &                  xp(3,list_sub(isub,ip))-min_phys(3,topoid))*dxi(3)) 
                    
                    IF(x3.LE.1.0_mk) THEN     
                       wx3(kkdec,isub,ip) =  1.0_mk - x3
                    END IF
                 END DO
              END IF

           END DO

        END DO

     CASE(ppm_param_rmsh_kernel_mp4)

        !-----------------------------------------------------------------
        ! M Prime Four
        !-----------------------------------------------------------------
        DO isub = 1,ppm_nsublist(topoid)
           !  loop over particles therein
           DO ip = 1,store_info(isub)
              isubl = ppm_isublist(isub,topoid)
              
              ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1,topoid) &
                   &         )*dxi(1))+2-istart(1,isubl)
              
              ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2,topoid) &
                   &         )*dxi(2))+2-istart(2,isubl)
              IF(ppm_dim.EQ.3) THEN
                 ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3,topoid) &
                   &         )*dxi(3))+2-istart(3,isubl)
              END IF
              
              iidec = 0
              jjdec = 0
              kkdec = 0

              DO ii = ip1-1,ip1+2
                 iidec = iidec + 1
                 x1 = ABS(REAL(ii-2+istart(1,isubl),mk)-( &
     &               xp(1,list_sub(isub,ip))-min_phys(1,topoid))*dxi(1))
                 IF(x1.LE.2.0_mk) THEN
                    IF(x1.LE.1.0_mk) THEN
                       wx1(iidec,isub,ip) = 1.0_mk - x1**2*(2.5_mk-1.5_mk*x1)
                    ELSE
                       wx1(iidec,isub,ip) = 2.0_MK + (-4.0_MK + &
     &                     (2.5_MK - 0.5_MK * x1)*x1)*x1
                    END IF
                 END IF
              END DO
              
              DO jj = ip2-1,ip2+2
                 jjdec = jjdec + 1
                 x2 = ABS(REAL(jj-2+istart(2,isubl),mk)-( &
     &               xp(2,list_sub(isub,ip))-min_phys(2,topoid))*dxi(2))
                 
                 IF(x2.LE.2.0_mk) THEN
                    IF(x2.LE.1.0_mk) THEN
                       wx2(jjdec,isub,ip) = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                    ELSE
                       wx2(jjdec,isub,ip) = 2.0_MK + (-4.0_MK + &
     &                     (2.5_MK - 0.5_MK * x2)*x2)*x2
                    END IF
                 END IF
                 
              END DO
              
              IF(ppm_dim.EQ.3) THEN
                 
                 DO kk    = ip3 - 1, ip3+2
                    
                    kkdec = kkdec + 1
                    x3 = ABS(REAL(kk-2+istart(3,isubl),mk)- (&
     &                  xp(3,list_sub(isub,ip))-min_phys(3,topoid))*dxi(3))
                    IF(x3.LE.2.0_MK) THEN
                       IF(x3.LE.1.0_MK) THEN     
                          wx3(kkdec,isub,ip) =  1.0_MK - x3**2*(2.5_MK - &
                               & 1.5_MK*x3)
                       ELSE
                          wx3(kkdec,isub,ip) =  2.0_MK + (-4.0_MK + &
                               & (2.5_MK - 0.5_MK*x3)*x3)*x3
                       END IF
                    END IF
                 END DO
              END IF
           END DO
        END DO
     END SELECT

     CALL ppm_util_time(tim1e)

     IF (ppm_debug.gt.0) THEN
        WRITE(msg,*) 'naked calc of weights took [ms]: ',1000.0*(tim1e-tim1s)
        CALL ppm_write(ppm_rank,'ppm_rmsh_comp_weights',msg,info)
     END IF

      !-------------------------------------------------------------------------
      ! Deallocation of the arrays....
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      lda(1) = 0
      CALL ppm_alloc(ilist1,lda,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
    &        'ppm_rmsh_comp_weights', &
    &        'pb in ilist1 deallocation',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,lda,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &       'ppm_rmsh_comp_weights',    &
     &       'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Now why on earth would you want to deallocate them?
      !-------------------------------------------------------------------------
#ifdef __COMPLETELY_LOCO
      IF (internal_weights) THEN
         CALL ppm_alloc(wx1,lda,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc, &
     &           'ppm_rmsh_comp_weights',   &
     &           'pb in weight deallocation',__LINE__,info)
            GOTO 9999
         END IF
         CALL ppm_alloc(wx2,lda,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc, &
     &           'ppm_rmsh_comp_weights',   &
     &           'pb in weight deallocation',__LINE__,info)
            GOTO 9999
         END IF
         IF(ppm_dim.EQ.3) THEN
            CALL ppm_alloc(wx3,lda,iopt,info)
            IF (info.NE.0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_dealloc, &
     &             'ppm_rmsh_comp_weights',    &
     &             'pb in weight deallocation',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF
#endif
      IF(.NOT.internal_weights) THEN
         NULLIFY(wx1_s,wx1_d)
         NULLIFY(wx2_s,wx2_d)
         NULLIFY(wx3_s,wx3_d)
      END IF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_rmsh_comp_weights',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_rmsh_comp_weights_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_rmsh_comp_weights_d
#endif

     
