      !------------------------------------------------------------------------!
      !     Subroutine   :                 ppm_interp_p2m
      !------------------------------------------------------------------------!
      !     
      !     Purpose      : This routine does particle to mesh
      !                    interpolation. Vector cases for lda.LT.5 are
      !                    explicitly unrolled for the 3D version. All 3D
      !                    versions are explicitly unrolled over the
      !                    kernel, the 2D versions are not.
      !      
      !     Input        : xp(:,:)            (F) particle positions
      !                    Np                 (I) number of particles.
      !                    lda                (I) leading dimension of up
      !                    up([:,]:)          (F) particle weights
      !                    topo_id            (I) topology identifier (user
      !                                           numbering)
      !                                           of target
      !                    mesh_id            (I) id of the mesh (user)
      !                    kernel             (I) Choice of the kernel used to
      !                                           compute the weights.
      !                    ghostsize(:)       (I) ghost size
      !
      !     Output       : info               (I) return status
      !                    field_up(:..:)     (F) field from which to interp
      !      
      !     Remarks      : 
      !     
      !     References   :
      !     
      !     Revisions    :
      !------------------------------------------------------------------------!
      !     $Log: ppm_interp_p2m.f,v $
      !     Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !     CBL version of the PPM library
      !
      !     Revision 1.26  2006/09/04 18:34:48  pchatela
      !     Fixes and cleanups to make it compile for
      !     release-1-0
      !
      !     Revision 1.25  2006/06/17 18:48:14  michaebe
      !     memory leak. Fixed.
      !
      !     Revision 1.24  2006/04/11 16:29:26  pchatela
      !     Added the witch-hat.
      !     Not thoroughly tested...
      !
      !     Revision 1.23  2006/02/03 09:34:02  ivos
      !     Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !     local subs in topo_store. Several mapping routines however need the
      !     info about all (global) subs.
      !     Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !     occurrences.
      !
      !     Revision 1.22  2005/08/22 11:41:54  michaebe
      !     removed include of mpif.h
      !
      !     Revision 1.21  2005/08/12 11:49:48  michaebe
      !     np=0 deadlock fix.
      !
      !     Revision 1.20  2005/06/30 18:14:16  michaebe
      !     Fixed a bug in the bugfix below
      !
      !     Revision 1.19  2005/06/28 16:59:47  ivos
      !     Slight code optimizations (save 3 additions and 1 INT per point),
      !     includes ifort bug fix as reported by Simone.
      !
      !     Revision 1.18  2005/06/23 21:02:53  ivos
      !     Added unrolled scalar version and unrolled (kernel) general lda.
      !     All 3D versions are unrolled now. Changed to use ppm_check_topoid
      !     and ppm_check_meshid.
      !
      !     Revision 1.17  2005/06/23 18:25:47  ivos
      !     Added overloaded version for lda=1 and lda=2 in 3D. Added
      !     default version for all non-overloaded lda. Added vector version
      !     for 2D.
      !
      !     Revision 1.16  2005/06/06 23:54:25  michaebe
      !     updated header and changed interface(!!)
      !
      !     Revision 1.15  2005/05/29 20:14:56  michaebe
      !     some fixes for the 2D case
      !
      !     Revision 1.13  2004/12/07 18:12:05  michaebe
      !     floating point arithmetic mania resolved min_phys -> min_sub..
      !
      !     Revision 1.12  2004/12/06 13:53:12  ivos
      !     Resolved CVS conflict.
      !
      !     Revision 1.11  2004/12/06 12:44:36  michaebe
      !     made the 2d case be empty
      !
      !     Revision 1.10  2004/12/06 12:34:57  michaebe
      !     added remarks and removed some unused variables
      !
      !     Revision 1.9  2004/12/06 11:50:13  michaebe
      !     changed up and xp from pointer to intent(in)
      !
      !     Revision 1.8  2004/12/06 11:22:05  michaebe
      !     modified the header
      !
      !     Revision 1.7  2004/11/10 09:11:51  michaebe
      !     added lda=3 and 4
      !
      !     Revision 1.6  2004/11/09 10:05:10  ivos
      !     CVS fix.
      !
      !     Revision 1.4  2004/11/03 08:12:58  ivos
      !     bugfix: ghostsize should be INTENT(IN) and not POINTER.
      !
      !     Revision 1.3  2004/11/02 16:33:52  michaebe
      !     added scalar case
      !
      !     Revision 1.2  2004/11/02 16:03:05  michaebe
      !     set insanity default
      !
      !     Revision 1.1  2004/11/02 12:52:03  michaebe
      !     inimp
      !
      !------------------------------------------------------------------------!
      !     Parallel Particle Mesh Library (PPM)
      !     Institute of Computational Science
      !     ETH Zentrum, Hirschengraben 84
      !     CH-8092 Zurich, Switzerland
      !------------------------------------------------------------------------!


#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_interp_p2m_ss_2d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_interp_p2m_ds_2d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_interp_p2m_sv_2d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_interp_p2m_dv_2d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_interp_p2m_ss_3d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_interp_p2m_ds_3d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_interp_p2m_sv_3d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_interp_p2m_dv_3d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif
#endif       

      !------------------------------------------------------------------------!
      !  INCLUDES
      !------------------------------------------------------------------------!

      !------------------------------------------------------------------------!
      !  Modules 
      !------------------------------------------------------------------------!
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_data
      USE ppm_module_data_rmsh
      USE ppm_module_data_mesh
      USE ppm_module_write
      USE ppm_module_map_field_ghost
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      IMPLICIT NONE
            
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
     !-------------------------------------------------------------------------!
     ! Arguments     
     !-------------------------------------------------------------------------!
#if   __MODE == __SCA
     REAL(MK) , DIMENSION(:)         , INTENT(IN)     :: up
#if   __DIME == __2D
     REAL(MK) , DIMENSION(:,:,:    ) , POINTER        :: field_up
#elif __DIME == __3D
     REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#endif
#elif __MODE == __VEC
     INTEGER                         , INTENT(in)     :: lda
     REAL(MK) , DIMENSION(:,:)       , INTENT(IN)        :: up
#if   __DIME == __2D
     REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#elif __DIME == __3D
     REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_up
#endif     
#endif     
     REAL(MK), DIMENSION(:,:)       , INTENT(IN)    :: xp
     INTEGER , DIMENSION(:  )       , INTENT(IN)    :: ghostsize
     INTEGER                        , INTENT(IN   ) :: Np
     INTEGER                        , INTENT(IN   ) :: topo_id,mesh_id
     INTEGER                        , INTENT(IN   ) :: kernel
     INTEGER                        , INTENT(  OUT) :: info
     !-------------------------------------------------------------------------!
     ! Local variables 
     !-------------------------------------------------------------------------!
     INTEGER,  DIMENSION(:,:)     , POINTER :: istart, ndata
     INTEGER,  DIMENSION(:)       , POINTER :: ilist1,ilist2
     REAL(mk), DIMENSION(:,:)     , POINTER :: min_phys,max_phys
     REAL(mk), DIMENSION(ppm_dim)           :: dxi,dx
     REAL(mk)                               :: dxx,dxy,dxz,dxxi,dxyi,dxzi
     REAL(mk), DIMENSION(ppm_dim)           :: len_phys
     REAL(mk)                               :: x1,x2,x3
     INTEGER                                :: kernel_support
     INTEGER,  DIMENSION(ppm_dim+2)         :: ldu,ldl
     INTEGER,  DIMENSION(ppm_dim)           :: Nc
     INTEGER                                :: i,j,k,ii,jj,kk,iidec,maptype
     INTEGER                                :: jjdec,nb_sub,npart,ipart
     INTEGER                                :: kkdec,ip1,nlist1
     INTEGER                                :: ip2,ip3,iface
     INTEGER                                :: isub,ifrom,ito,ip,dim,iopt,isubl
     INTEGER                                :: max_partnumber,idom,nlist2,idoml
     INTEGER, DIMENSION(ppm_dim)            :: Nm
     INTEGER                                :: nsubs
     INTEGER, DIMENSION(6)                  :: bcdef
     INTEGER                                :: topoid, meshid, iq
     LOGICAL                                :: internal_weights,lok
     ! aliases
     REAL(mk), DIMENSION(:,:,:),    POINTER :: min_sub, max_sub
     REAL(mk)                               :: myeps
     REAL(mk)                               :: tim1s, tim1e
     REAL(mk)                               :: xp1,xp2,xp3
     REAL(mk)                               :: wx1,wx2,wx3
     REAL(mk), DIMENSION(ppm_dim)           :: x0
     REAL(mk)                               :: x01,x02,x03
     CHARACTER(len=256)                     :: msg
     !-------------------------------------------------------------------------!
     !  Variables for unrolled versions
     !-------------------------------------------------------------------------!
     REAL(mk) :: x10,x11,x12,x13,x20,x21,x22,x23,x30,x31,x32,x33
     REAL(mk) :: a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33
     INTEGER  :: ip10,ip11,ip12,ip13,ip20,ip21,ip22,ip23,ip30,ip31,ip32,ip33,ldn
     REAL(mk) :: a10a20a30
     REAL(mk) :: a10a20a31
     REAL(mk) :: a10a20a32
     REAL(mk) :: a10a20a33
     REAL(mk) :: a10a21a30
     REAL(mk) :: a10a21a31
     REAL(mk) :: a10a21a32
     REAL(mk) :: a10a21a33
     REAL(mk) :: a10a22a30
     REAL(mk) :: a10a22a31
     REAL(mk) :: a10a22a32
     REAL(mk) :: a10a22a33
     REAL(mk) :: a10a23a30
     REAL(mk) :: a10a23a31
     REAL(mk) :: a10a23a32
     REAL(mk) :: a10a23a33
     REAL(mk) :: a11a20a30
     REAL(mk) :: a11a20a31
     REAL(mk) :: a11a20a32
     REAL(mk) :: a11a20a33
     REAL(mk) :: a11a21a30
     REAL(mk) :: a11a21a31
     REAL(mk) :: a11a21a32
     REAL(mk) :: a11a21a33
     REAL(mk) :: a11a22a30
     REAL(mk) :: a11a22a31
     REAL(mk) :: a11a22a32
     REAL(mk) :: a11a22a33
     REAL(mk) :: a11a23a30
     REAL(mk) :: a11a23a31
     REAL(mk) :: a11a23a32
     REAL(mk) :: a11a23a33
     REAL(mk) :: a12a20a30
     REAL(mk) :: a12a20a31
     REAL(mk) :: a12a20a32
     REAL(mk) :: a12a20a33
     REAL(mk) :: a12a21a30
     REAL(mk) :: a12a21a31
     REAL(mk) :: a12a21a32
     REAL(mk) :: a12a21a33
     REAL(mk) :: a12a22a30
     REAL(mk) :: a12a22a31
     REAL(mk) :: a12a22a32
     REAL(mk) :: a12a22a33
     REAL(mk) :: a12a23a30
     REAL(mk) :: a12a23a31
     REAL(mk) :: a12a23a32
     REAL(mk) :: a12a23a33
     REAL(mk) :: a13a20a30
     REAL(mk) :: a13a20a31
     REAL(mk) :: a13a20a32
     REAL(mk) :: a13a20a33
     REAL(mk) :: a13a21a30
     REAL(mk) :: a13a21a31
     REAL(mk) :: a13a21a32
     REAL(mk) :: a13a21a33
     REAL(mk) :: a13a22a30
     REAL(mk) :: a13a22a31
     REAL(mk) :: a13a22a32
     REAL(mk) :: a13a22a33
     REAL(mk) :: a13a23a30
     REAL(mk) :: a13a23a31
     REAL(mk) :: a13a23a32
     REAL(mk) :: a13a23a33

     
     !-------------------------------------------------------------------------!
     !  Initialise 
     !-------------------------------------------------------------------------!

     !-------------------------------------------------------------------------!
     !  This is a hack! Somehow having trouble with constructors in the
     !  ppm_module_data_rmsh module
     !-------------------------------------------------------------------------!
     ppm_rmsh_kernelsize = (/1,2,2,4/)
     
     CALL substart('ppm_interp_p2m',t0,info)

     dim = ppm_dim
     internal_weights = .FALSE.

     !-------------------------------------------------------------------------!
     !  Check arguments
     !-------------------------------------------------------------------------!
     IF (ppm_debug .GT. 0) THEN 
        IF (.NOT. ppm_initialized) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_ppm_noinit,'ppm_interp_p2m',  &
     &               'Please call ppm_init first!',__LINE__,info)
            GOTO 9999
        ENDIF
        IF (Np .GT. 0) THEN
           IF (SIZE(xp,2) .LT. Np) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'not enough particles contained in xp',__LINE__,info)
              GOTO 9999
           ENDIF
           IF (SIZE(xp,1) .LT.dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'leading dimension of xp insufficient',__LINE__,info)
              GOTO 9999
           ENDIF
        ENDIF
        IF (Np .LE. 0) THEN
           IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'particles must be specified',__LINE__,info)
              GOTO 9999
           END IF
           GOTO 9999
        END IF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'wrong kernel definition',__LINE__,info)
           GOTO 9999
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) & 
     &               .OR.(kernel_support.EQ.6))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'wrong kernel support',__LINE__,info)
           GOTO 9999
        END IF

     END IF
     
     !-------------------------------------------------------------------------!
     !  Check meshid and topoid validity
     !-------------------------------------------------------------------------!
     IF (ppm_debug .GT. 0) THEN
        CALL ppm_check_topoid(ppm_param_id_user,topo_id,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                 'topo_id is invalid!',__LINE__,info)
           GOTO 9999
        ENDIF
     END IF

     !-------------------------------------------------------------------------!
     !  Get the internal topoid and check the meshid
     !-------------------------------------------------------------------------!
     topoid = ppm_internal_topoid(topo_id)
     IF(ppm_debug.GT.0) THEN
        CALL ppm_check_meshid(ppm_param_id_user,mesh_id,topoid,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                 'mesh_id is invalid!',__LINE__,info)
           GOTO 9999
        ENDIF
     END IF

     !-------------------------------------------------------------------------!
     !  Get the internal meshid
     !-------------------------------------------------------------------------!
     meshid = ppm_meshid(topoid)%internal(mesh_id)
     
     !-------------------------------------------------------------------------!
     !  Get istart
     !-------------------------------------------------------------------------!
     istart => ppm_cart_mesh(meshid,topoid)%istart
     
     !-------------------------------------------------------------------------!
     !  Assignment of the useful arrays/scalar
     !-------------------------------------------------------------------------!
     Nm(1:dim) = ppm_cart_mesh(meshid,topoid)%Nm
     bcdef(1:(2*dim)) = ppm_bcdef(1:(2*dim),topoid)
     nsubs = ppm_nsublist(topoid)
    
     !-------------------------------------------------------------------------!
     !  If there is nothing to do, do nearly nothing
     !-------------------------------------------------------------------------!
     IF(Np.EQ.0) GOTO 9998
     
     !-------------------------------------------------------------------------!
     !  Alloc memory for particle lists
     !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
     !  or a finalize routine.
     !-------------------------------------------------------------------------!
     iopt   = ppm_param_alloc_fit
     ldu(1) = Np
     CALL ppm_alloc(ilist1,ldu,iopt,info)
     IF (info .NE. 0) THEN
       info = ppm_error_fatal
       CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                    'particle list 1 ILIST1',__LINE__,info)
       GOTO 9999
     ENDIF
     CALL ppm_alloc(ilist2,ldu,iopt,info)
     IF (info .NE. 0) THEN
       info = ppm_error_fatal
       CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                    'particle list 2 ILIST2',__LINE__,info)
       GOTO 9999
     ENDIF

     iopt   = ppm_param_alloc_fit
     ldu(1) = nsubs
     CALL ppm_alloc(store_info,ldu,iopt,info)
     IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'store_info allocation : problem',__LINE__,info)
        GOTO 9999
     ENDIF

     !-------------------------------------------------------------------------!
     !  Mesh spacing
     !-------------------------------------------------------------------------!
#if   __KIND == __SINGLE_PRECISION
     min_phys => ppm_min_physs
     max_phys => ppm_max_physs
#elif __KIND == __DOUBLE_PRECISION
     min_phys => ppm_min_physd
     max_phys => ppm_max_physd
#endif
     
     DO i = 1,dim
        Nc(i)       = Nm(i) - 1
        len_phys(i) = max_phys(i,topoid) - min_phys(i,topoid)
        dx(i)       = len_phys(i)/REAL(Nc(i),MK)
     ENDDO
     dxx = len_phys(1)/REAL(nc(1),mk)
     dxy = len_phys(2)/REAL(nc(2),mk)
     IF(dim.EQ.3) dxz = len_phys(3)/REAL(nc(3),mk)
     dxxi = 1.0_mk/dxx
     dxyi = 1.0_mk/dxy
     IF(dim.EQ.3) dxzi = 1.0_mk/dxz


     
     !-------------------------------------------------------------------------!
     !  Initialize the particle list
     !-------------------------------------------------------------------------!
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

     !-------------------------------------------------------------------------!
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !-------------------------------------------------------------------------!
     DO idom = ppm_nsublist(topoid),1,-1   
        idoml = ppm_isublist(idom,topoid)
        !----------------------------------------------------------------------!
        !  Loop over the remaining particles 
        !----------------------------------------------------------------------!
        nlist2 = 0
        npart = 0
        DO i=1,nlist1
           ipart = ilist1(i)

           !-------------------------------------------------------------------!
           !  If the particle is inside the current subdomain, assign it
           !-------------------------------------------------------------------!
           IF(ppm_dim.EQ.3) THEN
              !----------------------------------------------------------------!
              ! the particle is in the closure of the subdomain
              !----------------------------------------------------------------!
              IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
     &                   xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                   xp(3,ipart).GE.min_sub(3,idoml,topoid) .AND. &
     &                   xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                   xp(2,ipart).LE.max_sub(2,idoml,topoid) .AND. &
     &                   xp(3,ipart).LE.max_sub(3,idoml,topoid) ) ) THEN

                 IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                      (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                      (xp(3,ipart).LT.max_sub(3,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(6,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

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
     &                   xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                   xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                   xp(2,ipart).LE.max_sub(2,idoml,topoid) ) ) THEN

                 IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                      (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

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
        !----------------------------------------------------------------------!
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !----------------------------------------------------------------------!
        IF (nlist2.NE.nlist1) THEN
           nlist1 = nlist2
           DO i=1,nlist1
              ilist1(i) = ilist2(i)
           ENDDO
        ENDIF

        !----------------------------------------------------------------------!
        !  Exit if the list is empty
        !----------------------------------------------------------------------!
        IF (nlist1.EQ.0) EXIT
     END DO
     !-------------------------------------------------------------------------!
     !  Check that we sold all the particles 
     !-------------------------------------------------------------------------!
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_interp_p2m',  &
     &                    'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF
     
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
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'problem in internal allocation',__LINE__,info)
        GOTO 9999
     END IF
     
     list_sub=0
     
     !-------------------------------------------------------------------------!
     !  Initialize the particle list
     !-------------------------------------------------------------------------!
     nlist1     = 0
     DO ipart=1,Np
        nlist1         = nlist1 + 1
        ilist1(nlist1) = ipart
     ENDDO
     
     !-------------------------------------------------------------------------!
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !-------------------------------------------------------------------------!
     DO idom = ppm_nsublist(topoid),1,-1
        idoml = ppm_isublist(idom,topoid)
        !----------------------------------------------------------------------!
        !  loop over the remaining particles 
        !----------------------------------------------------------------------!
        nlist2 = 0
        npart = 0

        DO i=1,nlist1
           ipart = ilist1(i)
           !-------------------------------------------------------------------!
           !  If the particle is inside the current subdomain, assign it
           !-------------------------------------------------------------------!
           IF(ppm_dim.EQ.3) THEN
              IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
     &                   xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                   xp(3,ipart).GE.min_sub(3,idoml,topoid) .AND. &
     &                   xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                   xp(2,ipart).LE.max_sub(2,idoml,topoid) .AND. &
     &                   xp(3,ipart).LE.max_sub(3,idoml,topoid) ) ) THEN
                 
                 IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                      (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                      (xp(3,ipart).LT.max_sub(3,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(6,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
                    
                    
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
     &                   xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                   xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                   xp(2,ipart).LE.max_sub(2,idoml,topoid) ) ) THEN
                                    
                 IF(   (xp(1,ipart).LT.max_sub(1,idom,topoid) .OR.  &
     &                       (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                      (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                       (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                       bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                    
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
        !----------------------------------------------------------------------!
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !----------------------------------------------------------------------!
        IF (nlist2.NE.nlist1) THEN
           nlist1 = nlist2
           DO i=1,nlist1
              ilist1(i) = ilist2(i)
           ENDDO
        ENDIF
        
        !----------------------------------------------------------------------!
        !  Exit if the list is empty
        !----------------------------------------------------------------------!
        IF (nlist1.EQ.0) EXIT
     END DO

     !-------------------------------------------------------------------------!
     !  Check that we sold all the particles 
     !-------------------------------------------------------------------------!
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_interp_p2m_3d',  &
     &                    'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF
     
     !-------------------------------------------------------------------------!
     !  Allocate and alias the weights if we need them.
     !-------------------------------------------------------------------------!
     max_partnumber = 0
     DO idom = 1,ppm_nsublist(topoid)
        IF(store_info(idom).GE.max_partnumber) THEN
           max_partnumber = store_info(idom)  
        END IF
     END DO

9998 CONTINUE
#if __DIME == __3D
     !-------------------------------------------------------------------------!
     !  --- 3D ---
     !-------------------------------------------------------------------------!
     !  loop over subs
     ndata => ppm_cart_mesh(meshid,topoid)%nnodes

     DO isub = 1,ppm_nsublist(topoid)
        isubl = ppm_isublist(isub,topoid)
        
        DO k=1-ghostsize(3),ndata(3,isubl)+ghostsize(3)

           DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)

              DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
#if __MODE == __VEC
                 DO ldn=1,lda
                    field_up(ldn,i,j,k,isub) = 0.0_mk
                 END DO
#else      
                 field_up(i,j,k,isub) = 0.0_mk
#endif
              END DO
           END DO
        END DO
     END DO
     IF(np.EQ.0) GOTO 9997
     
     SELECT CASE(kernel)
        
     CASE(ppm_param_rmsh_kernel_mp4)
        
        !----------------------------------------------------------------------!
        ! M Prime Four
        !----------------------------------------------------------------------!
        DO isub = 1,ppm_nsublist(topoid)
           
#if __MODE == __VEC
           !-------------------------------------------------------------------!
           !  Unrolled versions for 4-vectors
           !-------------------------------------------------------------------!
           IF(lda.EQ.4) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)

                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01)
                 ip20 = INT(x02)
                 ip30 = INT(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1
                 
                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33


                 field_up(1,ip10,ip20,ip30,isub) = field_up(1,ip10,ip20,ip30,isub)+ &
     &                       a10*a20*a30*up(1,iq)
                 field_up(1,ip10,ip20,ip31,isub) = field_up(1,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(1,iq)
                 field_up(1,ip10,ip20,ip32,isub) = field_up(1,ip10,ip20,ip32,isub) + &
     &                       a10*a20*a32*up(1,iq)
                 field_up(1,ip10,ip20,ip33,isub) = field_up(1,ip10,ip20,ip33,isub) + &
     &                       a10*a20*a33*up(1,iq)
                 field_up(1,ip10,ip21,ip30,isub) = field_up(1,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(1,iq)
                 field_up(1,ip10,ip21,ip31,isub) = field_up(1,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(1,iq)
                 field_up(1,ip10,ip21,ip32,isub) = field_up(1,ip10,ip21,ip32,isub) + &
     &                       a10*a21*a32*up(1,iq)
                 field_up(1,ip10,ip21,ip33,isub) = field_up(1,ip10,ip21,ip33,isub) + &
     &                       a10*a21*a33*up(1,iq)
                 field_up(1,ip10,ip22,ip30,isub) = field_up(1,ip10,ip22,ip30,isub) + &
     &                       a10*a22*a30*up(1,iq)
                 field_up(1,ip10,ip22,ip31,isub) = field_up(1,ip10,ip22,ip31,isub) + &
     &                       a10*a22*a31*up(1,iq)
                 field_up(1,ip10,ip22,ip32,isub) = field_up(1,ip10,ip22,ip32,isub) + &
     &                       a10*a22*a32*up(1,iq)
                 field_up(1,ip10,ip22,ip33,isub) = field_up(1,ip10,ip22,ip33,isub) + &
     &                       a10*a22*a33*up(1,iq)
                 field_up(1,ip10,ip23,ip30,isub) = field_up(1,ip10,ip23,ip30,isub) + &
     &                       a10*a23*a30*up(1,iq)
                 field_up(1,ip10,ip23,ip31,isub) = field_up(1,ip10,ip23,ip31,isub) + &
     &                       a10*a23*a31*up(1,iq)
                 field_up(1,ip10,ip23,ip32,isub) = field_up(1,ip10,ip23,ip32,isub) + &
     &                       a10*a23*a32*up(1,iq)
                 field_up(1,ip10,ip23,ip33,isub) = field_up(1,ip10,ip23,ip33,isub) + &
     &                       a10*a23*a33*up(1,iq)
                 field_up(1,ip11,ip20,ip30,isub) = field_up(1,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(1,iq)
                 field_up(1,ip11,ip20,ip31,isub) = field_up(1,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(1,iq)
                 field_up(1,ip11,ip20,ip32,isub) = field_up(1,ip11,ip20,ip32,isub) + &
     &                       a11*a20*a32*up(1,iq)
                 field_up(1,ip11,ip20,ip33,isub) = field_up(1,ip11,ip20,ip33,isub) + &
     &                       a11*a20*a33*up(1,iq)
                 field_up(1,ip11,ip21,ip30,isub) = field_up(1,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(1,iq)
                 field_up(1,ip11,ip21,ip31,isub) = field_up(1,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(1,iq)
                 field_up(1,ip11,ip21,ip32,isub) = field_up(1,ip11,ip21,ip32,isub) + &
     &                       a11*a21*a32*up(1,iq)
                 field_up(1,ip11,ip21,ip33,isub) = field_up(1,ip11,ip21,ip33,isub) + &
     &                       a11*a21*a33*up(1,iq)
                 field_up(1,ip11,ip22,ip30,isub) = field_up(1,ip11,ip22,ip30,isub) + &
     &                       a11*a22*a30*up(1,iq)
                 field_up(1,ip11,ip22,ip31,isub) = field_up(1,ip11,ip22,ip31,isub) + &
     &                       a11*a22*a31*up(1,iq)
                 field_up(1,ip11,ip22,ip32,isub) = field_up(1,ip11,ip22,ip32,isub) + &
     &                       a11*a22*a32*up(1,iq)
                 field_up(1,ip11,ip22,ip33,isub) = field_up(1,ip11,ip22,ip33,isub) + &
     &                       a11*a22*a33*up(1,iq)
                 field_up(1,ip11,ip23,ip30,isub) = field_up(1,ip11,ip23,ip30,isub) + &
     &                       a11*a23*a30*up(1,iq)
                 field_up(1,ip11,ip23,ip31,isub) = field_up(1,ip11,ip23,ip31,isub) + &
     &                       a11*a23*a31*up(1,iq)
                 field_up(1,ip11,ip23,ip32,isub) = field_up(1,ip11,ip23,ip32,isub) + &
     &                       a11*a23*a32*up(1,iq)
                 field_up(1,ip11,ip23,ip33,isub) = field_up(1,ip11,ip23,ip33,isub) + &
     &                       a11*a23*a33*up(1,iq)
                 field_up(1,ip12,ip20,ip30,isub) = field_up(1,ip12,ip20,ip30,isub) + &
     &                       a12*a20*a30*up(1,iq)
                 field_up(1,ip12,ip20,ip31,isub) = field_up(1,ip12,ip20,ip31,isub) + &
     &                       a12*a20*a31*up(1,iq)
                 field_up(1,ip12,ip20,ip32,isub) = field_up(1,ip12,ip20,ip32,isub) + &
     &                       a12*a20*a32*up(1,iq)
                 field_up(1,ip12,ip20,ip33,isub) = field_up(1,ip12,ip20,ip33,isub) + &
     &                       a12*a20*a33*up(1,iq)
                 field_up(1,ip12,ip21,ip30,isub) = field_up(1,ip12,ip21,ip30,isub) + &
     &                       a12*a21*a30*up(1,iq)
                 field_up(1,ip12,ip21,ip31,isub) = field_up(1,ip12,ip21,ip31,isub) + &
     &                       a12*a21*a31*up(1,iq)
                 field_up(1,ip12,ip21,ip32,isub) = field_up(1,ip12,ip21,ip32,isub) + &
     &                       a12*a21*a32*up(1,iq)
                 field_up(1,ip12,ip21,ip33,isub) = field_up(1,ip12,ip21,ip33,isub) + &
     &                       a12*a21*a33*up(1,iq)
                 field_up(1,ip12,ip22,ip30,isub) = field_up(1,ip12,ip22,ip30,isub) + &
     &                       a12*a22*a30*up(1,iq)
                 field_up(1,ip12,ip22,ip31,isub) = field_up(1,ip12,ip22,ip31,isub) + &
     &                       a12*a22*a31*up(1,iq)
                 field_up(1,ip12,ip22,ip32,isub) = field_up(1,ip12,ip22,ip32,isub) + &
     &                       a12*a22*a32*up(1,iq)
                 field_up(1,ip12,ip22,ip33,isub) = field_up(1,ip12,ip22,ip33,isub) + &
     &                       a12*a22*a33*up(1,iq)
                 field_up(1,ip12,ip23,ip30,isub) = field_up(1,ip12,ip23,ip30,isub) + &
     &                       a12*a23*a30*up(1,iq)
                 field_up(1,ip12,ip23,ip31,isub) = field_up(1,ip12,ip23,ip31,isub) + &
     &                       a12*a23*a31*up(1,iq)
                 field_up(1,ip12,ip23,ip32,isub) = field_up(1,ip12,ip23,ip32,isub) + &
     &                       a12*a23*a32*up(1,iq)
                 field_up(1,ip12,ip23,ip33,isub) = field_up(1,ip12,ip23,ip33,isub) + &
     &                       a12*a23*a33*up(1,iq)
                 field_up(1,ip13,ip20,ip30,isub) = field_up(1,ip13,ip20,ip30,isub) + &
     &                       a13*a20*a30*up(1,iq)
                 field_up(1,ip13,ip20,ip31,isub) = field_up(1,ip13,ip20,ip31,isub) + &
     &                       a13*a20*a31*up(1,iq)
                 field_up(1,ip13,ip20,ip32,isub) = field_up(1,ip13,ip20,ip32,isub) + &
     &                       a13*a20*a32*up(1,iq)
                 field_up(1,ip13,ip20,ip33,isub) = field_up(1,ip13,ip20,ip33,isub) + &
     &                       a13*a20*a33*up(1,iq)
                 field_up(1,ip13,ip21,ip30,isub) = field_up(1,ip13,ip21,ip30,isub) + &
     &                       a13*a21*a30*up(1,iq)
                 field_up(1,ip13,ip21,ip31,isub) = field_up(1,ip13,ip21,ip31,isub) + &
     &                       a13*a21*a31*up(1,iq)
                 field_up(1,ip13,ip21,ip32,isub) = field_up(1,ip13,ip21,ip32,isub) + &
     &                       a13*a21*a32*up(1,iq)
                 field_up(1,ip13,ip21,ip33,isub) = field_up(1,ip13,ip21,ip33,isub) + &
     &                       a13*a21*a33*up(1,iq)
                 field_up(1,ip13,ip22,ip30,isub) = field_up(1,ip13,ip22,ip30,isub) + &
     &                       a13*a22*a30*up(1,iq)
                 field_up(1,ip13,ip22,ip31,isub) = field_up(1,ip13,ip22,ip31,isub) + &
     &                       a13*a22*a31*up(1,iq)
                 field_up(1,ip13,ip22,ip32,isub) = field_up(1,ip13,ip22,ip32,isub) + &
     &                       a13*a22*a32*up(1,iq)
                 field_up(1,ip13,ip22,ip33,isub) = field_up(1,ip13,ip22,ip33,isub) + &
     &                       a13*a22*a33*up(1,iq)
                 field_up(1,ip13,ip23,ip30,isub) = field_up(1,ip13,ip23,ip30,isub) + &
     &                       a13*a23*a30*up(1,iq)
                 field_up(1,ip13,ip23,ip31,isub) = field_up(1,ip13,ip23,ip31,isub) + &
     &                       a13*a23*a31*up(1,iq)
                 field_up(1,ip13,ip23,ip32,isub) = field_up(1,ip13,ip23,ip32,isub) + &
     &                       a13*a23*a32*up(1,iq)
                 field_up(1,ip13,ip23,ip33,isub) = field_up(1,ip13,ip23,ip33,isub) + &
     &                       a13*a23*a33*up(1,iq)


                 field_up(2,ip10,ip20,ip30,isub) = field_up(2,ip10,ip20,ip30,isub) + &
     &                       a10*a20*a30*up(2,iq)
                 field_up(2,ip10,ip20,ip31,isub) = field_up(2,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(2,iq)
                 field_up(2,ip10,ip20,ip32,isub) = field_up(2,ip10,ip20,ip32,isub) + &
     &                       a10*a20*a32*up(2,iq)
                 field_up(2,ip10,ip20,ip33,isub) = field_up(2,ip10,ip20,ip33,isub) + &
     &                       a10*a20*a33*up(2,iq)
                 field_up(2,ip10,ip21,ip30,isub) = field_up(2,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(2,iq)
                 field_up(2,ip10,ip21,ip31,isub) = field_up(2,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(2,iq)
                 field_up(2,ip10,ip21,ip32,isub) = field_up(2,ip10,ip21,ip32,isub) + &
     &                       a10*a21*a32*up(2,iq)
                 field_up(2,ip10,ip21,ip33,isub) = field_up(2,ip10,ip21,ip33,isub) + &
     &                       a10*a21*a33*up(2,iq)
                 field_up(2,ip10,ip22,ip30,isub) = field_up(2,ip10,ip22,ip30,isub) + &
     &                       a10*a22*a30*up(2,iq)
                 field_up(2,ip10,ip22,ip31,isub) = field_up(2,ip10,ip22,ip31,isub) + &
     &                       a10*a22*a31*up(2,iq)
                 field_up(2,ip10,ip22,ip32,isub) = field_up(2,ip10,ip22,ip32,isub) + &
     &                       a10*a22*a32*up(2,iq)
                 field_up(2,ip10,ip22,ip33,isub) = field_up(2,ip10,ip22,ip33,isub) + &
     &                       a10*a22*a33*up(2,iq)
                 field_up(2,ip10,ip23,ip30,isub) = field_up(2,ip10,ip23,ip30,isub) + &
     &                       a10*a23*a30*up(2,iq)
                 field_up(2,ip10,ip23,ip31,isub) = field_up(2,ip10,ip23,ip31,isub) + &
     &                       a10*a23*a31*up(2,iq)
                 field_up(2,ip10,ip23,ip32,isub) = field_up(2,ip10,ip23,ip32,isub) + &
     &                       a10*a23*a32*up(2,iq)
                 field_up(2,ip10,ip23,ip33,isub) = field_up(2,ip10,ip23,ip33,isub) + &
     &                       a10*a23*a33*up(2,iq)
                 field_up(2,ip11,ip20,ip30,isub) = field_up(2,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(2,iq)
                 field_up(2,ip11,ip20,ip31,isub) = field_up(2,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(2,iq)
                 field_up(2,ip11,ip20,ip32,isub) = field_up(2,ip11,ip20,ip32,isub) + &
     &                       a11*a20*a32*up(2,iq)
                 field_up(2,ip11,ip20,ip33,isub) = field_up(2,ip11,ip20,ip33,isub) + &
     &                       a11*a20*a33*up(2,iq)
                 field_up(2,ip11,ip21,ip30,isub) = field_up(2,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(2,iq)
                 field_up(2,ip11,ip21,ip31,isub) = field_up(2,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(2,iq)
                 field_up(2,ip11,ip21,ip32,isub) = field_up(2,ip11,ip21,ip32,isub) + &
     &                       a11*a21*a32*up(2,iq)
                 field_up(2,ip11,ip21,ip33,isub) = field_up(2,ip11,ip21,ip33,isub) + &
     &                       a11*a21*a33*up(2,iq)
                 field_up(2,ip11,ip22,ip30,isub) = field_up(2,ip11,ip22,ip30,isub) + &
     &                       a11*a22*a30*up(2,iq)
                 field_up(2,ip11,ip22,ip31,isub) = field_up(2,ip11,ip22,ip31,isub) + &
     &                       a11*a22*a31*up(2,iq)
                 field_up(2,ip11,ip22,ip32,isub) = field_up(2,ip11,ip22,ip32,isub) + &
     &                       a11*a22*a32*up(2,iq)
                 field_up(2,ip11,ip22,ip33,isub) = field_up(2,ip11,ip22,ip33,isub) + &
     &                       a11*a22*a33*up(2,iq)
                 field_up(2,ip11,ip23,ip30,isub) = field_up(2,ip11,ip23,ip30,isub) + &
     &                       a11*a23*a30*up(2,iq)
                 field_up(2,ip11,ip23,ip31,isub) = field_up(2,ip11,ip23,ip31,isub) + &
     &                       a11*a23*a31*up(2,iq)
                 field_up(2,ip11,ip23,ip32,isub) = field_up(2,ip11,ip23,ip32,isub) + &
     &                       a11*a23*a32*up(2,iq)
                 field_up(2,ip11,ip23,ip33,isub) = field_up(2,ip11,ip23,ip33,isub) + &
     &                       a11*a23*a33*up(2,iq)
                 field_up(2,ip12,ip20,ip30,isub) = field_up(2,ip12,ip20,ip30,isub) + &
     &                       a12*a20*a30*up(2,iq)
                 field_up(2,ip12,ip20,ip31,isub) = field_up(2,ip12,ip20,ip31,isub) + &
     &                       a12*a20*a31*up(2,iq)
                 field_up(2,ip12,ip20,ip32,isub) = field_up(2,ip12,ip20,ip32,isub) + &
     &                       a12*a20*a32*up(2,iq)
                 field_up(2,ip12,ip20,ip33,isub) = field_up(2,ip12,ip20,ip33,isub) + &
     &                       a12*a20*a33*up(2,iq)
                 field_up(2,ip12,ip21,ip30,isub) = field_up(2,ip12,ip21,ip30,isub) + &
     &                       a12*a21*a30*up(2,iq)
                 field_up(2,ip12,ip21,ip31,isub) = field_up(2,ip12,ip21,ip31,isub) + &
     &                       a12*a21*a31*up(2,iq)
                 field_up(2,ip12,ip21,ip32,isub) = field_up(2,ip12,ip21,ip32,isub) + &
     &                       a12*a21*a32*up(2,iq)
                 field_up(2,ip12,ip21,ip33,isub) = field_up(2,ip12,ip21,ip33,isub) + &
     &                       a12*a21*a33*up(2,iq)
                 field_up(2,ip12,ip22,ip30,isub) = field_up(2,ip12,ip22,ip30,isub) + &
     &                       a12*a22*a30*up(2,iq)
                 field_up(2,ip12,ip22,ip31,isub) = field_up(2,ip12,ip22,ip31,isub) + &
     &                       a12*a22*a31*up(2,iq)
                 field_up(2,ip12,ip22,ip32,isub) = field_up(2,ip12,ip22,ip32,isub) + &
     &                       a12*a22*a32*up(2,iq)
                 field_up(2,ip12,ip22,ip33,isub) = field_up(2,ip12,ip22,ip33,isub) + &
     &                       a12*a22*a33*up(2,iq)
                 field_up(2,ip12,ip23,ip30,isub) = field_up(2,ip12,ip23,ip30,isub) + &
     &                       a12*a23*a30*up(2,iq)
                 field_up(2,ip12,ip23,ip31,isub) = field_up(2,ip12,ip23,ip31,isub) + &
     &                       a12*a23*a31*up(2,iq)
                 field_up(2,ip12,ip23,ip32,isub) = field_up(2,ip12,ip23,ip32,isub) + &
     &                       a12*a23*a32*up(2,iq)
                 field_up(2,ip12,ip23,ip33,isub) = field_up(2,ip12,ip23,ip33,isub) + &
     &                       a12*a23*a33*up(2,iq)
                 field_up(2,ip13,ip20,ip30,isub) = field_up(2,ip13,ip20,ip30,isub) + &
     &                       a13*a20*a30*up(2,iq)
                 field_up(2,ip13,ip20,ip31,isub) = field_up(2,ip13,ip20,ip31,isub) + &
     &                       a13*a20*a31*up(2,iq)
                 field_up(2,ip13,ip20,ip32,isub) = field_up(2,ip13,ip20,ip32,isub) + &
     &                       a13*a20*a32*up(2,iq)
                 field_up(2,ip13,ip20,ip33,isub) = field_up(2,ip13,ip20,ip33,isub) + &
     &                       a13*a20*a33*up(2,iq)
                 field_up(2,ip13,ip21,ip30,isub) = field_up(2,ip13,ip21,ip30,isub) + &
     &                       a13*a21*a30*up(2,iq)
                 field_up(2,ip13,ip21,ip31,isub) = field_up(2,ip13,ip21,ip31,isub) + &
     &                       a13*a21*a31*up(2,iq)
                 field_up(2,ip13,ip21,ip32,isub) = field_up(2,ip13,ip21,ip32,isub) + &
     &                       a13*a21*a32*up(2,iq)
                 field_up(2,ip13,ip21,ip33,isub) = field_up(2,ip13,ip21,ip33,isub) + &
     &                       a13*a21*a33*up(2,iq)
                 field_up(2,ip13,ip22,ip30,isub) = field_up(2,ip13,ip22,ip30,isub) + &
     &                       a13*a22*a30*up(2,iq)
                 field_up(2,ip13,ip22,ip31,isub) = field_up(2,ip13,ip22,ip31,isub) + &
     &                       a13*a22*a31*up(2,iq)
                 field_up(2,ip13,ip22,ip32,isub) = field_up(2,ip13,ip22,ip32,isub) + &
     &                       a13*a22*a32*up(2,iq)
                 field_up(2,ip13,ip22,ip33,isub) = field_up(2,ip13,ip22,ip33,isub) + &
     &                       a13*a22*a33*up(2,iq)
                 field_up(2,ip13,ip23,ip30,isub) = field_up(2,ip13,ip23,ip30,isub) + &
     &                       a13*a23*a30*up(2,iq)
                 field_up(2,ip13,ip23,ip31,isub) = field_up(2,ip13,ip23,ip31,isub) + &
     &                       a13*a23*a31*up(2,iq)
                 field_up(2,ip13,ip23,ip32,isub) = field_up(2,ip13,ip23,ip32,isub) + &
     &                       a13*a23*a32*up(2,iq)
                 field_up(2,ip13,ip23,ip33,isub) = field_up(2,ip13,ip23,ip33,isub) + &
     &                       a13*a23*a33*up(2,iq)

                 field_up(3,ip10,ip20,ip30,isub) = field_up(3,ip10,ip20,ip30,isub) + &
     &                       a10*a20*a30*up(3,iq)
                 field_up(3,ip10,ip20,ip31,isub) = field_up(3,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(3,iq)
                 field_up(3,ip10,ip20,ip32,isub) = field_up(3,ip10,ip20,ip32,isub) + &
     &                       a10*a20*a32*up(3,iq)
                 field_up(3,ip10,ip20,ip33,isub) = field_up(3,ip10,ip20,ip33,isub) + &
     &                       a10*a20*a33*up(3,iq)
                 field_up(3,ip10,ip21,ip30,isub) = field_up(3,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(3,iq)
                 field_up(3,ip10,ip21,ip31,isub) = field_up(3,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(3,iq)
                 field_up(3,ip10,ip21,ip32,isub) = field_up(3,ip10,ip21,ip32,isub) + &
     &                       a10*a21*a32*up(3,iq)
                 field_up(3,ip10,ip21,ip33,isub) = field_up(3,ip10,ip21,ip33,isub) + &
     &                       a10*a21*a33*up(3,iq)
                 field_up(3,ip10,ip22,ip30,isub) = field_up(3,ip10,ip22,ip30,isub) + &
     &                       a10*a22*a30*up(3,iq)
                 field_up(3,ip10,ip22,ip31,isub) = field_up(3,ip10,ip22,ip31,isub) + &
     &                       a10*a22*a31*up(3,iq)
                 field_up(3,ip10,ip22,ip32,isub) = field_up(3,ip10,ip22,ip32,isub) + &
     &                       a10*a22*a32*up(3,iq)
                 field_up(3,ip10,ip22,ip33,isub) = field_up(3,ip10,ip22,ip33,isub) + &
     &                       a10*a22*a33*up(3,iq)
                 field_up(3,ip10,ip23,ip30,isub) = field_up(3,ip10,ip23,ip30,isub) + &
     &                       a10*a23*a30*up(3,iq)
                 field_up(3,ip10,ip23,ip31,isub) = field_up(3,ip10,ip23,ip31,isub) + &
     &                       a10*a23*a31*up(3,iq)
                 field_up(3,ip10,ip23,ip32,isub) = field_up(3,ip10,ip23,ip32,isub) + &
     &                       a10*a23*a32*up(3,iq)
                 field_up(3,ip10,ip23,ip33,isub) = field_up(3,ip10,ip23,ip33,isub) + &
     &                       a10*a23*a33*up(3,iq)
                 field_up(3,ip11,ip20,ip30,isub) = field_up(3,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(3,iq)
                 field_up(3,ip11,ip20,ip31,isub) = field_up(3,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(3,iq)
                 field_up(3,ip11,ip20,ip32,isub) = field_up(3,ip11,ip20,ip32,isub) + &
     &                       a11*a20*a32*up(3,iq)
                 field_up(3,ip11,ip20,ip33,isub) = field_up(3,ip11,ip20,ip33,isub) + &
     &                       a11*a20*a33*up(3,iq)
                 field_up(3,ip11,ip21,ip30,isub) = field_up(3,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(3,iq)
                 field_up(3,ip11,ip21,ip31,isub) = field_up(3,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(3,iq)
                 field_up(3,ip11,ip21,ip32,isub) = field_up(3,ip11,ip21,ip32,isub) + &
     &                       a11*a21*a32*up(3,iq)
                 field_up(3,ip11,ip21,ip33,isub) = field_up(3,ip11,ip21,ip33,isub) + &
     &                       a11*a21*a33*up(3,iq)
                 field_up(3,ip11,ip22,ip30,isub) = field_up(3,ip11,ip22,ip30,isub) + &
     &                       a11*a22*a30*up(3,iq)
                 field_up(3,ip11,ip22,ip31,isub) = field_up(3,ip11,ip22,ip31,isub) + &
     &                       a11*a22*a31*up(3,iq)
                 field_up(3,ip11,ip22,ip32,isub) = field_up(3,ip11,ip22,ip32,isub) + &
     &                       a11*a22*a32*up(3,iq)
                 field_up(3,ip11,ip22,ip33,isub) = field_up(3,ip11,ip22,ip33,isub) + &
     &                       a11*a22*a33*up(3,iq)
                 field_up(3,ip11,ip23,ip30,isub) = field_up(3,ip11,ip23,ip30,isub) + &
     &                       a11*a23*a30*up(3,iq)
                 field_up(3,ip11,ip23,ip31,isub) = field_up(3,ip11,ip23,ip31,isub) + &
     &                       a11*a23*a31*up(3,iq)
                 field_up(3,ip11,ip23,ip32,isub) = field_up(3,ip11,ip23,ip32,isub) + &
     &                       a11*a23*a32*up(3,iq)
                 field_up(3,ip11,ip23,ip33,isub) = field_up(3,ip11,ip23,ip33,isub) + &
     &                       a11*a23*a33*up(3,iq)
                 field_up(3,ip12,ip20,ip30,isub) = field_up(3,ip12,ip20,ip30,isub) + &
     &                       a12*a20*a30*up(3,iq)
                 field_up(3,ip12,ip20,ip31,isub) = field_up(3,ip12,ip20,ip31,isub) + &
     &                       a12*a20*a31*up(3,iq)
                 field_up(3,ip12,ip20,ip32,isub) = field_up(3,ip12,ip20,ip32,isub) + &
     &                       a12*a20*a32*up(3,iq)
                 field_up(3,ip12,ip20,ip33,isub) = field_up(3,ip12,ip20,ip33,isub) + &
     &                       a12*a20*a33*up(3,iq)
                 field_up(3,ip12,ip21,ip30,isub) = field_up(3,ip12,ip21,ip30,isub) + &
     &                       a12*a21*a30*up(3,iq)
                 field_up(3,ip12,ip21,ip31,isub) = field_up(3,ip12,ip21,ip31,isub) + &
     &                       a12*a21*a31*up(3,iq)
                 field_up(3,ip12,ip21,ip32,isub) = field_up(3,ip12,ip21,ip32,isub) + &
     &                       a12*a21*a32*up(3,iq)
                 field_up(3,ip12,ip21,ip33,isub) = field_up(3,ip12,ip21,ip33,isub) + &
     &                       a12*a21*a33*up(3,iq)
                 field_up(3,ip12,ip22,ip30,isub) = field_up(3,ip12,ip22,ip30,isub) + &
     &                       a12*a22*a30*up(3,iq)
                 field_up(3,ip12,ip22,ip31,isub) = field_up(3,ip12,ip22,ip31,isub) + &
     &                       a12*a22*a31*up(3,iq)
                 field_up(3,ip12,ip22,ip32,isub) = field_up(3,ip12,ip22,ip32,isub) + &
     &                       a12*a22*a32*up(3,iq)
                 field_up(3,ip12,ip22,ip33,isub) = field_up(3,ip12,ip22,ip33,isub) + &
     &                       a12*a22*a33*up(3,iq)
                 field_up(3,ip12,ip23,ip30,isub) = field_up(3,ip12,ip23,ip30,isub) + &
     &                       a12*a23*a30*up(3,iq)
                 field_up(3,ip12,ip23,ip31,isub) = field_up(3,ip12,ip23,ip31,isub) + &
     &                       a12*a23*a31*up(3,iq)
                 field_up(3,ip12,ip23,ip32,isub) = field_up(3,ip12,ip23,ip32,isub) + &
     &                       a12*a23*a32*up(3,iq)
                 field_up(3,ip12,ip23,ip33,isub) = field_up(3,ip12,ip23,ip33,isub) + &
     &                       a12*a23*a33*up(3,iq)
                 field_up(3,ip13,ip20,ip30,isub) = field_up(3,ip13,ip20,ip30,isub) + &
     &                       a13*a20*a30*up(3,iq)
                 field_up(3,ip13,ip20,ip31,isub) = field_up(3,ip13,ip20,ip31,isub) + &
     &                       a13*a20*a31*up(3,iq)
                 field_up(3,ip13,ip20,ip32,isub) = field_up(3,ip13,ip20,ip32,isub) + &
     &                       a13*a20*a32*up(3,iq)
                 field_up(3,ip13,ip20,ip33,isub) = field_up(3,ip13,ip20,ip33,isub) + &
     &                       a13*a20*a33*up(3,iq)
                 field_up(3,ip13,ip21,ip30,isub) = field_up(3,ip13,ip21,ip30,isub) + &
     &                       a13*a21*a30*up(3,iq)
                 field_up(3,ip13,ip21,ip31,isub) = field_up(3,ip13,ip21,ip31,isub) + &
     &                       a13*a21*a31*up(3,iq)
                 field_up(3,ip13,ip21,ip32,isub) = field_up(3,ip13,ip21,ip32,isub) + &
     &                       a13*a21*a32*up(3,iq)
                 field_up(3,ip13,ip21,ip33,isub) = field_up(3,ip13,ip21,ip33,isub) + &
     &                       a13*a21*a33*up(3,iq)
                 field_up(3,ip13,ip22,ip30,isub) = field_up(3,ip13,ip22,ip30,isub) + &
     &                       a13*a22*a30*up(3,iq)
                 field_up(3,ip13,ip22,ip31,isub) = field_up(3,ip13,ip22,ip31,isub) + &
     &                       a13*a22*a31*up(3,iq)
                 field_up(3,ip13,ip22,ip32,isub) = field_up(3,ip13,ip22,ip32,isub) + &
     &                       a13*a22*a32*up(3,iq)
                 field_up(3,ip13,ip22,ip33,isub) = field_up(3,ip13,ip22,ip33,isub) + &
     &                       a13*a22*a33*up(3,iq)
                 field_up(3,ip13,ip23,ip30,isub) = field_up(3,ip13,ip23,ip30,isub) + &
     &                       a13*a23*a30*up(3,iq)
                 field_up(3,ip13,ip23,ip31,isub) = field_up(3,ip13,ip23,ip31,isub) + &
     &                       a13*a23*a31*up(3,iq)
                 field_up(3,ip13,ip23,ip32,isub) = field_up(3,ip13,ip23,ip32,isub) + &
     &                       a13*a23*a32*up(3,iq)
                 field_up(3,ip13,ip23,ip33,isub) = field_up(3,ip13,ip23,ip33,isub) + &
     &                       a13*a23*a33*up(3,iq)

                 field_up(4,ip10,ip20,ip30,isub) = field_up(4,ip10,ip20,ip30,isub) + &
     &                       a10*a20*a30*up(4,iq)
                 field_up(4,ip10,ip20,ip31,isub) = field_up(4,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(4,iq)
                 field_up(4,ip10,ip20,ip32,isub) = field_up(4,ip10,ip20,ip32,isub) + &
     &                       a10*a20*a32*up(4,iq)
                 field_up(4,ip10,ip20,ip33,isub) = field_up(4,ip10,ip20,ip33,isub) + &
     &                       a10*a20*a33*up(4,iq)
                 field_up(4,ip10,ip21,ip30,isub) = field_up(4,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(4,iq)
                 field_up(4,ip10,ip21,ip31,isub) = field_up(4,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(4,iq)
                 field_up(4,ip10,ip21,ip32,isub) = field_up(4,ip10,ip21,ip32,isub) + &
     &                       a10*a21*a32*up(4,iq)
                 field_up(4,ip10,ip21,ip33,isub) = field_up(4,ip10,ip21,ip33,isub) + &
     &                       a10*a21*a33*up(4,iq)
                 field_up(4,ip10,ip22,ip30,isub) = field_up(4,ip10,ip22,ip30,isub) + &
     &                       a10*a22*a30*up(4,iq)
                 field_up(4,ip10,ip22,ip31,isub) = field_up(4,ip10,ip22,ip31,isub) + &
     &                       a10*a22*a31*up(4,iq)
                 field_up(4,ip10,ip22,ip32,isub) = field_up(4,ip10,ip22,ip32,isub) + &
     &                       a10*a22*a32*up(4,iq)
                 field_up(4,ip10,ip22,ip33,isub) = field_up(4,ip10,ip22,ip33,isub) + &
     &                       a10*a22*a33*up(4,iq)
                 field_up(4,ip10,ip23,ip30,isub) = field_up(4,ip10,ip23,ip30,isub) + &
     &                       a10*a23*a30*up(4,iq)
                 field_up(4,ip10,ip23,ip31,isub) = field_up(4,ip10,ip23,ip31,isub) + &
     &                       a10*a23*a31*up(4,iq)
                 field_up(4,ip10,ip23,ip32,isub) = field_up(4,ip10,ip23,ip32,isub) + &
     &                       a10*a23*a32*up(4,iq)
                 field_up(4,ip10,ip23,ip33,isub) = field_up(4,ip10,ip23,ip33,isub) + &
     &                       a10*a23*a33*up(4,iq)
                 field_up(4,ip11,ip20,ip30,isub) = field_up(4,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(4,iq)
                 field_up(4,ip11,ip20,ip31,isub) = field_up(4,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(4,iq)
                 field_up(4,ip11,ip20,ip32,isub) = field_up(4,ip11,ip20,ip32,isub) + &
     &                       a11*a20*a32*up(4,iq)
                 field_up(4,ip11,ip20,ip33,isub) = field_up(4,ip11,ip20,ip33,isub) + &
     &                       a11*a20*a33*up(4,iq)
                 field_up(4,ip11,ip21,ip30,isub) = field_up(4,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(4,iq)
                 field_up(4,ip11,ip21,ip31,isub) = field_up(4,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(4,iq)
                 field_up(4,ip11,ip21,ip32,isub) = field_up(4,ip11,ip21,ip32,isub) + &
     &                       a11*a21*a32*up(4,iq)
                 field_up(4,ip11,ip21,ip33,isub) = field_up(4,ip11,ip21,ip33,isub) + &
     &                       a11*a21*a33*up(4,iq)
                 field_up(4,ip11,ip22,ip30,isub) = field_up(4,ip11,ip22,ip30,isub) + &
     &                       a11*a22*a30*up(4,iq)
                 field_up(4,ip11,ip22,ip31,isub) = field_up(4,ip11,ip22,ip31,isub) + &
     &                       a11*a22*a31*up(4,iq)
                 field_up(4,ip11,ip22,ip32,isub) = field_up(4,ip11,ip22,ip32,isub) + &
     &                       a11*a22*a32*up(4,iq)
                 field_up(4,ip11,ip22,ip33,isub) = field_up(4,ip11,ip22,ip33,isub) + &
     &                       a11*a22*a33*up(4,iq)
                 field_up(4,ip11,ip23,ip30,isub) = field_up(4,ip11,ip23,ip30,isub) + &
     &                       a11*a23*a30*up(4,iq)
                 field_up(4,ip11,ip23,ip31,isub) = field_up(4,ip11,ip23,ip31,isub) + &
     &                       a11*a23*a31*up(4,iq)
                 field_up(4,ip11,ip23,ip32,isub) = field_up(4,ip11,ip23,ip32,isub) + &
     &                       a11*a23*a32*up(4,iq)
                 field_up(4,ip11,ip23,ip33,isub) = field_up(4,ip11,ip23,ip33,isub) + &
     &                       a11*a23*a33*up(4,iq)
                 field_up(4,ip12,ip20,ip30,isub) = field_up(4,ip12,ip20,ip30,isub) + &
     &                       a12*a20*a30*up(4,iq)
                 field_up(4,ip12,ip20,ip31,isub) = field_up(4,ip12,ip20,ip31,isub) + &
     &                       a12*a20*a31*up(4,iq)
                 field_up(4,ip12,ip20,ip32,isub) = field_up(4,ip12,ip20,ip32,isub) + &
     &                       a12*a20*a32*up(4,iq)
                 field_up(4,ip12,ip20,ip33,isub) = field_up(4,ip12,ip20,ip33,isub) + &
     &                       a12*a20*a33*up(4,iq)
                 field_up(4,ip12,ip21,ip30,isub) = field_up(4,ip12,ip21,ip30,isub) + &
     &                       a12*a21*a30*up(4,iq)
                 field_up(4,ip12,ip21,ip31,isub) = field_up(4,ip12,ip21,ip31,isub) + &
     &                       a12*a21*a31*up(4,iq)
                 field_up(4,ip12,ip21,ip32,isub) = field_up(4,ip12,ip21,ip32,isub) + &
     &                       a12*a21*a32*up(4,iq)
                 field_up(4,ip12,ip21,ip33,isub) = field_up(4,ip12,ip21,ip33,isub) + &
     &                       a12*a21*a33*up(4,iq)
                 field_up(4,ip12,ip22,ip30,isub) = field_up(4,ip12,ip22,ip30,isub) + &
     &                       a12*a22*a30*up(4,iq)
                 field_up(4,ip12,ip22,ip31,isub) = field_up(4,ip12,ip22,ip31,isub) + &
     &                       a12*a22*a31*up(4,iq)
                 field_up(4,ip12,ip22,ip32,isub) = field_up(4,ip12,ip22,ip32,isub) + &
     &                       a12*a22*a32*up(4,iq)
                 field_up(4,ip12,ip22,ip33,isub) = field_up(4,ip12,ip22,ip33,isub) + &
     &                       a12*a22*a33*up(4,iq)
                 field_up(4,ip12,ip23,ip30,isub) = field_up(4,ip12,ip23,ip30,isub) + &
     &                       a12*a23*a30*up(4,iq)
                 field_up(4,ip12,ip23,ip31,isub) = field_up(4,ip12,ip23,ip31,isub) + &
     &                       a12*a23*a31*up(4,iq)
                 field_up(4,ip12,ip23,ip32,isub) = field_up(4,ip12,ip23,ip32,isub) + &
     &                       a12*a23*a32*up(4,iq)
                 field_up(4,ip12,ip23,ip33,isub) = field_up(4,ip12,ip23,ip33,isub) + &
     &                       a12*a23*a33*up(4,iq)
                 field_up(4,ip13,ip20,ip30,isub) = field_up(4,ip13,ip20,ip30,isub) + &
     &                       a13*a20*a30*up(4,iq)
                 field_up(4,ip13,ip20,ip31,isub) = field_up(4,ip13,ip20,ip31,isub) + &
     &                       a13*a20*a31*up(4,iq)
                 field_up(4,ip13,ip20,ip32,isub) = field_up(4,ip13,ip20,ip32,isub) + &
     &                       a13*a20*a32*up(4,iq)
                 field_up(4,ip13,ip20,ip33,isub) = field_up(4,ip13,ip20,ip33,isub) + &
     &                       a13*a20*a33*up(4,iq)
                 field_up(4,ip13,ip21,ip30,isub) = field_up(4,ip13,ip21,ip30,isub) + &
     &                       a13*a21*a30*up(4,iq)
                 field_up(4,ip13,ip21,ip31,isub) = field_up(4,ip13,ip21,ip31,isub) + &
     &                       a13*a21*a31*up(4,iq)
                 field_up(4,ip13,ip21,ip32,isub) = field_up(4,ip13,ip21,ip32,isub) + &
     &                       a13*a21*a32*up(4,iq)
                 field_up(4,ip13,ip21,ip33,isub) = field_up(4,ip13,ip21,ip33,isub) + &
     &                       a13*a21*a33*up(4,iq)
                 field_up(4,ip13,ip22,ip30,isub) = field_up(4,ip13,ip22,ip30,isub) + &
     &                       a13*a22*a30*up(4,iq)
                 field_up(4,ip13,ip22,ip31,isub) = field_up(4,ip13,ip22,ip31,isub) + &
     &                       a13*a22*a31*up(4,iq)
                 field_up(4,ip13,ip22,ip32,isub) = field_up(4,ip13,ip22,ip32,isub) + &
     &                       a13*a22*a32*up(4,iq)
                 field_up(4,ip13,ip22,ip33,isub) = field_up(4,ip13,ip22,ip33,isub) + &
     &                       a13*a22*a33*up(4,iq)
                 field_up(4,ip13,ip23,ip30,isub) = field_up(4,ip13,ip23,ip30,isub) + &
     &                       a13*a23*a30*up(4,iq)
                 field_up(4,ip13,ip23,ip31,isub) = field_up(4,ip13,ip23,ip31,isub) + &
     &                       a13*a23*a31*up(4,iq)
                 field_up(4,ip13,ip23,ip32,isub) = field_up(4,ip13,ip23,ip32,isub) + &
     &                       a13*a23*a32*up(4,iq)
                 field_up(4,ip13,ip23,ip33,isub) = field_up(4,ip13,ip23,ip33,isub) + &
     &                       a13*a23*a33*up(4,iq)


              END DO
              !----------------------------------------------------------------!
              !  Unrolled versions for 3-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 3) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)

                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01)
                 ip20 = INT(x02)
                 ip30 = INT(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33


                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

                 field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(1,iq)
                 field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(2,iq)
                 field_up(3,ip10,ip20,ip30,isub)=field_up(3,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(3,iq)
                 field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(1,iq)
                 field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(2,iq)
                 field_up(3,ip11,ip20,ip30,isub)=field_up(3,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(3,iq)
                 field_up(1,ip12,ip20,ip30,isub)=field_up(1,ip12,ip20,ip30,isub)+&
     &                      a12a20a30*up(1,iq)
                 field_up(2,ip12,ip20,ip30,isub)=field_up(2,ip12,ip20,ip30,isub)+&
     &                      a12a20a30*up(2,iq)
                 field_up(3,ip12,ip20,ip30,isub)=field_up(3,ip12,ip20,ip30,isub)+&
     &                      a12a20a30*up(3,iq)
                 field_up(1,ip13,ip20,ip30,isub)=field_up(1,ip13,ip20,ip30,isub)+&
     &                      a13a20a30*up(1,iq)
                 field_up(2,ip13,ip20,ip30,isub)=field_up(2,ip13,ip20,ip30,isub)+&
     &                      a13a20a30*up(2,iq)
                 field_up(3,ip13,ip20,ip30,isub)=field_up(3,ip13,ip20,ip30,isub)+&
     &                      a13a20a30*up(3,iq)
                 field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(1,iq)
                 field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(2,iq)
                 field_up(3,ip10,ip21,ip30,isub)=field_up(3,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(3,iq)
                 field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(1,iq)
                 field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(2,iq)
                 field_up(3,ip11,ip21,ip30,isub)=field_up(3,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(3,iq)
                 field_up(1,ip12,ip21,ip30,isub)=field_up(1,ip12,ip21,ip30,isub)+&
     &                      a12a21a30*up(1,iq)
                 field_up(2,ip12,ip21,ip30,isub)=field_up(2,ip12,ip21,ip30,isub)+&
     &                      a12a21a30*up(2,iq)
                 field_up(3,ip12,ip21,ip30,isub)=field_up(3,ip12,ip21,ip30,isub)+&
     &                      a12a21a30*up(3,iq)
                 field_up(1,ip13,ip21,ip30,isub)=field_up(1,ip13,ip21,ip30,isub)+&
     &                      a13a21a30*up(1,iq)
                 field_up(2,ip13,ip21,ip30,isub)=field_up(2,ip13,ip21,ip30,isub)+&
     &                      a13a21a30*up(2,iq)
                 field_up(3,ip13,ip21,ip30,isub)=field_up(3,ip13,ip21,ip30,isub)+&
     &                      a13a21a30*up(3,iq)
                 field_up(1,ip10,ip22,ip30,isub)=field_up(1,ip10,ip22,ip30,isub)+&
     &                      a10a22a30*up(1,iq)
                 field_up(2,ip10,ip22,ip30,isub)=field_up(2,ip10,ip22,ip30,isub)+&
     &                      a10a22a30*up(2,iq)
                 field_up(3,ip10,ip22,ip30,isub)=field_up(3,ip10,ip22,ip30,isub)+&
     &                      a10a22a30*up(3,iq)
                 field_up(1,ip11,ip22,ip30,isub)=field_up(1,ip11,ip22,ip30,isub)+&
     &                      a11a22a30*up(1,iq)
                 field_up(2,ip11,ip22,ip30,isub)=field_up(2,ip11,ip22,ip30,isub)+&
     &                      a11a22a30*up(2,iq)
                 field_up(3,ip11,ip22,ip30,isub)=field_up(3,ip11,ip22,ip30,isub)+&
     &                      a11a22a30*up(3,iq)
                 field_up(1,ip12,ip22,ip30,isub)=field_up(1,ip12,ip22,ip30,isub)+&
     &                      a12a22a30*up(1,iq)
                 field_up(2,ip12,ip22,ip30,isub)=field_up(2,ip12,ip22,ip30,isub)+&
     &                      a12a22a30*up(2,iq)
                 field_up(3,ip12,ip22,ip30,isub)=field_up(3,ip12,ip22,ip30,isub)+&
     &                      a12a22a30*up(3,iq)
                 field_up(1,ip13,ip22,ip30,isub)=field_up(1,ip13,ip22,ip30,isub)+&
     &                      a13a22a30*up(1,iq)
                 field_up(2,ip13,ip22,ip30,isub)=field_up(2,ip13,ip22,ip30,isub)+&
     &                      a13a22a30*up(2,iq)
                 field_up(3,ip13,ip22,ip30,isub)=field_up(3,ip13,ip22,ip30,isub)+&
     &                      a13a22a30*up(3,iq)
                 field_up(1,ip10,ip23,ip30,isub)=field_up(1,ip10,ip23,ip30,isub)+&
     &                      a10a23a30*up(1,iq)
                 field_up(2,ip10,ip23,ip30,isub)=field_up(2,ip10,ip23,ip30,isub)+&
     &                      a10a23a30*up(2,iq)
                 field_up(3,ip10,ip23,ip30,isub)=field_up(3,ip10,ip23,ip30,isub)+&
     &                      a10a23a30*up(3,iq)
                 field_up(1,ip11,ip23,ip30,isub)=field_up(1,ip11,ip23,ip30,isub)+&
     &                      a11a23a30*up(1,iq)
                 field_up(2,ip11,ip23,ip30,isub)=field_up(2,ip11,ip23,ip30,isub)+&
     &                      a11a23a30*up(2,iq)
                 field_up(3,ip11,ip23,ip30,isub)=field_up(3,ip11,ip23,ip30,isub)+&
     &                      a11a23a30*up(3,iq)
                 field_up(1,ip12,ip23,ip30,isub)=field_up(1,ip12,ip23,ip30,isub)+&
     &                      a12a23a30*up(1,iq)
                 field_up(2,ip12,ip23,ip30,isub)=field_up(2,ip12,ip23,ip30,isub)+&
     &                      a12a23a30*up(2,iq)
                 field_up(3,ip12,ip23,ip30,isub)=field_up(3,ip12,ip23,ip30,isub)+&
     &                      a12a23a30*up(3,iq)
                 field_up(1,ip13,ip23,ip30,isub)=field_up(1,ip13,ip23,ip30,isub)+&
     &                      a13a23a30*up(1,iq)
                 field_up(2,ip13,ip23,ip30,isub)=field_up(2,ip13,ip23,ip30,isub)+&
     &                      a13a23a30*up(2,iq)
                 field_up(3,ip13,ip23,ip30,isub)=field_up(3,ip13,ip23,ip30,isub)+&
     &                      a13a23a30*up(3,iq)
                 field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(1,iq)
                 field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(2,iq)
                 field_up(3,ip10,ip20,ip31,isub)=field_up(3,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(3,iq)
                 field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(1,iq)
                 field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(2,iq)
                 field_up(3,ip11,ip20,ip31,isub)=field_up(3,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(3,iq)
                 field_up(1,ip12,ip20,ip31,isub)=field_up(1,ip12,ip20,ip31,isub)+&
     &                      a12a20a31*up(1,iq)
                 field_up(2,ip12,ip20,ip31,isub)=field_up(2,ip12,ip20,ip31,isub)+&
     &                      a12a20a31*up(2,iq)
                 field_up(3,ip12,ip20,ip31,isub)=field_up(3,ip12,ip20,ip31,isub)+&
     &                      a12a20a31*up(3,iq)
                 field_up(1,ip13,ip20,ip31,isub)=field_up(1,ip13,ip20,ip31,isub)+&
     &                      a13a20a31*up(1,iq)
                 field_up(2,ip13,ip20,ip31,isub)=field_up(2,ip13,ip20,ip31,isub)+&
     &                      a13a20a31*up(2,iq)
                 field_up(3,ip13,ip20,ip31,isub)=field_up(3,ip13,ip20,ip31,isub)+&
     &                      a13a20a31*up(3,iq)
                 field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(1,iq)
                 field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(2,iq)
                 field_up(3,ip10,ip21,ip31,isub)=field_up(3,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(3,iq)
                 field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(1,iq)
                 field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(2,iq)
                 field_up(3,ip11,ip21,ip31,isub)=field_up(3,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(3,iq)
                 field_up(1,ip12,ip21,ip31,isub)=field_up(1,ip12,ip21,ip31,isub)+&
     &                      a12a21a31*up(1,iq)
                 field_up(2,ip12,ip21,ip31,isub)=field_up(2,ip12,ip21,ip31,isub)+&
     &                      a12a21a31*up(2,iq)
                 field_up(3,ip12,ip21,ip31,isub)=field_up(3,ip12,ip21,ip31,isub)+&
     &                      a12a21a31*up(3,iq)
                 field_up(1,ip13,ip21,ip31,isub)=field_up(1,ip13,ip21,ip31,isub)+&
     &                      a13a21a31*up(1,iq)
                 field_up(2,ip13,ip21,ip31,isub)=field_up(2,ip13,ip21,ip31,isub)+&
     &                      a13a21a31*up(2,iq)
                 field_up(3,ip13,ip21,ip31,isub)=field_up(3,ip13,ip21,ip31,isub)+&
     &                      a13a21a31*up(3,iq)
                 field_up(1,ip10,ip22,ip31,isub)=field_up(1,ip10,ip22,ip31,isub)+&
     &                      a10a22a31*up(1,iq)
                 field_up(2,ip10,ip22,ip31,isub)=field_up(2,ip10,ip22,ip31,isub)+&
     &                      a10a22a31*up(2,iq)
                 field_up(3,ip10,ip22,ip31,isub)=field_up(3,ip10,ip22,ip31,isub)+&
     &                      a10a22a31*up(3,iq)
                 field_up(1,ip11,ip22,ip31,isub)=field_up(1,ip11,ip22,ip31,isub)+&
     &                      a11a22a31*up(1,iq)
                 field_up(2,ip11,ip22,ip31,isub)=field_up(2,ip11,ip22,ip31,isub)+&
     &                      a11a22a31*up(2,iq)
                 field_up(3,ip11,ip22,ip31,isub)=field_up(3,ip11,ip22,ip31,isub)+&
     &                      a11a22a31*up(3,iq)
                 field_up(1,ip12,ip22,ip31,isub)=field_up(1,ip12,ip22,ip31,isub)+&
     &                      a12a22a31*up(1,iq)
                 field_up(2,ip12,ip22,ip31,isub)=field_up(2,ip12,ip22,ip31,isub)+&
     &                      a12a22a31*up(2,iq)
                 field_up(3,ip12,ip22,ip31,isub)=field_up(3,ip12,ip22,ip31,isub)+&
     &                      a12a22a31*up(3,iq)
                 field_up(1,ip13,ip22,ip31,isub)=field_up(1,ip13,ip22,ip31,isub)+&
     &                      a13a22a31*up(1,iq)
                 field_up(2,ip13,ip22,ip31,isub)=field_up(2,ip13,ip22,ip31,isub)+&
     &                      a13a22a31*up(2,iq)
                 field_up(3,ip13,ip22,ip31,isub)=field_up(3,ip13,ip22,ip31,isub)+&
     &                      a13a22a31*up(3,iq)
                 field_up(1,ip10,ip23,ip31,isub)=field_up(1,ip10,ip23,ip31,isub)+&
     &                      a10a23a31*up(1,iq)
                 field_up(2,ip10,ip23,ip31,isub)=field_up(2,ip10,ip23,ip31,isub)+&
     &                      a10a23a31*up(2,iq)
                 field_up(3,ip10,ip23,ip31,isub)=field_up(3,ip10,ip23,ip31,isub)+&
     &                      a10a23a31*up(3,iq)
                 field_up(1,ip11,ip23,ip31,isub)=field_up(1,ip11,ip23,ip31,isub)+&
     &                      a11a23a31*up(1,iq)
                 field_up(2,ip11,ip23,ip31,isub)=field_up(2,ip11,ip23,ip31,isub)+&
     &                      a11a23a31*up(2,iq)
                 field_up(3,ip11,ip23,ip31,isub)=field_up(3,ip11,ip23,ip31,isub)+&
     &                      a11a23a31*up(3,iq)
                 field_up(1,ip12,ip23,ip31,isub)=field_up(1,ip12,ip23,ip31,isub)+&
     &                      a12a23a31*up(1,iq)
                 field_up(2,ip12,ip23,ip31,isub)=field_up(2,ip12,ip23,ip31,isub)+&
     &                      a12a23a31*up(2,iq)
                 field_up(3,ip12,ip23,ip31,isub)=field_up(3,ip12,ip23,ip31,isub)+&
     &                      a12a23a31*up(3,iq)
                 field_up(1,ip13,ip23,ip31,isub)=field_up(1,ip13,ip23,ip31,isub)+&
     &                      a13a23a31*up(1,iq)
                 field_up(2,ip13,ip23,ip31,isub)=field_up(2,ip13,ip23,ip31,isub)+&
     &                      a13a23a31*up(2,iq)
                 field_up(3,ip13,ip23,ip31,isub)=field_up(3,ip13,ip23,ip31,isub)+&
     &                      a13a23a31*up(3,iq)
                 field_up(1,ip10,ip20,ip32,isub)=field_up(1,ip10,ip20,ip32,isub)+&
     &                      a10a20a32*up(1,iq)
                 field_up(2,ip10,ip20,ip32,isub)=field_up(2,ip10,ip20,ip32,isub)+&
     &                      a10a20a32*up(2,iq)
                 field_up(3,ip10,ip20,ip32,isub)=field_up(3,ip10,ip20,ip32,isub)+&
     &                      a10a20a32*up(3,iq)
                 field_up(1,ip11,ip20,ip32,isub)=field_up(1,ip11,ip20,ip32,isub)+&
     &                      a11a20a32*up(1,iq)
                 field_up(2,ip11,ip20,ip32,isub)=field_up(2,ip11,ip20,ip32,isub)+&
     &                      a11a20a32*up(2,iq)
                 field_up(3,ip11,ip20,ip32,isub)=field_up(3,ip11,ip20,ip32,isub)+&
     &                      a11a20a32*up(3,iq)
                 field_up(1,ip12,ip20,ip32,isub)=field_up(1,ip12,ip20,ip32,isub)+&
     &                      a12a20a32*up(1,iq)
                 field_up(2,ip12,ip20,ip32,isub)=field_up(2,ip12,ip20,ip32,isub)+&
     &                      a12a20a32*up(2,iq)
                 field_up(3,ip12,ip20,ip32,isub)=field_up(3,ip12,ip20,ip32,isub)+&
     &                      a12a20a32*up(3,iq)
                 field_up(1,ip13,ip20,ip32,isub)=field_up(1,ip13,ip20,ip32,isub)+&
     &                      a13a20a32*up(1,iq)
                 field_up(2,ip13,ip20,ip32,isub)=field_up(2,ip13,ip20,ip32,isub)+&
     &                      a13a20a32*up(2,iq)
                 field_up(3,ip13,ip20,ip32,isub)=field_up(3,ip13,ip20,ip32,isub)+&
     &                      a13a20a32*up(3,iq)
                 field_up(1,ip10,ip21,ip32,isub)=field_up(1,ip10,ip21,ip32,isub)+&
     &                      a10a21a32*up(1,iq)
                 field_up(2,ip10,ip21,ip32,isub)=field_up(2,ip10,ip21,ip32,isub)+&
     &                      a10a21a32*up(2,iq)
                 field_up(3,ip10,ip21,ip32,isub)=field_up(3,ip10,ip21,ip32,isub)+&
     &                      a10a21a32*up(3,iq)
                 field_up(1,ip11,ip21,ip32,isub)=field_up(1,ip11,ip21,ip32,isub)+&
     &                      a11a21a32*up(1,iq)
                 field_up(2,ip11,ip21,ip32,isub)=field_up(2,ip11,ip21,ip32,isub)+&
     &                      a11a21a32*up(2,iq)
                 field_up(3,ip11,ip21,ip32,isub)=field_up(3,ip11,ip21,ip32,isub)+&
     &                      a11a21a32*up(3,iq)
                 field_up(1,ip12,ip21,ip32,isub)=field_up(1,ip12,ip21,ip32,isub)+&
     &                      a12a21a32*up(1,iq)
                 field_up(2,ip12,ip21,ip32,isub)=field_up(2,ip12,ip21,ip32,isub)+&
     &                      a12a21a32*up(2,iq)
                 field_up(3,ip12,ip21,ip32,isub)=field_up(3,ip12,ip21,ip32,isub)+&
     &                      a12a21a32*up(3,iq)
                 field_up(1,ip13,ip21,ip32,isub)=field_up(1,ip13,ip21,ip32,isub)+&
     &                      a13a21a32*up(1,iq)
                 field_up(2,ip13,ip21,ip32,isub)=field_up(2,ip13,ip21,ip32,isub)+&
     &                      a13a21a32*up(2,iq)
                 field_up(3,ip13,ip21,ip32,isub)=field_up(3,ip13,ip21,ip32,isub)+&
     &                      a13a21a32*up(3,iq)
                 field_up(1,ip10,ip22,ip32,isub)=field_up(1,ip10,ip22,ip32,isub)+&
     &                      a10a22a32*up(1,iq)
                 field_up(2,ip10,ip22,ip32,isub)=field_up(2,ip10,ip22,ip32,isub)+&
     &                      a10a22a32*up(2,iq)
                 field_up(3,ip10,ip22,ip32,isub)=field_up(3,ip10,ip22,ip32,isub)+&
     &                      a10a22a32*up(3,iq)
                 field_up(1,ip11,ip22,ip32,isub)=field_up(1,ip11,ip22,ip32,isub)+&
     &                      a11a22a32*up(1,iq)
                 field_up(2,ip11,ip22,ip32,isub)=field_up(2,ip11,ip22,ip32,isub)+&
     &                      a11a22a32*up(2,iq)
                 field_up(3,ip11,ip22,ip32,isub)=field_up(3,ip11,ip22,ip32,isub)+&
     &                      a11a22a32*up(3,iq)
                 field_up(1,ip12,ip22,ip32,isub)=field_up(1,ip12,ip22,ip32,isub)+&
     &                      a12a22a32*up(1,iq)
                 field_up(2,ip12,ip22,ip32,isub)=field_up(2,ip12,ip22,ip32,isub)+&
     &                      a12a22a32*up(2,iq)
                 field_up(3,ip12,ip22,ip32,isub)=field_up(3,ip12,ip22,ip32,isub)+&
     &                      a12a22a32*up(3,iq)
                 field_up(1,ip13,ip22,ip32,isub)=field_up(1,ip13,ip22,ip32,isub)+&
     &                      a13a22a32*up(1,iq)
                 field_up(2,ip13,ip22,ip32,isub)=field_up(2,ip13,ip22,ip32,isub)+&
     &                      a13a22a32*up(2,iq)
                 field_up(3,ip13,ip22,ip32,isub)=field_up(3,ip13,ip22,ip32,isub)+&
     &                      a13a22a32*up(3,iq)
                 field_up(1,ip10,ip23,ip32,isub)=field_up(1,ip10,ip23,ip32,isub)+&
     &                      a10a23a32*up(1,iq)
                 field_up(2,ip10,ip23,ip32,isub)=field_up(2,ip10,ip23,ip32,isub)+&
     &                      a10a23a32*up(2,iq)
                 field_up(3,ip10,ip23,ip32,isub)=field_up(3,ip10,ip23,ip32,isub)+&
     &                      a10a23a32*up(3,iq)
                 field_up(1,ip11,ip23,ip32,isub)=field_up(1,ip11,ip23,ip32,isub)+&
     &                      a11a23a32*up(1,iq)
                 field_up(2,ip11,ip23,ip32,isub)=field_up(2,ip11,ip23,ip32,isub)+&
     &                      a11a23a32*up(2,iq)
                 field_up(3,ip11,ip23,ip32,isub)=field_up(3,ip11,ip23,ip32,isub)+&
     &                      a11a23a32*up(3,iq)
                 field_up(1,ip12,ip23,ip32,isub)=field_up(1,ip12,ip23,ip32,isub)+&
     &                      a12a23a32*up(1,iq)
                 field_up(2,ip12,ip23,ip32,isub)=field_up(2,ip12,ip23,ip32,isub)+&
     &                      a12a23a32*up(2,iq)
                 field_up(3,ip12,ip23,ip32,isub)=field_up(3,ip12,ip23,ip32,isub)+&
     &                      a12a23a32*up(3,iq)
                 field_up(1,ip13,ip23,ip32,isub)=field_up(1,ip13,ip23,ip32,isub)+&
     &                      a13a23a32*up(1,iq)
                 field_up(2,ip13,ip23,ip32,isub)=field_up(2,ip13,ip23,ip32,isub)+&
     &                      a13a23a32*up(2,iq)
                 field_up(3,ip13,ip23,ip32,isub)=field_up(3,ip13,ip23,ip32,isub)+&
     &                      a13a23a32*up(3,iq)
                 field_up(1,ip10,ip20,ip33,isub)=field_up(1,ip10,ip20,ip33,isub)+&
     &                      a10a20a33*up(1,iq)
                 field_up(2,ip10,ip20,ip33,isub)=field_up(2,ip10,ip20,ip33,isub)+&
     &                      a10a20a33*up(2,iq)
                 field_up(3,ip10,ip20,ip33,isub)=field_up(3,ip10,ip20,ip33,isub)+&
     &                      a10a20a33*up(3,iq)
                 field_up(1,ip11,ip20,ip33,isub)=field_up(1,ip11,ip20,ip33,isub)+&
     &                      a11a20a33*up(1,iq)
                 field_up(2,ip11,ip20,ip33,isub)=field_up(2,ip11,ip20,ip33,isub)+&
     &                      a11a20a33*up(2,iq)
                 field_up(3,ip11,ip20,ip33,isub)=field_up(3,ip11,ip20,ip33,isub)+&
     &                      a11a20a33*up(3,iq)
                 field_up(1,ip12,ip20,ip33,isub)=field_up(1,ip12,ip20,ip33,isub)+&
     &                      a12a20a33*up(1,iq)
                 field_up(2,ip12,ip20,ip33,isub)=field_up(2,ip12,ip20,ip33,isub)+&
     &                      a12a20a33*up(2,iq)
                 field_up(3,ip12,ip20,ip33,isub)=field_up(3,ip12,ip20,ip33,isub)+&
     &                      a12a20a33*up(3,iq)
                 field_up(1,ip13,ip20,ip33,isub)=field_up(1,ip13,ip20,ip33,isub)+&
     &                      a13a20a33*up(1,iq)
                 field_up(2,ip13,ip20,ip33,isub)=field_up(2,ip13,ip20,ip33,isub)+&
     &                      a13a20a33*up(2,iq)
                 field_up(3,ip13,ip20,ip33,isub)=field_up(3,ip13,ip20,ip33,isub)+&
     &                      a13a20a33*up(3,iq)
                 field_up(1,ip10,ip21,ip33,isub)=field_up(1,ip10,ip21,ip33,isub)+&
     &                      a10a21a33*up(1,iq)
                 field_up(2,ip10,ip21,ip33,isub)=field_up(2,ip10,ip21,ip33,isub)+&
     &                      a10a21a33*up(2,iq)
                 field_up(3,ip10,ip21,ip33,isub)=field_up(3,ip10,ip21,ip33,isub)+&
     &                      a10a21a33*up(3,iq)
                 field_up(1,ip11,ip21,ip33,isub)=field_up(1,ip11,ip21,ip33,isub)+&
     &                      a11a21a33*up(1,iq)
                 field_up(2,ip11,ip21,ip33,isub)=field_up(2,ip11,ip21,ip33,isub)+&
     &                      a11a21a33*up(2,iq)
                 field_up(3,ip11,ip21,ip33,isub)=field_up(3,ip11,ip21,ip33,isub)+&
     &                      a11a21a33*up(3,iq)
                 field_up(1,ip12,ip21,ip33,isub)=field_up(1,ip12,ip21,ip33,isub)+&
     &                      a12a21a33*up(1,iq)
                 field_up(2,ip12,ip21,ip33,isub)=field_up(2,ip12,ip21,ip33,isub)+&
     &                      a12a21a33*up(2,iq)
                 field_up(3,ip12,ip21,ip33,isub)=field_up(3,ip12,ip21,ip33,isub)+&
     &                      a12a21a33*up(3,iq)
                 field_up(1,ip13,ip21,ip33,isub)=field_up(1,ip13,ip21,ip33,isub)+&
     &                      a13a21a33*up(1,iq)
                 field_up(2,ip13,ip21,ip33,isub)=field_up(2,ip13,ip21,ip33,isub)+&
     &                      a13a21a33*up(2,iq)
                 field_up(3,ip13,ip21,ip33,isub)=field_up(3,ip13,ip21,ip33,isub)+&
     &                      a13a21a33*up(3,iq)
                 field_up(1,ip10,ip22,ip33,isub)=field_up(1,ip10,ip22,ip33,isub)+&
     &                      a10a22a33*up(1,iq)
                 field_up(2,ip10,ip22,ip33,isub)=field_up(2,ip10,ip22,ip33,isub)+&
     &                      a10a22a33*up(2,iq)
                 field_up(3,ip10,ip22,ip33,isub)=field_up(3,ip10,ip22,ip33,isub)+&
     &                      a10a22a33*up(3,iq)
                 field_up(1,ip11,ip22,ip33,isub)=field_up(1,ip11,ip22,ip33,isub)+&
     &                      a11a22a33*up(1,iq)
                 field_up(2,ip11,ip22,ip33,isub)=field_up(2,ip11,ip22,ip33,isub)+&
     &                      a11a22a33*up(2,iq)
                 field_up(3,ip11,ip22,ip33,isub)=field_up(3,ip11,ip22,ip33,isub)+&
     &                      a11a22a33*up(3,iq)
                 field_up(1,ip12,ip22,ip33,isub)=field_up(1,ip12,ip22,ip33,isub)+&
     &                      a12a22a33*up(1,iq)
                 field_up(2,ip12,ip22,ip33,isub)=field_up(2,ip12,ip22,ip33,isub)+&
     &                      a12a22a33*up(2,iq)
                 field_up(3,ip12,ip22,ip33,isub)=field_up(3,ip12,ip22,ip33,isub)+&
     &                      a12a22a33*up(3,iq)
                 field_up(1,ip13,ip22,ip33,isub)=field_up(1,ip13,ip22,ip33,isub)+&
     &                      a13a22a33*up(1,iq)
                 field_up(2,ip13,ip22,ip33,isub)=field_up(2,ip13,ip22,ip33,isub)+&
     &                      a13a22a33*up(2,iq)
                 field_up(3,ip13,ip22,ip33,isub)=field_up(3,ip13,ip22,ip33,isub)+&
     &                      a13a22a33*up(3,iq)
                 field_up(1,ip10,ip23,ip33,isub)=field_up(1,ip10,ip23,ip33,isub)+&
     &                      a10a23a33*up(1,iq)
                 field_up(2,ip10,ip23,ip33,isub)=field_up(2,ip10,ip23,ip33,isub)+&
     &                      a10a23a33*up(2,iq)
                 field_up(3,ip10,ip23,ip33,isub)=field_up(3,ip10,ip23,ip33,isub)+&
     &                      a10a23a33*up(3,iq)
                 field_up(1,ip11,ip23,ip33,isub)=field_up(1,ip11,ip23,ip33,isub)+&
     &                      a11a23a33*up(1,iq)
                 field_up(2,ip11,ip23,ip33,isub)=field_up(2,ip11,ip23,ip33,isub)+&
     &                      a11a23a33*up(2,iq)
                 field_up(3,ip11,ip23,ip33,isub)=field_up(3,ip11,ip23,ip33,isub)+&
     &                      a11a23a33*up(3,iq)
                 field_up(1,ip12,ip23,ip33,isub)=field_up(1,ip12,ip23,ip33,isub)+&
     &                      a12a23a33*up(1,iq)
                 field_up(2,ip12,ip23,ip33,isub)=field_up(2,ip12,ip23,ip33,isub)+&
     &                      a12a23a33*up(2,iq)
                 field_up(3,ip12,ip23,ip33,isub)=field_up(3,ip12,ip23,ip33,isub)+&
     &                      a12a23a33*up(3,iq)
                 field_up(1,ip13,ip23,ip33,isub)=field_up(1,ip13,ip23,ip33,isub)+&
     &                      a13a23a33*up(1,iq)
                 field_up(2,ip13,ip23,ip33,isub)=field_up(2,ip13,ip23,ip33,isub)+&
     &                      a13a23a33*up(2,iq)
                 field_up(3,ip13,ip23,ip33,isub)=field_up(3,ip13,ip23,ip33,isub)+&
     &                      a13a23a33*up(3,iq)
              END DO
              !----------------------------------------------------------------!
              !  Unrolled versions for 2-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 2) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01)
                 ip20 = INT(x02)
                 ip30 = INT(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

                 field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(1,iq)
                 field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(2,iq)
                 field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(1,iq)
                 field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(2,iq)
                 field_up(1,ip12,ip20,ip30,isub)=field_up(1,ip12,ip20,ip30,isub)+&
     &                      a12a20a30*up(1,iq)
                 field_up(2,ip12,ip20,ip30,isub)=field_up(2,ip12,ip20,ip30,isub)+&
     &                      a12a20a30*up(2,iq)
                 field_up(1,ip13,ip20,ip30,isub)=field_up(1,ip13,ip20,ip30,isub)+&
     &                      a13a20a30*up(1,iq)
                 field_up(2,ip13,ip20,ip30,isub)=field_up(2,ip13,ip20,ip30,isub)+&
     &                      a13a20a30*up(2,iq)
                 field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(1,iq)
                 field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(2,iq)
                 field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(1,iq)
                 field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(2,iq)
                 field_up(1,ip12,ip21,ip30,isub)=field_up(1,ip12,ip21,ip30,isub)+&
     &                      a12a21a30*up(1,iq)
                 field_up(2,ip12,ip21,ip30,isub)=field_up(2,ip12,ip21,ip30,isub)+&
     &                      a12a21a30*up(2,iq)
                 field_up(1,ip13,ip21,ip30,isub)=field_up(1,ip13,ip21,ip30,isub)+&
     &                      a13a21a30*up(1,iq)
                 field_up(2,ip13,ip21,ip30,isub)=field_up(2,ip13,ip21,ip30,isub)+&
     &                      a13a21a30*up(2,iq)
                 field_up(1,ip10,ip22,ip30,isub)=field_up(1,ip10,ip22,ip30,isub)+&
     &                      a10a22a30*up(1,iq)
                 field_up(2,ip10,ip22,ip30,isub)=field_up(2,ip10,ip22,ip30,isub)+&
     &                      a10a22a30*up(2,iq)
                 field_up(1,ip11,ip22,ip30,isub)=field_up(1,ip11,ip22,ip30,isub)+&
     &                      a11a22a30*up(1,iq)
                 field_up(2,ip11,ip22,ip30,isub)=field_up(2,ip11,ip22,ip30,isub)+&
     &                      a11a22a30*up(2,iq)
                 field_up(1,ip12,ip22,ip30,isub)=field_up(1,ip12,ip22,ip30,isub)+&
     &                      a12a22a30*up(1,iq)
                 field_up(2,ip12,ip22,ip30,isub)=field_up(2,ip12,ip22,ip30,isub)+&
     &                      a12a22a30*up(2,iq)
                 field_up(1,ip13,ip22,ip30,isub)=field_up(1,ip13,ip22,ip30,isub)+&
     &                      a13a22a30*up(1,iq)
                 field_up(2,ip13,ip22,ip30,isub)=field_up(2,ip13,ip22,ip30,isub)+&
     &                      a13a22a30*up(2,iq)
                 field_up(1,ip10,ip23,ip30,isub)=field_up(1,ip10,ip23,ip30,isub)+&
     &                      a10a23a30*up(1,iq)
                 field_up(2,ip10,ip23,ip30,isub)=field_up(2,ip10,ip23,ip30,isub)+&
     &                      a10a23a30*up(2,iq)
                 field_up(1,ip11,ip23,ip30,isub)=field_up(1,ip11,ip23,ip30,isub)+&
     &                      a11a23a30*up(1,iq)
                 field_up(2,ip11,ip23,ip30,isub)=field_up(2,ip11,ip23,ip30,isub)+&
     &                      a11a23a30*up(2,iq)
                 field_up(1,ip12,ip23,ip30,isub)=field_up(1,ip12,ip23,ip30,isub)+&
     &                      a12a23a30*up(1,iq)
                 field_up(2,ip12,ip23,ip30,isub)=field_up(2,ip12,ip23,ip30,isub)+&
     &                      a12a23a30*up(2,iq)
                 field_up(1,ip13,ip23,ip30,isub)=field_up(1,ip13,ip23,ip30,isub)+&
     &                      a13a23a30*up(1,iq)
                 field_up(2,ip13,ip23,ip30,isub)=field_up(2,ip13,ip23,ip30,isub)+&
     &                      a13a23a30*up(2,iq)
                 field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(1,iq)
                 field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(2,iq)
                 field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(1,iq)
                 field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(2,iq)
                 field_up(1,ip12,ip20,ip31,isub)=field_up(1,ip12,ip20,ip31,isub)+&
     &                      a12a20a31*up(1,iq)
                 field_up(2,ip12,ip20,ip31,isub)=field_up(2,ip12,ip20,ip31,isub)+&
     &                      a12a20a31*up(2,iq)
                 field_up(1,ip13,ip20,ip31,isub)=field_up(1,ip13,ip20,ip31,isub)+&
     &                      a13a20a31*up(1,iq)
                 field_up(2,ip13,ip20,ip31,isub)=field_up(2,ip13,ip20,ip31,isub)+&
     &                      a13a20a31*up(2,iq)
                 field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(1,iq)
                 field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(2,iq)
                 field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(1,iq)
                 field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(2,iq)
                 field_up(1,ip12,ip21,ip31,isub)=field_up(1,ip12,ip21,ip31,isub)+&
     &                      a12a21a31*up(1,iq)
                 field_up(2,ip12,ip21,ip31,isub)=field_up(2,ip12,ip21,ip31,isub)+&
     &                      a12a21a31*up(2,iq)
                 field_up(1,ip13,ip21,ip31,isub)=field_up(1,ip13,ip21,ip31,isub)+&
     &                      a13a21a31*up(1,iq)
                 field_up(2,ip13,ip21,ip31,isub)=field_up(2,ip13,ip21,ip31,isub)+&
     &                      a13a21a31*up(2,iq)
                 field_up(1,ip10,ip22,ip31,isub)=field_up(1,ip10,ip22,ip31,isub)+&
     &                      a10a22a31*up(1,iq)
                 field_up(2,ip10,ip22,ip31,isub)=field_up(2,ip10,ip22,ip31,isub)+&
     &                      a10a22a31*up(2,iq)
                 field_up(1,ip11,ip22,ip31,isub)=field_up(1,ip11,ip22,ip31,isub)+&
     &                      a11a22a31*up(1,iq)
                 field_up(2,ip11,ip22,ip31,isub)=field_up(2,ip11,ip22,ip31,isub)+&
     &                      a11a22a31*up(2,iq)
                 field_up(1,ip12,ip22,ip31,isub)=field_up(1,ip12,ip22,ip31,isub)+&
     &                      a12a22a31*up(1,iq)
                 field_up(2,ip12,ip22,ip31,isub)=field_up(2,ip12,ip22,ip31,isub)+&
     &                      a12a22a31*up(2,iq)
                 field_up(1,ip13,ip22,ip31,isub)=field_up(1,ip13,ip22,ip31,isub)+&
     &                      a13a22a31*up(1,iq)
                 field_up(2,ip13,ip22,ip31,isub)=field_up(2,ip13,ip22,ip31,isub)+&
     &                      a13a22a31*up(2,iq)
                 field_up(1,ip10,ip23,ip31,isub)=field_up(1,ip10,ip23,ip31,isub)+&
     &                      a10a23a31*up(1,iq)
                 field_up(2,ip10,ip23,ip31,isub)=field_up(2,ip10,ip23,ip31,isub)+&
     &                      a10a23a31*up(2,iq)
                 field_up(1,ip11,ip23,ip31,isub)=field_up(1,ip11,ip23,ip31,isub)+&
     &                      a11a23a31*up(1,iq)
                 field_up(2,ip11,ip23,ip31,isub)=field_up(2,ip11,ip23,ip31,isub)+&
     &                      a11a23a31*up(2,iq)
                 field_up(1,ip12,ip23,ip31,isub)=field_up(1,ip12,ip23,ip31,isub)+&
     &                      a12a23a31*up(1,iq)
                 field_up(2,ip12,ip23,ip31,isub)=field_up(2,ip12,ip23,ip31,isub)+&
     &                      a12a23a31*up(2,iq)
                 field_up(1,ip13,ip23,ip31,isub)=field_up(1,ip13,ip23,ip31,isub)+&
     &                      a13a23a31*up(1,iq)
                 field_up(2,ip13,ip23,ip31,isub)=field_up(2,ip13,ip23,ip31,isub)+&
     &                      a13a23a31*up(2,iq)
                 field_up(1,ip10,ip20,ip32,isub)=field_up(1,ip10,ip20,ip32,isub)+&
     &                      a10a20a32*up(1,iq)
                 field_up(2,ip10,ip20,ip32,isub)=field_up(2,ip10,ip20,ip32,isub)+&
     &                      a10a20a32*up(2,iq)
                 field_up(1,ip11,ip20,ip32,isub)=field_up(1,ip11,ip20,ip32,isub)+&
     &                      a11a20a32*up(1,iq)
                 field_up(2,ip11,ip20,ip32,isub)=field_up(2,ip11,ip20,ip32,isub)+&
     &                      a11a20a32*up(2,iq)
                 field_up(1,ip12,ip20,ip32,isub)=field_up(1,ip12,ip20,ip32,isub)+&
     &                      a12a20a32*up(1,iq)
                 field_up(2,ip12,ip20,ip32,isub)=field_up(2,ip12,ip20,ip32,isub)+&
     &                      a12a20a32*up(2,iq)
                 field_up(1,ip13,ip20,ip32,isub)=field_up(1,ip13,ip20,ip32,isub)+&
     &                      a13a20a32*up(1,iq)
                 field_up(2,ip13,ip20,ip32,isub)=field_up(2,ip13,ip20,ip32,isub)+&
     &                      a13a20a32*up(2,iq)
                 field_up(1,ip10,ip21,ip32,isub)=field_up(1,ip10,ip21,ip32,isub)+&
     &                      a10a21a32*up(1,iq)
                 field_up(2,ip10,ip21,ip32,isub)=field_up(2,ip10,ip21,ip32,isub)+&
     &                      a10a21a32*up(2,iq)
                 field_up(1,ip11,ip21,ip32,isub)=field_up(1,ip11,ip21,ip32,isub)+&
     &                      a11a21a32*up(1,iq)
                 field_up(2,ip11,ip21,ip32,isub)=field_up(2,ip11,ip21,ip32,isub)+&
     &                      a11a21a32*up(2,iq)
                 field_up(1,ip12,ip21,ip32,isub)=field_up(1,ip12,ip21,ip32,isub)+&
     &                      a12a21a32*up(1,iq)
                 field_up(2,ip12,ip21,ip32,isub)=field_up(2,ip12,ip21,ip32,isub)+&
     &                      a12a21a32*up(2,iq)
                 field_up(1,ip13,ip21,ip32,isub)=field_up(1,ip13,ip21,ip32,isub)+&
     &                      a13a21a32*up(1,iq)
                 field_up(2,ip13,ip21,ip32,isub)=field_up(2,ip13,ip21,ip32,isub)+&
     &                      a13a21a32*up(2,iq)
                 field_up(1,ip10,ip22,ip32,isub)=field_up(1,ip10,ip22,ip32,isub)+&
     &                      a10a22a32*up(1,iq)
                 field_up(2,ip10,ip22,ip32,isub)=field_up(2,ip10,ip22,ip32,isub)+&
     &                      a10a22a32*up(2,iq)
                 field_up(1,ip11,ip22,ip32,isub)=field_up(1,ip11,ip22,ip32,isub)+&
     &                      a11a22a32*up(1,iq)
                 field_up(2,ip11,ip22,ip32,isub)=field_up(2,ip11,ip22,ip32,isub)+&
     &                      a11a22a32*up(2,iq)
                 field_up(1,ip12,ip22,ip32,isub)=field_up(1,ip12,ip22,ip32,isub)+&
     &                      a12a22a32*up(1,iq)
                 field_up(2,ip12,ip22,ip32,isub)=field_up(2,ip12,ip22,ip32,isub)+&
     &                      a12a22a32*up(2,iq)
                 field_up(1,ip13,ip22,ip32,isub)=field_up(1,ip13,ip22,ip32,isub)+&
     &                      a13a22a32*up(1,iq)
                 field_up(2,ip13,ip22,ip32,isub)=field_up(2,ip13,ip22,ip32,isub)+&
     &                      a13a22a32*up(2,iq)
                 field_up(1,ip10,ip23,ip32,isub)=field_up(1,ip10,ip23,ip32,isub)+&
     &                      a10a23a32*up(1,iq)
                 field_up(2,ip10,ip23,ip32,isub)=field_up(2,ip10,ip23,ip32,isub)+&
     &                      a10a23a32*up(2,iq)
                 field_up(1,ip11,ip23,ip32,isub)=field_up(1,ip11,ip23,ip32,isub)+&
     &                      a11a23a32*up(1,iq)
                 field_up(2,ip11,ip23,ip32,isub)=field_up(2,ip11,ip23,ip32,isub)+&
     &                      a11a23a32*up(2,iq)
                 field_up(1,ip12,ip23,ip32,isub)=field_up(1,ip12,ip23,ip32,isub)+&
     &                      a12a23a32*up(1,iq)
                 field_up(2,ip12,ip23,ip32,isub)=field_up(2,ip12,ip23,ip32,isub)+&
     &                      a12a23a32*up(2,iq)
                 field_up(1,ip13,ip23,ip32,isub)=field_up(1,ip13,ip23,ip32,isub)+&
     &                      a13a23a32*up(1,iq)
                 field_up(2,ip13,ip23,ip32,isub)=field_up(2,ip13,ip23,ip32,isub)+&
     &                      a13a23a32*up(2,iq)
                 field_up(1,ip10,ip20,ip33,isub)=field_up(1,ip10,ip20,ip33,isub)+&
     &                      a10a20a33*up(1,iq)
                 field_up(2,ip10,ip20,ip33,isub)=field_up(2,ip10,ip20,ip33,isub)+&
     &                      a10a20a33*up(2,iq)
                 field_up(1,ip11,ip20,ip33,isub)=field_up(1,ip11,ip20,ip33,isub)+&
     &                      a11a20a33*up(1,iq)
                 field_up(2,ip11,ip20,ip33,isub)=field_up(2,ip11,ip20,ip33,isub)+&
     &                      a11a20a33*up(2,iq)
                 field_up(1,ip12,ip20,ip33,isub)=field_up(1,ip12,ip20,ip33,isub)+&
     &                      a12a20a33*up(1,iq)
                 field_up(2,ip12,ip20,ip33,isub)=field_up(2,ip12,ip20,ip33,isub)+&
     &                      a12a20a33*up(2,iq)
                 field_up(1,ip13,ip20,ip33,isub)=field_up(1,ip13,ip20,ip33,isub)+&
     &                      a13a20a33*up(1,iq)
                 field_up(2,ip13,ip20,ip33,isub)=field_up(2,ip13,ip20,ip33,isub)+&
     &                      a13a20a33*up(2,iq)
                 field_up(1,ip10,ip21,ip33,isub)=field_up(1,ip10,ip21,ip33,isub)+&
     &                      a10a21a33*up(1,iq)
                 field_up(2,ip10,ip21,ip33,isub)=field_up(2,ip10,ip21,ip33,isub)+&
     &                      a10a21a33*up(2,iq)
                 field_up(1,ip11,ip21,ip33,isub)=field_up(1,ip11,ip21,ip33,isub)+&
     &                      a11a21a33*up(1,iq)
                 field_up(2,ip11,ip21,ip33,isub)=field_up(2,ip11,ip21,ip33,isub)+&
     &                      a11a21a33*up(2,iq)
                 field_up(1,ip12,ip21,ip33,isub)=field_up(1,ip12,ip21,ip33,isub)+&
     &                      a12a21a33*up(1,iq)
                 field_up(2,ip12,ip21,ip33,isub)=field_up(2,ip12,ip21,ip33,isub)+&
     &                      a12a21a33*up(2,iq)
                 field_up(1,ip13,ip21,ip33,isub)=field_up(1,ip13,ip21,ip33,isub)+&
     &                      a13a21a33*up(1,iq)
                 field_up(2,ip13,ip21,ip33,isub)=field_up(2,ip13,ip21,ip33,isub)+&
     &                      a13a21a33*up(2,iq)
                 field_up(1,ip10,ip22,ip33,isub)=field_up(1,ip10,ip22,ip33,isub)+&
     &                      a10a22a33*up(1,iq)
                 field_up(2,ip10,ip22,ip33,isub)=field_up(2,ip10,ip22,ip33,isub)+&
     &                      a10a22a33*up(2,iq)
                 field_up(1,ip11,ip22,ip33,isub)=field_up(1,ip11,ip22,ip33,isub)+&
     &                      a11a22a33*up(1,iq)
                 field_up(2,ip11,ip22,ip33,isub)=field_up(2,ip11,ip22,ip33,isub)+&
     &                      a11a22a33*up(2,iq)
                 field_up(1,ip12,ip22,ip33,isub)=field_up(1,ip12,ip22,ip33,isub)+&
     &                      a12a22a33*up(1,iq)
                 field_up(2,ip12,ip22,ip33,isub)=field_up(2,ip12,ip22,ip33,isub)+&
     &                      a12a22a33*up(2,iq)
                 field_up(1,ip13,ip22,ip33,isub)=field_up(1,ip13,ip22,ip33,isub)+&
     &                      a13a22a33*up(1,iq)
                 field_up(2,ip13,ip22,ip33,isub)=field_up(2,ip13,ip22,ip33,isub)+&
     &                      a13a22a33*up(2,iq)
                 field_up(1,ip10,ip23,ip33,isub)=field_up(1,ip10,ip23,ip33,isub)+&
     &                      a10a23a33*up(1,iq)
                 field_up(2,ip10,ip23,ip33,isub)=field_up(2,ip10,ip23,ip33,isub)+&
     &                      a10a23a33*up(2,iq)
                 field_up(1,ip11,ip23,ip33,isub)=field_up(1,ip11,ip23,ip33,isub)+&
     &                      a11a23a33*up(1,iq)
                 field_up(2,ip11,ip23,ip33,isub)=field_up(2,ip11,ip23,ip33,isub)+&
     &                      a11a23a33*up(2,iq)
                 field_up(1,ip12,ip23,ip33,isub)=field_up(1,ip12,ip23,ip33,isub)+&
     &                      a12a23a33*up(1,iq)
                 field_up(2,ip12,ip23,ip33,isub)=field_up(2,ip12,ip23,ip33,isub)+&
     &                      a12a23a33*up(2,iq)
                 field_up(1,ip13,ip23,ip33,isub)=field_up(1,ip13,ip23,ip33,isub)+&
     &                      a13a23a33*up(1,iq)
                 field_up(2,ip13,ip23,ip33,isub)=field_up(2,ip13,ip23,ip33,isub)+&
     &                      a13a23a33*up(2,iq)

              ENDDO
              !----------------------------------------------------------------!
              !  Unrolled versions for 1-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 1) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01)
                 ip20 = INT(x02)
                 ip30 = INT(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33


                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

                 field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(1,iq)
                 field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(1,iq)
                 field_up(1,ip12,ip20,ip30,isub)=field_up(1,ip12,ip20,ip30,isub)+&
     &                      a12a20a30*up(1,iq)
                 field_up(1,ip13,ip20,ip30,isub)=field_up(1,ip13,ip20,ip30,isub)+&
     &                      a13a20a30*up(1,iq)
                 field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(1,iq)
                 field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(1,iq)
                 field_up(1,ip12,ip21,ip30,isub)=field_up(1,ip12,ip21,ip30,isub)+&
     &                      a12a21a30*up(1,iq)
                 field_up(1,ip13,ip21,ip30,isub)=field_up(1,ip13,ip21,ip30,isub)+&
     &                      a13a21a30*up(1,iq)
                 field_up(1,ip10,ip22,ip30,isub)=field_up(1,ip10,ip22,ip30,isub)+&
     &                      a10a22a30*up(1,iq)
                 field_up(1,ip11,ip22,ip30,isub)=field_up(1,ip11,ip22,ip30,isub)+&
     &                      a11a22a30*up(1,iq)
                 field_up(1,ip12,ip22,ip30,isub)=field_up(1,ip12,ip22,ip30,isub)+&
     &                      a12a22a30*up(1,iq)
                 field_up(1,ip13,ip22,ip30,isub)=field_up(1,ip13,ip22,ip30,isub)+&
     &                      a13a22a30*up(1,iq)
                 field_up(1,ip10,ip23,ip30,isub)=field_up(1,ip10,ip23,ip30,isub)+&
     &                      a10a23a30*up(1,iq)
                 field_up(1,ip11,ip23,ip30,isub)=field_up(1,ip11,ip23,ip30,isub)+&
     &                      a11a23a30*up(1,iq)
                 field_up(1,ip12,ip23,ip30,isub)=field_up(1,ip12,ip23,ip30,isub)+&
     &                      a12a23a30*up(1,iq)
                 field_up(1,ip13,ip23,ip30,isub)=field_up(1,ip13,ip23,ip30,isub)+&
     &                      a13a23a30*up(1,iq)
                 field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(1,iq)
                 field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(1,iq)
                 field_up(1,ip12,ip20,ip31,isub)=field_up(1,ip12,ip20,ip31,isub)+&
     &                      a12a20a31*up(1,iq)
                 field_up(1,ip13,ip20,ip31,isub)=field_up(1,ip13,ip20,ip31,isub)+&
     &                      a13a20a31*up(1,iq)
                 field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(1,iq)
                 field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(1,iq)
                 field_up(1,ip12,ip21,ip31,isub)=field_up(1,ip12,ip21,ip31,isub)+&
     &                      a12a21a31*up(1,iq)
                 field_up(1,ip13,ip21,ip31,isub)=field_up(1,ip13,ip21,ip31,isub)+&
     &                      a13a21a31*up(1,iq)
                 field_up(1,ip10,ip22,ip31,isub)=field_up(1,ip10,ip22,ip31,isub)+&
     &                      a10a22a31*up(1,iq)
                 field_up(1,ip11,ip22,ip31,isub)=field_up(1,ip11,ip22,ip31,isub)+&
     &                      a11a22a31*up(1,iq)
                 field_up(1,ip12,ip22,ip31,isub)=field_up(1,ip12,ip22,ip31,isub)+&
     &                      a12a22a31*up(1,iq)
                 field_up(1,ip13,ip22,ip31,isub)=field_up(1,ip13,ip22,ip31,isub)+&
     &                      a13a22a31*up(1,iq)
                 field_up(1,ip10,ip23,ip31,isub)=field_up(1,ip10,ip23,ip31,isub)+&
     &                      a10a23a31*up(1,iq)
                 field_up(1,ip11,ip23,ip31,isub)=field_up(1,ip11,ip23,ip31,isub)+&
     &                      a11a23a31*up(1,iq)
                 field_up(1,ip12,ip23,ip31,isub)=field_up(1,ip12,ip23,ip31,isub)+&
     &                      a12a23a31*up(1,iq)
                 field_up(1,ip13,ip23,ip31,isub)=field_up(1,ip13,ip23,ip31,isub)+&
     &                      a13a23a31*up(1,iq)
                 field_up(1,ip10,ip20,ip32,isub)=field_up(1,ip10,ip20,ip32,isub)+&
     &                      a10a20a32*up(1,iq)
                 field_up(1,ip11,ip20,ip32,isub)=field_up(1,ip11,ip20,ip32,isub)+&
     &                      a11a20a32*up(1,iq)
                 field_up(1,ip12,ip20,ip32,isub)=field_up(1,ip12,ip20,ip32,isub)+&
     &                      a12a20a32*up(1,iq)
                 field_up(1,ip13,ip20,ip32,isub)=field_up(1,ip13,ip20,ip32,isub)+&
     &                      a13a20a32*up(1,iq)
                 field_up(1,ip10,ip21,ip32,isub)=field_up(1,ip10,ip21,ip32,isub)+&
     &                      a10a21a32*up(1,iq)
                 field_up(1,ip11,ip21,ip32,isub)=field_up(1,ip11,ip21,ip32,isub)+&
     &                      a11a21a32*up(1,iq)
                 field_up(1,ip12,ip21,ip32,isub)=field_up(1,ip12,ip21,ip32,isub)+&
     &                      a12a21a32*up(1,iq)
                 field_up(1,ip13,ip21,ip32,isub)=field_up(1,ip13,ip21,ip32,isub)+&
     &                      a13a21a32*up(1,iq)
                 field_up(1,ip10,ip22,ip32,isub)=field_up(1,ip10,ip22,ip32,isub)+&
     &                      a10a22a32*up(1,iq)
                 field_up(1,ip11,ip22,ip32,isub)=field_up(1,ip11,ip22,ip32,isub)+&
     &                      a11a22a32*up(1,iq)
                 field_up(1,ip12,ip22,ip32,isub)=field_up(1,ip12,ip22,ip32,isub)+&
     &                      a12a22a32*up(1,iq)
                 field_up(1,ip13,ip22,ip32,isub)=field_up(1,ip13,ip22,ip32,isub)+&
     &                      a13a22a32*up(1,iq)
                 field_up(1,ip10,ip23,ip32,isub)=field_up(1,ip10,ip23,ip32,isub)+&
     &                      a10a23a32*up(1,iq)
                 field_up(1,ip11,ip23,ip32,isub)=field_up(1,ip11,ip23,ip32,isub)+&
     &                      a11a23a32*up(1,iq)
                 field_up(1,ip12,ip23,ip32,isub)=field_up(1,ip12,ip23,ip32,isub)+&
     &                      a12a23a32*up(1,iq)
                 field_up(1,ip13,ip23,ip32,isub)=field_up(1,ip13,ip23,ip32,isub)+&
     &                      a13a23a32*up(1,iq)
                 field_up(1,ip10,ip20,ip33,isub)=field_up(1,ip10,ip20,ip33,isub)+&
     &                      a10a20a33*up(1,iq)
                 field_up(1,ip11,ip20,ip33,isub)=field_up(1,ip11,ip20,ip33,isub)+&
     &                      a11a20a33*up(1,iq)
                 field_up(1,ip12,ip20,ip33,isub)=field_up(1,ip12,ip20,ip33,isub)+&
     &                      a12a20a33*up(1,iq)
                 field_up(1,ip13,ip20,ip33,isub)=field_up(1,ip13,ip20,ip33,isub)+&
     &                      a13a20a33*up(1,iq)
                 field_up(1,ip10,ip21,ip33,isub)=field_up(1,ip10,ip21,ip33,isub)+&
     &                      a10a21a33*up(1,iq)
                 field_up(1,ip11,ip21,ip33,isub)=field_up(1,ip11,ip21,ip33,isub)+&
     &                      a11a21a33*up(1,iq)
                 field_up(1,ip12,ip21,ip33,isub)=field_up(1,ip12,ip21,ip33,isub)+&
     &                      a12a21a33*up(1,iq)
                 field_up(1,ip13,ip21,ip33,isub)=field_up(1,ip13,ip21,ip33,isub)+&
     &                      a13a21a33*up(1,iq)
                 field_up(1,ip10,ip22,ip33,isub)=field_up(1,ip10,ip22,ip33,isub)+&
     &                      a10a22a33*up(1,iq)
                 field_up(1,ip11,ip22,ip33,isub)=field_up(1,ip11,ip22,ip33,isub)+&
     &                      a11a22a33*up(1,iq)
                 field_up(1,ip12,ip22,ip33,isub)=field_up(1,ip12,ip22,ip33,isub)+&
     &                      a12a22a33*up(1,iq)
                 field_up(1,ip13,ip22,ip33,isub)=field_up(1,ip13,ip22,ip33,isub)+&
     &                      a13a22a33*up(1,iq)
                 field_up(1,ip10,ip23,ip33,isub)=field_up(1,ip10,ip23,ip33,isub)+&
     &                      a10a23a33*up(1,iq)
                 field_up(1,ip11,ip23,ip33,isub)=field_up(1,ip11,ip23,ip33,isub)+&
     &                      a11a23a33*up(1,iq)
                 field_up(1,ip12,ip23,ip33,isub)=field_up(1,ip12,ip23,ip33,isub)+&
     &                      a12a23a33*up(1,iq)
                 field_up(1,ip13,ip23,ip33,isub)=field_up(1,ip13,ip23,ip33,isub)+&
     &                      a13a23a33*up(1,iq)

              ENDDO
              !----------------------------------------------------------------!
              !  All other lda are NOT UNROLLED. Vectorization will be over lda!
              !----------------------------------------------------------------!
           ELSE
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01)
                 ip20 = INT(x02)
                 ip30 = INT(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

                 DO ldn=1,lda

                    field_up(ldn,ip10,ip20,ip30,isub)=field_up(ldn,ip10,ip20,ip30,isub)+&
     &                         a10a20a30*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip30,isub)=field_up(ldn,ip11,ip20,ip30,isub)+&
     &                         a11a20a30*up(ldn,iq)
                    field_up(ldn,ip12,ip20,ip30,isub)=field_up(ldn,ip12,ip20,ip30,isub)+&
     &                         a12a20a30*up(ldn,iq)
                    field_up(ldn,ip13,ip20,ip30,isub)=field_up(ldn,ip13,ip20,ip30,isub)+&
     &                         a13a20a30*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip30,isub)=field_up(ldn,ip10,ip21,ip30,isub)+&
     &                         a10a21a30*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip30,isub)=field_up(ldn,ip11,ip21,ip30,isub)+&
     &                         a11a21a30*up(ldn,iq)
                    field_up(ldn,ip12,ip21,ip30,isub)=field_up(ldn,ip12,ip21,ip30,isub)+&
     &                         a12a21a30*up(ldn,iq)
                    field_up(ldn,ip13,ip21,ip30,isub)=field_up(ldn,ip13,ip21,ip30,isub)+&
     &                         a13a21a30*up(ldn,iq)
                    field_up(ldn,ip10,ip22,ip30,isub)=field_up(ldn,ip10,ip22,ip30,isub)+&
     &                         a10a22a30*up(ldn,iq)
                    field_up(ldn,ip11,ip22,ip30,isub)=field_up(ldn,ip11,ip22,ip30,isub)+&
     &                         a11a22a30*up(ldn,iq)
                    field_up(ldn,ip12,ip22,ip30,isub)=field_up(ldn,ip12,ip22,ip30,isub)+&
     &                         a12a22a30*up(ldn,iq)
                    field_up(ldn,ip13,ip22,ip30,isub)=field_up(ldn,ip13,ip22,ip30,isub)+&
     &                         a13a22a30*up(ldn,iq)
                    field_up(ldn,ip10,ip23,ip30,isub)=field_up(ldn,ip10,ip23,ip30,isub)+&
     &                         a10a23a30*up(ldn,iq)
                    field_up(ldn,ip11,ip23,ip30,isub)=field_up(ldn,ip11,ip23,ip30,isub)+&
     &                         a11a23a30*up(ldn,iq)
                    field_up(ldn,ip12,ip23,ip30,isub)=field_up(ldn,ip12,ip23,ip30,isub)+&
     &                         a12a23a30*up(ldn,iq)
                    field_up(ldn,ip13,ip23,ip30,isub)=field_up(ldn,ip13,ip23,ip30,isub)+&
     &                         a13a23a30*up(ldn,iq)
                    field_up(ldn,ip10,ip20,ip31,isub)=field_up(ldn,ip10,ip20,ip31,isub)+&
     &                         a10a20a31*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip31,isub)=field_up(ldn,ip11,ip20,ip31,isub)+&
     &                         a11a20a31*up(ldn,iq)
                    field_up(ldn,ip12,ip20,ip31,isub)=field_up(ldn,ip12,ip20,ip31,isub)+&
     &                         a12a20a31*up(ldn,iq)
                    field_up(ldn,ip13,ip20,ip31,isub)=field_up(ldn,ip13,ip20,ip31,isub)+&
     &                         a13a20a31*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip31,isub)=field_up(ldn,ip10,ip21,ip31,isub)+&
     &                         a10a21a31*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip31,isub)=field_up(ldn,ip11,ip21,ip31,isub)+&
     &                         a11a21a31*up(ldn,iq)
                    field_up(ldn,ip12,ip21,ip31,isub)=field_up(ldn,ip12,ip21,ip31,isub)+&
     &                         a12a21a31*up(ldn,iq)
                    field_up(ldn,ip13,ip21,ip31,isub)=field_up(ldn,ip13,ip21,ip31,isub)+&
     &                         a13a21a31*up(ldn,iq)
                    field_up(ldn,ip10,ip22,ip31,isub)=field_up(ldn,ip10,ip22,ip31,isub)+&
     &                         a10a22a31*up(ldn,iq)
                    field_up(ldn,ip11,ip22,ip31,isub)=field_up(ldn,ip11,ip22,ip31,isub)+&
     &                         a11a22a31*up(ldn,iq)
                    field_up(ldn,ip12,ip22,ip31,isub)=field_up(ldn,ip12,ip22,ip31,isub)+&
     &                         a12a22a31*up(ldn,iq)
                    field_up(ldn,ip13,ip22,ip31,isub)=field_up(ldn,ip13,ip22,ip31,isub)+&
     &                         a13a22a31*up(ldn,iq)
                    field_up(ldn,ip10,ip23,ip31,isub)=field_up(ldn,ip10,ip23,ip31,isub)+&
     &                         a10a23a31*up(ldn,iq)
                    field_up(ldn,ip11,ip23,ip31,isub)=field_up(ldn,ip11,ip23,ip31,isub)+&
     &                         a11a23a31*up(ldn,iq)
                    field_up(ldn,ip12,ip23,ip31,isub)=field_up(ldn,ip12,ip23,ip31,isub)+&
     &                         a12a23a31*up(ldn,iq)
                    field_up(ldn,ip13,ip23,ip31,isub)=field_up(ldn,ip13,ip23,ip31,isub)+&
     &                         a13a23a31*up(ldn,iq)
                    field_up(ldn,ip10,ip20,ip32,isub)=field_up(ldn,ip10,ip20,ip32,isub)+&
     &                         a10a20a32*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip32,isub)=field_up(ldn,ip11,ip20,ip32,isub)+&
     &                         a11a20a32*up(ldn,iq)
                    field_up(ldn,ip12,ip20,ip32,isub)=field_up(ldn,ip12,ip20,ip32,isub)+&
     &                         a12a20a32*up(ldn,iq)
                    field_up(ldn,ip13,ip20,ip32,isub)=field_up(ldn,ip13,ip20,ip32,isub)+&
     &                         a13a20a32*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip32,isub)=field_up(ldn,ip10,ip21,ip32,isub)+&
     &                         a10a21a32*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip32,isub)=field_up(ldn,ip11,ip21,ip32,isub)+&
     &                         a11a21a32*up(ldn,iq)
                    field_up(ldn,ip12,ip21,ip32,isub)=field_up(ldn,ip12,ip21,ip32,isub)+&
     &                         a12a21a32*up(ldn,iq)
                    field_up(ldn,ip13,ip21,ip32,isub)=field_up(ldn,ip13,ip21,ip32,isub)+&
     &                         a13a21a32*up(ldn,iq)
                    field_up(ldn,ip10,ip22,ip32,isub)=field_up(ldn,ip10,ip22,ip32,isub)+&
     &                         a10a22a32*up(ldn,iq)
                    field_up(ldn,ip11,ip22,ip32,isub)=field_up(ldn,ip11,ip22,ip32,isub)+&
     &                         a11a22a32*up(ldn,iq)
                    field_up(ldn,ip12,ip22,ip32,isub)=field_up(ldn,ip12,ip22,ip32,isub)+&
     &                         a12a22a32*up(ldn,iq)
                    field_up(ldn,ip13,ip22,ip32,isub)=field_up(ldn,ip13,ip22,ip32,isub)+&
     &                         a13a22a32*up(ldn,iq)
                    field_up(ldn,ip10,ip23,ip32,isub)=field_up(ldn,ip10,ip23,ip32,isub)+&
     &                         a10a23a32*up(ldn,iq)
                    field_up(ldn,ip11,ip23,ip32,isub)=field_up(ldn,ip11,ip23,ip32,isub)+&
     &                         a11a23a32*up(ldn,iq)
                    field_up(ldn,ip12,ip23,ip32,isub)=field_up(ldn,ip12,ip23,ip32,isub)+&
     &                         a12a23a32*up(ldn,iq)
                    field_up(ldn,ip13,ip23,ip32,isub)=field_up(ldn,ip13,ip23,ip32,isub)+&
     &                         a13a23a32*up(ldn,iq)
                    field_up(ldn,ip10,ip20,ip33,isub)=field_up(ldn,ip10,ip20,ip33,isub)+&
     &                         a10a20a33*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip33,isub)=field_up(ldn,ip11,ip20,ip33,isub)+&
     &                         a11a20a33*up(ldn,iq)
                    field_up(ldn,ip12,ip20,ip33,isub)=field_up(ldn,ip12,ip20,ip33,isub)+&
     &                         a12a20a33*up(ldn,iq)
                    field_up(ldn,ip13,ip20,ip33,isub)=field_up(ldn,ip13,ip20,ip33,isub)+&
     &                         a13a20a33*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip33,isub)=field_up(ldn,ip10,ip21,ip33,isub)+&
     &                         a10a21a33*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip33,isub)=field_up(ldn,ip11,ip21,ip33,isub)+&
     &                         a11a21a33*up(ldn,iq)
                    field_up(ldn,ip12,ip21,ip33,isub)=field_up(ldn,ip12,ip21,ip33,isub)+&
     &                         a12a21a33*up(ldn,iq)
                    field_up(ldn,ip13,ip21,ip33,isub)=field_up(ldn,ip13,ip21,ip33,isub)+&
     &                         a13a21a33*up(ldn,iq)
                    field_up(ldn,ip10,ip22,ip33,isub)=field_up(ldn,ip10,ip22,ip33,isub)+&
     &                         a10a22a33*up(ldn,iq)
                    field_up(ldn,ip11,ip22,ip33,isub)=field_up(ldn,ip11,ip22,ip33,isub)+&
     &                         a11a22a33*up(ldn,iq)
                    field_up(ldn,ip12,ip22,ip33,isub)=field_up(ldn,ip12,ip22,ip33,isub)+&
     &                         a12a22a33*up(ldn,iq)
                    field_up(ldn,ip13,ip22,ip33,isub)=field_up(ldn,ip13,ip22,ip33,isub)+&
     &                         a13a22a33*up(ldn,iq)
                    field_up(ldn,ip10,ip23,ip33,isub)=field_up(ldn,ip10,ip23,ip33,isub)+&
     &                         a10a23a33*up(ldn,iq)
                    field_up(ldn,ip11,ip23,ip33,isub)=field_up(ldn,ip11,ip23,ip33,isub)+&
     &                         a11a23a33*up(ldn,iq)
                    field_up(ldn,ip12,ip23,ip33,isub)=field_up(ldn,ip12,ip23,ip33,isub)+&
     &                         a12a23a33*up(ldn,iq)
                    field_up(ldn,ip13,ip23,ip33,isub)=field_up(ldn,ip13,ip23,ip33,isub)+&
     &                         a13a23a33*up(ldn,iq)

                 ENDDO    ! lda
              ENDDO        ! iq
           END IF          ! unrolled lda cases
#elif __MODE == __SCA
           isubl = ppm_isublist(isub,topoid)
           DO ip = 1,store_info(isub)
              
              iq    = list_sub(isub,ip)

              x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
              x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
              x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi
              
              ip10 = INT(x01)
              ip20 = INT(x02)
              ip30 = INT(x03)
               
              ip11 = ip10 + 1
              ip21 = ip20 + 1
              ip31 = ip30 + 1

              ip12 = ip11 + 1
              ip22 = ip21 + 1
              ip32 = ip31 + 1

              ip13 = ip11 + 2
              ip23 = ip21 + 2
              ip33 = ip31 + 2

              xp1 = x01-REAL(ip10,mk)
              xp2 = x02-REAL(ip20,mk)
              xp3 = x03-REAL(ip30,mk)

              x10 = xp1 + 1.0_mk
              x11 = x10 - 1.0_mk
              x12 = x10 - 2.0_mk
              x13 = x10 - 3.0_mk

              x20 = xp2 + 1.0_mk
              x21 = x20 - 1.0_mk
              x22 = x20 - 2.0_mk
              x23 = x20 - 3.0_mk

              x30 = xp3 + 1.0_mk
              x31 = x30 - 1.0_mk
              x32 = x30 - 2.0_mk
              x33 = x30 - 3.0_mk

              a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
              a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
              a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

              a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
              a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
              a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

              a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
              a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
              a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

              a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
              a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
              a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

              a10a20a30 = a10*a20*a30
              a10a20a31 = a10*a20*a31
              a10a20a32 = a10*a20*a32
              a10a20a33 = a10*a20*a33
              a10a21a30 = a10*a21*a30
              a10a21a31 = a10*a21*a31
              a10a21a32 = a10*a21*a32
              a10a21a33 = a10*a21*a33
              a10a22a30 = a10*a22*a30
              a10a22a31 = a10*a22*a31
              a10a22a32 = a10*a22*a32
              a10a22a33 = a10*a22*a33
              a10a23a30 = a10*a23*a30
              a10a23a31 = a10*a23*a31
              a10a23a32 = a10*a23*a32
              a10a23a33 = a10*a23*a33
              a11a20a30 = a11*a20*a30
              a11a20a31 = a11*a20*a31
              a11a20a32 = a11*a20*a32
              a11a20a33 = a11*a20*a33
              a11a21a30 = a11*a21*a30
              a11a21a31 = a11*a21*a31
              a11a21a32 = a11*a21*a32
              a11a21a33 = a11*a21*a33
              a11a22a30 = a11*a22*a30
              a11a22a31 = a11*a22*a31
              a11a22a32 = a11*a22*a32
              a11a22a33 = a11*a22*a33
              a11a23a30 = a11*a23*a30
              a11a23a31 = a11*a23*a31
              a11a23a32 = a11*a23*a32
              a11a23a33 = a11*a23*a33
              a12a20a30 = a12*a20*a30
              a12a20a31 = a12*a20*a31
              a12a20a32 = a12*a20*a32
              a12a20a33 = a12*a20*a33
              a12a21a30 = a12*a21*a30
              a12a21a31 = a12*a21*a31
              a12a21a32 = a12*a21*a32
              a12a21a33 = a12*a21*a33
              a12a22a30 = a12*a22*a30
              a12a22a31 = a12*a22*a31
              a12a22a32 = a12*a22*a32
              a12a22a33 = a12*a22*a33
              a12a23a30 = a12*a23*a30
              a12a23a31 = a12*a23*a31
              a12a23a32 = a12*a23*a32
              a12a23a33 = a12*a23*a33
              a13a20a30 = a13*a20*a30
              a13a20a31 = a13*a20*a31
              a13a20a32 = a13*a20*a32
              a13a20a33 = a13*a20*a33
              a13a21a30 = a13*a21*a30
              a13a21a31 = a13*a21*a31
              a13a21a32 = a13*a21*a32
              a13a21a33 = a13*a21*a33
              a13a22a30 = a13*a22*a30
              a13a22a31 = a13*a22*a31
              a13a22a32 = a13*a22*a32
              a13a22a33 = a13*a22*a33
              a13a23a30 = a13*a23*a30
              a13a23a31 = a13*a23*a31
              a13a23a32 = a13*a23*a32
              a13a23a33 = a13*a23*a33

              field_up(ip10,ip20,ip30,isub)=field_up(ip10,ip20,ip30,isub)+&
     &                   a10a20a30*up(iq)
              field_up(ip11,ip20,ip30,isub)=field_up(ip11,ip20,ip30,isub)+&
     &                   a11a20a30*up(iq)
              field_up(ip12,ip20,ip30,isub)=field_up(ip12,ip20,ip30,isub)+&
     &                   a12a20a30*up(iq)
              field_up(ip13,ip20,ip30,isub)=field_up(ip13,ip20,ip30,isub)+&
     &                   a13a20a30*up(iq)
              field_up(ip10,ip21,ip30,isub)=field_up(ip10,ip21,ip30,isub)+&
     &                   a10a21a30*up(iq)
              field_up(ip11,ip21,ip30,isub)=field_up(ip11,ip21,ip30,isub)+&
     &                   a11a21a30*up(iq)
              field_up(ip12,ip21,ip30,isub)=field_up(ip12,ip21,ip30,isub)+&
     &                   a12a21a30*up(iq)
              field_up(ip13,ip21,ip30,isub)=field_up(ip13,ip21,ip30,isub)+&
     &                   a13a21a30*up(iq)
              field_up(ip10,ip22,ip30,isub)=field_up(ip10,ip22,ip30,isub)+&
     &                   a10a22a30*up(iq)
              field_up(ip11,ip22,ip30,isub)=field_up(ip11,ip22,ip30,isub)+&
     &                   a11a22a30*up(iq)
              field_up(ip12,ip22,ip30,isub)=field_up(ip12,ip22,ip30,isub)+&
     &                   a12a22a30*up(iq)
              field_up(ip13,ip22,ip30,isub)=field_up(ip13,ip22,ip30,isub)+&
     &                   a13a22a30*up(iq)
              field_up(ip10,ip23,ip30,isub)=field_up(ip10,ip23,ip30,isub)+&
     &                   a10a23a30*up(iq)
              field_up(ip11,ip23,ip30,isub)=field_up(ip11,ip23,ip30,isub)+&
     &                   a11a23a30*up(iq)
              field_up(ip12,ip23,ip30,isub)=field_up(ip12,ip23,ip30,isub)+&
     &                   a12a23a30*up(iq)
              field_up(ip13,ip23,ip30,isub)=field_up(ip13,ip23,ip30,isub)+&
     &                   a13a23a30*up(iq)
              field_up(ip10,ip20,ip31,isub)=field_up(ip10,ip20,ip31,isub)+&
     &                   a10a20a31*up(iq)
              field_up(ip11,ip20,ip31,isub)=field_up(ip11,ip20,ip31,isub)+&
     &                   a11a20a31*up(iq)
              field_up(ip12,ip20,ip31,isub)=field_up(ip12,ip20,ip31,isub)+&
     &                   a12a20a31*up(iq)
              field_up(ip13,ip20,ip31,isub)=field_up(ip13,ip20,ip31,isub)+&
     &                   a13a20a31*up(iq)
              field_up(ip10,ip21,ip31,isub)=field_up(ip10,ip21,ip31,isub)+&
     &                   a10a21a31*up(iq)
              field_up(ip11,ip21,ip31,isub)=field_up(ip11,ip21,ip31,isub)+&
     &                   a11a21a31*up(iq)
              field_up(ip12,ip21,ip31,isub)=field_up(ip12,ip21,ip31,isub)+&
     &                   a12a21a31*up(iq)
              field_up(ip13,ip21,ip31,isub)=field_up(ip13,ip21,ip31,isub)+&
     &                   a13a21a31*up(iq)
              field_up(ip10,ip22,ip31,isub)=field_up(ip10,ip22,ip31,isub)+&
     &                   a10a22a31*up(iq)
              field_up(ip11,ip22,ip31,isub)=field_up(ip11,ip22,ip31,isub)+&
     &                   a11a22a31*up(iq)
              field_up(ip12,ip22,ip31,isub)=field_up(ip12,ip22,ip31,isub)+&
     &                   a12a22a31*up(iq)
              field_up(ip13,ip22,ip31,isub)=field_up(ip13,ip22,ip31,isub)+&
     &                   a13a22a31*up(iq)
              field_up(ip10,ip23,ip31,isub)=field_up(ip10,ip23,ip31,isub)+&
     &                   a10a23a31*up(iq)
              field_up(ip11,ip23,ip31,isub)=field_up(ip11,ip23,ip31,isub)+&
     &                   a11a23a31*up(iq)
              field_up(ip12,ip23,ip31,isub)=field_up(ip12,ip23,ip31,isub)+&
     &                   a12a23a31*up(iq)
              field_up(ip13,ip23,ip31,isub)=field_up(ip13,ip23,ip31,isub)+&
     &                   a13a23a31*up(iq)
              field_up(ip10,ip20,ip32,isub)=field_up(ip10,ip20,ip32,isub)+&
     &                   a10a20a32*up(iq)
              field_up(ip11,ip20,ip32,isub)=field_up(ip11,ip20,ip32,isub)+&
     &                   a11a20a32*up(iq)
              field_up(ip12,ip20,ip32,isub)=field_up(ip12,ip20,ip32,isub)+&
     &                   a12a20a32*up(iq)
              field_up(ip13,ip20,ip32,isub)=field_up(ip13,ip20,ip32,isub)+&
     &                   a13a20a32*up(iq)
              field_up(ip10,ip21,ip32,isub)=field_up(ip10,ip21,ip32,isub)+&
     &                   a10a21a32*up(iq)
              field_up(ip11,ip21,ip32,isub)=field_up(ip11,ip21,ip32,isub)+&
     &                   a11a21a32*up(iq)
              field_up(ip12,ip21,ip32,isub)=field_up(ip12,ip21,ip32,isub)+&
     &                   a12a21a32*up(iq)
              field_up(ip13,ip21,ip32,isub)=field_up(ip13,ip21,ip32,isub)+&
     &                   a13a21a32*up(iq)
              field_up(ip10,ip22,ip32,isub)=field_up(ip10,ip22,ip32,isub)+&
     &                   a10a22a32*up(iq)
              field_up(ip11,ip22,ip32,isub)=field_up(ip11,ip22,ip32,isub)+&
     &                   a11a22a32*up(iq)
              field_up(ip12,ip22,ip32,isub)=field_up(ip12,ip22,ip32,isub)+&
     &                   a12a22a32*up(iq)
              field_up(ip13,ip22,ip32,isub)=field_up(ip13,ip22,ip32,isub)+&
     &                   a13a22a32*up(iq)
              field_up(ip10,ip23,ip32,isub)=field_up(ip10,ip23,ip32,isub)+&
     &                   a10a23a32*up(iq)
              field_up(ip11,ip23,ip32,isub)=field_up(ip11,ip23,ip32,isub)+&
     &                   a11a23a32*up(iq)
              field_up(ip12,ip23,ip32,isub)=field_up(ip12,ip23,ip32,isub)+&
     &                   a12a23a32*up(iq)
              field_up(ip13,ip23,ip32,isub)=field_up(ip13,ip23,ip32,isub)+&
     &                   a13a23a32*up(iq)
              field_up(ip10,ip20,ip33,isub)=field_up(ip10,ip20,ip33,isub)+&
     &                   a10a20a33*up(iq)
              field_up(ip11,ip20,ip33,isub)=field_up(ip11,ip20,ip33,isub)+&
     &                   a11a20a33*up(iq)
              field_up(ip12,ip20,ip33,isub)=field_up(ip12,ip20,ip33,isub)+&
     &                   a12a20a33*up(iq)
              field_up(ip13,ip20,ip33,isub)=field_up(ip13,ip20,ip33,isub)+&
     &                   a13a20a33*up(iq)
              field_up(ip10,ip21,ip33,isub)=field_up(ip10,ip21,ip33,isub)+&
     &                   a10a21a33*up(iq)
              field_up(ip11,ip21,ip33,isub)=field_up(ip11,ip21,ip33,isub)+&
     &                   a11a21a33*up(iq)
              field_up(ip12,ip21,ip33,isub)=field_up(ip12,ip21,ip33,isub)+&
     &                   a12a21a33*up(iq)
              field_up(ip13,ip21,ip33,isub)=field_up(ip13,ip21,ip33,isub)+&
     &                   a13a21a33*up(iq)
              field_up(ip10,ip22,ip33,isub)=field_up(ip10,ip22,ip33,isub)+&
     &                   a10a22a33*up(iq)
              field_up(ip11,ip22,ip33,isub)=field_up(ip11,ip22,ip33,isub)+&
     &                   a11a22a33*up(iq)
              field_up(ip12,ip22,ip33,isub)=field_up(ip12,ip22,ip33,isub)+&
     &                   a12a22a33*up(iq)
              field_up(ip13,ip22,ip33,isub)=field_up(ip13,ip22,ip33,isub)+&
     &                   a13a22a33*up(iq)
              field_up(ip10,ip23,ip33,isub)=field_up(ip10,ip23,ip33,isub)+&
     &                   a10a23a33*up(iq)
              field_up(ip11,ip23,ip33,isub)=field_up(ip11,ip23,ip33,isub)+&
     &                   a11a23a33*up(iq)
              field_up(ip12,ip23,ip33,isub)=field_up(ip12,ip23,ip33,isub)+&
     &                   a12a23a33*up(iq)
              field_up(ip13,ip23,ip33,isub)=field_up(ip13,ip23,ip33,isub)+&
     &                   a13a23a33*up(iq)
           ENDDO        ! iq
#endif     
        END DO              ! loop over subs
     CASE(ppm_param_rmsh_kernel_bsp2)
        !----------------------------------------------------------------------!
        !  B-spline 2 (Witch hat)
        !----------------------------------------------------------------------!
        DO isub = 1,ppm_nsublist(topoid)
           
#if __MODE == __VEC
           !-------------------------------------------------------------------!
           !  Unrolled versions for 4-vectors
           !-------------------------------------------------------------------!
           IF(lda.EQ.4) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)                              

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31

                 field_up(1,ip10,ip20,ip30,isub) = field_up(1,ip10,ip20,ip30,isub) + &
     &                       a10*a20*a30*up(1,iq)
                 field_up(1,ip10,ip20,ip31,isub) = field_up(1,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(1,iq)
                 field_up(1,ip10,ip21,ip30,isub) = field_up(1,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(1,iq)
                 field_up(1,ip10,ip21,ip31,isub) = field_up(1,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(1,iq)
                 field_up(1,ip11,ip20,ip30,isub) = field_up(1,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(1,iq)
                 field_up(1,ip11,ip20,ip31,isub) = field_up(1,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(1,iq)
                 field_up(1,ip11,ip21,ip30,isub) = field_up(1,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(1,iq)
                 field_up(1,ip11,ip21,ip31,isub) = field_up(1,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(1,iq)

                 field_up(2,ip10,ip20,ip30,isub) = field_up(2,ip10,ip20,ip30,isub) + &
     &                       a10*a20*a30*up(2,iq)
                 field_up(2,ip10,ip20,ip31,isub) = field_up(2,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(2,iq)
                 field_up(2,ip10,ip21,ip30,isub) = field_up(2,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(2,iq)
                 field_up(2,ip10,ip21,ip31,isub) = field_up(2,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(2,iq)
                 field_up(2,ip11,ip20,ip30,isub) = field_up(2,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(2,iq)
                 field_up(2,ip11,ip20,ip31,isub) = field_up(2,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(2,iq)
                 field_up(2,ip11,ip21,ip30,isub) = field_up(2,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(2,iq)
                 field_up(2,ip11,ip21,ip31,isub) = field_up(2,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(2,iq)

                 field_up(3,ip10,ip20,ip30,isub) = field_up(3,ip10,ip20,ip30,isub) + &
     &                       a10*a20*a30*up(3,iq)
                 field_up(3,ip10,ip20,ip31,isub) = field_up(3,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(3,iq)
                 field_up(3,ip10,ip21,ip30,isub) = field_up(3,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(3,iq)
                 field_up(3,ip10,ip21,ip31,isub) = field_up(3,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(3,iq)
                 field_up(3,ip11,ip20,ip30,isub) = field_up(3,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(3,iq)
                 field_up(3,ip11,ip20,ip31,isub) = field_up(3,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(3,iq)
                 field_up(3,ip11,ip21,ip30,isub) = field_up(3,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(3,iq)
                 field_up(3,ip11,ip21,ip31,isub) = field_up(3,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(3,iq)

                 field_up(4,ip10,ip20,ip30,isub) = field_up(4,ip10,ip20,ip30,isub) + &
     &                       a10*a20*a30*up(4,iq)
                 field_up(4,ip10,ip20,ip31,isub) = field_up(4,ip10,ip20,ip31,isub) + &
     &                       a10*a20*a31*up(4,iq)
                 field_up(4,ip10,ip21,ip30,isub) = field_up(4,ip10,ip21,ip30,isub) + &
     &                       a10*a21*a30*up(4,iq)
                 field_up(4,ip10,ip21,ip31,isub) = field_up(4,ip10,ip21,ip31,isub) + &
     &                       a10*a21*a31*up(4,iq)
                 field_up(4,ip11,ip20,ip30,isub) = field_up(4,ip11,ip20,ip30,isub) + &
     &                       a11*a20*a30*up(4,iq)
                 field_up(4,ip11,ip20,ip31,isub) = field_up(4,ip11,ip20,ip31,isub) + &
     &                       a11*a20*a31*up(4,iq)
                 field_up(4,ip11,ip21,ip30,isub) = field_up(4,ip11,ip21,ip30,isub) + &
     &                       a11*a21*a30*up(4,iq)
                 field_up(4,ip11,ip21,ip31,isub) = field_up(4,ip11,ip21,ip31,isub) + &
     &                       a11*a21*a31*up(4,iq)

              END DO
              !----------------------------------------------------------------!
              !  Unrolled versions for 3-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 3) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)                              

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31

                 
                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                                
                 field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(1,iq)
                 field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(2,iq)
                 field_up(3,ip10,ip20,ip30,isub)=field_up(3,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(3,iq)
                 field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(1,iq)
                 field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(2,iq)
                 field_up(3,ip11,ip20,ip30,isub)=field_up(3,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(3,iq)
                 field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(1,iq)
                 field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(2,iq)
                 field_up(3,ip10,ip21,ip30,isub)=field_up(3,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(3,iq)
                 field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(1,iq)
                 field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(2,iq)
                 field_up(3,ip11,ip21,ip30,isub)=field_up(3,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(3,iq)
                 field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(1,iq)
                 field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(2,iq)
                 field_up(3,ip10,ip20,ip31,isub)=field_up(3,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(3,iq)
                 field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(1,iq)
                 field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(2,iq)
                 field_up(3,ip11,ip20,ip31,isub)=field_up(3,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(3,iq)
                 field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(1,iq)
                 field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(2,iq)
                 field_up(3,ip10,ip21,ip31,isub)=field_up(3,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(3,iq)
                 field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(1,iq)
                 field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(2,iq)
                 field_up(3,ip11,ip21,ip31,isub)=field_up(3,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(3,iq)
              END DO
              !----------------------------------------------------------------!
              !  Unrolled versions for 2-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 2) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)
                 
                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi
                 
                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)                              

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31
                 
                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 
                 field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(1,iq)
                 field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(2,iq)
                 field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(1,iq)
                 field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(2,iq)
                 field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(1,iq)
                 field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(2,iq)
                 field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(1,iq)
                 field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(2,iq)
                 field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(1,iq)
                 field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(2,iq)
                 field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(1,iq)
                 field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(2,iq)
                 field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(1,iq)
                 field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(2,iq)
                 field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(1,iq)
                 field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(2,iq)
              ENDDO
              !----------------------------------------------------------------!
              !  Unrolled versions for 1-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 1) THEN
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)
                 
                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi
                 
                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1
                 
                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1
                 
                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)                              
                 
                 x10 = xp1
                 x11 = x10 - 1.0_mk
                 
                 x20 = xp2
                 x21 = x20 - 1.0_mk
                 
                 x30 = xp3
                 x31 = x30 - 1.0_mk
                 
                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30
                 
                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31
                 
                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 
                 field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
     &                      a10a20a30*up(1,iq)
                 field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
     &                      a11a20a30*up(1,iq)
                 field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
     &                      a10a21a30*up(1,iq)
                 field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
     &                      a11a21a30*up(1,iq)
                 field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
     &                      a10a20a31*up(1,iq)
                 field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
     &                      a11a20a31*up(1,iq)
                 field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
     &                      a10a21a31*up(1,iq)
                 field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
     &                      a11a21a31*up(1,iq)

              ENDDO
              !----------------------------------------------------------------!
              !  All other lda are NOT UNROLLED. Vectorization will be over lda!
              !----------------------------------------------------------------!
           ELSE
              isubl = ppm_isublist(isub,topoid)
              DO ip = 1,store_info(isub)
                 
                 iq    = list_sub(isub,ip)
                 
                 x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)                              

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 DO ldn=1,lda
                    
                    field_up(ldn,ip10,ip20,ip30,isub)=field_up(ldn,ip10,ip20,ip30,isub)+&
     &                         a10a20a30*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip30,isub)=field_up(ldn,ip11,ip20,ip30,isub)+&
     &                         a11a20a30*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip30,isub)=field_up(ldn,ip10,ip21,ip30,isub)+&
     &                         a10a21a30*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip30,isub)=field_up(ldn,ip11,ip21,ip30,isub)+&
     &                         a11a21a30*up(ldn,iq)
                    field_up(ldn,ip10,ip20,ip31,isub)=field_up(ldn,ip10,ip20,ip31,isub)+&
     &                         a10a20a31*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip31,isub)=field_up(ldn,ip11,ip20,ip31,isub)+&
     &                         a11a20a31*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip31,isub)=field_up(ldn,ip10,ip21,ip31,isub)+&
     &                         a10a21a31*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip31,isub)=field_up(ldn,ip11,ip21,ip31,isub)+&
     &                         a11a21a31*up(ldn,iq)
                 ENDDO    ! lda
              ENDDO        ! iq
           END IF          ! unrolled lda cases
#elif __MODE == __SCA
           isubl = ppm_isublist(isub,topoid)
           DO ip = 1,store_info(isub)
              
              iq    = list_sub(isub,ip)
              
              x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
              x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
              x03 = (xp(3,iq)-min_sub(3,isubl,topoid))*dxzi

              ip10 = INT(x01) + 1
              ip20 = INT(x02) + 1
              ip30 = INT(x03) + 1

              ip11 = ip10 + 1
              ip21 = ip20 + 1
              ip31 = ip30 + 1

              xp1 = x01-REAL(ip10-1,mk)
              xp2 = x02-REAL(ip20-1,mk)
              xp3 = x03-REAL(ip30-1,mk)                         

              x10 = xp1
              x11 = x10 - 1.0_mk

              x20 = xp2
              x21 = x20 - 1.0_mk

              x30 = xp3
              x31 = x30 - 1.0_mk

              a10 = 1.0_mk - x10
              a20 = 1.0_mk - x20
              a30 = 1.0_mk - x30

              a11 = 1.0_mk + x11
              a21 = 1.0_mk + x21
              a31 = 1.0_mk + x31
              
              a10a20a30 = a10*a20*a30
              a10a20a31 = a10*a20*a31
              a10a21a30 = a10*a21*a30
              a10a21a31 = a10*a21*a31
              a11a20a30 = a11*a20*a30
              a11a20a31 = a11*a20*a31
              a11a21a30 = a11*a21*a30
              a11a21a31 = a11*a21*a31

              
              field_up(ip10,ip20,ip30,isub)=field_up(ip10,ip20,ip30,isub)+&
     &                   a10a20a30*up(iq)
              field_up(ip11,ip20,ip30,isub)=field_up(ip11,ip20,ip30,isub)+&
     &                   a11a20a30*up(iq)
              field_up(ip10,ip21,ip30,isub)=field_up(ip10,ip21,ip30,isub)+&
     &                   a10a21a30*up(iq)
              field_up(ip11,ip21,ip30,isub)=field_up(ip11,ip21,ip30,isub)+&
     &                   a11a21a30*up(iq)
              field_up(ip10,ip20,ip31,isub)=field_up(ip10,ip20,ip31,isub)+&
     &                   a10a20a31*up(iq)
              field_up(ip11,ip20,ip31,isub)=field_up(ip11,ip20,ip31,isub)+&
     &                   a11a20a31*up(iq)
              field_up(ip10,ip21,ip31,isub)=field_up(ip10,ip21,ip31,isub)+&
     &                   a10a21a31*up(iq)
              field_up(ip11,ip21,ip31,isub)=field_up(ip11,ip21,ip31,isub)+&
     &                   a11a21a31*up(iq)
           ENDDO        ! iq
#endif
        END DO              ! loop over subs
     CASE DEFAULT
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',    &
     &                    'Only Mp4 and Bsp2 are available. Use ppm_rmsh_remesh for other kernels.', &
     &                    __LINE__,info)
     END SELECT         ! kernel type
#else
     !-------------------------------------------------------------------------!
     !  --- 2D --- 
     !-------------------------------------------------------------------------!
     !  loop over subs
     ndata => ppm_cart_mesh(meshid,topoid)%nnodes
     
     DO isub = 1,ppm_nsublist(topoid)
        isubl = ppm_isublist(isub,topoid)
        
        DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
           
           DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
#if __MODE == __VEC
              DO ldn=1,lda
                 field_up(ldn,i,j,isub) = 0.0_mk
              END DO
#else      
              field_up(i,j,isub) = 0.0_mk
#endif
           END DO
        END DO
     END DO
     IF(np.EQ.0) GOTO 9997
     SELECT CASE(kernel)
     CASE(ppm_param_rmsh_kernel_mp4)
        
        !----------------------------------------------------------------------!
        ! M Prime Four
        !----------------------------------------------------------------------!
        DO isub = 1,ppm_nsublist(topoid)
           
           DO ip = 1,store_info(isub)
              
              isubl = ppm_isublist(isub,topoid)
              iq    = list_sub(isub,ip)
              
              x01 = (xp(1,iq)-min_sub(1,isubl,topoid))*dxxi
              x02 = (xp(2,iq)-min_sub(2,isubl,topoid))*dxyi
              
              ip1 = INT(x01)+1
              ip2 = INT(x02)+1
              
              xp1 = x01-AINT(x01)
              xp2 = x02-AINT(x02)
              
              DO jj = -1,2
                 
                 x2 = ABS(xp2 - REAL(jj,mk))
                 
                 IF(x2.LT.1.0_mk) THEN
                    wx2 = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                 ELSE
                    wx2 = 2.0_mk + (-4.0_mk + &
     &                         (2.5_mk - 0.5_mk * x2)*x2)*x2
                 END IF
                 
                 DO ii    = - 1,2
                    
                    x1 = ABS(xp1 - REAL(ii,mk))
                    
                    IF(x1.LT.1.0_MK) THEN
                       wx1 =  1.0_mk - x1**2*(2.5_mk - &
     &                             1.5_mk*x1)
                    ELSE
                       wx1 =  2.0_mk + (-4.0_mk + &
     &                             (2.5_mk - 0.5_mk*x1)*x1)*x1
                    END IF
#if __MODE == __SCA                    
                    field_up(ii+ip1,jj+ip2,isub) &
     &                         = field_up(ii+ip1,jj+ip2,isub) &
     &                                          + wx1*wx2*up(iq)
#else
                    DO ldn=1,lda
                       field_up(ldn,ii+ip1,jj+ip2,isub) &
     &                            = field_up(ldn,ii+ip1,jj+ip2,isub) &
     &                                             + wx1*wx2*up(ldn,iq)
                    ENDDO
#endif                    
                 END DO
              END DO
           END DO
        END DO          ! loop over subs
     CASE DEFAULT
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',    &
     &           'Only Mp4 is available. Use ppm_rmsh_remesh for other kernels.', &
     &           __LINE__,info)
     END SELECT
#endif

9997 CONTINUE
     !-------------------------------------------------------------------------!
     !  Now map the ghosts in order to get consistent values at the border of
     !  the subdomains.
     !-------------------------------------------------------------------------!
     maptype = ppm_param_map_init
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#endif     
     IF (info .NE. 0) GOTO 9999
     maptype = ppm_param_map_ghost_put
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#endif     

     IF (info .NE. 0) GOTO 9999
     maptype = ppm_param_map_push
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#endif     

     IF (info .NE. 0) GOTO 9999
     maptype = ppm_param_map_send
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#endif     

     IF (info .NE. 0) GOTO 9999
     maptype = ppm_param_map_pop
#if   __MODE == __SCA     
     CALL ppm_map_field_ghost(field_up,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#elif __MODE == __VEC     
     CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,ghostsize,maptype, &
     &           info)
#endif     

     IF (info .NE. 0) GOTO 9999
     IF(np.EQ.0) GOTO 9999
     !-------------------------------------------------------------------------!
     ! Deallocation of the arrays....
     !-------------------------------------------------------------------------!
     iopt = ppm_param_dealloc
     ldu(1) = 0
     CALL ppm_alloc(ilist1,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m', &
     &              'pb in ilist1 deallocation',__LINE__,info)
        GOTO 9999
     END IF
     CALL ppm_alloc(ilist2,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m',  &
     &              'pb in ilist2 deallocation',__LINE__,info)
        GOTO 9999
     END IF
     CALL ppm_alloc(store_info,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m',  &
     &              'pb in ilist2 deallocation',__LINE__,info)
        GOTO 9999
     END IF
     CALL ppm_alloc(list_sub,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m',  &
     &              'pb in ilist2 deallocation',__LINE__,info)
        GOTO 9999
     END IF

     !-------------------------------------------------------------------------!
     !  Return 
     !-------------------------------------------------------------------------!
9999 CONTINUE
     CALL substop('ppm_interp_p2m',t0,info)
     RETURN
     
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_ss_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_sv_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_ss_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_sv_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_interp_p2m_dv_3d
#endif
#endif
#endif       

