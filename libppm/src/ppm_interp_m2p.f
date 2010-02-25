      !-------------------------------------------------------------------------
      !     Subroutine   :                   ppm_interp_m2p
      !-------------------------------------------------------------------------
      !     
      !     Purpose      : This subroutine carries out mesh to particle 
      !                    interpolation.
      !                    The 3D versions for lda.LT.5 are explicitly unrolled. 
      !                    All 3D version are unrolled over the kernel support.
      !                    The 2D versions are not.
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
      !                    field_up(:..:)     (F) field onto which to interp
      !      
      !     Remarks      : 
      !     
      !     References   :
      !     
      !     Revisions    :
      !-------------------------------------------------------------------------
      !     $Log: ppm_interp_m2p.f,v $
      !     Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !     CBL version of the PPM library
      !
      !     Revision 1.25  2006/10/20 14:07:45  ivos
      !     2d_vec and 3d_sca cases were commented out and there was a bug in the
      !     argument declaration in 2d_vec. Fixed and put them back in.
      !
      !     Revision 1.24  2006/09/04 18:34:48  pchatela
      !     Fixes and cleanups to make it compile for
      !     release-1-0
      !
      !     Revision 1.23  2006/06/17 18:45:24  michaebe
      !     memory leak. Fixed.
      !
      !     Revision 1.22  2006/02/03 09:34:03  ivos
      !     Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !     local subs in topo_store. Several mapping routines however need the
      !     info about all (global) subs.
      !     Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !     occurrences.
      !
      !     Revision 1.21  2005/08/22 11:42:56  michaebe
      !     removed include of mpif.h
      !
      !     Revision 1.20  2005/07/25 00:29:58  ivos
      !     removed unused USE ppm_module_map_field.
      !
      !     Revision 1.19  2005/06/30 18:16:49  michaebe
      !     fixed a bug in the 1.14 1.15 bugfix
      !
      !     Revision 1.18  2005/06/30 08:56:19  pchatela
      !     Fixed off-by-one index bug in Bsp2 (classic mistake by a C guy...).
      !
      !     Revision 1.17  2005/06/29 09:21:24  pchatela
      !     Screwed up the log in previous version, then fixed it...
      !
      !     Revision 1.16  2005/06/29 09:10:03  pchatela
      !     Added B-Spline 2 with Simone s fix. It could be optimized, by 
      !     skipping a couple of variable definitions in the BSp2 scheme.
      !
      !     Revision 1.15  2005/06/28 16:59:47  ivos
      !     Slight code optimizations (save 3 additions and 1 INT per point),
      !     includes ifort bug fix as reported by Simone.
      !
      !     Revision 1.14  2005/06/23 21:03:56  ivos
      !     Added unrolled versions for all vector cases. Added 3D scalar case 
      !     and 2D vector case. All 3D cases are unrolled over the kernel.
      !     Changed to use ppm_check_topoid and ppm_check_meshid.
      !
      !     Revision 1.13  2005/06/06 23:54:00  michaebe
      !     mandatory post-commit cosmetics
      !
      !     Revision 1.12  2005/06/06 23:49:44  michaebe
      !     updated header and changed interface (!!)
      !
      !     Revision 1.11  2005/05/29 20:14:57  michaebe
      !     some fixes for the 2D case
      !
      !     Revision 1.10  2005/02/16 02:12:01  michaebe
      !     added inclusion of ppm_define.h
      !
      !     Revision 1.9  2005/02/16 01:24:52  michaebe
      !     added sxf90 pragmas
      !
      !     Revision 1.8  2005/02/15 23:16:00  michaebe
      !     took if(lda) out of the loop
      !
      !     Revision 1.7  2005/02/15 20:24:23  michaebe
      !     added nsublist=1 optimization
      !
      !     Revision 1.6  2005/01/17 11:17:28  michaebe
      !     further optimization
      !
      !     Revision 1.5  2004/12/07 18:11:39  michaebe
      !     floating point arithmetic mania resolved min_phys -> min_sub..
      !
      !     Revision 1.4  2004/11/10 09:12:07  michaebe
      !     replaced floor by aint
      !
      !     Revision 1.3  2004/11/09 10:06:13  ivos
      !     CVS fix.
      !
      !     Revision 1.1  2004/11/02 12:51:46  michaebe
      !     inimp
      !
      !-------------------------------------------------------------------------
      !     Parallel Particle Mesh Library (PPM)
      !     Institute of Computational Science
      !     ETH Zentrum, Hirschengraben 84
      !     CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------


#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_interp_m2p_ss_2d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_interp_m2p_ds_2d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_interp_m2p_sv_2d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_interp_m2p_dv_2d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_interp_m2p_ss_3d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_interp_m2p_ds_3d(xp,Np,up,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_interp_m2p_sv_3d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_interp_m2p_dv_3d(xp,Np,up,lda,topo_id,mesh_id,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif
#endif       

      !-------------------------------------------------------------------------
      !  INCLUDES
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_data
      USE ppm_module_data_rmsh
      USE ppm_module_data_mesh
      USE ppm_module_write
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      IMPLICIT NONE
            
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      ! Arguments     
      !-------------------------------------------------------------------------
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER        :: up
      REAL(MK)                                         :: tup
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#endif
#elif __MODE == __VEC
      INTEGER                         , INTENT(in)     :: lda
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: up
      REAL(MK) , DIMENSION(lda)                        :: tup
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_up
#endif     
#endif     
      REAL(MK), DIMENSION(:,:)       , POINTER       :: xp
      INTEGER , DIMENSION(:  )       , INTENT(in)    :: ghostsize
      INTEGER                        , INTENT(IN   ) :: Np
      INTEGER                        , INTENT(IN   ) :: topo_id,mesh_id
      INTEGER                        , INTENT(IN   ) :: kernel
      INTEGER                        , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      ! Local variables 
      !-------------------------------------------------------------------------
      INTEGER,  DIMENSION(:,:)     , POINTER :: istart
      INTEGER,  DIMENSION(:)       , POINTER :: ilist1,ilist2
      REAL(MK), DIMENSION(:,:)     , POINTER :: min_phys,max_phys
      REAL   ,  DIMENSION(ppm_dim)           :: dxi,dx
      REAL   ,  DIMENSION(ppm_dim)           :: len_phys
      REAL(MK)                               :: x1,x2,x3,epsilon
      INTEGER                                :: kernel_support
      INTEGER,  DIMENSION(ppm_dim+2)         :: ldu,ldl
      INTEGER,  DIMENSION(ppm_dim)           :: Nc
      INTEGER                                :: i,j,k,ii,jj,kk,iidec,maptype
      INTEGER                                :: jjdec,nb_sub,npart,ipart
      INTEGER                                :: kkdec,ip1,nbpt_z,nlist1
      INTEGER                                :: ip2,ip3,nbpt_x,nbpt_y,iface
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
      CHARACTER(len=256)                     :: msg
      !-------------------------------------------------------------------------
      !  Variables for unrolled versions
      !-------------------------------------------------------------------------
      REAL(mk) :: x10,x11,x12,x13,x20,x21,x22,x23,x30,x31,x32,x33
      REAL(mk) :: a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33
      INTEGER  :: ip10,ip11,ip12,ip13,ip20,ip21,ip22,ip23,ip30,ip31,ip32,ip33,ldn

#if   __DIME == __2D
      REAL(mk) :: a10a20
      REAL(mk) :: a10a21
      REAL(mk) :: a11a20
      REAL(mk) :: a11a21
#elif __DIME == __3D
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
#endif


     
      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_interp_m2p',t0,info)

      !-------------------------------------------------------------------------
      !  Some compilers have problems initializing constants that are declared
      ! in the module therefore we do it here:
      !-------------------------------------------------------------------------
      ppm_rmsh_kernelsize = (/1,2,2,4/)
     
      dim = ppm_dim
      internal_weights = .FALSE.

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN 
         IF (.NOT. ppm_initialized) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_ppm_noinit,'ppm_interp_m2p',  &
     &               'Please call ppm_init first!',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Np.GT.0) THEN
            IF (SIZE(xp,2).LT.Np) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                    'not enough particles contained in xp',__LINE__,info)
               GOTO 9999
            ENDIF
            IF (SIZE(xp,1).LT.dim) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                    'leading dimension of xp insufficient',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
         IF (Np.LE.0) THEN
            IF (Np.LT.0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                       'particles must be specified',__LINE__,info)
               GOTO 9999
            END IF
            GOTO 9999
         END IF
         IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                       'wrong kernel definition',__LINE__,info)
            GOTO 9999
         END IF
         kernel_support = ppm_rmsh_kernelsize(kernel)*2
         IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) & 
     &                 .OR.(kernel_support.EQ.6))) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                       'wrong kernel support',__LINE__,info)
            GOTO 9999
         END IF
      END IF
     
      !-------------------------------------------------------------------------
      !  Check meshid and topoid validity
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_check_topoid(ppm_param_id_user,topo_id,lok,info)
         IF (.NOT.lok) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                   'topo_id is invalid!',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF

      !-------------------------------------------------------------------------
      !  If theres nothing to do, do nothing
      !-------------------------------------------------------------------------
      IF(Np.EQ.0) GOTO 9999
     
      !-------------------------------------------------------------------------
      !  Get the internal topoid
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)
      IF(ppm_debug.GT.0) THEN
         CALL ppm_check_meshid(ppm_param_id_user,mesh_id,topoid,lok,info)
         IF (.NOT.lok) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                   'mesh_id is invalid!',__LINE__,info)
            GOTO 9999
         END IF
      END IF

      !-------------------------------------------------------------------------
      !  Get the internal meshid
      !-------------------------------------------------------------------------
      meshid = ppm_meshid(topoid)%internal(mesh_id)
     
      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------
      istart => ppm_cart_mesh(meshid,topoid)%istart
     
      !-------------------------------------------------------------------------
      !  Assignment of the useful arrays/scalar
      !-------------------------------------------------------------------------
      Nm(1:dim) = ppm_cart_mesh(meshid,topoid)%Nm
      bcdef(1:(2*dim)) = ppm_bcdef(1:(2*dim),topoid)
      nsubs = ppm_nsublist(topoid)
    
      !-------------------------------------------------------------------------
      !  Allocate memory for particle lists
      !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
      !  or a finalize routine.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Np
      IF(nsubs.NE.1) CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                       'particle list 1 ILIST1',__LINE__,info)
         GOTO 9999
      END IF
      IF(nsubs.NE.1) CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                       'particle list 2 ILIST2',__LINE__,info)
         GOTO 9999
      END IF

      iopt   = ppm_param_alloc_fit
      ldu(1) = nsubs
      IF(nsubs.NE.1) CALL ppm_alloc(store_info,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                       'store_info allocation : problem',__LINE__,info)
         GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Mesh spacing
      !-------------------------------------------------------------------------
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
      END DO
      dxi = 1.0_mk/dx

      !-------------------------------------------------------------------------
      !  Initialize the particle list
      !-------------------------------------------------------------------------
      IF(nsubs.NE.1) THEN
         nlist1     = 0
         store_info = 0
         DO ipart=1,Np
            nlist1         = nlist1 + 1
            ilist1(nlist1) = ipart
         END DO
      END IF
#if   __KIND == __SINGLE_PRECISION
      myeps = ppm_myepss
      min_sub => ppm_min_subs
      max_sub => ppm_max_subs
#elif __KIND == __DOUBLE_PRECISION
      myeps = ppm_myepsd
      min_sub => ppm_min_subd
      max_sub => ppm_max_subd
#endif

      !-------------------------------------------------------------------------
      !  Loop over the subdomains (since the first domains are most likely
      !  to be empty, we look backwards to reduce the number of elements in
      !  nlist2 as fast as possible)
      !-------------------------------------------------------------------------
      IF(nsubs.NE.1) THEN
         DO idom = ppm_nsublist(topoid),1,-1   
            idoml = ppm_isublist(idom,topoid)
            !-------------------------------------------------------------------
            !  Loop over the remaining particles 
            !-------------------------------------------------------------------
            nlist2 = 0
            npart = 0
            DO i=1,nlist1
               ipart = ilist1(i)

               !----------------------------------------------------------------
               !  If the particle is inside the current subdomain, assign it
               !----------------------------------------------------------------
               IF(ppm_dim.EQ.3) THEN
                  !-------------------------------------------------------------
                  !  The particle is in the closure of the subdomain
                  !-------------------------------------------------------------
                  IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
     &                       xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                       xp(3,ipart).GE.min_sub(3,idoml,topoid) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml,topoid) .AND. &
     &                       xp(3,ipart).LE.max_sub(3,idoml,topoid) ) ) THEN

                     IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(3,ipart).LT.max_sub(3,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(6,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

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
     &                       xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml,topoid) ) ) THEN                    
                     IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                        npart = npart + 1
                        store_info(idom) = npart
                     ELSE
                        nlist2         = nlist2 + 1
                        ilist2(nlist2) = ipart
                     END IF
                  ELSE
                     nlist2         = nlist2 + 1
                     ilist2(nlist2) = ipart
                  END IF
               END IF
            END DO ! end loop over remaining parts in ilist1
            !-------------------------------------------------------------------
            !  Copy the lists (well, only if nlist2 changed - decreased)
            !-------------------------------------------------------------------
            IF (nlist2.NE.nlist1) THEN
               nlist1 = nlist2
               DO i=1,nlist1
                  ilist1(i) = ilist2(i)
               END DO
            ENDIF

            !-------------------------------------------------------------------
            !  Exit if the list is empty
            !-------------------------------------------------------------------
            IF (nlist1.EQ.0) EXIT
         END DO ! end loop over subs
         !----------------------------------------------------------------------
         !  Check that we sold all the particles 
         !----------------------------------------------------------------------
         IF (nlist2.GT.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_part_unass,'ppm_interp_m2p',  &
     &                        'MAJOR PROBLEM',__LINE__,info)
            GOTO 9999
         END IF

         !----------------------------------------------------------------------
         !  Whats the maximum number of particles per subdomain
         !----------------------------------------------------------------------
         max_partnumber = 0
         DO idom=1,ppm_nsublist(topoid)
            IF(store_info(idom).GE.max_partnumber) THEN
               max_partnumber = store_info(idom)
            END IF
         END DO
         iopt   = ppm_param_alloc_fit
         ldu(1) = ppm_nsublist(topoid)
         ldu(2) = max_partnumber

         !----------------------------------------------------------------------
         !  Allocate particle list
         !----------------------------------------------------------------------
         CALL ppm_alloc(list_sub,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                         'problem in internal allocation',__LINE__,info)
            GOTO 9999
         END IF

         list_sub = 0
         !----------------------------------------------------------------------
         !  Initialize the particle list
         !----------------------------------------------------------------------
         nlist1     = 0
         DO ipart=1,Np
            nlist1         = nlist1 + 1
            ilist1(nlist1) = ipart
         ENDDO

         !----------------------------------------------------------------------
         !  Loop over the subdomains (since the first domains are most likely
         !  to be empty, we look backwards to reduce the number of elements in
         !  nlist2 as fast as possible)
         !----------------------------------------------------------------------
         DO idom = ppm_nsublist(topoid),1,-1
            idoml = ppm_isublist(idom,topoid)
            !-------------------------------------------------------------------
            !  loop over the remaining particles 
            !-------------------------------------------------------------------
            nlist2 = 0
            npart = 0
            DO i=1,nlist1
               ipart = ilist1(i)
               !----------------------------------------------------------------
               !  If the particle is inside the current subdomain, assign it
               !----------------------------------------------------------------
               IF(ppm_dim.EQ.3) THEN
                  IF( ( xp(1,ipart).GE.min_sub(1,idoml,topoid) .AND. &
     &                       xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                       xp(3,ipart).GE.min_sub(3,idoml,topoid) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml,topoid) .AND. &
     &                       xp(3,ipart).LE.max_sub(3,idoml,topoid) ) ) THEN

                     IF(   (xp(1,ipart).LT.max_sub(1,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(3,ipart).LT.max_sub(3,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(6,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

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
     &                       xp(2,ipart).GE.min_sub(2,idoml,topoid) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml,topoid) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml,topoid) ) ) THEN
                     IF(   (xp(1,ipart).LT.max_sub(1,idom,topoid) .OR.  &
     &                           (ppm_subs_bc(2,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml,topoid) .OR.  &
     &                           (ppm_subs_bc(4,idoml,topoid).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                        npart = npart + 1
                        list_sub(idom,npart) = ipart
                     ELSE
                        nlist2         = nlist2 + 1
                        ilist2(nlist2) = ipart
                     END IF
                  ELSE
                     nlist2         = nlist2 + 1
                     ilist2(nlist2) = ipart
                  END IF
               END IF
            END DO ! end loop over subs
            !-------------------------------------------------------------------
            !  Copy the lists (well, only if nlist2 changed - decreased)
            !-------------------------------------------------------------------
            IF (nlist2.NE.nlist1) THEN
               nlist1 = nlist2
               DO i=1,nlist1
                  ilist1(i) = ilist2(i)
               END DO
            END IF

            !-------------------------------------------------------------------
            !  Exit if the list is empty
            !-------------------------------------------------------------------
            IF (nlist1.EQ.0) EXIT
         END DO

         !----------------------------------------------------------------------
         !  Check that we sold all the particles 
         !----------------------------------------------------------------------
         IF (nlist2.GT.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_part_unass,'ppm_interp_m2p_3d',  &
     &                        'MAJOR PROBLEM',__LINE__,info)
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  Allocate and alias the weights if we need them.
         !----------------------------------------------------------------------
         max_partnumber = 0
         DO idom = 1,ppm_nsublist(topoid)
            IF(store_info(idom).GE.max_partnumber) THEN
               max_partnumber = store_info(idom)  
            END IF
         END DO
      ELSE ! now the case of one subdomain on this processor
         max_partnumber = np
         iopt = ppm_param_alloc_fit
         ldu(1) = 1
         CALL ppm_alloc(store_info,ldu,iopt,info)
         IF(info.NE.0) THEN
            CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                         'problem in internal allocation',__LINE__,info)
            GOTO 9999
         END IF
         store_info(1) = max_partnumber
         ldu(1) = 1
         ldu(2) = max_partnumber
         CALL ppm_alloc(list_sub,ldu,iopt,info)
         IF(info.NE.0) THEN
            CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                         'problem in internal allocation',__LINE__,info)
            GOTO 9999
         END IF
         DO i=1,max_partnumber
            list_sub(1,i) = i
         END DO
      END IF


      !-------------------------------------------------------------------------
      !  Beginning of the computation
      !-------------------------------------------------------------------------
      ! Particle quantity reset
#if   __MODE == __SCA
      DO ip=1,np
         up(ip) = 0.0_mk
      END DO
#elif __MODE == __VEC
      IF(lda.EQ.3) THEN
         DO ip=1,np
            up(1,ip) = 0.0_mk
            up(2,ip) = 0.0_mk
            up(3,ip) = 0.0_mk
         END DO
      ELSEIF(lda.EQ.5) THEN
         DO ip=1,np
            up(1,ip) = 0.0_mk
            up(2,ip) = 0.0_mk
            up(3,ip) = 0.0_mk
            up(4,ip) = 0.0_mk
            up(5,ip) = 0.0_mk
         END DO
      ELSE
         up = 0.0_mk
      END IF
#endif
      SELECT CASE(kernel)

      CASE(ppm_param_rmsh_kernel_mp4)

         !----------------------------------------------------------------------
         !  M Prime Four
         !----------------------------------------------------------------------
         !  loop over subs
         DO isub = 1,ppm_nsublist(topoid)
#if __DIME == __2D
            !-------------------------------------------------------------------
            !  --- 2D ---
            !-------------------------------------------------------------------
#if __MODE == __SCA
            isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)
               x0(1) = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
               x0(2) = (xp(2,iq)-min_phys(2,topoid))*dxi(2)

               ip1 = INT(x0(1))+2-istart(1,isubl)
               ip2 = INT(x0(2))+2-istart(2,isubl)

               xp1 = x0(1)-AINT(x0(1))
               xp2 = x0(2)-AINT(x0(2))

               DO jj = -1,2
                  x2 = ABS(xp2 - REAL(jj,mk))
                  IF(x2.LT.1.0_mk) THEN
                     wx2 = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                  ELSE
                     wx2 = 2.0_mk + (-4.0_mk + &
     &                          (2.5_mk - 0.5_mk * x2)*x2)*x2
                  END IF

                  DO ii    = - 1,2
                     x1 = ABS(xp1 - REAL(ii,mk))
                     IF(x1.LT.1.0_MK) THEN     
                        wx1 =  1.0_mk - x1**2*(2.5_mk - &
     &                              1.5_mk*x1)
                     ELSE
                        wx1 =  2.0_mk + (-4.0_mk + &
     &                              (2.5_mk - 0.5_mk*x1)*x1)*x1
                     END IF
                     up(iq) = up(iq) + wx1*wx2*field_up(ii+ip1,jj+ip2,isub)
                  END DO
               END DO
            END DO ! end loop over particles in the current subdomain
#elif __MODE == __VEC
            ! This will only vectorize over lda
            isubl = ppm_isublist(isub,topoid)
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)
               x0(1) = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
               x0(2) = (xp(2,iq)-min_phys(2,topoid))*dxi(2)

               ip1 = INT(x0(1))+2-istart(1,isubl)
               ip2 = INT(x0(2))+2-istart(2,isubl)

               xp1 = x0(1)-AINT(x0(1))
               xp2 = x0(2)-AINT(x0(2))

               DO jj = -1,2
                  x2 = ABS(xp2 - REAL(jj,mk))
                  IF(x2.LT.1.0_mk) THEN
                     wx2 = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                  ELSE
                     wx2 = 2.0_mk + (-4.0_mk + &
     &                          (2.5_mk - 0.5_mk * x2)*x2)*x2
                  END IF

                  DO ii    = - 1,2
                     x1 = ABS(xp1 - REAL(ii,mk))
                     IF(x1.LT.1.0_MK) THEN     
                        wx1 =  1.0_mk - x1**2*(2.5_mk - &
     &                              1.5_mk*x1)
                     ELSE
                        wx1 =  2.0_mk + (-4.0_mk + &
     &                              (2.5_mk - 0.5_mk*x1)*x1)*x1
                     END IF
                     DO ldn=1,lda
                        up(ldn,iq) = up(ldn,iq) + wx1*wx2*     &
     &                                                   field_up(ldn,ii+ip1,jj+ip2,isub)
                     END DO
                  END DO
               END DO
            END DO ! end loop over particles in the current subdomain
#endif           
#else
            !-------------------------------------------------------------------
            !  --- 3D ---
            !-------------------------------------------------------------------
#if   __MODE == __SCA
            isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)

               x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

               ip10 = INT(x0(1))
               ip20 = INT(x0(2))
               ip30 = INT(x0(3))

               ip11 = ip10 + 1
               ip21 = ip20 + 1
               ip31 = ip30 + 1

               ip12 = ip11 + 1
               ip22 = ip21 + 1
               ip32 = ip31 + 1

               ip13 = ip11 + 2
               ip23 = ip21 + 2
               ip33 = ip31 + 2

               xp1 = x0(1)-REAL(ip10,mk)
               xp2 = x0(2)-REAL(ip20,mk)
               xp3 = x0(3)-REAL(ip30,mk)

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

               up(iq) = up(iq) + &
     &                     a10a20a30*field_up(ip10,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a20a31*field_up(ip10,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a20a32*field_up(ip10,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a20a33*field_up(ip10,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a10a21a30*field_up(ip10,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a21a31*field_up(ip10,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a21a32*field_up(ip10,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a21a33*field_up(ip10,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a10a22a30*field_up(ip10,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a22a31*field_up(ip10,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a22a32*field_up(ip10,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a22a33*field_up(ip10,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a10a23a30*field_up(ip10,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a23a31*field_up(ip10,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a23a32*field_up(ip10,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a23a33*field_up(ip10,ip23,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a20a30*field_up(ip11,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a20a31*field_up(ip11,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a20a32*field_up(ip11,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a20a33*field_up(ip11,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a21a30*field_up(ip11,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a21a31*field_up(ip11,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a21a32*field_up(ip11,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a21a33*field_up(ip11,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a22a30*field_up(ip11,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a22a31*field_up(ip11,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a22a32*field_up(ip11,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a22a33*field_up(ip11,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a23a30*field_up(ip11,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a23a31*field_up(ip11,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a23a32*field_up(ip11,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a23a33*field_up(ip11,ip23,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a20a30*field_up(ip12,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a20a31*field_up(ip12,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a20a32*field_up(ip12,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a20a33*field_up(ip12,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a21a30*field_up(ip12,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a21a31*field_up(ip12,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a21a32*field_up(ip12,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a21a33*field_up(ip12,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a22a30*field_up(ip12,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a22a31*field_up(ip12,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a22a32*field_up(ip12,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a22a33*field_up(ip12,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a23a30*field_up(ip12,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a23a31*field_up(ip12,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a23a32*field_up(ip12,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a23a33*field_up(ip12,ip23,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a20a30*field_up(ip13,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a20a31*field_up(ip13,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a20a32*field_up(ip13,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a20a33*field_up(ip13,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a21a30*field_up(ip13,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a21a31*field_up(ip13,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a21a32*field_up(ip13,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a21a33*field_up(ip13,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a22a30*field_up(ip13,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a22a31*field_up(ip13,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a22a32*field_up(ip13,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a22a33*field_up(ip13,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a23a30*field_up(ip13,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a23a31*field_up(ip13,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a23a32*field_up(ip13,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a23a33*field_up(ip13,ip23,ip33,isub)
            END DO ! end loop over particles in the current subdomain
#elif __MODE == __VEC
            !-------------------------------------------------------------------
            !  Unrolled version for 1-vectors
            !-------------------------------------------------------------------
            IF(lda.EQ.1) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a32*field_up(1,ip10,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a33*field_up(1,ip10,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a32*field_up(1,ip10,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a33*field_up(1,ip10,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a30*field_up(1,ip10,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a31*field_up(1,ip10,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a32*field_up(1,ip10,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a33*field_up(1,ip10,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a30*field_up(1,ip10,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a31*field_up(1,ip10,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a32*field_up(1,ip10,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a33*field_up(1,ip10,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a32*field_up(1,ip11,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a33*field_up(1,ip11,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a32*field_up(1,ip11,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a33*field_up(1,ip11,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a30*field_up(1,ip11,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a31*field_up(1,ip11,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a32*field_up(1,ip11,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a33*field_up(1,ip11,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a30*field_up(1,ip11,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a31*field_up(1,ip11,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a32*field_up(1,ip11,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a33*field_up(1,ip11,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a30*field_up(1,ip12,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a31*field_up(1,ip12,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a32*field_up(1,ip12,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a33*field_up(1,ip12,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a30*field_up(1,ip12,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a31*field_up(1,ip12,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a32*field_up(1,ip12,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a33*field_up(1,ip12,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a30*field_up(1,ip12,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a31*field_up(1,ip12,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a32*field_up(1,ip12,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a33*field_up(1,ip12,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a30*field_up(1,ip12,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a31*field_up(1,ip12,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a32*field_up(1,ip12,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a33*field_up(1,ip12,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a30*field_up(1,ip13,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a31*field_up(1,ip13,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a32*field_up(1,ip13,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a33*field_up(1,ip13,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a30*field_up(1,ip13,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a31*field_up(1,ip13,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a32*field_up(1,ip13,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a33*field_up(1,ip13,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a30*field_up(1,ip13,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a31*field_up(1,ip13,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a32*field_up(1,ip13,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a33*field_up(1,ip13,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a30*field_up(1,ip13,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a31*field_up(1,ip13,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a32*field_up(1,ip13,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a33*field_up(1,ip13,ip23,ip33,isub)

               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 2-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.2) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a32*field_up(1,ip10,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a33*field_up(1,ip10,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a32*field_up(1,ip10,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a33*field_up(1,ip10,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a30*field_up(1,ip10,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a31*field_up(1,ip10,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a32*field_up(1,ip10,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a33*field_up(1,ip10,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a30*field_up(1,ip10,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a31*field_up(1,ip10,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a32*field_up(1,ip10,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a33*field_up(1,ip10,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a32*field_up(1,ip11,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a33*field_up(1,ip11,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a32*field_up(1,ip11,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a33*field_up(1,ip11,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a30*field_up(1,ip11,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a31*field_up(1,ip11,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a32*field_up(1,ip11,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a33*field_up(1,ip11,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a30*field_up(1,ip11,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a31*field_up(1,ip11,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a32*field_up(1,ip11,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a33*field_up(1,ip11,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a30*field_up(1,ip12,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a31*field_up(1,ip12,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a32*field_up(1,ip12,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a33*field_up(1,ip12,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a30*field_up(1,ip12,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a31*field_up(1,ip12,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a32*field_up(1,ip12,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a33*field_up(1,ip12,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a30*field_up(1,ip12,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a31*field_up(1,ip12,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a32*field_up(1,ip12,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a33*field_up(1,ip12,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a30*field_up(1,ip12,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a31*field_up(1,ip12,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a32*field_up(1,ip12,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a33*field_up(1,ip12,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a30*field_up(1,ip13,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a31*field_up(1,ip13,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a32*field_up(1,ip13,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a33*field_up(1,ip13,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a30*field_up(1,ip13,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a31*field_up(1,ip13,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a32*field_up(1,ip13,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a33*field_up(1,ip13,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a30*field_up(1,ip13,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a31*field_up(1,ip13,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a32*field_up(1,ip13,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a33*field_up(1,ip13,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a30*field_up(1,ip13,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a31*field_up(1,ip13,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a32*field_up(1,ip13,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a33*field_up(1,ip13,ip23,ip33,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a32*field_up(2,ip10,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a33*field_up(2,ip10,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a32*field_up(2,ip10,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a33*field_up(2,ip10,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a30*field_up(2,ip10,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a31*field_up(2,ip10,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a32*field_up(2,ip10,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a33*field_up(2,ip10,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a30*field_up(2,ip10,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a31*field_up(2,ip10,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a32*field_up(2,ip10,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a33*field_up(2,ip10,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a32*field_up(2,ip11,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a33*field_up(2,ip11,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a32*field_up(2,ip11,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a33*field_up(2,ip11,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a30*field_up(2,ip11,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a31*field_up(2,ip11,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a32*field_up(2,ip11,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a33*field_up(2,ip11,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a30*field_up(2,ip11,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a31*field_up(2,ip11,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a32*field_up(2,ip11,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a33*field_up(2,ip11,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a30*field_up(2,ip12,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a31*field_up(2,ip12,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a32*field_up(2,ip12,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a33*field_up(2,ip12,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a30*field_up(2,ip12,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a31*field_up(2,ip12,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a32*field_up(2,ip12,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a33*field_up(2,ip12,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a30*field_up(2,ip12,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a31*field_up(2,ip12,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a32*field_up(2,ip12,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a33*field_up(2,ip12,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a30*field_up(2,ip12,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a31*field_up(2,ip12,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a32*field_up(2,ip12,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a33*field_up(2,ip12,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a30*field_up(2,ip13,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a31*field_up(2,ip13,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a32*field_up(2,ip13,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a33*field_up(2,ip13,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a30*field_up(2,ip13,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a31*field_up(2,ip13,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a32*field_up(2,ip13,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a33*field_up(2,ip13,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a30*field_up(2,ip13,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a31*field_up(2,ip13,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a32*field_up(2,ip13,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a33*field_up(2,ip13,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a30*field_up(2,ip13,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a31*field_up(2,ip13,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a32*field_up(2,ip13,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a33*field_up(2,ip13,ip23,ip33,isub)


               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 3-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.3) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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


                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a32*field_up(1,ip10,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a33*field_up(1,ip10,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a32*field_up(1,ip10,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a33*field_up(1,ip10,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a30*field_up(1,ip10,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a31*field_up(1,ip10,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a32*field_up(1,ip10,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a33*field_up(1,ip10,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a30*field_up(1,ip10,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a31*field_up(1,ip10,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a32*field_up(1,ip10,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a33*field_up(1,ip10,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a32*field_up(1,ip11,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a33*field_up(1,ip11,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a32*field_up(1,ip11,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a33*field_up(1,ip11,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a30*field_up(1,ip11,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a31*field_up(1,ip11,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a32*field_up(1,ip11,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a33*field_up(1,ip11,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a30*field_up(1,ip11,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a31*field_up(1,ip11,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a32*field_up(1,ip11,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a33*field_up(1,ip11,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a30*field_up(1,ip12,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a31*field_up(1,ip12,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a32*field_up(1,ip12,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a33*field_up(1,ip12,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a30*field_up(1,ip12,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a31*field_up(1,ip12,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a32*field_up(1,ip12,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a33*field_up(1,ip12,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a30*field_up(1,ip12,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a31*field_up(1,ip12,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a32*field_up(1,ip12,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a33*field_up(1,ip12,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a30*field_up(1,ip12,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a31*field_up(1,ip12,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a32*field_up(1,ip12,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a33*field_up(1,ip12,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a30*field_up(1,ip13,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a31*field_up(1,ip13,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a32*field_up(1,ip13,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a33*field_up(1,ip13,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a30*field_up(1,ip13,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a31*field_up(1,ip13,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a32*field_up(1,ip13,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a33*field_up(1,ip13,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a30*field_up(1,ip13,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a31*field_up(1,ip13,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a32*field_up(1,ip13,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a33*field_up(1,ip13,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a30*field_up(1,ip13,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a31*field_up(1,ip13,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a32*field_up(1,ip13,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a33*field_up(1,ip13,ip23,ip33,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a32*field_up(2,ip10,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a33*field_up(2,ip10,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a32*field_up(2,ip10,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a33*field_up(2,ip10,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a30*field_up(2,ip10,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a31*field_up(2,ip10,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a32*field_up(2,ip10,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a33*field_up(2,ip10,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a30*field_up(2,ip10,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a31*field_up(2,ip10,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a32*field_up(2,ip10,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a33*field_up(2,ip10,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a32*field_up(2,ip11,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a33*field_up(2,ip11,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a32*field_up(2,ip11,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a33*field_up(2,ip11,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a30*field_up(2,ip11,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a31*field_up(2,ip11,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a32*field_up(2,ip11,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a33*field_up(2,ip11,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a30*field_up(2,ip11,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a31*field_up(2,ip11,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a32*field_up(2,ip11,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a33*field_up(2,ip11,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a30*field_up(2,ip12,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a31*field_up(2,ip12,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a32*field_up(2,ip12,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a33*field_up(2,ip12,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a30*field_up(2,ip12,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a31*field_up(2,ip12,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a32*field_up(2,ip12,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a33*field_up(2,ip12,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a30*field_up(2,ip12,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a31*field_up(2,ip12,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a32*field_up(2,ip12,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a33*field_up(2,ip12,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a30*field_up(2,ip12,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a31*field_up(2,ip12,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a32*field_up(2,ip12,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a33*field_up(2,ip12,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a30*field_up(2,ip13,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a31*field_up(2,ip13,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a32*field_up(2,ip13,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a33*field_up(2,ip13,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a30*field_up(2,ip13,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a31*field_up(2,ip13,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a32*field_up(2,ip13,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a33*field_up(2,ip13,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a30*field_up(2,ip13,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a31*field_up(2,ip13,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a32*field_up(2,ip13,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a33*field_up(2,ip13,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a30*field_up(2,ip13,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a31*field_up(2,ip13,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a32*field_up(2,ip13,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a33*field_up(2,ip13,ip23,ip33,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a20a30*field_up(3,ip10,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a31*field_up(3,ip10,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a32*field_up(3,ip10,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a33*field_up(3,ip10,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a30*field_up(3,ip10,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a31*field_up(3,ip10,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a32*field_up(3,ip10,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a33*field_up(3,ip10,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a30*field_up(3,ip10,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a31*field_up(3,ip10,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a32*field_up(3,ip10,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a33*field_up(3,ip10,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a30*field_up(3,ip10,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a31*field_up(3,ip10,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a32*field_up(3,ip10,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a33*field_up(3,ip10,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a30*field_up(3,ip11,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a31*field_up(3,ip11,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a32*field_up(3,ip11,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a33*field_up(3,ip11,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a30*field_up(3,ip11,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a31*field_up(3,ip11,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a32*field_up(3,ip11,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a33*field_up(3,ip11,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a30*field_up(3,ip11,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a31*field_up(3,ip11,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a32*field_up(3,ip11,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a33*field_up(3,ip11,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a30*field_up(3,ip11,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a31*field_up(3,ip11,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a32*field_up(3,ip11,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a33*field_up(3,ip11,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a30*field_up(3,ip12,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a31*field_up(3,ip12,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a32*field_up(3,ip12,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a33*field_up(3,ip12,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a30*field_up(3,ip12,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a31*field_up(3,ip12,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a32*field_up(3,ip12,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a33*field_up(3,ip12,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a30*field_up(3,ip12,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a31*field_up(3,ip12,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a32*field_up(3,ip12,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a33*field_up(3,ip12,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a30*field_up(3,ip12,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a31*field_up(3,ip12,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a32*field_up(3,ip12,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a33*field_up(3,ip12,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a30*field_up(3,ip13,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a31*field_up(3,ip13,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a32*field_up(3,ip13,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a33*field_up(3,ip13,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a30*field_up(3,ip13,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a31*field_up(3,ip13,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a32*field_up(3,ip13,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a33*field_up(3,ip13,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a30*field_up(3,ip13,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a31*field_up(3,ip13,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a32*field_up(3,ip13,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a33*field_up(3,ip13,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a30*field_up(3,ip13,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a31*field_up(3,ip13,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a32*field_up(3,ip13,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a33*field_up(3,ip13,ip23,ip33,isub)

               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 4-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.4) THEN
               isubl = ppm_isublist(isub,topoid)                 
#ifdef __SXF90
!CDIR NODEP
#endif              
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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


                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a32*field_up(1,ip10,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a33*field_up(1,ip10,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a32*field_up(1,ip10,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a33*field_up(1,ip10,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a30*field_up(1,ip10,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a31*field_up(1,ip10,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a32*field_up(1,ip10,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a33*field_up(1,ip10,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a30*field_up(1,ip10,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a31*field_up(1,ip10,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a32*field_up(1,ip10,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a33*field_up(1,ip10,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a32*field_up(1,ip11,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a33*field_up(1,ip11,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a32*field_up(1,ip11,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a33*field_up(1,ip11,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a30*field_up(1,ip11,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a31*field_up(1,ip11,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a32*field_up(1,ip11,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a33*field_up(1,ip11,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a30*field_up(1,ip11,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a31*field_up(1,ip11,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a32*field_up(1,ip11,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a33*field_up(1,ip11,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a30*field_up(1,ip12,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a31*field_up(1,ip12,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a32*field_up(1,ip12,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a33*field_up(1,ip12,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a30*field_up(1,ip12,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a31*field_up(1,ip12,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a32*field_up(1,ip12,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a33*field_up(1,ip12,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a30*field_up(1,ip12,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a31*field_up(1,ip12,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a32*field_up(1,ip12,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a33*field_up(1,ip12,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a30*field_up(1,ip12,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a31*field_up(1,ip12,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a32*field_up(1,ip12,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a33*field_up(1,ip12,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a30*field_up(1,ip13,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a31*field_up(1,ip13,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a32*field_up(1,ip13,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a33*field_up(1,ip13,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a30*field_up(1,ip13,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a31*field_up(1,ip13,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a32*field_up(1,ip13,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a33*field_up(1,ip13,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a30*field_up(1,ip13,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a31*field_up(1,ip13,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a32*field_up(1,ip13,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a33*field_up(1,ip13,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a30*field_up(1,ip13,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a31*field_up(1,ip13,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a32*field_up(1,ip13,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a33*field_up(1,ip13,ip23,ip33,isub)


                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a32*field_up(2,ip10,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a33*field_up(2,ip10,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a32*field_up(2,ip10,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a33*field_up(2,ip10,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a30*field_up(2,ip10,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a31*field_up(2,ip10,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a32*field_up(2,ip10,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a33*field_up(2,ip10,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a30*field_up(2,ip10,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a31*field_up(2,ip10,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a32*field_up(2,ip10,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a33*field_up(2,ip10,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a32*field_up(2,ip11,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a33*field_up(2,ip11,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a32*field_up(2,ip11,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a33*field_up(2,ip11,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a30*field_up(2,ip11,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a31*field_up(2,ip11,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a32*field_up(2,ip11,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a33*field_up(2,ip11,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a30*field_up(2,ip11,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a31*field_up(2,ip11,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a32*field_up(2,ip11,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a33*field_up(2,ip11,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a30*field_up(2,ip12,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a31*field_up(2,ip12,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a32*field_up(2,ip12,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a33*field_up(2,ip12,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a30*field_up(2,ip12,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a31*field_up(2,ip12,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a32*field_up(2,ip12,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a33*field_up(2,ip12,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a30*field_up(2,ip12,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a31*field_up(2,ip12,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a32*field_up(2,ip12,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a33*field_up(2,ip12,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a30*field_up(2,ip12,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a31*field_up(2,ip12,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a32*field_up(2,ip12,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a33*field_up(2,ip12,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a30*field_up(2,ip13,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a31*field_up(2,ip13,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a32*field_up(2,ip13,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a33*field_up(2,ip13,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a30*field_up(2,ip13,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a31*field_up(2,ip13,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a32*field_up(2,ip13,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a33*field_up(2,ip13,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a30*field_up(2,ip13,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a31*field_up(2,ip13,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a32*field_up(2,ip13,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a33*field_up(2,ip13,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a30*field_up(2,ip13,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a31*field_up(2,ip13,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a32*field_up(2,ip13,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a33*field_up(2,ip13,ip23,ip33,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a20a30*field_up(3,ip10,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a31*field_up(3,ip10,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a32*field_up(3,ip10,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a33*field_up(3,ip10,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a30*field_up(3,ip10,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a31*field_up(3,ip10,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a32*field_up(3,ip10,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a33*field_up(3,ip10,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a30*field_up(3,ip10,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a31*field_up(3,ip10,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a32*field_up(3,ip10,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a33*field_up(3,ip10,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a30*field_up(3,ip10,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a31*field_up(3,ip10,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a32*field_up(3,ip10,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a33*field_up(3,ip10,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a30*field_up(3,ip11,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a31*field_up(3,ip11,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a32*field_up(3,ip11,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a33*field_up(3,ip11,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a30*field_up(3,ip11,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a31*field_up(3,ip11,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a32*field_up(3,ip11,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a33*field_up(3,ip11,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a30*field_up(3,ip11,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a31*field_up(3,ip11,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a32*field_up(3,ip11,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a33*field_up(3,ip11,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a30*field_up(3,ip11,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a31*field_up(3,ip11,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a32*field_up(3,ip11,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a33*field_up(3,ip11,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a30*field_up(3,ip12,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a31*field_up(3,ip12,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a32*field_up(3,ip12,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a33*field_up(3,ip12,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a30*field_up(3,ip12,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a31*field_up(3,ip12,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a32*field_up(3,ip12,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a33*field_up(3,ip12,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a30*field_up(3,ip12,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a31*field_up(3,ip12,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a32*field_up(3,ip12,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a33*field_up(3,ip12,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a30*field_up(3,ip12,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a31*field_up(3,ip12,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a32*field_up(3,ip12,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a33*field_up(3,ip12,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a30*field_up(3,ip13,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a31*field_up(3,ip13,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a32*field_up(3,ip13,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a33*field_up(3,ip13,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a30*field_up(3,ip13,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a31*field_up(3,ip13,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a32*field_up(3,ip13,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a33*field_up(3,ip13,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a30*field_up(3,ip13,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a31*field_up(3,ip13,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a32*field_up(3,ip13,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a33*field_up(3,ip13,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a30*field_up(3,ip13,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a31*field_up(3,ip13,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a32*field_up(3,ip13,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a33*field_up(3,ip13,ip23,ip33,isub)

                  up(4,iq) = up(4,iq) + &
     &                        a10a20a30*field_up(4,ip10,ip20,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a20a31*field_up(4,ip10,ip20,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a20a32*field_up(4,ip10,ip20,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a20a33*field_up(4,ip10,ip20,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a21a30*field_up(4,ip10,ip21,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a21a31*field_up(4,ip10,ip21,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a21a32*field_up(4,ip10,ip21,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a21a33*field_up(4,ip10,ip21,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a22a30*field_up(4,ip10,ip22,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a22a31*field_up(4,ip10,ip22,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a22a32*field_up(4,ip10,ip22,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a22a33*field_up(4,ip10,ip22,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a23a30*field_up(4,ip10,ip23,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a23a31*field_up(4,ip10,ip23,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a23a32*field_up(4,ip10,ip23,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a23a33*field_up(4,ip10,ip23,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a20a30*field_up(4,ip11,ip20,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a20a31*field_up(4,ip11,ip20,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a20a32*field_up(4,ip11,ip20,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a20a33*field_up(4,ip11,ip20,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a21a30*field_up(4,ip11,ip21,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a21a31*field_up(4,ip11,ip21,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a21a32*field_up(4,ip11,ip21,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a21a33*field_up(4,ip11,ip21,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a22a30*field_up(4,ip11,ip22,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a22a31*field_up(4,ip11,ip22,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a22a32*field_up(4,ip11,ip22,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a22a33*field_up(4,ip11,ip22,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a23a30*field_up(4,ip11,ip23,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a23a31*field_up(4,ip11,ip23,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a23a32*field_up(4,ip11,ip23,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a23a33*field_up(4,ip11,ip23,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a20a30*field_up(4,ip12,ip20,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a20a31*field_up(4,ip12,ip20,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a20a32*field_up(4,ip12,ip20,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a20a33*field_up(4,ip12,ip20,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a21a30*field_up(4,ip12,ip21,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a21a31*field_up(4,ip12,ip21,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a21a32*field_up(4,ip12,ip21,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a21a33*field_up(4,ip12,ip21,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a22a30*field_up(4,ip12,ip22,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a22a31*field_up(4,ip12,ip22,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a22a32*field_up(4,ip12,ip22,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a22a33*field_up(4,ip12,ip22,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a23a30*field_up(4,ip12,ip23,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a23a31*field_up(4,ip12,ip23,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a23a32*field_up(4,ip12,ip23,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a12a23a33*field_up(4,ip12,ip23,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a20a30*field_up(4,ip13,ip20,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a20a31*field_up(4,ip13,ip20,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a20a32*field_up(4,ip13,ip20,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a20a33*field_up(4,ip13,ip20,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a21a30*field_up(4,ip13,ip21,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a21a31*field_up(4,ip13,ip21,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a21a32*field_up(4,ip13,ip21,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a21a33*field_up(4,ip13,ip21,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a22a30*field_up(4,ip13,ip22,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a22a31*field_up(4,ip13,ip22,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a22a32*field_up(4,ip13,ip22,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a22a33*field_up(4,ip13,ip22,ip33,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a23a30*field_up(4,ip13,ip23,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a23a31*field_up(4,ip13,ip23,ip31,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a23a32*field_up(4,ip13,ip23,ip32,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a13a23a33*field_up(4,ip13,ip23,ip33,isub)

               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  All other lda are not unrolled. This will vectorize over lda!
               !----------------------------------------------------------------
            ELSE
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a30*field_up(ldn,ip10,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a31*field_up(ldn,ip10,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a32*field_up(ldn,ip10,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a33*field_up(ldn,ip10,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a30*field_up(ldn,ip10,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a31*field_up(ldn,ip10,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a32*field_up(ldn,ip10,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a33*field_up(ldn,ip10,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a30*field_up(ldn,ip10,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a31*field_up(ldn,ip10,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a32*field_up(ldn,ip10,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a33*field_up(ldn,ip10,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a30*field_up(ldn,ip10,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a31*field_up(ldn,ip10,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a32*field_up(ldn,ip10,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a33*field_up(ldn,ip10,ip23,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a30*field_up(ldn,ip11,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a31*field_up(ldn,ip11,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a32*field_up(ldn,ip11,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a33*field_up(ldn,ip11,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a30*field_up(ldn,ip11,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a31*field_up(ldn,ip11,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a32*field_up(ldn,ip11,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a33*field_up(ldn,ip11,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a30*field_up(ldn,ip11,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a31*field_up(ldn,ip11,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a32*field_up(ldn,ip11,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a33*field_up(ldn,ip11,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a30*field_up(ldn,ip11,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a31*field_up(ldn,ip11,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a32*field_up(ldn,ip11,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a33*field_up(ldn,ip11,ip23,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a30*field_up(ldn,ip12,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a31*field_up(ldn,ip12,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a32*field_up(ldn,ip12,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a33*field_up(ldn,ip12,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a30*field_up(ldn,ip12,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a31*field_up(ldn,ip12,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a32*field_up(ldn,ip12,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a33*field_up(ldn,ip12,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a30*field_up(ldn,ip12,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a31*field_up(ldn,ip12,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a32*field_up(ldn,ip12,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a33*field_up(ldn,ip12,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a30*field_up(ldn,ip12,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a31*field_up(ldn,ip12,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a32*field_up(ldn,ip12,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a33*field_up(ldn,ip12,ip23,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a30*field_up(ldn,ip13,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a31*field_up(ldn,ip13,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a32*field_up(ldn,ip13,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a33*field_up(ldn,ip13,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a30*field_up(ldn,ip13,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a31*field_up(ldn,ip13,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a32*field_up(ldn,ip13,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a33*field_up(ldn,ip13,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a30*field_up(ldn,ip13,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a31*field_up(ldn,ip13,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a32*field_up(ldn,ip13,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a33*field_up(ldn,ip13,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a30*field_up(ldn,ip13,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a31*field_up(ldn,ip13,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a32*field_up(ldn,ip13,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a33*field_up(ldn,ip13,ip23,ip33,isub)
                  END DO ! lda
               END DO ! end loop over particles in the current subdomain
            END IF ! unrolled lda cases 
#endif              
#endif              
         END DO ! isub
         !----------------------------------------------------------------------
         !  B-spline 2 (Witch hat)
         !----------------------------------------------------------------------
      CASE(ppm_param_rmsh_kernel_bsp2)
         DO isub = 1,ppm_nsublist(topoid)
#if  __DIME == __2D
            !-------------------------------------------------------------------
            !  --- 2D ---
            !-------------------------------------------------------------------
#if  __MODE == __SCA
            isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)
               
               x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)

               ip10 = INT(x0(1)) + 1
               ip20 = INT(x0(2)) + 1

               ip11 = ip10 + 1
               ip21 = ip20 + 1

               xp1 = x0(1)-REAL(ip10-1,mk)
               xp2 = x0(2)-REAL(ip20-1,mk)

               x10 = xp1
               x11 = x10 - 1.0_mk

               x20 = xp2
               x21 = x20 - 1.0_mk

               a10 = 1.0_mk - x10
               a20 = 1.0_mk - x20

               a11 = 1.0_mk + x11
               a21 = 1.0_mk + x21

               a10a20 = a10*a20
               a10a21 = a10*a21

               a11a20 = a11*a20
               a11a21 = a11*a21

               up(iq) = up(iq) + &
     &                     a10a20*field_up(ip10,ip20,isub)
               up(iq) = up(iq) + &
     &                     a10a21*field_up(ip10,ip21,isub)
               up(iq) = up(iq) + &
     &                     a11a20*field_up(ip11,ip20,isub)
               up(iq) = up(iq) + &
     &                     a11a21*field_up(ip11,ip21,isub)
            END DO  ! end loop over particles in the current subdomain
#elif __MODE == __VEC
            !-------------------------------------------------------------------
            !  Unrolled version for 1-vectors
            !-------------------------------------------------------------------
            IF(lda.EQ.1) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)

                  x10 = xp1
                  x11 = x10 - 1.0_mk

                  x20 = xp2
                  x21 = x20 - 1.0_mk

                  a10 = 1.0_mk - x10
                  a20 = 1.0_mk - x20

                  a11 = 1.0_mk + x11
                  a21 = 1.0_mk + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21

                  up(1,iq) = up(1,iq) + &
     &                        a10a20*field_up(1,ip10,ip20,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21*field_up(1,ip10,ip21,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20*field_up(1,ip11,ip20,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21*field_up(1,ip11,ip21,isub)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 2-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.2) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)

                  x10 = xp1
                  x11 = x10 - 1.0_mk

                  x20 = xp2
                  x21 = x20 - 1.0_mk

                  a10 = 1.0_mk - x10
                  a20 = 1.0_mk - x20

                  a11 = 1.0_mk + x11
                  a21 = 1.0_mk + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21

                  up(1,iq) = up(1,iq) + &
     &                        a10a20*field_up(1,ip10,ip20,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21*field_up(1,ip10,ip21,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20*field_up(1,ip11,ip20,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21*field_up(1,ip11,ip21,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20*field_up(2,ip10,ip20,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21*field_up(2,ip10,ip21,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20*field_up(2,ip11,ip20,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21*field_up(2,ip11,ip21,isub)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 3-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.3) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)

                  x10 = xp1
                  x11 = x10 - 1.0_mk

                  x20 = xp2
                  x21 = x20 - 1.0_mk

                  a10 = 1.0_mk - x10
                  a20 = 1.0_mk - x20

                  a11 = 1.0_mk + x11
                  a21 = 1.0_mk + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21

                  up(1,iq) = up(1,iq) + &
     &                        a10a20*field_up(1,ip10,ip20,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21*field_up(1,ip10,ip21,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20*field_up(1,ip11,ip20,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21*field_up(1,ip11,ip21,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20*field_up(2,ip10,ip20,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21*field_up(2,ip10,ip21,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20*field_up(2,ip11,ip20,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21*field_up(2,ip11,ip21,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a20*field_up(3,ip10,ip20,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21*field_up(3,ip10,ip21,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20*field_up(3,ip11,ip20,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21*field_up(3,ip11,ip21,isub)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  All other lda are not unrolled. This will vectorize over lda!
               !----------------------------------------------------------------
            ELSE
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)

                  x10 = xp1
                  x11 = x10 - 1.0_mk

                  x20 = xp2
                  x21 = x20 - 1.0_mk

                  a10 = 1.0_mk - x10
                  a20 = 1.0_mk - x20

                  a11 = 1.0_mk + x11
                  a21 = 1.0_mk + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21
                  DO ldn=1,lda
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20*field_up(ldn,ip10,ip20,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21*field_up(ldn,ip10,ip21,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20*field_up(ldn,ip11,ip20,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21*field_up(ldn,ip11,ip21,isub)
                  END DO   ! ldn
               END DO ! end loop over particles in the current subdomain
            END IF
#endif
#elif __DIME == __3D
            !-------------------------------------------------------------------
            !  --- 3D ---
            !-------------------------------------------------------------------
#if   __MODE == __SCA
            isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)
               x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)
               
               ip10 = INT(x0(1)) + 1
               ip20 = INT(x0(2)) + 1
               ip30 = INT(x0(3)) + 1

               ip11 = ip10 + 1
               ip21 = ip20 + 1
               ip31 = ip30 + 1

               xp1 = x0(1)-REAL(ip10-1,mk)
               xp2 = x0(2)-REAL(ip20-1,mk)
               xp3 = x0(3)-REAL(ip30-1,mk)

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

               up(iq) = up(iq) + &
     &                     a10a20a30*field_up(ip10,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a20a31*field_up(ip10,ip20,ip31,isub)

               up(iq) = up(iq) + &
     &                     a10a21a30*field_up(ip10,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a21a31*field_up(ip10,ip21,ip31,isub)

               up(iq) = up(iq) + &
     &                     a11a20a30*field_up(ip11,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a20a31*field_up(ip11,ip20,ip31,isub)

               up(iq) = up(iq) + &
     &                     a11a21a30*field_up(ip11,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a21a31*field_up(ip11,ip21,ip31,isub)
            END DO  ! iq
#elif __MODE == __VEC
            !-------------------------------------------------------------------
            !  Unrolled version for 1-vectors
            !-------------------------------------------------------------------
            IF(lda.EQ.1) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)
                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1
                  ip30 = INT(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)
                  xp3 = x0(3)-REAL(ip30-1,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 2-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.2) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif    
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)
                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1
                  ip30 = INT(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)
                  xp3 = x0(3)-REAL(ip30-1,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 3-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.3) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)
                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1
                  ip30 = INT(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)
                  xp3 = x0(3)-REAL(ip30-1,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a20a30*field_up(3,ip10,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a31*field_up(3,ip10,ip20,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a21a30*field_up(3,ip10,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a31*field_up(3,ip10,ip21,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a11a20a30*field_up(3,ip11,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a31*field_up(3,ip11,ip20,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a11a21a30*field_up(3,ip11,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a31*field_up(3,ip11,ip21,ip31,isub)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 4-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.4) THEN              
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif              
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)
                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)

                  ip10 = INT(x0(1)) + 1
                  ip20 = INT(x0(2)) + 1
                  ip30 = INT(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)
                  xp3 = x0(3)-REAL(ip30-1,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)

                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a20a30*field_up(3,ip10,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a31*field_up(3,ip10,ip20,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a21a30*field_up(3,ip10,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a31*field_up(3,ip10,ip21,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a11a20a30*field_up(3,ip11,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a31*field_up(3,ip11,ip20,ip31,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a11a21a30*field_up(3,ip11,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a31*field_up(3,ip11,ip21,ip31,isub)

                  up(4,iq) = up(4,iq) + &
     &                        a10a20a30*field_up(4,ip10,ip20,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a20a31*field_up(4,ip10,ip20,ip31,isub)

                  up(4,iq) = up(4,iq) + &
     &                        a10a21a30*field_up(4,ip10,ip21,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a10a21a31*field_up(4,ip10,ip21,ip31,isub)

                  up(4,iq) = up(4,iq) + &
     &                        a11a20a30*field_up(4,ip11,ip20,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a20a31*field_up(4,ip11,ip20,ip31,isub)

                  up(4,iq) = up(4,iq) + &
     &                        a11a21a30*field_up(4,ip11,ip21,ip30,isub)
                  up(4,iq) = up(4,iq) + &
     &                        a11a21a31*field_up(4,ip11,ip21,ip31,isub)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  All other lda are not unrolled. This will vectorize over lda!
               !----------------------------------------------------------------
            ELSE
               isubl = ppm_isublist(isub,topoid)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)
                  x0(1) = (xp(1,iq)-min_sub(1,isubl,topoid))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl,topoid))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl,topoid))*dxi(3)
                  
                  ip10 = INT(x0(1)) + 1 
                  ip20 = INT(x0(2)) + 1
                  ip30 = INT(x0(3)) + 1
                  
                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1
                  
                  xp1 = x0(1)-REAL(ip10-1,mk)
                  xp2 = x0(2)-REAL(ip20-1,mk)
                  xp3 = x0(3)-REAL(ip30-1,mk)
                  
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
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a30*field_up(ldn,ip10,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a31*field_up(ldn,ip10,ip20,ip31,isub)
                     
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a30*field_up(ldn,ip10,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a31*field_up(ldn,ip10,ip21,ip31,isub)
                     
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a30*field_up(ldn,ip11,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a31*field_up(ldn,ip11,ip20,ip31,isub)
                     
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a30*field_up(ldn,ip11,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a31*field_up(ldn,ip11,ip21,ip31,isub)
                  END DO   ! ldn
               END DO ! end loop over particles in the current subdomain
            END IF ! lda unroll
#endif
#endif
         END DO ! isub
      CASE DEFAULT
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',    &
     &                     'Only Mp4 and BSp2 are available. Use  ppm_rmsh_remesh for other kernels.', &
     &                     __LINE__,info)
      END SELECT ! kernel type

      !-------------------------------------------------------------------------
      !  Dont deallocate if something went wrong, as the arrays may not even be
      !  associated
      !-------------------------------------------------------------------------
      IF (info.NE.0) GOTO 9999
      !-------------------------------------------------------------------------
      !  Deallocation of the arrays....
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      ldu(1) = 0
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &               'ppm_interp_m2p', &
     &               'pb in ilist1 deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &               'ppm_interp_m2p',  &
     &               'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(store_info,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &               'ppm_interp_p2m',  &
     &               'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(list_sub,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &               'ppm_interp_p2m',  &
     &               'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      END IF
      
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_interp_m2p',t0,info)
      RETURN
     
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_interp_m2p_dv_3d
#endif
#endif
#endif       


     
