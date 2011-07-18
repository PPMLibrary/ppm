      !-------------------------------------------------------------------------
      !     Subroutine   :                 ppm_interp_p2m_renorm
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


#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_renorm_ss_2d(topoid,meshid,xp,Np,up,kernel, &
           & ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_renorm_ds_2d(topoid,meshid,xp,Np,up,kernel, &
           & ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_renorm_sv_2d(topoid,meshid,xp,Np,up,lda,kernel,&
           & ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_renorm_dv_2d(topoid,meshid,xp,Np,up,lda,kernel,&
           & ghostsize,field_up,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_renorm_ss_3d(topoid,meshid,xp,Np,up,kernel, &
           & ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_renorm_ds_3d(topoid,meshid,xp,Np,up,kernel, &
           & ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_renorm_sv_3d(topoid,meshid,xp,Np,up,lda,kernel,&
           & ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_renorm_dv_3d(topoid,meshid,xp,Np,up,lda,kernel,&
           & ghostsize,field_up,info)
#endif
#endif
#endif       
      !!! This routine does particle to mesh interpolation. Vector cases for
      !!! lda < 5 are explicitly unrolled for the 3D version. All 3D
      !!! versions are explicitly unrolled over the kernel, the 2D versions
      !!! are not.
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
      USE ppm_module_data_rmsh
      USE ppm_module_data_mesh
      USE ppm_module_write
      USE ppm_module_map
      USE ppm_module_check_id
      IMPLICIT NONE
             
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !--------------------------------------------------------------------------
      ! Arguments     
      !--------------------------------------------------------------------------
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , INTENT(IN)   :: up
      !!! particle weights from which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
#elif __MODE == __VEC
      INTEGER                         , INTENT(in)   :: lda
      !!! leading dimension of up
      REAL(MK) , DIMENSION(:,:)       , INTENT(IN)   :: up
      !!! particle weights from which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
#endif
      REAL(MK), DIMENSION(:,:)       , INTENT(IN)    :: xp
      !!! particle positions
      INTEGER , DIMENSION(:  )       , INTENT(IN)    :: ghostsize
      !!! The size (width) of the ghost layer
      INTEGER                        , INTENT(IN   ) :: Np
      !!! number of particles
      INTEGER                        , INTENT(IN   ) :: meshid
      !!! id of the mesh (user)
      INTEGER                        , INTENT(IN   ) :: topoid
      !!! topology identifier (user numbering) of target
      INTEGER                        , INTENT(IN   ) :: kernel
      !!! Choice of the kernel used to compute the weights.
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !--------------------------------------------------------------------------
      ! Local variables
      !--------------------------------------------------------------------------
      INTEGER,  DIMENSION(:,:)     , POINTER :: istart   => NULL()
      INTEGER,  DIMENSION(:,:)     , POINTER :: ndata    => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist1   => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist2   => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: min_phys => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: max_phys => NULL()
      REAL(MK),  DIMENSION(ppm_dim)           :: dxi,dx
      REAL(MK),  DIMENSION(ppm_dim)           :: len_phys
      REAL(MK)                               :: x1,x2,x3,epsilon
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
      INTEGER                                :: iq
      LOGICAL                                :: internal_weights,lok
      ! aliases
      REAL(mk), DIMENSION(:,:),      POINTER :: min_sub => NULL()
      REAL(mk), DIMENSION(:,:),      POINTER :: max_sub => NULL()
      REAL(MK)                               :: myeps
      REAL(MK)                               :: tim1s, tim1e
      REAL(MK)                               :: xp1,xp2,xp3
      REAL(MK)                               :: wx1,wx2,wx3
      REAL(MK), DIMENSION(ppm_dim)           :: x0
      CHARACTER(len=256)                     :: msg
      TYPE(ppm_t_equi_mesh), POINTER         :: p_mesh => NULL()
      TYPE(ppm_t_topo)     , POINTER         :: topo   => NULL()
      !-----------------------------------------------------
      !                          Renormalization
      !-----------------------------------------------------
#if __DIME == __2D
#error two dimensional renormalization is not yet implemented
      REAL(MK), DIMENSION(:,:,:),   POINTER  :: field_reno => NULL()
#elif __DIME == __3D
      REAL(MK), DIMENSION(:,:,:,:), POINTER  :: field_reno => NULL()
#endif
      
      !-------------------------------------------------------------------------
      !  Variables for unrolled versions
      !-------------------------------------------------------------------------
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

      
      !--------------------------------------------------------------------------
      !  Initialise 
      !--------------------------------------------------------------------------

      !--------------------------------------------------------------------------
      !  This is a hack! Somehow having trouble with constructors in the
      !  ppm_module_data_rmsh module
      !--------------------------------------------------------------------------
      ppm_rmsh_kernelsize = (/1,2,2,4/)
      
      CALL substart('ppm_interp_p2m_renorm',t0,info)

      dim = ppm_dim
      internal_weights = .FALSE.

      !--------------------------------------------------------------------------
      !  Check arguments
      !--------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN 
        CALL check
        IF (info .NE. 0) GOTO 9999
      END IF
      



      topo => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  Get the meshid
      !-------------------------------------------------------------------------
      p_mesh => topo%mesh(meshid)
      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------
      istart => p_mesh%istart

      !--------------------------------------------------------------------------
      !  Assignment of the useful arrays/scalar
      !--------------------------------------------------------------------------
      Nm(1:dim) = p_mesh%Nm
      bcdef(1:(2*dim)) = topo%bcdef(1:(2*dim))
      nsubs = topo%nsublist
     
      !--------------------------------------------------------------------------
      !  If there is nothing to do, do nearly nothing
      !--------------------------------------------------------------------------
      IF(Np.EQ.0) GOTO 9998
      
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
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_renorm_3d',     &
             &        'particle list 1 ILIST1',__LINE__,info)
        GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_renorm_3d',     &
             &        'particle list 2 ILIST2',__LINE__,info)
        GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = nsubs
      CALL ppm_alloc(store_info,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_renorm_3d',     &
              &        'store_info allocation : problem',__LINE__,info)
         GOTO 9999
      ENDIF

      !--------------------------------------------------------------------------
      !  Mesh spacing
      !--------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_phys => topo%min_physs
      max_phys => topo%max_physs
#elif __KIND == __DOUBLE_PRECISION
      min_phys => topo%min_physd
      max_phys => topo%max_physd
#endif
      
      DO i = 1,dim
         Nc(i)       = Nm(i) - 1
         len_phys(i) = max_phys(i) - min_phys(i)
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
      min_sub => topo%min_subs
      max_sub => topo%max_subs
#elif __KIND == __DOUBLE_PRECISION
      myeps = ppm_myepsd
      min_sub => topo%min_subd
      max_sub => topo%max_subd
#endif

      !--------------------------------------------------------------------------
      !  Loop over the subdomains (since the first domains are most likely
      !  to be empty, we look backwards to reduce the number of elements in
      !  nlist2 as fast as possible)
      !--------------------------------------------------------------------------
      DO idom = topo%nsublist,1,-1
         idoml = topo%isublist(idom)
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
               IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                    &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                    &xp(3,ipart).GE.min_sub(3,idoml) .AND. &
                    &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                    &xp(2,ipart).LE.max_sub(2,idoml) .AND. &
                    &xp(3,ipart).LE.max_sub(3,idoml) ) ) THEN

                  IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
                       & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                       & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                       &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                       & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
                       & bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
                       &(xp(3,ipart).LT.max_sub(3,idoml) .OR.  &
                       & (topo%subs_bc(6,idoml).EQ.1   .AND.    &
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
               
               IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                    &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                    &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                    &xp(2,ipart).LE.max_sub(2,idoml) ) ) THEN
                  
                  IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
                       & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                       & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                       &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                       & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
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
         CALL ppm_error(ppm_err_part_unass,'ppm_interp_p2m_renorm',  &
              &       'MAJOR PROBLEM',__LINE__,info)
         GOTO 9999
      ENDIF
      
      !--------------------------------------------------------------------------
      ! select the particle for with several domain computation (french)
      !--------------------------------------------------------------------------
      
      max_partnumber = 0
      DO idom=1,topo%nsublist
         IF(store_info(idom).GE.max_partnumber) THEN
            max_partnumber = store_info(idom)
         END IF
      END DO
      iopt   = ppm_param_alloc_fit
      ldu(1) = topo%nsublist
      ldu(2) = max_partnumber
      
      CALL ppm_alloc(list_sub,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_renorm_3d',     &
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
      DO idom = topo%nsublist,1,-1
         idoml = topo%isublist(idom)
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
               IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                    &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                    &xp(3,ipart).GE.min_sub(3,idoml) .AND. &
                    &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                    &xp(2,ipart).LE.max_sub(2,idoml) .AND. &
                    &xp(3,ipart).LE.max_sub(3,idoml) ) ) THEN
                  
                  IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
                       & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                       & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                       &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                       & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
                       & bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
                       &(xp(3,ipart).LT.max_sub(3,idoml) .OR.  &
                       & (topo%subs_bc(6,idoml).EQ.1   .AND.    &
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

               IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                    &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                    &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                    &xp(2,ipart).LE.max_sub(2,idoml) ) ) THEN
                                     
                  IF(   (xp(1,ipart).LT.max_sub(1,idom) .OR.  &
                       & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                       & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                       &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                       & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
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
         CALL ppm_error(ppm_err_part_unass,'ppm_interp_p2m_renorm_3d',  &
              &       'MAJOR PROBLEM',__LINE__,info)
         GOTO 9999
      ENDIF
      
      !--------------------------------------------------------------------------
      !  Allocate and alias the weights if we need them.
      !--------------------------------------------------------------------------
      max_partnumber = 0
      DO idom = 1,topo%nsublist
         IF(store_info(idom).GE.max_partnumber) THEN
            max_partnumber = store_info(idom)  
         END IF
      END DO

9998 CONTINUE
#if __DIME == __3D
      !--------------------------------------------------------------------------
      !  --- 3D ---
      !--------------------------------------------------------------------------
      !  loop over subs
      ndata => p_mesh%nnodes
      ALLOCATE(field_reno((1-ghostsize(1)):(ndata(1,1)+ghostsize(1)),&
     &                    (1-ghostsize(2)):(ndata(2,1)+ghostsize(2)),&
     &                    (1-ghostsize(3)):(ndata(3,1)+ghostsize(3)),&
     &                    topo%nsublist))
      
      DO isub = 1,topo%nsublist
         isubl = topo%isublist(isub)
         
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
                  field_reno(i,j,k,isub) = 0.0_mk
               END DO
            END DO
         END DO
      END DO
      IF(np.EQ.0) GOTO 9997
      
      SELECT CASE(kernel)
         
      CASE(ppm_param_rmsh_kernel_mp4)
         
         !-----------------------------------------------------------------
         ! M Prime Four
         !-----------------------------------------------------------------
         DO isub = 1,topo%nsublist

#if __MODE == __VEC

        !------------------------------------------------------------------------
        !  Unrolled versions for 3-vectors
        !------------------------------------------------------------------------
         IF (lda .EQ. 3) THEN
            DO ip = 1,store_info(isub)
               
               isubl = topo%isublist(isub)
               iq    = list_sub(isub,ip)


               !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
               !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
               !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

               x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
               !ip11 = INT(x0(1))+2-istart(1,isubl)
               !ip21 = INT(x0(2))+2-istart(2,isubl)
               !ip31 = INT(x0(3))+2-istart(3,isubl)

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

               field_reno(ip10,ip20,ip30,isub) = &
                    &   field_reno(ip10,ip20,ip30,isub) + &
                    &   a10a20a30
               field_reno(ip11,ip20,ip30,isub) = &
                    &   field_reno(ip11,ip20,ip30,isub) + &
                    &   a11a20a30
               field_reno(ip12,ip20,ip30,isub) = &
                    &   field_reno(ip12,ip20,ip30,isub) + &
                    &   a12a20a30
               field_reno(ip13,ip20,ip30,isub) = &
                    &   field_reno(ip13,ip20,ip30,isub) + &
                    &   a13a20a30
               field_reno(ip10,ip21,ip30,isub) = &
                    &   field_reno(ip10,ip21,ip30,isub) + &
                    &   a10a21a30
               field_reno(ip11,ip21,ip30,isub) = &
                    &   field_reno(ip11,ip21,ip30,isub) + &
                    &   a11a21a30
               field_reno(ip12,ip21,ip30,isub) = &
                    &   field_reno(ip12,ip21,ip30,isub) + &
                    &   a12a21a30
               field_reno(ip13,ip21,ip30,isub) = &
                    &   field_reno(ip13,ip21,ip30,isub) + &
                    &   a13a21a30
               field_reno(ip10,ip22,ip30,isub) = &
                    &   field_reno(ip10,ip22,ip30,isub) + &
                    &   a10a22a30
               field_reno(ip11,ip22,ip30,isub) = &
                    &   field_reno(ip11,ip22,ip30,isub) + &
                    &   a11a22a30
               field_reno(ip12,ip22,ip30,isub) = &
                    &   field_reno(ip12,ip22,ip30,isub) + &
                    &   a12a22a30
               field_reno(ip13,ip22,ip30,isub) = &
                    &   field_reno(ip13,ip22,ip30,isub) + &
                    &   a13a22a30
               field_reno(ip10,ip23,ip30,isub) = &
                    &   field_reno(ip10,ip23,ip30,isub) + &
                    &   a10a23a30
               field_reno(ip11,ip23,ip30,isub) = &
                    &   field_reno(ip11,ip23,ip30,isub) + &
                    &   a11a23a30
               field_reno(ip12,ip23,ip30,isub) = &
                    &   field_reno(ip12,ip23,ip30,isub) + &
                    &   a12a23a30
               field_reno(ip13,ip23,ip30,isub) = &
                    &   field_reno(ip13,ip23,ip30,isub) + &
                    &   a13a23a30
               field_reno(ip10,ip20,ip31,isub) = &
                    &   field_reno(ip10,ip20,ip31,isub) + &
                    &   a10a20a31
               field_reno(ip11,ip20,ip31,isub) = &
                    &   field_reno(ip11,ip20,ip31,isub) + &
                    &   a11a20a31
               field_reno(ip12,ip20,ip31,isub) = &
                    &   field_reno(ip12,ip20,ip31,isub) + &
                    &   a12a20a31
               field_reno(ip13,ip20,ip31,isub) = &
                    &   field_reno(ip13,ip20,ip31,isub) + &
                    &   a13a20a31
               field_reno(ip10,ip21,ip31,isub) = &
                    &   field_reno(ip10,ip21,ip31,isub) + &
                    &   a10a21a31
               field_reno(ip11,ip21,ip31,isub) = &
                    &   field_reno(ip11,ip21,ip31,isub) + &
                    &   a11a21a31
               field_reno(ip12,ip21,ip31,isub) = &
                    &   field_reno(ip12,ip21,ip31,isub) + &
                    &   a12a21a31
               field_reno(ip13,ip21,ip31,isub) = &
                    &   field_reno(ip13,ip21,ip31,isub) + &
                    &   a13a21a31
               field_reno(ip10,ip22,ip31,isub) = &
                    &   field_reno(ip10,ip22,ip31,isub) + &
                    &   a10a22a31
               field_reno(ip11,ip22,ip31,isub) = &
                    &   field_reno(ip11,ip22,ip31,isub) + &
                    &   a11a22a31
               field_reno(ip12,ip22,ip31,isub) = &
                    &   field_reno(ip12,ip22,ip31,isub) + &
                    &   a12a22a31
               field_reno(ip13,ip22,ip31,isub) = &
                    &   field_reno(ip13,ip22,ip31,isub) + &
                    &   a13a22a31
               field_reno(ip10,ip23,ip31,isub) = &
                    &   field_reno(ip10,ip23,ip31,isub) + &
                    &   a10a23a31
               field_reno(ip11,ip23,ip31,isub) = &
                    &   field_reno(ip11,ip23,ip31,isub) + &
                    &   a11a23a31
               field_reno(ip12,ip23,ip31,isub) = &
                    &   field_reno(ip12,ip23,ip31,isub) + &
                    &   a12a23a31
               field_reno(ip13,ip23,ip31,isub) = &
                    &   field_reno(ip13,ip23,ip31,isub) + &
                    &   a13a23a31
               field_reno(ip10,ip20,ip32,isub) = &
                    &   field_reno(ip10,ip20,ip32,isub) + &
                    &   a10a20a32
               field_reno(ip11,ip20,ip32,isub) = &
                    &   field_reno(ip11,ip20,ip32,isub) + &
                    &   a11a20a32
               field_reno(ip12,ip20,ip32,isub) = &
                    &   field_reno(ip12,ip20,ip32,isub) + &
                    &   a12a20a32
               field_reno(ip13,ip20,ip32,isub) = &
                    &   field_reno(ip13,ip20,ip32,isub) + &
                    &   a13a20a32
               field_reno(ip10,ip21,ip32,isub) = &
                    &   field_reno(ip10,ip21,ip32,isub) + &
                    &   a10a21a32
               field_reno(ip11,ip21,ip32,isub) = &
                    &   field_reno(ip11,ip21,ip32,isub) + &
                    &   a11a21a32
               field_reno(ip12,ip21,ip32,isub) = &
                    &   field_reno(ip12,ip21,ip32,isub) + &
                    &   a12a21a32
               field_reno(ip13,ip21,ip32,isub) = &
                    &   field_reno(ip13,ip21,ip32,isub) + &
                    &   a13a21a32
               field_reno(ip10,ip22,ip32,isub) = &
                    &   field_reno(ip10,ip22,ip32,isub) + &
                    &   a10a22a32
               field_reno(ip11,ip22,ip32,isub) = &
                    &   field_reno(ip11,ip22,ip32,isub) + &
                    &   a11a22a32
               field_reno(ip12,ip22,ip32,isub) = &
                    &   field_reno(ip12,ip22,ip32,isub) + &
                    &   a12a22a32
               field_reno(ip13,ip22,ip32,isub) = &
                    &   field_reno(ip13,ip22,ip32,isub) + &
                    &   a13a22a32
               field_reno(ip10,ip23,ip32,isub) = &
                    &   field_reno(ip10,ip23,ip32,isub) + &
                    &   a10a23a32
               field_reno(ip11,ip23,ip32,isub) = &
                    &   field_reno(ip11,ip23,ip32,isub) + &
                    &   a11a23a32
               field_reno(ip12,ip23,ip32,isub) = &
                    &   field_reno(ip12,ip23,ip32,isub) + &
                    &   a12a23a32
               field_reno(ip13,ip23,ip32,isub) = &
                    &   field_reno(ip13,ip23,ip32,isub) + &
                    &   a13a23a32
               field_reno(ip10,ip20,ip33,isub) = &
                    &   field_reno(ip10,ip20,ip33,isub) + &
                    &   a10a20a33
               field_reno(ip11,ip20,ip33,isub) = &
                    &   field_reno(ip11,ip20,ip33,isub) + &
                    &   a11a20a33
               field_reno(ip12,ip20,ip33,isub) = &
                    &   field_reno(ip12,ip20,ip33,isub) + &
                    &   a12a20a33
               field_reno(ip13,ip20,ip33,isub) = &
                    &   field_reno(ip13,ip20,ip33,isub) + &
                    &   a13a20a33
               field_reno(ip10,ip21,ip33,isub) = &
                    &   field_reno(ip10,ip21,ip33,isub) + &
                    &   a10a21a33
               field_reno(ip11,ip21,ip33,isub) = &
                    &   field_reno(ip11,ip21,ip33,isub) + &
                    &   a11a21a33
               field_reno(ip12,ip21,ip33,isub) = &
                    &   field_reno(ip12,ip21,ip33,isub) + &
                    &   a12a21a33
               field_reno(ip13,ip21,ip33,isub) = &
                    &   field_reno(ip13,ip21,ip33,isub) + &
                    &   a13a21a33
               field_reno(ip10,ip22,ip33,isub) = &
                    &   field_reno(ip10,ip22,ip33,isub) + &
                    &   a10a22a33
               field_reno(ip11,ip22,ip33,isub) = &
                    &   field_reno(ip11,ip22,ip33,isub) + &
                    &   a11a22a33
               field_reno(ip12,ip22,ip33,isub) = &
                    &   field_reno(ip12,ip22,ip33,isub) + &
                    &   a12a22a33
               field_reno(ip13,ip22,ip33,isub) = &
                    &   field_reno(ip13,ip22,ip33,isub) + &
                    &   a13a22a33
               field_reno(ip10,ip23,ip33,isub) = &
                    &   field_reno(ip10,ip23,ip33,isub) + &
                    &   a10a23a33
               field_reno(ip11,ip23,ip33,isub) = &
                    &   field_reno(ip11,ip23,ip33,isub) + &
                    &   a11a23a33
               field_reno(ip12,ip23,ip33,isub) = &
                    &   field_reno(ip12,ip23,ip33,isub) + &
                    &   a12a23a33
               field_reno(ip13,ip23,ip33,isub) = &
                    &   field_reno(ip13,ip23,ip33,isub) + &
                    &   a13a23a33

               
#ifdef __NOMICROINSTRUCTIONS
field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
&a10a20a30*up(1,iq)
field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
&a10a20a30*up(2,iq)
field_up(3,ip10,ip20,ip30,isub)=field_up(3,ip10,ip20,ip30,isub)+&
&a10a20a30*up(3,iq)
field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
&a11a20a30*up(1,iq)
field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
&a11a20a30*up(2,iq)
field_up(3,ip11,ip20,ip30,isub)=field_up(3,ip11,ip20,ip30,isub)+&
&a11a20a30*up(3,iq)
field_up(1,ip12,ip20,ip30,isub)=field_up(1,ip12,ip20,ip30,isub)+&
&a12a20a30*up(1,iq)
field_up(2,ip12,ip20,ip30,isub)=field_up(2,ip12,ip20,ip30,isub)+&
&a12a20a30*up(2,iq)
field_up(3,ip12,ip20,ip30,isub)=field_up(3,ip12,ip20,ip30,isub)+&
&a12a20a30*up(3,iq)
field_up(1,ip13,ip20,ip30,isub)=field_up(1,ip13,ip20,ip30,isub)+&
&a13a20a30*up(1,iq)
field_up(2,ip13,ip20,ip30,isub)=field_up(2,ip13,ip20,ip30,isub)+&
&a13a20a30*up(2,iq)
field_up(3,ip13,ip20,ip30,isub)=field_up(3,ip13,ip20,ip30,isub)+&
&a13a20a30*up(3,iq)
field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
&a10a21a30*up(1,iq)
field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
&a10a21a30*up(2,iq)
field_up(3,ip10,ip21,ip30,isub)=field_up(3,ip10,ip21,ip30,isub)+&
&a10a21a30*up(3,iq)
field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
&a11a21a30*up(1,iq)
field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
&a11a21a30*up(2,iq)
field_up(3,ip11,ip21,ip30,isub)=field_up(3,ip11,ip21,ip30,isub)+&
&a11a21a30*up(3,iq)
field_up(1,ip12,ip21,ip30,isub)=field_up(1,ip12,ip21,ip30,isub)+&
&a12a21a30*up(1,iq)
field_up(2,ip12,ip21,ip30,isub)=field_up(2,ip12,ip21,ip30,isub)+&
&a12a21a30*up(2,iq)
field_up(3,ip12,ip21,ip30,isub)=field_up(3,ip12,ip21,ip30,isub)+&
&a12a21a30*up(3,iq)
field_up(1,ip13,ip21,ip30,isub)=field_up(1,ip13,ip21,ip30,isub)+&
&a13a21a30*up(1,iq)
field_up(2,ip13,ip21,ip30,isub)=field_up(2,ip13,ip21,ip30,isub)+&
&a13a21a30*up(2,iq)
field_up(3,ip13,ip21,ip30,isub)=field_up(3,ip13,ip21,ip30,isub)+&
&a13a21a30*up(3,iq)
field_up(1,ip10,ip22,ip30,isub)=field_up(1,ip10,ip22,ip30,isub)+&
&a10a22a30*up(1,iq)
field_up(2,ip10,ip22,ip30,isub)=field_up(2,ip10,ip22,ip30,isub)+&
&a10a22a30*up(2,iq)
field_up(3,ip10,ip22,ip30,isub)=field_up(3,ip10,ip22,ip30,isub)+&
&a10a22a30*up(3,iq)
field_up(1,ip11,ip22,ip30,isub)=field_up(1,ip11,ip22,ip30,isub)+&
&a11a22a30*up(1,iq)
field_up(2,ip11,ip22,ip30,isub)=field_up(2,ip11,ip22,ip30,isub)+&
&a11a22a30*up(2,iq)
field_up(3,ip11,ip22,ip30,isub)=field_up(3,ip11,ip22,ip30,isub)+&
&a11a22a30*up(3,iq)
field_up(1,ip12,ip22,ip30,isub)=field_up(1,ip12,ip22,ip30,isub)+&
&a12a22a30*up(1,iq)
field_up(2,ip12,ip22,ip30,isub)=field_up(2,ip12,ip22,ip30,isub)+&
&a12a22a30*up(2,iq)
field_up(3,ip12,ip22,ip30,isub)=field_up(3,ip12,ip22,ip30,isub)+&
&a12a22a30*up(3,iq)
field_up(1,ip13,ip22,ip30,isub)=field_up(1,ip13,ip22,ip30,isub)+&
&a13a22a30*up(1,iq)
field_up(2,ip13,ip22,ip30,isub)=field_up(2,ip13,ip22,ip30,isub)+&
&a13a22a30*up(2,iq)
field_up(3,ip13,ip22,ip30,isub)=field_up(3,ip13,ip22,ip30,isub)+&
&a13a22a30*up(3,iq)
field_up(1,ip10,ip23,ip30,isub)=field_up(1,ip10,ip23,ip30,isub)+&
&a10a23a30*up(1,iq)
field_up(2,ip10,ip23,ip30,isub)=field_up(2,ip10,ip23,ip30,isub)+&
&a10a23a30*up(2,iq)
field_up(3,ip10,ip23,ip30,isub)=field_up(3,ip10,ip23,ip30,isub)+&
&a10a23a30*up(3,iq)
field_up(1,ip11,ip23,ip30,isub)=field_up(1,ip11,ip23,ip30,isub)+&
&a11a23a30*up(1,iq)
field_up(2,ip11,ip23,ip30,isub)=field_up(2,ip11,ip23,ip30,isub)+&
&a11a23a30*up(2,iq)
field_up(3,ip11,ip23,ip30,isub)=field_up(3,ip11,ip23,ip30,isub)+&
&a11a23a30*up(3,iq)
field_up(1,ip12,ip23,ip30,isub)=field_up(1,ip12,ip23,ip30,isub)+&
&a12a23a30*up(1,iq)
field_up(2,ip12,ip23,ip30,isub)=field_up(2,ip12,ip23,ip30,isub)+&
&a12a23a30*up(2,iq)
field_up(3,ip12,ip23,ip30,isub)=field_up(3,ip12,ip23,ip30,isub)+&
&a12a23a30*up(3,iq)
field_up(1,ip13,ip23,ip30,isub)=field_up(1,ip13,ip23,ip30,isub)+&
&a13a23a30*up(1,iq)
field_up(2,ip13,ip23,ip30,isub)=field_up(2,ip13,ip23,ip30,isub)+&
&a13a23a30*up(2,iq)
field_up(3,ip13,ip23,ip30,isub)=field_up(3,ip13,ip23,ip30,isub)+&
&a13a23a30*up(3,iq)
field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
&a10a20a31*up(1,iq)
field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
&a10a20a31*up(2,iq)
field_up(3,ip10,ip20,ip31,isub)=field_up(3,ip10,ip20,ip31,isub)+&
&a10a20a31*up(3,iq)
field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
&a11a20a31*up(1,iq)
field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
&a11a20a31*up(2,iq)
field_up(3,ip11,ip20,ip31,isub)=field_up(3,ip11,ip20,ip31,isub)+&
&a11a20a31*up(3,iq)
field_up(1,ip12,ip20,ip31,isub)=field_up(1,ip12,ip20,ip31,isub)+&
&a12a20a31*up(1,iq)
field_up(2,ip12,ip20,ip31,isub)=field_up(2,ip12,ip20,ip31,isub)+&
&a12a20a31*up(2,iq)
field_up(3,ip12,ip20,ip31,isub)=field_up(3,ip12,ip20,ip31,isub)+&
&a12a20a31*up(3,iq)
field_up(1,ip13,ip20,ip31,isub)=field_up(1,ip13,ip20,ip31,isub)+&
&a13a20a31*up(1,iq)
field_up(2,ip13,ip20,ip31,isub)=field_up(2,ip13,ip20,ip31,isub)+&
&a13a20a31*up(2,iq)
field_up(3,ip13,ip20,ip31,isub)=field_up(3,ip13,ip20,ip31,isub)+&
&a13a20a31*up(3,iq)
field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
&a10a21a31*up(1,iq)
field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
&a10a21a31*up(2,iq)
field_up(3,ip10,ip21,ip31,isub)=field_up(3,ip10,ip21,ip31,isub)+&
&a10a21a31*up(3,iq)
field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
&a11a21a31*up(1,iq)
field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
&a11a21a31*up(2,iq)
field_up(3,ip11,ip21,ip31,isub)=field_up(3,ip11,ip21,ip31,isub)+&
&a11a21a31*up(3,iq)
field_up(1,ip12,ip21,ip31,isub)=field_up(1,ip12,ip21,ip31,isub)+&
&a12a21a31*up(1,iq)
field_up(2,ip12,ip21,ip31,isub)=field_up(2,ip12,ip21,ip31,isub)+&
&a12a21a31*up(2,iq)
field_up(3,ip12,ip21,ip31,isub)=field_up(3,ip12,ip21,ip31,isub)+&
&a12a21a31*up(3,iq)
field_up(1,ip13,ip21,ip31,isub)=field_up(1,ip13,ip21,ip31,isub)+&
&a13a21a31*up(1,iq)
field_up(2,ip13,ip21,ip31,isub)=field_up(2,ip13,ip21,ip31,isub)+&
&a13a21a31*up(2,iq)
field_up(3,ip13,ip21,ip31,isub)=field_up(3,ip13,ip21,ip31,isub)+&
&a13a21a31*up(3,iq)
field_up(1,ip10,ip22,ip31,isub)=field_up(1,ip10,ip22,ip31,isub)+&
&a10a22a31*up(1,iq)
field_up(2,ip10,ip22,ip31,isub)=field_up(2,ip10,ip22,ip31,isub)+&
&a10a22a31*up(2,iq)
field_up(3,ip10,ip22,ip31,isub)=field_up(3,ip10,ip22,ip31,isub)+&
&a10a22a31*up(3,iq)
field_up(1,ip11,ip22,ip31,isub)=field_up(1,ip11,ip22,ip31,isub)+&
&a11a22a31*up(1,iq)
field_up(2,ip11,ip22,ip31,isub)=field_up(2,ip11,ip22,ip31,isub)+&
&a11a22a31*up(2,iq)
field_up(3,ip11,ip22,ip31,isub)=field_up(3,ip11,ip22,ip31,isub)+&
&a11a22a31*up(3,iq)
field_up(1,ip12,ip22,ip31,isub)=field_up(1,ip12,ip22,ip31,isub)+&
&a12a22a31*up(1,iq)
field_up(2,ip12,ip22,ip31,isub)=field_up(2,ip12,ip22,ip31,isub)+&
&a12a22a31*up(2,iq)
field_up(3,ip12,ip22,ip31,isub)=field_up(3,ip12,ip22,ip31,isub)+&
&a12a22a31*up(3,iq)
field_up(1,ip13,ip22,ip31,isub)=field_up(1,ip13,ip22,ip31,isub)+&
&a13a22a31*up(1,iq)
field_up(2,ip13,ip22,ip31,isub)=field_up(2,ip13,ip22,ip31,isub)+&
&a13a22a31*up(2,iq)
field_up(3,ip13,ip22,ip31,isub)=field_up(3,ip13,ip22,ip31,isub)+&
&a13a22a31*up(3,iq)
field_up(1,ip10,ip23,ip31,isub)=field_up(1,ip10,ip23,ip31,isub)+&
&a10a23a31*up(1,iq)
field_up(2,ip10,ip23,ip31,isub)=field_up(2,ip10,ip23,ip31,isub)+&
&a10a23a31*up(2,iq)
field_up(3,ip10,ip23,ip31,isub)=field_up(3,ip10,ip23,ip31,isub)+&
&a10a23a31*up(3,iq)
field_up(1,ip11,ip23,ip31,isub)=field_up(1,ip11,ip23,ip31,isub)+&
&a11a23a31*up(1,iq)
field_up(2,ip11,ip23,ip31,isub)=field_up(2,ip11,ip23,ip31,isub)+&
&a11a23a31*up(2,iq)
field_up(3,ip11,ip23,ip31,isub)=field_up(3,ip11,ip23,ip31,isub)+&
&a11a23a31*up(3,iq)
field_up(1,ip12,ip23,ip31,isub)=field_up(1,ip12,ip23,ip31,isub)+&
&a12a23a31*up(1,iq)
field_up(2,ip12,ip23,ip31,isub)=field_up(2,ip12,ip23,ip31,isub)+&
&a12a23a31*up(2,iq)
field_up(3,ip12,ip23,ip31,isub)=field_up(3,ip12,ip23,ip31,isub)+&
&a12a23a31*up(3,iq)
field_up(1,ip13,ip23,ip31,isub)=field_up(1,ip13,ip23,ip31,isub)+&
&a13a23a31*up(1,iq)
field_up(2,ip13,ip23,ip31,isub)=field_up(2,ip13,ip23,ip31,isub)+&
&a13a23a31*up(2,iq)
field_up(3,ip13,ip23,ip31,isub)=field_up(3,ip13,ip23,ip31,isub)+&
&a13a23a31*up(3,iq)
field_up(1,ip10,ip20,ip32,isub)=field_up(1,ip10,ip20,ip32,isub)+&
&a10a20a32*up(1,iq)
field_up(2,ip10,ip20,ip32,isub)=field_up(2,ip10,ip20,ip32,isub)+&
&a10a20a32*up(2,iq)
field_up(3,ip10,ip20,ip32,isub)=field_up(3,ip10,ip20,ip32,isub)+&
&a10a20a32*up(3,iq)
field_up(1,ip11,ip20,ip32,isub)=field_up(1,ip11,ip20,ip32,isub)+&
&a11a20a32*up(1,iq)
field_up(2,ip11,ip20,ip32,isub)=field_up(2,ip11,ip20,ip32,isub)+&
&a11a20a32*up(2,iq)
field_up(3,ip11,ip20,ip32,isub)=field_up(3,ip11,ip20,ip32,isub)+&
&a11a20a32*up(3,iq)
field_up(1,ip12,ip20,ip32,isub)=field_up(1,ip12,ip20,ip32,isub)+&
&a12a20a32*up(1,iq)
field_up(2,ip12,ip20,ip32,isub)=field_up(2,ip12,ip20,ip32,isub)+&
&a12a20a32*up(2,iq)
field_up(3,ip12,ip20,ip32,isub)=field_up(3,ip12,ip20,ip32,isub)+&
&a12a20a32*up(3,iq)
field_up(1,ip13,ip20,ip32,isub)=field_up(1,ip13,ip20,ip32,isub)+&
&a13a20a32*up(1,iq)
field_up(2,ip13,ip20,ip32,isub)=field_up(2,ip13,ip20,ip32,isub)+&
&a13a20a32*up(2,iq)
field_up(3,ip13,ip20,ip32,isub)=field_up(3,ip13,ip20,ip32,isub)+&
&a13a20a32*up(3,iq)
field_up(1,ip10,ip21,ip32,isub)=field_up(1,ip10,ip21,ip32,isub)+&
&a10a21a32*up(1,iq)
field_up(2,ip10,ip21,ip32,isub)=field_up(2,ip10,ip21,ip32,isub)+&
&a10a21a32*up(2,iq)
field_up(3,ip10,ip21,ip32,isub)=field_up(3,ip10,ip21,ip32,isub)+&
&a10a21a32*up(3,iq)
field_up(1,ip11,ip21,ip32,isub)=field_up(1,ip11,ip21,ip32,isub)+&
&a11a21a32*up(1,iq)
field_up(2,ip11,ip21,ip32,isub)=field_up(2,ip11,ip21,ip32,isub)+&
&a11a21a32*up(2,iq)
field_up(3,ip11,ip21,ip32,isub)=field_up(3,ip11,ip21,ip32,isub)+&
&a11a21a32*up(3,iq)
field_up(1,ip12,ip21,ip32,isub)=field_up(1,ip12,ip21,ip32,isub)+&
&a12a21a32*up(1,iq)
field_up(2,ip12,ip21,ip32,isub)=field_up(2,ip12,ip21,ip32,isub)+&
&a12a21a32*up(2,iq)
field_up(3,ip12,ip21,ip32,isub)=field_up(3,ip12,ip21,ip32,isub)+&
&a12a21a32*up(3,iq)
field_up(1,ip13,ip21,ip32,isub)=field_up(1,ip13,ip21,ip32,isub)+&
&a13a21a32*up(1,iq)
field_up(2,ip13,ip21,ip32,isub)=field_up(2,ip13,ip21,ip32,isub)+&
&a13a21a32*up(2,iq)
field_up(3,ip13,ip21,ip32,isub)=field_up(3,ip13,ip21,ip32,isub)+&
&a13a21a32*up(3,iq)
field_up(1,ip10,ip22,ip32,isub)=field_up(1,ip10,ip22,ip32,isub)+&
&a10a22a32*up(1,iq)
field_up(2,ip10,ip22,ip32,isub)=field_up(2,ip10,ip22,ip32,isub)+&
&a10a22a32*up(2,iq)
field_up(3,ip10,ip22,ip32,isub)=field_up(3,ip10,ip22,ip32,isub)+&
&a10a22a32*up(3,iq)
field_up(1,ip11,ip22,ip32,isub)=field_up(1,ip11,ip22,ip32,isub)+&
&a11a22a32*up(1,iq)
field_up(2,ip11,ip22,ip32,isub)=field_up(2,ip11,ip22,ip32,isub)+&
&a11a22a32*up(2,iq)
field_up(3,ip11,ip22,ip32,isub)=field_up(3,ip11,ip22,ip32,isub)+&
&a11a22a32*up(3,iq)
field_up(1,ip12,ip22,ip32,isub)=field_up(1,ip12,ip22,ip32,isub)+&
&a12a22a32*up(1,iq)
field_up(2,ip12,ip22,ip32,isub)=field_up(2,ip12,ip22,ip32,isub)+&
&a12a22a32*up(2,iq)
field_up(3,ip12,ip22,ip32,isub)=field_up(3,ip12,ip22,ip32,isub)+&
&a12a22a32*up(3,iq)
field_up(1,ip13,ip22,ip32,isub)=field_up(1,ip13,ip22,ip32,isub)+&
&a13a22a32*up(1,iq)
field_up(2,ip13,ip22,ip32,isub)=field_up(2,ip13,ip22,ip32,isub)+&
&a13a22a32*up(2,iq)
field_up(3,ip13,ip22,ip32,isub)=field_up(3,ip13,ip22,ip32,isub)+&
&a13a22a32*up(3,iq)
field_up(1,ip10,ip23,ip32,isub)=field_up(1,ip10,ip23,ip32,isub)+&
&a10a23a32*up(1,iq)
field_up(2,ip10,ip23,ip32,isub)=field_up(2,ip10,ip23,ip32,isub)+&
&a10a23a32*up(2,iq)
field_up(3,ip10,ip23,ip32,isub)=field_up(3,ip10,ip23,ip32,isub)+&
&a10a23a32*up(3,iq)
field_up(1,ip11,ip23,ip32,isub)=field_up(1,ip11,ip23,ip32,isub)+&
&a11a23a32*up(1,iq)
field_up(2,ip11,ip23,ip32,isub)=field_up(2,ip11,ip23,ip32,isub)+&
&a11a23a32*up(2,iq)
field_up(3,ip11,ip23,ip32,isub)=field_up(3,ip11,ip23,ip32,isub)+&
&a11a23a32*up(3,iq)
field_up(1,ip12,ip23,ip32,isub)=field_up(1,ip12,ip23,ip32,isub)+&
&a12a23a32*up(1,iq)
field_up(2,ip12,ip23,ip32,isub)=field_up(2,ip12,ip23,ip32,isub)+&
&a12a23a32*up(2,iq)
field_up(3,ip12,ip23,ip32,isub)=field_up(3,ip12,ip23,ip32,isub)+&
&a12a23a32*up(3,iq)
field_up(1,ip13,ip23,ip32,isub)=field_up(1,ip13,ip23,ip32,isub)+&
&a13a23a32*up(1,iq)
field_up(2,ip13,ip23,ip32,isub)=field_up(2,ip13,ip23,ip32,isub)+&
&a13a23a32*up(2,iq)
field_up(3,ip13,ip23,ip32,isub)=field_up(3,ip13,ip23,ip32,isub)+&
&a13a23a32*up(3,iq)
field_up(1,ip10,ip20,ip33,isub)=field_up(1,ip10,ip20,ip33,isub)+&
&a10a20a33*up(1,iq)
field_up(2,ip10,ip20,ip33,isub)=field_up(2,ip10,ip20,ip33,isub)+&
&a10a20a33*up(2,iq)
field_up(3,ip10,ip20,ip33,isub)=field_up(3,ip10,ip20,ip33,isub)+&
&a10a20a33*up(3,iq)
field_up(1,ip11,ip20,ip33,isub)=field_up(1,ip11,ip20,ip33,isub)+&
&a11a20a33*up(1,iq)
field_up(2,ip11,ip20,ip33,isub)=field_up(2,ip11,ip20,ip33,isub)+&
&a11a20a33*up(2,iq)
field_up(3,ip11,ip20,ip33,isub)=field_up(3,ip11,ip20,ip33,isub)+&
&a11a20a33*up(3,iq)
field_up(1,ip12,ip20,ip33,isub)=field_up(1,ip12,ip20,ip33,isub)+&
&a12a20a33*up(1,iq)
field_up(2,ip12,ip20,ip33,isub)=field_up(2,ip12,ip20,ip33,isub)+&
&a12a20a33*up(2,iq)
field_up(3,ip12,ip20,ip33,isub)=field_up(3,ip12,ip20,ip33,isub)+&
&a12a20a33*up(3,iq)
field_up(1,ip13,ip20,ip33,isub)=field_up(1,ip13,ip20,ip33,isub)+&
&a13a20a33*up(1,iq)
field_up(2,ip13,ip20,ip33,isub)=field_up(2,ip13,ip20,ip33,isub)+&
&a13a20a33*up(2,iq)
field_up(3,ip13,ip20,ip33,isub)=field_up(3,ip13,ip20,ip33,isub)+&
&a13a20a33*up(3,iq)
field_up(1,ip10,ip21,ip33,isub)=field_up(1,ip10,ip21,ip33,isub)+&
&a10a21a33*up(1,iq)
field_up(2,ip10,ip21,ip33,isub)=field_up(2,ip10,ip21,ip33,isub)+&
&a10a21a33*up(2,iq)
field_up(3,ip10,ip21,ip33,isub)=field_up(3,ip10,ip21,ip33,isub)+&
&a10a21a33*up(3,iq)
field_up(1,ip11,ip21,ip33,isub)=field_up(1,ip11,ip21,ip33,isub)+&
&a11a21a33*up(1,iq)
field_up(2,ip11,ip21,ip33,isub)=field_up(2,ip11,ip21,ip33,isub)+&
&a11a21a33*up(2,iq)
field_up(3,ip11,ip21,ip33,isub)=field_up(3,ip11,ip21,ip33,isub)+&
&a11a21a33*up(3,iq)
field_up(1,ip12,ip21,ip33,isub)=field_up(1,ip12,ip21,ip33,isub)+&
&a12a21a33*up(1,iq)
field_up(2,ip12,ip21,ip33,isub)=field_up(2,ip12,ip21,ip33,isub)+&
&a12a21a33*up(2,iq)
field_up(3,ip12,ip21,ip33,isub)=field_up(3,ip12,ip21,ip33,isub)+&
&a12a21a33*up(3,iq)
field_up(1,ip13,ip21,ip33,isub)=field_up(1,ip13,ip21,ip33,isub)+&
&a13a21a33*up(1,iq)
field_up(2,ip13,ip21,ip33,isub)=field_up(2,ip13,ip21,ip33,isub)+&
&a13a21a33*up(2,iq)
field_up(3,ip13,ip21,ip33,isub)=field_up(3,ip13,ip21,ip33,isub)+&
&a13a21a33*up(3,iq)
field_up(1,ip10,ip22,ip33,isub)=field_up(1,ip10,ip22,ip33,isub)+&
&a10a22a33*up(1,iq)
field_up(2,ip10,ip22,ip33,isub)=field_up(2,ip10,ip22,ip33,isub)+&
&a10a22a33*up(2,iq)
field_up(3,ip10,ip22,ip33,isub)=field_up(3,ip10,ip22,ip33,isub)+&
&a10a22a33*up(3,iq)
field_up(1,ip11,ip22,ip33,isub)=field_up(1,ip11,ip22,ip33,isub)+&
&a11a22a33*up(1,iq)
field_up(2,ip11,ip22,ip33,isub)=field_up(2,ip11,ip22,ip33,isub)+&
&a11a22a33*up(2,iq)
field_up(3,ip11,ip22,ip33,isub)=field_up(3,ip11,ip22,ip33,isub)+&
&a11a22a33*up(3,iq)
field_up(1,ip12,ip22,ip33,isub)=field_up(1,ip12,ip22,ip33,isub)+&
&a12a22a33*up(1,iq)
field_up(2,ip12,ip22,ip33,isub)=field_up(2,ip12,ip22,ip33,isub)+&
&a12a22a33*up(2,iq)
field_up(3,ip12,ip22,ip33,isub)=field_up(3,ip12,ip22,ip33,isub)+&
&a12a22a33*up(3,iq)
field_up(1,ip13,ip22,ip33,isub)=field_up(1,ip13,ip22,ip33,isub)+&
&a13a22a33*up(1,iq)
field_up(2,ip13,ip22,ip33,isub)=field_up(2,ip13,ip22,ip33,isub)+&
&a13a22a33*up(2,iq)
field_up(3,ip13,ip22,ip33,isub)=field_up(3,ip13,ip22,ip33,isub)+&
&a13a22a33*up(3,iq)
field_up(1,ip10,ip23,ip33,isub)=field_up(1,ip10,ip23,ip33,isub)+&
&a10a23a33*up(1,iq)
field_up(2,ip10,ip23,ip33,isub)=field_up(2,ip10,ip23,ip33,isub)+&
&a10a23a33*up(2,iq)
field_up(3,ip10,ip23,ip33,isub)=field_up(3,ip10,ip23,ip33,isub)+&
&a10a23a33*up(3,iq)
field_up(1,ip11,ip23,ip33,isub)=field_up(1,ip11,ip23,ip33,isub)+&
&a11a23a33*up(1,iq)
field_up(2,ip11,ip23,ip33,isub)=field_up(2,ip11,ip23,ip33,isub)+&
&a11a23a33*up(2,iq)
field_up(3,ip11,ip23,ip33,isub)=field_up(3,ip11,ip23,ip33,isub)+&
&a11a23a33*up(3,iq)
field_up(1,ip12,ip23,ip33,isub)=field_up(1,ip12,ip23,ip33,isub)+&
&a12a23a33*up(1,iq)
field_up(2,ip12,ip23,ip33,isub)=field_up(2,ip12,ip23,ip33,isub)+&
&a12a23a33*up(2,iq)
field_up(3,ip12,ip23,ip33,isub)=field_up(3,ip12,ip23,ip33,isub)+&
&a12a23a33*up(3,iq)
field_up(1,ip13,ip23,ip33,isub)=field_up(1,ip13,ip23,ip33,isub)+&
&a13a23a33*up(1,iq)
field_up(2,ip13,ip23,ip33,isub)=field_up(2,ip13,ip23,ip33,isub)+&
&a13a23a33*up(2,iq)
field_up(3,ip13,ip23,ip33,isub)=field_up(3,ip13,ip23,ip33,isub)+&
&a13a23a33*up(3,iq)
#else
field_up(1:3,ip10,ip20,ip30,isub)=field_up(1:3,ip10,ip20,ip30,isub)+&
&a10a20a30*up(1:3,iq)
field_up(1:3,ip11,ip20,ip30,isub)=field_up(1:3,ip11,ip20,ip30,isub)+&
&a11a20a30*up(1:3,iq)
field_up(1:3,ip12,ip20,ip30,isub)=field_up(1:3,ip12,ip20,ip30,isub)+&
&a12a20a30*up(1:3,iq)
field_up(1:3,ip13,ip20,ip30,isub)=field_up(1:3,ip13,ip20,ip30,isub)+&
&a13a20a30*up(1:3,iq)
field_up(1:3,ip10,ip21,ip30,isub)=field_up(1:3,ip10,ip21,ip30,isub)+&
&a10a21a30*up(1:3,iq)
field_up(1:3,ip11,ip21,ip30,isub)=field_up(1:3,ip11,ip21,ip30,isub)+&
&a11a21a30*up(1:3,iq)
field_up(1:3,ip12,ip21,ip30,isub)=field_up(1:3,ip12,ip21,ip30,isub)+&
&a12a21a30*up(1:3,iq)
field_up(1:3,ip13,ip21,ip30,isub)=field_up(1:3,ip13,ip21,ip30,isub)+&
&a13a21a30*up(1:3,iq)
field_up(1:3,ip10,ip22,ip30,isub)=field_up(1:3,ip10,ip22,ip30,isub)+&
&a10a22a30*up(1:3,iq)
field_up(1:3,ip11,ip22,ip30,isub)=field_up(1:3,ip11,ip22,ip30,isub)+&
&a11a22a30*up(1:3,iq)
field_up(1:3,ip12,ip22,ip30,isub)=field_up(1:3,ip12,ip22,ip30,isub)+&
&a12a22a30*up(1:3,iq)
field_up(1:3,ip13,ip22,ip30,isub)=field_up(1:3,ip13,ip22,ip30,isub)+&
&a13a22a30*up(1:3,iq)
field_up(1:3,ip10,ip23,ip30,isub)=field_up(1:3,ip10,ip23,ip30,isub)+&
&a10a23a30*up(1:3,iq)
field_up(1:3,ip11,ip23,ip30,isub)=field_up(1:3,ip11,ip23,ip30,isub)+&
&a11a23a30*up(1:3,iq)
field_up(1:3,ip12,ip23,ip30,isub)=field_up(1:3,ip12,ip23,ip30,isub)+&
&a12a23a30*up(1:3,iq)
field_up(1:3,ip13,ip23,ip30,isub)=field_up(1:3,ip13,ip23,ip30,isub)+&
&a13a23a30*up(1:3,iq)
field_up(1:3,ip10,ip20,ip31,isub)=field_up(1:3,ip10,ip20,ip31,isub)+&
&a10a20a31*up(1:3,iq)
field_up(1:3,ip11,ip20,ip31,isub)=field_up(1:3,ip11,ip20,ip31,isub)+&
&a11a20a31*up(1:3,iq)
field_up(1:3,ip12,ip20,ip31,isub)=field_up(1:3,ip12,ip20,ip31,isub)+&
&a12a20a31*up(1:3,iq)
field_up(1:3,ip13,ip20,ip31,isub)=field_up(1:3,ip13,ip20,ip31,isub)+&
&a13a20a31*up(1:3,iq)
field_up(1:3,ip10,ip21,ip31,isub)=field_up(1:3,ip10,ip21,ip31,isub)+&
&a10a21a31*up(1:3,iq)
field_up(1:3,ip11,ip21,ip31,isub)=field_up(1:3,ip11,ip21,ip31,isub)+&
&a11a21a31*up(1:3,iq)
field_up(1:3,ip12,ip21,ip31,isub)=field_up(1:3,ip12,ip21,ip31,isub)+&
&a12a21a31*up(1:3,iq)
field_up(1:3,ip13,ip21,ip31,isub)=field_up(1:3,ip13,ip21,ip31,isub)+&
&a13a21a31*up(1:3,iq)
field_up(1:3,ip10,ip22,ip31,isub)=field_up(1:3,ip10,ip22,ip31,isub)+&
&a10a22a31*up(1:3,iq)
field_up(1:3,ip11,ip22,ip31,isub)=field_up(1:3,ip11,ip22,ip31,isub)+&
&a11a22a31*up(1:3,iq)
field_up(1:3,ip12,ip22,ip31,isub)=field_up(1:3,ip12,ip22,ip31,isub)+&
&a12a22a31*up(1:3,iq)
field_up(1:3,ip13,ip22,ip31,isub)=field_up(1:3,ip13,ip22,ip31,isub)+&
&a13a22a31*up(1:3,iq)
field_up(1:3,ip10,ip23,ip31,isub)=field_up(1:3,ip10,ip23,ip31,isub)+&
&a10a23a31*up(1:3,iq)
field_up(1:3,ip11,ip23,ip31,isub)=field_up(1:3,ip11,ip23,ip31,isub)+&
&a11a23a31*up(1:3,iq)
field_up(1:3,ip12,ip23,ip31,isub)=field_up(1:3,ip12,ip23,ip31,isub)+&
&a12a23a31*up(1:3,iq)
field_up(1:3,ip13,ip23,ip31,isub)=field_up(1:3,ip13,ip23,ip31,isub)+&
&a13a23a31*up(1:3,iq)
field_up(1:3,ip10,ip20,ip32,isub)=field_up(1:3,ip10,ip20,ip32,isub)+&
&a10a20a32*up(1:3,iq)
field_up(1:3,ip11,ip20,ip32,isub)=field_up(1:3,ip11,ip20,ip32,isub)+&
&a11a20a32*up(1:3,iq)
field_up(1:3,ip12,ip20,ip32,isub)=field_up(1:3,ip12,ip20,ip32,isub)+&
&a12a20a32*up(1:3,iq)
field_up(1:3,ip13,ip20,ip32,isub)=field_up(1:3,ip13,ip20,ip32,isub)+&
&a13a20a32*up(1:3,iq)
field_up(1:3,ip10,ip21,ip32,isub)=field_up(1:3,ip10,ip21,ip32,isub)+&
&a10a21a32*up(1:3,iq)
field_up(1:3,ip11,ip21,ip32,isub)=field_up(1:3,ip11,ip21,ip32,isub)+&
&a11a21a32*up(1:3,iq)
field_up(1:3,ip12,ip21,ip32,isub)=field_up(1:3,ip12,ip21,ip32,isub)+&
&a12a21a32*up(1:3,iq)
field_up(1:3,ip13,ip21,ip32,isub)=field_up(1:3,ip13,ip21,ip32,isub)+&
&a13a21a32*up(1:3,iq)
field_up(1:3,ip10,ip22,ip32,isub)=field_up(1:3,ip10,ip22,ip32,isub)+&
&a10a22a32*up(1:3,iq)
field_up(1:3,ip11,ip22,ip32,isub)=field_up(1:3,ip11,ip22,ip32,isub)+&
&a11a22a32*up(1:3,iq)
field_up(1:3,ip12,ip22,ip32,isub)=field_up(1:3,ip12,ip22,ip32,isub)+&
&a12a22a32*up(1:3,iq)
field_up(1:3,ip13,ip22,ip32,isub)=field_up(1:3,ip13,ip22,ip32,isub)+&
&a13a22a32*up(1:3,iq)
field_up(1:3,ip10,ip23,ip32,isub)=field_up(1:3,ip10,ip23,ip32,isub)+&
&a10a23a32*up(1:3,iq)
field_up(1:3,ip11,ip23,ip32,isub)=field_up(1:3,ip11,ip23,ip32,isub)+&
&a11a23a32*up(1:3,iq)
field_up(1:3,ip12,ip23,ip32,isub)=field_up(1:3,ip12,ip23,ip32,isub)+&
&a12a23a32*up(1:3,iq)
field_up(1:3,ip13,ip23,ip32,isub)=field_up(1:3,ip13,ip23,ip32,isub)+&
&a13a23a32*up(1:3,iq)
field_up(1:3,ip10,ip20,ip33,isub)=field_up(1:3,ip10,ip20,ip33,isub)+&
&a10a20a33*up(1:3,iq)
field_up(1:3,ip11,ip20,ip33,isub)=field_up(1:3,ip11,ip20,ip33,isub)+&
&a11a20a33*up(1:3,iq)
field_up(1:3,ip12,ip20,ip33,isub)=field_up(1:3,ip12,ip20,ip33,isub)+&
&a12a20a33*up(1:3,iq)
field_up(1:3,ip13,ip20,ip33,isub)=field_up(1:3,ip13,ip20,ip33,isub)+&
&a13a20a33*up(1:3,iq)
field_up(1:3,ip10,ip21,ip33,isub)=field_up(1:3,ip10,ip21,ip33,isub)+&
&a10a21a33*up(1:3,iq)
field_up(1:3,ip11,ip21,ip33,isub)=field_up(1:3,ip11,ip21,ip33,isub)+&
&a11a21a33*up(1:3,iq)
field_up(1:3,ip12,ip21,ip33,isub)=field_up(1:3,ip12,ip21,ip33,isub)+&
&a12a21a33*up(1:3,iq)
field_up(1:3,ip13,ip21,ip33,isub)=field_up(1:3,ip13,ip21,ip33,isub)+&
&a13a21a33*up(1:3,iq)
field_up(1:3,ip10,ip22,ip33,isub)=field_up(1:3,ip10,ip22,ip33,isub)+&
&a10a22a33*up(1:3,iq)
field_up(1:3,ip11,ip22,ip33,isub)=field_up(1:3,ip11,ip22,ip33,isub)+&
&a11a22a33*up(1:3,iq)
field_up(1:3,ip12,ip22,ip33,isub)=field_up(1:3,ip12,ip22,ip33,isub)+&
&a12a22a33*up(1:3,iq)
field_up(1:3,ip13,ip22,ip33,isub)=field_up(1:3,ip13,ip22,ip33,isub)+&
&a13a22a33*up(1:3,iq)
field_up(1:3,ip10,ip23,ip33,isub)=field_up(1:3,ip10,ip23,ip33,isub)+&
&a10a23a33*up(1:3,iq)
field_up(1:3,ip11,ip23,ip33,isub)=field_up(1:3,ip11,ip23,ip33,isub)+&
&a11a23a33*up(1:3,iq)
field_up(1:3,ip12,ip23,ip33,isub)=field_up(1:3,ip12,ip23,ip33,isub)+&
&a12a23a33*up(1:3,iq)
field_up(1:3,ip13,ip23,ip33,isub)=field_up(1:3,ip13,ip23,ip33,isub)+&
&a13a23a33*up(1:3,iq)

#endif
            END DO
        !------------------------------------------------------------------------
        !  Unrolled versions for 2-vectors
        !------------------------------------------------------------------------
         ELSEIF (lda .EQ. 2) THEN
            DO ip = 1,store_info(isub)

               isubl = topo%isublist(isub)
               iq    = list_sub(isub,ip)


               !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
               !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
               !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

               x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

               !ip11 = INT(x0(1))+2-istart(1,isubl)
               !ip21 = INT(x0(2))+2-istart(2,isubl)
               !ip31 = INT(x0(3))+2-istart(3,isubl)

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

               PRINT*,'###################################################'
               PRINT*,'RENORMALIZATION NOT IMPLEMENTED FOR LDA = 2'
               PRINT*,'###################################################'
               STOP

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

#ifdef __NOMICROINSTRUCTIONS
field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
&a10a20a30*up(1,iq)
field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
&a10a20a30*up(2,iq)
field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
&a11a20a30*up(1,iq)
field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
&a11a20a30*up(2,iq)
field_up(1,ip12,ip20,ip30,isub)=field_up(1,ip12,ip20,ip30,isub)+&
&a12a20a30*up(1,iq)
field_up(2,ip12,ip20,ip30,isub)=field_up(2,ip12,ip20,ip30,isub)+&
&a12a20a30*up(2,iq)
field_up(1,ip13,ip20,ip30,isub)=field_up(1,ip13,ip20,ip30,isub)+&
&a13a20a30*up(1,iq)
field_up(2,ip13,ip20,ip30,isub)=field_up(2,ip13,ip20,ip30,isub)+&
&a13a20a30*up(2,iq)
field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
&a10a21a30*up(1,iq)
field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
&a10a21a30*up(2,iq)
field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
&a11a21a30*up(1,iq)
field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
&a11a21a30*up(2,iq)
field_up(1,ip12,ip21,ip30,isub)=field_up(1,ip12,ip21,ip30,isub)+&
&a12a21a30*up(1,iq)
field_up(2,ip12,ip21,ip30,isub)=field_up(2,ip12,ip21,ip30,isub)+&
&a12a21a30*up(2,iq)
field_up(1,ip13,ip21,ip30,isub)=field_up(1,ip13,ip21,ip30,isub)+&
&a13a21a30*up(1,iq)
field_up(2,ip13,ip21,ip30,isub)=field_up(2,ip13,ip21,ip30,isub)+&
&a13a21a30*up(2,iq)
field_up(1,ip10,ip22,ip30,isub)=field_up(1,ip10,ip22,ip30,isub)+&
&a10a22a30*up(1,iq)
field_up(2,ip10,ip22,ip30,isub)=field_up(2,ip10,ip22,ip30,isub)+&
&a10a22a30*up(2,iq)
field_up(1,ip11,ip22,ip30,isub)=field_up(1,ip11,ip22,ip30,isub)+&
&a11a22a30*up(1,iq)
field_up(2,ip11,ip22,ip30,isub)=field_up(2,ip11,ip22,ip30,isub)+&
&a11a22a30*up(2,iq)
field_up(1,ip12,ip22,ip30,isub)=field_up(1,ip12,ip22,ip30,isub)+&
&a12a22a30*up(1,iq)
field_up(2,ip12,ip22,ip30,isub)=field_up(2,ip12,ip22,ip30,isub)+&
&a12a22a30*up(2,iq)
field_up(1,ip13,ip22,ip30,isub)=field_up(1,ip13,ip22,ip30,isub)+&
&a13a22a30*up(1,iq)
field_up(2,ip13,ip22,ip30,isub)=field_up(2,ip13,ip22,ip30,isub)+&
&a13a22a30*up(2,iq)
field_up(1,ip10,ip23,ip30,isub)=field_up(1,ip10,ip23,ip30,isub)+&
&a10a23a30*up(1,iq)
field_up(2,ip10,ip23,ip30,isub)=field_up(2,ip10,ip23,ip30,isub)+&
&a10a23a30*up(2,iq)
field_up(1,ip11,ip23,ip30,isub)=field_up(1,ip11,ip23,ip30,isub)+&
&a11a23a30*up(1,iq)
field_up(2,ip11,ip23,ip30,isub)=field_up(2,ip11,ip23,ip30,isub)+&
&a11a23a30*up(2,iq)
field_up(1,ip12,ip23,ip30,isub)=field_up(1,ip12,ip23,ip30,isub)+&
&a12a23a30*up(1,iq)
field_up(2,ip12,ip23,ip30,isub)=field_up(2,ip12,ip23,ip30,isub)+&
&a12a23a30*up(2,iq)
field_up(1,ip13,ip23,ip30,isub)=field_up(1,ip13,ip23,ip30,isub)+&
&a13a23a30*up(1,iq)
field_up(2,ip13,ip23,ip30,isub)=field_up(2,ip13,ip23,ip30,isub)+&
&a13a23a30*up(2,iq)
field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
&a10a20a31*up(1,iq)
field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
&a10a20a31*up(2,iq)
field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
&a11a20a31*up(1,iq)
field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
&a11a20a31*up(2,iq)
field_up(1,ip12,ip20,ip31,isub)=field_up(1,ip12,ip20,ip31,isub)+&
&a12a20a31*up(1,iq)
field_up(2,ip12,ip20,ip31,isub)=field_up(2,ip12,ip20,ip31,isub)+&
&a12a20a31*up(2,iq)
field_up(1,ip13,ip20,ip31,isub)=field_up(1,ip13,ip20,ip31,isub)+&
&a13a20a31*up(1,iq)
field_up(2,ip13,ip20,ip31,isub)=field_up(2,ip13,ip20,ip31,isub)+&
&a13a20a31*up(2,iq)
field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
&a10a21a31*up(1,iq)
field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
&a10a21a31*up(2,iq)
field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
&a11a21a31*up(1,iq)
field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
&a11a21a31*up(2,iq)
field_up(1,ip12,ip21,ip31,isub)=field_up(1,ip12,ip21,ip31,isub)+&
&a12a21a31*up(1,iq)
field_up(2,ip12,ip21,ip31,isub)=field_up(2,ip12,ip21,ip31,isub)+&
&a12a21a31*up(2,iq)
field_up(1,ip13,ip21,ip31,isub)=field_up(1,ip13,ip21,ip31,isub)+&
&a13a21a31*up(1,iq)
field_up(2,ip13,ip21,ip31,isub)=field_up(2,ip13,ip21,ip31,isub)+&
&a13a21a31*up(2,iq)
field_up(1,ip10,ip22,ip31,isub)=field_up(1,ip10,ip22,ip31,isub)+&
&a10a22a31*up(1,iq)
field_up(2,ip10,ip22,ip31,isub)=field_up(2,ip10,ip22,ip31,isub)+&
&a10a22a31*up(2,iq)
field_up(1,ip11,ip22,ip31,isub)=field_up(1,ip11,ip22,ip31,isub)+&
&a11a22a31*up(1,iq)
field_up(2,ip11,ip22,ip31,isub)=field_up(2,ip11,ip22,ip31,isub)+&
&a11a22a31*up(2,iq)
field_up(1,ip12,ip22,ip31,isub)=field_up(1,ip12,ip22,ip31,isub)+&
&a12a22a31*up(1,iq)
field_up(2,ip12,ip22,ip31,isub)=field_up(2,ip12,ip22,ip31,isub)+&
&a12a22a31*up(2,iq)
field_up(1,ip13,ip22,ip31,isub)=field_up(1,ip13,ip22,ip31,isub)+&
&a13a22a31*up(1,iq)
field_up(2,ip13,ip22,ip31,isub)=field_up(2,ip13,ip22,ip31,isub)+&
&a13a22a31*up(2,iq)
field_up(1,ip10,ip23,ip31,isub)=field_up(1,ip10,ip23,ip31,isub)+&
&a10a23a31*up(1,iq)
field_up(2,ip10,ip23,ip31,isub)=field_up(2,ip10,ip23,ip31,isub)+&
&a10a23a31*up(2,iq)
field_up(1,ip11,ip23,ip31,isub)=field_up(1,ip11,ip23,ip31,isub)+&
&a11a23a31*up(1,iq)
field_up(2,ip11,ip23,ip31,isub)=field_up(2,ip11,ip23,ip31,isub)+&
&a11a23a31*up(2,iq)
field_up(1,ip12,ip23,ip31,isub)=field_up(1,ip12,ip23,ip31,isub)+&
&a12a23a31*up(1,iq)
field_up(2,ip12,ip23,ip31,isub)=field_up(2,ip12,ip23,ip31,isub)+&
&a12a23a31*up(2,iq)
field_up(1,ip13,ip23,ip31,isub)=field_up(1,ip13,ip23,ip31,isub)+&
&a13a23a31*up(1,iq)
field_up(2,ip13,ip23,ip31,isub)=field_up(2,ip13,ip23,ip31,isub)+&
&a13a23a31*up(2,iq)
field_up(1,ip10,ip20,ip32,isub)=field_up(1,ip10,ip20,ip32,isub)+&
&a10a20a32*up(1,iq)
field_up(2,ip10,ip20,ip32,isub)=field_up(2,ip10,ip20,ip32,isub)+&
&a10a20a32*up(2,iq)
field_up(1,ip11,ip20,ip32,isub)=field_up(1,ip11,ip20,ip32,isub)+&
&a11a20a32*up(1,iq)
field_up(2,ip11,ip20,ip32,isub)=field_up(2,ip11,ip20,ip32,isub)+&
&a11a20a32*up(2,iq)
field_up(1,ip12,ip20,ip32,isub)=field_up(1,ip12,ip20,ip32,isub)+&
&a12a20a32*up(1,iq)
field_up(2,ip12,ip20,ip32,isub)=field_up(2,ip12,ip20,ip32,isub)+&
&a12a20a32*up(2,iq)
field_up(1,ip13,ip20,ip32,isub)=field_up(1,ip13,ip20,ip32,isub)+&
&a13a20a32*up(1,iq)
field_up(2,ip13,ip20,ip32,isub)=field_up(2,ip13,ip20,ip32,isub)+&
&a13a20a32*up(2,iq)
field_up(1,ip10,ip21,ip32,isub)=field_up(1,ip10,ip21,ip32,isub)+&
&a10a21a32*up(1,iq)
field_up(2,ip10,ip21,ip32,isub)=field_up(2,ip10,ip21,ip32,isub)+&
&a10a21a32*up(2,iq)
field_up(1,ip11,ip21,ip32,isub)=field_up(1,ip11,ip21,ip32,isub)+&
&a11a21a32*up(1,iq)
field_up(2,ip11,ip21,ip32,isub)=field_up(2,ip11,ip21,ip32,isub)+&
&a11a21a32*up(2,iq)
field_up(1,ip12,ip21,ip32,isub)=field_up(1,ip12,ip21,ip32,isub)+&
&a12a21a32*up(1,iq)
field_up(2,ip12,ip21,ip32,isub)=field_up(2,ip12,ip21,ip32,isub)+&
&a12a21a32*up(2,iq)
field_up(1,ip13,ip21,ip32,isub)=field_up(1,ip13,ip21,ip32,isub)+&
&a13a21a32*up(1,iq)
field_up(2,ip13,ip21,ip32,isub)=field_up(2,ip13,ip21,ip32,isub)+&
&a13a21a32*up(2,iq)
field_up(1,ip10,ip22,ip32,isub)=field_up(1,ip10,ip22,ip32,isub)+&
&a10a22a32*up(1,iq)
field_up(2,ip10,ip22,ip32,isub)=field_up(2,ip10,ip22,ip32,isub)+&
&a10a22a32*up(2,iq)
field_up(1,ip11,ip22,ip32,isub)=field_up(1,ip11,ip22,ip32,isub)+&
&a11a22a32*up(1,iq)
field_up(2,ip11,ip22,ip32,isub)=field_up(2,ip11,ip22,ip32,isub)+&
&a11a22a32*up(2,iq)
field_up(1,ip12,ip22,ip32,isub)=field_up(1,ip12,ip22,ip32,isub)+&
&a12a22a32*up(1,iq)
field_up(2,ip12,ip22,ip32,isub)=field_up(2,ip12,ip22,ip32,isub)+&
&a12a22a32*up(2,iq)
field_up(1,ip13,ip22,ip32,isub)=field_up(1,ip13,ip22,ip32,isub)+&
&a13a22a32*up(1,iq)
field_up(2,ip13,ip22,ip32,isub)=field_up(2,ip13,ip22,ip32,isub)+&
&a13a22a32*up(2,iq)
field_up(1,ip10,ip23,ip32,isub)=field_up(1,ip10,ip23,ip32,isub)+&
&a10a23a32*up(1,iq)
field_up(2,ip10,ip23,ip32,isub)=field_up(2,ip10,ip23,ip32,isub)+&
&a10a23a32*up(2,iq)
field_up(1,ip11,ip23,ip32,isub)=field_up(1,ip11,ip23,ip32,isub)+&
&a11a23a32*up(1,iq)
field_up(2,ip11,ip23,ip32,isub)=field_up(2,ip11,ip23,ip32,isub)+&
&a11a23a32*up(2,iq)
field_up(1,ip12,ip23,ip32,isub)=field_up(1,ip12,ip23,ip32,isub)+&
&a12a23a32*up(1,iq)
field_up(2,ip12,ip23,ip32,isub)=field_up(2,ip12,ip23,ip32,isub)+&
&a12a23a32*up(2,iq)
field_up(1,ip13,ip23,ip32,isub)=field_up(1,ip13,ip23,ip32,isub)+&
&a13a23a32*up(1,iq)
field_up(2,ip13,ip23,ip32,isub)=field_up(2,ip13,ip23,ip32,isub)+&
&a13a23a32*up(2,iq)
field_up(1,ip10,ip20,ip33,isub)=field_up(1,ip10,ip20,ip33,isub)+&
&a10a20a33*up(1,iq)
field_up(2,ip10,ip20,ip33,isub)=field_up(2,ip10,ip20,ip33,isub)+&
&a10a20a33*up(2,iq)
field_up(1,ip11,ip20,ip33,isub)=field_up(1,ip11,ip20,ip33,isub)+&
&a11a20a33*up(1,iq)
field_up(2,ip11,ip20,ip33,isub)=field_up(2,ip11,ip20,ip33,isub)+&
&a11a20a33*up(2,iq)
field_up(1,ip12,ip20,ip33,isub)=field_up(1,ip12,ip20,ip33,isub)+&
&a12a20a33*up(1,iq)
field_up(2,ip12,ip20,ip33,isub)=field_up(2,ip12,ip20,ip33,isub)+&
&a12a20a33*up(2,iq)
field_up(1,ip13,ip20,ip33,isub)=field_up(1,ip13,ip20,ip33,isub)+&
&a13a20a33*up(1,iq)
field_up(2,ip13,ip20,ip33,isub)=field_up(2,ip13,ip20,ip33,isub)+&
&a13a20a33*up(2,iq)
field_up(1,ip10,ip21,ip33,isub)=field_up(1,ip10,ip21,ip33,isub)+&
&a10a21a33*up(1,iq)
field_up(2,ip10,ip21,ip33,isub)=field_up(2,ip10,ip21,ip33,isub)+&
&a10a21a33*up(2,iq)
field_up(1,ip11,ip21,ip33,isub)=field_up(1,ip11,ip21,ip33,isub)+&
&a11a21a33*up(1,iq)
field_up(2,ip11,ip21,ip33,isub)=field_up(2,ip11,ip21,ip33,isub)+&
&a11a21a33*up(2,iq)
field_up(1,ip12,ip21,ip33,isub)=field_up(1,ip12,ip21,ip33,isub)+&
&a12a21a33*up(1,iq)
field_up(2,ip12,ip21,ip33,isub)=field_up(2,ip12,ip21,ip33,isub)+&
&a12a21a33*up(2,iq)
field_up(1,ip13,ip21,ip33,isub)=field_up(1,ip13,ip21,ip33,isub)+&
&a13a21a33*up(1,iq)
field_up(2,ip13,ip21,ip33,isub)=field_up(2,ip13,ip21,ip33,isub)+&
&a13a21a33*up(2,iq)
field_up(1,ip10,ip22,ip33,isub)=field_up(1,ip10,ip22,ip33,isub)+&
&a10a22a33*up(1,iq)
field_up(2,ip10,ip22,ip33,isub)=field_up(2,ip10,ip22,ip33,isub)+&
&a10a22a33*up(2,iq)
field_up(1,ip11,ip22,ip33,isub)=field_up(1,ip11,ip22,ip33,isub)+&
&a11a22a33*up(1,iq)
field_up(2,ip11,ip22,ip33,isub)=field_up(2,ip11,ip22,ip33,isub)+&
&a11a22a33*up(2,iq)
field_up(1,ip12,ip22,ip33,isub)=field_up(1,ip12,ip22,ip33,isub)+&
&a12a22a33*up(1,iq)
field_up(2,ip12,ip22,ip33,isub)=field_up(2,ip12,ip22,ip33,isub)+&
&a12a22a33*up(2,iq)
field_up(1,ip13,ip22,ip33,isub)=field_up(1,ip13,ip22,ip33,isub)+&
&a13a22a33*up(1,iq)
field_up(2,ip13,ip22,ip33,isub)=field_up(2,ip13,ip22,ip33,isub)+&
&a13a22a33*up(2,iq)
field_up(1,ip10,ip23,ip33,isub)=field_up(1,ip10,ip23,ip33,isub)+&
&a10a23a33*up(1,iq)
field_up(2,ip10,ip23,ip33,isub)=field_up(2,ip10,ip23,ip33,isub)+&
&a10a23a33*up(2,iq)
field_up(1,ip11,ip23,ip33,isub)=field_up(1,ip11,ip23,ip33,isub)+&
&a11a23a33*up(1,iq)
field_up(2,ip11,ip23,ip33,isub)=field_up(2,ip11,ip23,ip33,isub)+&
&a11a23a33*up(2,iq)
field_up(1,ip12,ip23,ip33,isub)=field_up(1,ip12,ip23,ip33,isub)+&
&a12a23a33*up(1,iq)
field_up(2,ip12,ip23,ip33,isub)=field_up(2,ip12,ip23,ip33,isub)+&
&a12a23a33*up(2,iq)
field_up(1,ip13,ip23,ip33,isub)=field_up(1,ip13,ip23,ip33,isub)+&
&a13a23a33*up(1,iq)
field_up(2,ip13,ip23,ip33,isub)=field_up(2,ip13,ip23,ip33,isub)+&
&a13a23a33*up(2,iq)
#else
field_up(1:2,ip10,ip20,ip30,isub)=field_up(1:2,ip10,ip20,ip30,isub)+&
&a10a20a30*up(1:2,iq)
field_up(1:2,ip11,ip20,ip30,isub)=field_up(1:2,ip11,ip20,ip30,isub)+&
&a11a20a30*up(1:2,iq)
field_up(1:2,ip12,ip20,ip30,isub)=field_up(1:2,ip12,ip20,ip30,isub)+&
&a12a20a30*up(1:2,iq)
field_up(1:2,ip13,ip20,ip30,isub)=field_up(1:2,ip13,ip20,ip30,isub)+&
&a13a20a30*up(1:2,iq)
field_up(1:2,ip10,ip21,ip30,isub)=field_up(1:2,ip10,ip21,ip30,isub)+&
&a10a21a30*up(1:2,iq)
field_up(1:2,ip11,ip21,ip30,isub)=field_up(1:2,ip11,ip21,ip30,isub)+&
&a11a21a30*up(1:2,iq)
field_up(1:2,ip12,ip21,ip30,isub)=field_up(1:2,ip12,ip21,ip30,isub)+&
&a12a21a30*up(1:2,iq)
field_up(1:2,ip13,ip21,ip30,isub)=field_up(1:2,ip13,ip21,ip30,isub)+&
&a13a21a30*up(1:2,iq)
field_up(1:2,ip10,ip22,ip30,isub)=field_up(1:2,ip10,ip22,ip30,isub)+&
&a10a22a30*up(1:2,iq)
field_up(1:2,ip11,ip22,ip30,isub)=field_up(1:2,ip11,ip22,ip30,isub)+&
&a11a22a30*up(1:2,iq)
field_up(1:2,ip12,ip22,ip30,isub)=field_up(1:2,ip12,ip22,ip30,isub)+&
&a12a22a30*up(1:2,iq)
field_up(1:2,ip13,ip22,ip30,isub)=field_up(1:2,ip13,ip22,ip30,isub)+&
&a13a22a30*up(1:2,iq)
field_up(1:2,ip10,ip23,ip30,isub)=field_up(1:2,ip10,ip23,ip30,isub)+&
&a10a23a30*up(1:2,iq)
field_up(1:2,ip11,ip23,ip30,isub)=field_up(1:2,ip11,ip23,ip30,isub)+&
&a11a23a30*up(1:2,iq)
field_up(1:2,ip12,ip23,ip30,isub)=field_up(1:2,ip12,ip23,ip30,isub)+&
&a12a23a30*up(1:2,iq)
field_up(1:2,ip13,ip23,ip30,isub)=field_up(1:2,ip13,ip23,ip30,isub)+&
&a13a23a30*up(1:2,iq)
field_up(1:2,ip10,ip20,ip31,isub)=field_up(1:2,ip10,ip20,ip31,isub)+&
&a10a20a31*up(1:2,iq)
field_up(1:2,ip11,ip20,ip31,isub)=field_up(1:2,ip11,ip20,ip31,isub)+&
&a11a20a31*up(1:2,iq)
field_up(1:2,ip12,ip20,ip31,isub)=field_up(1:2,ip12,ip20,ip31,isub)+&
&a12a20a31*up(1:2,iq)
field_up(1:2,ip13,ip20,ip31,isub)=field_up(1:2,ip13,ip20,ip31,isub)+&
&a13a20a31*up(1:2,iq)
field_up(1:2,ip10,ip21,ip31,isub)=field_up(1:2,ip10,ip21,ip31,isub)+&
&a10a21a31*up(1:2,iq)
field_up(1:2,ip11,ip21,ip31,isub)=field_up(1:2,ip11,ip21,ip31,isub)+&
&a11a21a31*up(1:2,iq)
field_up(1:2,ip12,ip21,ip31,isub)=field_up(1:2,ip12,ip21,ip31,isub)+&
&a12a21a31*up(1:2,iq)
field_up(1:2,ip13,ip21,ip31,isub)=field_up(1:2,ip13,ip21,ip31,isub)+&
&a13a21a31*up(1:2,iq)
field_up(1:2,ip10,ip22,ip31,isub)=field_up(1:2,ip10,ip22,ip31,isub)+&
&a10a22a31*up(1:2,iq)
field_up(1:2,ip11,ip22,ip31,isub)=field_up(1:2,ip11,ip22,ip31,isub)+&
&a11a22a31*up(1:2,iq)
field_up(1:2,ip12,ip22,ip31,isub)=field_up(1:2,ip12,ip22,ip31,isub)+&
&a12a22a31*up(1:2,iq)
field_up(1:2,ip13,ip22,ip31,isub)=field_up(1:2,ip13,ip22,ip31,isub)+&
&a13a22a31*up(1:2,iq)
field_up(1:2,ip10,ip23,ip31,isub)=field_up(1:2,ip10,ip23,ip31,isub)+&
&a10a23a31*up(1:2,iq)
field_up(1:2,ip11,ip23,ip31,isub)=field_up(1:2,ip11,ip23,ip31,isub)+&
&a11a23a31*up(1:2,iq)
field_up(1:2,ip12,ip23,ip31,isub)=field_up(1:2,ip12,ip23,ip31,isub)+&
&a12a23a31*up(1:2,iq)
field_up(1:2,ip13,ip23,ip31,isub)=field_up(1:2,ip13,ip23,ip31,isub)+&
&a13a23a31*up(1:2,iq)
field_up(1:2,ip10,ip20,ip32,isub)=field_up(1:2,ip10,ip20,ip32,isub)+&
&a10a20a32*up(1:2,iq)
field_up(1:2,ip11,ip20,ip32,isub)=field_up(1:2,ip11,ip20,ip32,isub)+&
&a11a20a32*up(1:2,iq)
field_up(1:2,ip12,ip20,ip32,isub)=field_up(1:2,ip12,ip20,ip32,isub)+&
&a12a20a32*up(1:2,iq)
field_up(1:2,ip13,ip20,ip32,isub)=field_up(1:2,ip13,ip20,ip32,isub)+&
&a13a20a32*up(1:2,iq)
field_up(1:2,ip10,ip21,ip32,isub)=field_up(1:2,ip10,ip21,ip32,isub)+&
&a10a21a32*up(1:2,iq)
field_up(1:2,ip11,ip21,ip32,isub)=field_up(1:2,ip11,ip21,ip32,isub)+&
&a11a21a32*up(1:2,iq)
field_up(1:2,ip12,ip21,ip32,isub)=field_up(1:2,ip12,ip21,ip32,isub)+&
&a12a21a32*up(1:2,iq)
field_up(1:2,ip13,ip21,ip32,isub)=field_up(1:2,ip13,ip21,ip32,isub)+&
&a13a21a32*up(1:2,iq)
field_up(1:2,ip10,ip22,ip32,isub)=field_up(1:2,ip10,ip22,ip32,isub)+&
&a10a22a32*up(1:2,iq)
field_up(1:2,ip11,ip22,ip32,isub)=field_up(1:2,ip11,ip22,ip32,isub)+&
&a11a22a32*up(1:2,iq)
field_up(1:2,ip12,ip22,ip32,isub)=field_up(1:2,ip12,ip22,ip32,isub)+&
&a12a22a32*up(1:2,iq)
field_up(1:2,ip13,ip22,ip32,isub)=field_up(1:2,ip13,ip22,ip32,isub)+&
&a13a22a32*up(1:2,iq)
field_up(1:2,ip10,ip23,ip32,isub)=field_up(1:2,ip10,ip23,ip32,isub)+&
&a10a23a32*up(1:2,iq)
field_up(1:2,ip11,ip23,ip32,isub)=field_up(1:2,ip11,ip23,ip32,isub)+&
&a11a23a32*up(1:2,iq)
field_up(1:2,ip12,ip23,ip32,isub)=field_up(1:2,ip12,ip23,ip32,isub)+&
&a12a23a32*up(1:2,iq)
field_up(1:2,ip13,ip23,ip32,isub)=field_up(1:2,ip13,ip23,ip32,isub)+&
&a13a23a32*up(1:2,iq)
field_up(1:2,ip10,ip20,ip33,isub)=field_up(1:2,ip10,ip20,ip33,isub)+&
&a10a20a33*up(1:2,iq)
field_up(1:2,ip11,ip20,ip33,isub)=field_up(1:2,ip11,ip20,ip33,isub)+&
&a11a20a33*up(1:2,iq)
field_up(1:2,ip12,ip20,ip33,isub)=field_up(1:2,ip12,ip20,ip33,isub)+&
&a12a20a33*up(1:2,iq)
field_up(1:2,ip13,ip20,ip33,isub)=field_up(1:2,ip13,ip20,ip33,isub)+&
&a13a20a33*up(1:2,iq)
field_up(1:2,ip10,ip21,ip33,isub)=field_up(1:2,ip10,ip21,ip33,isub)+&
&a10a21a33*up(1:2,iq)
field_up(1:2,ip11,ip21,ip33,isub)=field_up(1:2,ip11,ip21,ip33,isub)+&
&a11a21a33*up(1:2,iq)
field_up(1:2,ip12,ip21,ip33,isub)=field_up(1:2,ip12,ip21,ip33,isub)+&
&a12a21a33*up(1:2,iq)
field_up(1:2,ip13,ip21,ip33,isub)=field_up(1:2,ip13,ip21,ip33,isub)+&
&a13a21a33*up(1:2,iq)
field_up(1:2,ip10,ip22,ip33,isub)=field_up(1:2,ip10,ip22,ip33,isub)+&
&a10a22a33*up(1:2,iq)
field_up(1:2,ip11,ip22,ip33,isub)=field_up(1:2,ip11,ip22,ip33,isub)+&
&a11a22a33*up(1:2,iq)
field_up(1:2,ip12,ip22,ip33,isub)=field_up(1:2,ip12,ip22,ip33,isub)+&
&a12a22a33*up(1:2,iq)
field_up(1:2,ip13,ip22,ip33,isub)=field_up(1:2,ip13,ip22,ip33,isub)+&
&a13a22a33*up(1:2,iq)
field_up(1:2,ip10,ip23,ip33,isub)=field_up(1:2,ip10,ip23,ip33,isub)+&
&a10a23a33*up(1:2,iq)
field_up(1:2,ip11,ip23,ip33,isub)=field_up(1:2,ip11,ip23,ip33,isub)+&
&a11a23a33*up(1:2,iq)
field_up(1:2,ip12,ip23,ip33,isub)=field_up(1:2,ip12,ip23,ip33,isub)+&
&a12a23a33*up(1:2,iq)
field_up(1:2,ip13,ip23,ip33,isub)=field_up(1:2,ip13,ip23,ip33,isub)+&
&a13a23a33*up(1:2,iq)

#endif
           ENDDO
        !------------------------------------------------------------------------
        !  Unrolled versions for 1-vectors
        !------------------------------------------------------------------------
        ELSEIF (lda .EQ. 1) THEN
            DO ip = 1,store_info(isub)

               isubl = topo%isublist(isub)
               iq    = list_sub(isub,ip)

               !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
               !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
               !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

               x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

               !ip11 = INT(x0(1))+2-istart(1,isubl)
               !ip21 = INT(x0(2))+2-istart(2,isubl)
               !ip31 = INT(x0(3))+2-istart(3,isubl)

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
               PRINT*,'###################################################'
               PRINT*,'RENORMALIZATION NOT IMPLEMENTED FOR LDA = 1'
               PRINT*,'###################################################'
               STOP

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
&a10a20a30*up(1,iq)
field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
&a11a20a30*up(1,iq)
field_up(1,ip12,ip20,ip30,isub)=field_up(1,ip12,ip20,ip30,isub)+&
&a12a20a30*up(1,iq)
field_up(1,ip13,ip20,ip30,isub)=field_up(1,ip13,ip20,ip30,isub)+&
&a13a20a30*up(1,iq)
field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
&a10a21a30*up(1,iq)
field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
&a11a21a30*up(1,iq)
field_up(1,ip12,ip21,ip30,isub)=field_up(1,ip12,ip21,ip30,isub)+&
&a12a21a30*up(1,iq)
field_up(1,ip13,ip21,ip30,isub)=field_up(1,ip13,ip21,ip30,isub)+&
&a13a21a30*up(1,iq)
field_up(1,ip10,ip22,ip30,isub)=field_up(1,ip10,ip22,ip30,isub)+&
&a10a22a30*up(1,iq)
field_up(1,ip11,ip22,ip30,isub)=field_up(1,ip11,ip22,ip30,isub)+&
&a11a22a30*up(1,iq)
field_up(1,ip12,ip22,ip30,isub)=field_up(1,ip12,ip22,ip30,isub)+&
&a12a22a30*up(1,iq)
field_up(1,ip13,ip22,ip30,isub)=field_up(1,ip13,ip22,ip30,isub)+&
&a13a22a30*up(1,iq)
field_up(1,ip10,ip23,ip30,isub)=field_up(1,ip10,ip23,ip30,isub)+&
&a10a23a30*up(1,iq)
field_up(1,ip11,ip23,ip30,isub)=field_up(1,ip11,ip23,ip30,isub)+&
&a11a23a30*up(1,iq)
field_up(1,ip12,ip23,ip30,isub)=field_up(1,ip12,ip23,ip30,isub)+&
&a12a23a30*up(1,iq)
field_up(1,ip13,ip23,ip30,isub)=field_up(1,ip13,ip23,ip30,isub)+&
&a13a23a30*up(1,iq)
field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
&a10a20a31*up(1,iq)
field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
&a11a20a31*up(1,iq)
field_up(1,ip12,ip20,ip31,isub)=field_up(1,ip12,ip20,ip31,isub)+&
&a12a20a31*up(1,iq)
field_up(1,ip13,ip20,ip31,isub)=field_up(1,ip13,ip20,ip31,isub)+&
&a13a20a31*up(1,iq)
field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
&a10a21a31*up(1,iq)
field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
&a11a21a31*up(1,iq)
field_up(1,ip12,ip21,ip31,isub)=field_up(1,ip12,ip21,ip31,isub)+&
&a12a21a31*up(1,iq)
field_up(1,ip13,ip21,ip31,isub)=field_up(1,ip13,ip21,ip31,isub)+&
&a13a21a31*up(1,iq)
field_up(1,ip10,ip22,ip31,isub)=field_up(1,ip10,ip22,ip31,isub)+&
&a10a22a31*up(1,iq)
field_up(1,ip11,ip22,ip31,isub)=field_up(1,ip11,ip22,ip31,isub)+&
&a11a22a31*up(1,iq)
field_up(1,ip12,ip22,ip31,isub)=field_up(1,ip12,ip22,ip31,isub)+&
&a12a22a31*up(1,iq)
field_up(1,ip13,ip22,ip31,isub)=field_up(1,ip13,ip22,ip31,isub)+&
&a13a22a31*up(1,iq)
field_up(1,ip10,ip23,ip31,isub)=field_up(1,ip10,ip23,ip31,isub)+&
&a10a23a31*up(1,iq)
field_up(1,ip11,ip23,ip31,isub)=field_up(1,ip11,ip23,ip31,isub)+&
&a11a23a31*up(1,iq)
field_up(1,ip12,ip23,ip31,isub)=field_up(1,ip12,ip23,ip31,isub)+&
&a12a23a31*up(1,iq)
field_up(1,ip13,ip23,ip31,isub)=field_up(1,ip13,ip23,ip31,isub)+&
&a13a23a31*up(1,iq)
field_up(1,ip10,ip20,ip32,isub)=field_up(1,ip10,ip20,ip32,isub)+&
&a10a20a32*up(1,iq)
field_up(1,ip11,ip20,ip32,isub)=field_up(1,ip11,ip20,ip32,isub)+&
&a11a20a32*up(1,iq)
field_up(1,ip12,ip20,ip32,isub)=field_up(1,ip12,ip20,ip32,isub)+&
&a12a20a32*up(1,iq)
field_up(1,ip13,ip20,ip32,isub)=field_up(1,ip13,ip20,ip32,isub)+&
&a13a20a32*up(1,iq)
field_up(1,ip10,ip21,ip32,isub)=field_up(1,ip10,ip21,ip32,isub)+&
&a10a21a32*up(1,iq)
field_up(1,ip11,ip21,ip32,isub)=field_up(1,ip11,ip21,ip32,isub)+&
&a11a21a32*up(1,iq)
field_up(1,ip12,ip21,ip32,isub)=field_up(1,ip12,ip21,ip32,isub)+&
&a12a21a32*up(1,iq)
field_up(1,ip13,ip21,ip32,isub)=field_up(1,ip13,ip21,ip32,isub)+&
&a13a21a32*up(1,iq)
field_up(1,ip10,ip22,ip32,isub)=field_up(1,ip10,ip22,ip32,isub)+&
&a10a22a32*up(1,iq)
field_up(1,ip11,ip22,ip32,isub)=field_up(1,ip11,ip22,ip32,isub)+&
&a11a22a32*up(1,iq)
field_up(1,ip12,ip22,ip32,isub)=field_up(1,ip12,ip22,ip32,isub)+&
&a12a22a32*up(1,iq)
field_up(1,ip13,ip22,ip32,isub)=field_up(1,ip13,ip22,ip32,isub)+&
&a13a22a32*up(1,iq)
field_up(1,ip10,ip23,ip32,isub)=field_up(1,ip10,ip23,ip32,isub)+&
&a10a23a32*up(1,iq)
field_up(1,ip11,ip23,ip32,isub)=field_up(1,ip11,ip23,ip32,isub)+&
&a11a23a32*up(1,iq)
field_up(1,ip12,ip23,ip32,isub)=field_up(1,ip12,ip23,ip32,isub)+&
&a12a23a32*up(1,iq)
field_up(1,ip13,ip23,ip32,isub)=field_up(1,ip13,ip23,ip32,isub)+&
&a13a23a32*up(1,iq)
field_up(1,ip10,ip20,ip33,isub)=field_up(1,ip10,ip20,ip33,isub)+&
&a10a20a33*up(1,iq)
field_up(1,ip11,ip20,ip33,isub)=field_up(1,ip11,ip20,ip33,isub)+&
&a11a20a33*up(1,iq)
field_up(1,ip12,ip20,ip33,isub)=field_up(1,ip12,ip20,ip33,isub)+&
&a12a20a33*up(1,iq)
field_up(1,ip13,ip20,ip33,isub)=field_up(1,ip13,ip20,ip33,isub)+&
&a13a20a33*up(1,iq)
field_up(1,ip10,ip21,ip33,isub)=field_up(1,ip10,ip21,ip33,isub)+&
&a10a21a33*up(1,iq)
field_up(1,ip11,ip21,ip33,isub)=field_up(1,ip11,ip21,ip33,isub)+&
&a11a21a33*up(1,iq)
field_up(1,ip12,ip21,ip33,isub)=field_up(1,ip12,ip21,ip33,isub)+&
&a12a21a33*up(1,iq)
field_up(1,ip13,ip21,ip33,isub)=field_up(1,ip13,ip21,ip33,isub)+&
&a13a21a33*up(1,iq)
field_up(1,ip10,ip22,ip33,isub)=field_up(1,ip10,ip22,ip33,isub)+&
&a10a22a33*up(1,iq)
field_up(1,ip11,ip22,ip33,isub)=field_up(1,ip11,ip22,ip33,isub)+&
&a11a22a33*up(1,iq)
field_up(1,ip12,ip22,ip33,isub)=field_up(1,ip12,ip22,ip33,isub)+&
&a12a22a33*up(1,iq)
field_up(1,ip13,ip22,ip33,isub)=field_up(1,ip13,ip22,ip33,isub)+&
&a13a22a33*up(1,iq)
field_up(1,ip10,ip23,ip33,isub)=field_up(1,ip10,ip23,ip33,isub)+&
&a10a23a33*up(1,iq)
field_up(1,ip11,ip23,ip33,isub)=field_up(1,ip11,ip23,ip33,isub)+&
&a11a23a33*up(1,iq)
field_up(1,ip12,ip23,ip33,isub)=field_up(1,ip12,ip23,ip33,isub)+&
&a12a23a33*up(1,iq)
field_up(1,ip13,ip23,ip33,isub)=field_up(1,ip13,ip23,ip33,isub)+&
&a13a23a33*up(1,iq)

           ENDDO
        !------------------------------------------------------------------------
        !  All other lda are NOT UNROLLED. Vectorization will be over lda!
        !------------------------------------------------------------------------
        ELSE
            DO ip = 1,store_info(isub)

               isubl = topo%isublist(isub)
               iq    = list_sub(isub,ip)

               !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
               !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
               !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

               x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

               !ip11 = INT(x0(1))+2-istart(1,isubl)
               !ip21 = INT(x0(2))+2-istart(2,isubl)
               !ip31 = INT(x0(3))+2-istart(3,isubl)

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

field_up(ldn,ip10,ip20,ip30,isub)=field_up(ldn,ip10,ip20,ip30,isub)+&
&a10a20a30*up(ldn,iq)
field_up(ldn,ip11,ip20,ip30,isub)=field_up(ldn,ip11,ip20,ip30,isub)+&
&a11a20a30*up(ldn,iq)
field_up(ldn,ip12,ip20,ip30,isub)=field_up(ldn,ip12,ip20,ip30,isub)+&
&a12a20a30*up(ldn,iq)
field_up(ldn,ip13,ip20,ip30,isub)=field_up(ldn,ip13,ip20,ip30,isub)+&
&a13a20a30*up(ldn,iq)
field_up(ldn,ip10,ip21,ip30,isub)=field_up(ldn,ip10,ip21,ip30,isub)+&
&a10a21a30*up(ldn,iq)
field_up(ldn,ip11,ip21,ip30,isub)=field_up(ldn,ip11,ip21,ip30,isub)+&
&a11a21a30*up(ldn,iq)
field_up(ldn,ip12,ip21,ip30,isub)=field_up(ldn,ip12,ip21,ip30,isub)+&
&a12a21a30*up(ldn,iq)
field_up(ldn,ip13,ip21,ip30,isub)=field_up(ldn,ip13,ip21,ip30,isub)+&
&a13a21a30*up(ldn,iq)
field_up(ldn,ip10,ip22,ip30,isub)=field_up(ldn,ip10,ip22,ip30,isub)+&
&a10a22a30*up(ldn,iq)
field_up(ldn,ip11,ip22,ip30,isub)=field_up(ldn,ip11,ip22,ip30,isub)+&
&a11a22a30*up(ldn,iq)
field_up(ldn,ip12,ip22,ip30,isub)=field_up(ldn,ip12,ip22,ip30,isub)+&
&a12a22a30*up(ldn,iq)
field_up(ldn,ip13,ip22,ip30,isub)=field_up(ldn,ip13,ip22,ip30,isub)+&
&a13a22a30*up(ldn,iq)
field_up(ldn,ip10,ip23,ip30,isub)=field_up(ldn,ip10,ip23,ip30,isub)+&
&a10a23a30*up(ldn,iq)
field_up(ldn,ip11,ip23,ip30,isub)=field_up(ldn,ip11,ip23,ip30,isub)+&
&a11a23a30*up(ldn,iq)
field_up(ldn,ip12,ip23,ip30,isub)=field_up(ldn,ip12,ip23,ip30,isub)+&
&a12a23a30*up(ldn,iq)
field_up(ldn,ip13,ip23,ip30,isub)=field_up(ldn,ip13,ip23,ip30,isub)+&
&a13a23a30*up(ldn,iq)
field_up(ldn,ip10,ip20,ip31,isub)=field_up(ldn,ip10,ip20,ip31,isub)+&
&a10a20a31*up(ldn,iq)
field_up(ldn,ip11,ip20,ip31,isub)=field_up(ldn,ip11,ip20,ip31,isub)+&
&a11a20a31*up(ldn,iq)
field_up(ldn,ip12,ip20,ip31,isub)=field_up(ldn,ip12,ip20,ip31,isub)+&
&a12a20a31*up(ldn,iq)
field_up(ldn,ip13,ip20,ip31,isub)=field_up(ldn,ip13,ip20,ip31,isub)+&
&a13a20a31*up(ldn,iq)
field_up(ldn,ip10,ip21,ip31,isub)=field_up(ldn,ip10,ip21,ip31,isub)+&
&a10a21a31*up(ldn,iq)
field_up(ldn,ip11,ip21,ip31,isub)=field_up(ldn,ip11,ip21,ip31,isub)+&
&a11a21a31*up(ldn,iq)
field_up(ldn,ip12,ip21,ip31,isub)=field_up(ldn,ip12,ip21,ip31,isub)+&
&a12a21a31*up(ldn,iq)
field_up(ldn,ip13,ip21,ip31,isub)=field_up(ldn,ip13,ip21,ip31,isub)+&
&a13a21a31*up(ldn,iq)
field_up(ldn,ip10,ip22,ip31,isub)=field_up(ldn,ip10,ip22,ip31,isub)+&
&a10a22a31*up(ldn,iq)
field_up(ldn,ip11,ip22,ip31,isub)=field_up(ldn,ip11,ip22,ip31,isub)+&
&a11a22a31*up(ldn,iq)
field_up(ldn,ip12,ip22,ip31,isub)=field_up(ldn,ip12,ip22,ip31,isub)+&
&a12a22a31*up(ldn,iq)
field_up(ldn,ip13,ip22,ip31,isub)=field_up(ldn,ip13,ip22,ip31,isub)+&
&a13a22a31*up(ldn,iq)
field_up(ldn,ip10,ip23,ip31,isub)=field_up(ldn,ip10,ip23,ip31,isub)+&
&a10a23a31*up(ldn,iq)
field_up(ldn,ip11,ip23,ip31,isub)=field_up(ldn,ip11,ip23,ip31,isub)+&
&a11a23a31*up(ldn,iq)
field_up(ldn,ip12,ip23,ip31,isub)=field_up(ldn,ip12,ip23,ip31,isub)+&
&a12a23a31*up(ldn,iq)
field_up(ldn,ip13,ip23,ip31,isub)=field_up(ldn,ip13,ip23,ip31,isub)+&
&a13a23a31*up(ldn,iq)
field_up(ldn,ip10,ip20,ip32,isub)=field_up(ldn,ip10,ip20,ip32,isub)+&
&a10a20a32*up(ldn,iq)
field_up(ldn,ip11,ip20,ip32,isub)=field_up(ldn,ip11,ip20,ip32,isub)+&
&a11a20a32*up(ldn,iq)
field_up(ldn,ip12,ip20,ip32,isub)=field_up(ldn,ip12,ip20,ip32,isub)+&
&a12a20a32*up(ldn,iq)
field_up(ldn,ip13,ip20,ip32,isub)=field_up(ldn,ip13,ip20,ip32,isub)+&
&a13a20a32*up(ldn,iq)
field_up(ldn,ip10,ip21,ip32,isub)=field_up(ldn,ip10,ip21,ip32,isub)+&
&a10a21a32*up(ldn,iq)
field_up(ldn,ip11,ip21,ip32,isub)=field_up(ldn,ip11,ip21,ip32,isub)+&
&a11a21a32*up(ldn,iq)
field_up(ldn,ip12,ip21,ip32,isub)=field_up(ldn,ip12,ip21,ip32,isub)+&
&a12a21a32*up(ldn,iq)
field_up(ldn,ip13,ip21,ip32,isub)=field_up(ldn,ip13,ip21,ip32,isub)+&
&a13a21a32*up(ldn,iq)
field_up(ldn,ip10,ip22,ip32,isub)=field_up(ldn,ip10,ip22,ip32,isub)+&
&a10a22a32*up(ldn,iq)
field_up(ldn,ip11,ip22,ip32,isub)=field_up(ldn,ip11,ip22,ip32,isub)+&
&a11a22a32*up(ldn,iq)
field_up(ldn,ip12,ip22,ip32,isub)=field_up(ldn,ip12,ip22,ip32,isub)+&
&a12a22a32*up(ldn,iq)
field_up(ldn,ip13,ip22,ip32,isub)=field_up(ldn,ip13,ip22,ip32,isub)+&
&a13a22a32*up(ldn,iq)
field_up(ldn,ip10,ip23,ip32,isub)=field_up(ldn,ip10,ip23,ip32,isub)+&
&a10a23a32*up(ldn,iq)
field_up(ldn,ip11,ip23,ip32,isub)=field_up(ldn,ip11,ip23,ip32,isub)+&
&a11a23a32*up(ldn,iq)
field_up(ldn,ip12,ip23,ip32,isub)=field_up(ldn,ip12,ip23,ip32,isub)+&
&a12a23a32*up(ldn,iq)
field_up(ldn,ip13,ip23,ip32,isub)=field_up(ldn,ip13,ip23,ip32,isub)+&
&a13a23a32*up(ldn,iq)
field_up(ldn,ip10,ip20,ip33,isub)=field_up(ldn,ip10,ip20,ip33,isub)+&
&a10a20a33*up(ldn,iq)
field_up(ldn,ip11,ip20,ip33,isub)=field_up(ldn,ip11,ip20,ip33,isub)+&
&a11a20a33*up(ldn,iq)
field_up(ldn,ip12,ip20,ip33,isub)=field_up(ldn,ip12,ip20,ip33,isub)+&
&a12a20a33*up(ldn,iq)
field_up(ldn,ip13,ip20,ip33,isub)=field_up(ldn,ip13,ip20,ip33,isub)+&
&a13a20a33*up(ldn,iq)
field_up(ldn,ip10,ip21,ip33,isub)=field_up(ldn,ip10,ip21,ip33,isub)+&
&a10a21a33*up(ldn,iq)
field_up(ldn,ip11,ip21,ip33,isub)=field_up(ldn,ip11,ip21,ip33,isub)+&
&a11a21a33*up(ldn,iq)
field_up(ldn,ip12,ip21,ip33,isub)=field_up(ldn,ip12,ip21,ip33,isub)+&
&a12a21a33*up(ldn,iq)
field_up(ldn,ip13,ip21,ip33,isub)=field_up(ldn,ip13,ip21,ip33,isub)+&
&a13a21a33*up(ldn,iq)
field_up(ldn,ip10,ip22,ip33,isub)=field_up(ldn,ip10,ip22,ip33,isub)+&
&a10a22a33*up(ldn,iq)
field_up(ldn,ip11,ip22,ip33,isub)=field_up(ldn,ip11,ip22,ip33,isub)+&
&a11a22a33*up(ldn,iq)
field_up(ldn,ip12,ip22,ip33,isub)=field_up(ldn,ip12,ip22,ip33,isub)+&
&a12a22a33*up(ldn,iq)
field_up(ldn,ip13,ip22,ip33,isub)=field_up(ldn,ip13,ip22,ip33,isub)+&
&a13a22a33*up(ldn,iq)
field_up(ldn,ip10,ip23,ip33,isub)=field_up(ldn,ip10,ip23,ip33,isub)+&
&a10a23a33*up(ldn,iq)
field_up(ldn,ip11,ip23,ip33,isub)=field_up(ldn,ip11,ip23,ip33,isub)+&
&a11a23a33*up(ldn,iq)
field_up(ldn,ip12,ip23,ip33,isub)=field_up(ldn,ip12,ip23,ip33,isub)+&
&a12a23a33*up(ldn,iq)
field_up(ldn,ip13,ip23,ip33,isub)=field_up(ldn,ip13,ip23,ip33,isub)+&
&a13a23a33*up(ldn,iq)

               ENDDO    ! lda
           ENDDO        ! iq
        END IF          ! unrolled lda cases
#elif __MODE == __SCA
            DO ip = 1,store_info(isub)

               isubl = topo%isublist(isub)
               iq    = list_sub(isub,ip)

               !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
               !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
               !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

               x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
               !ip11 = INT(x0(1))+2-istart(1,isubl)
               !ip21 = INT(x0(2))+2-istart(2,isubl)
               !ip31 = INT(x0(3))+2-istart(3,isubl)

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

field_up(ip10,ip20,ip30,isub)=field_up(ip10,ip20,ip30,isub)+&
&a10a20a30*up(iq)
field_up(ip11,ip20,ip30,isub)=field_up(ip11,ip20,ip30,isub)+&
&a11a20a30*up(iq)
field_up(ip12,ip20,ip30,isub)=field_up(ip12,ip20,ip30,isub)+&
&a12a20a30*up(iq)
field_up(ip13,ip20,ip30,isub)=field_up(ip13,ip20,ip30,isub)+&
&a13a20a30*up(iq)
field_up(ip10,ip21,ip30,isub)=field_up(ip10,ip21,ip30,isub)+&
&a10a21a30*up(iq)
field_up(ip11,ip21,ip30,isub)=field_up(ip11,ip21,ip30,isub)+&
&a11a21a30*up(iq)
field_up(ip12,ip21,ip30,isub)=field_up(ip12,ip21,ip30,isub)+&
&a12a21a30*up(iq)
field_up(ip13,ip21,ip30,isub)=field_up(ip13,ip21,ip30,isub)+&
&a13a21a30*up(iq)
field_up(ip10,ip22,ip30,isub)=field_up(ip10,ip22,ip30,isub)+&
&a10a22a30*up(iq)
field_up(ip11,ip22,ip30,isub)=field_up(ip11,ip22,ip30,isub)+&
&a11a22a30*up(iq)
field_up(ip12,ip22,ip30,isub)=field_up(ip12,ip22,ip30,isub)+&
&a12a22a30*up(iq)
field_up(ip13,ip22,ip30,isub)=field_up(ip13,ip22,ip30,isub)+&
&a13a22a30*up(iq)
field_up(ip10,ip23,ip30,isub)=field_up(ip10,ip23,ip30,isub)+&
&a10a23a30*up(iq)
field_up(ip11,ip23,ip30,isub)=field_up(ip11,ip23,ip30,isub)+&
&a11a23a30*up(iq)
field_up(ip12,ip23,ip30,isub)=field_up(ip12,ip23,ip30,isub)+&
&a12a23a30*up(iq)
field_up(ip13,ip23,ip30,isub)=field_up(ip13,ip23,ip30,isub)+&
&a13a23a30*up(iq)
field_up(ip10,ip20,ip31,isub)=field_up(ip10,ip20,ip31,isub)+&
&a10a20a31*up(iq)
field_up(ip11,ip20,ip31,isub)=field_up(ip11,ip20,ip31,isub)+&
&a11a20a31*up(iq)
field_up(ip12,ip20,ip31,isub)=field_up(ip12,ip20,ip31,isub)+&
&a12a20a31*up(iq)
field_up(ip13,ip20,ip31,isub)=field_up(ip13,ip20,ip31,isub)+&
&a13a20a31*up(iq)
field_up(ip10,ip21,ip31,isub)=field_up(ip10,ip21,ip31,isub)+&
&a10a21a31*up(iq)
field_up(ip11,ip21,ip31,isub)=field_up(ip11,ip21,ip31,isub)+&
&a11a21a31*up(iq)
field_up(ip12,ip21,ip31,isub)=field_up(ip12,ip21,ip31,isub)+&
&a12a21a31*up(iq)
field_up(ip13,ip21,ip31,isub)=field_up(ip13,ip21,ip31,isub)+&
&a13a21a31*up(iq)
field_up(ip10,ip22,ip31,isub)=field_up(ip10,ip22,ip31,isub)+&
&a10a22a31*up(iq)
field_up(ip11,ip22,ip31,isub)=field_up(ip11,ip22,ip31,isub)+&
&a11a22a31*up(iq)
field_up(ip12,ip22,ip31,isub)=field_up(ip12,ip22,ip31,isub)+&
&a12a22a31*up(iq)
field_up(ip13,ip22,ip31,isub)=field_up(ip13,ip22,ip31,isub)+&
&a13a22a31*up(iq)
field_up(ip10,ip23,ip31,isub)=field_up(ip10,ip23,ip31,isub)+&
&a10a23a31*up(iq)
field_up(ip11,ip23,ip31,isub)=field_up(ip11,ip23,ip31,isub)+&
&a11a23a31*up(iq)
field_up(ip12,ip23,ip31,isub)=field_up(ip12,ip23,ip31,isub)+&
&a12a23a31*up(iq)
field_up(ip13,ip23,ip31,isub)=field_up(ip13,ip23,ip31,isub)+&
&a13a23a31*up(iq)
field_up(ip10,ip20,ip32,isub)=field_up(ip10,ip20,ip32,isub)+&
&a10a20a32*up(iq)
field_up(ip11,ip20,ip32,isub)=field_up(ip11,ip20,ip32,isub)+&
&a11a20a32*up(iq)
field_up(ip12,ip20,ip32,isub)=field_up(ip12,ip20,ip32,isub)+&
&a12a20a32*up(iq)
field_up(ip13,ip20,ip32,isub)=field_up(ip13,ip20,ip32,isub)+&
&a13a20a32*up(iq)
field_up(ip10,ip21,ip32,isub)=field_up(ip10,ip21,ip32,isub)+&
&a10a21a32*up(iq)
field_up(ip11,ip21,ip32,isub)=field_up(ip11,ip21,ip32,isub)+&
&a11a21a32*up(iq)
field_up(ip12,ip21,ip32,isub)=field_up(ip12,ip21,ip32,isub)+&
&a12a21a32*up(iq)
field_up(ip13,ip21,ip32,isub)=field_up(ip13,ip21,ip32,isub)+&
&a13a21a32*up(iq)
field_up(ip10,ip22,ip32,isub)=field_up(ip10,ip22,ip32,isub)+&
&a10a22a32*up(iq)
field_up(ip11,ip22,ip32,isub)=field_up(ip11,ip22,ip32,isub)+&
&a11a22a32*up(iq)
field_up(ip12,ip22,ip32,isub)=field_up(ip12,ip22,ip32,isub)+&
&a12a22a32*up(iq)
field_up(ip13,ip22,ip32,isub)=field_up(ip13,ip22,ip32,isub)+&
&a13a22a32*up(iq)
field_up(ip10,ip23,ip32,isub)=field_up(ip10,ip23,ip32,isub)+&
&a10a23a32*up(iq)
field_up(ip11,ip23,ip32,isub)=field_up(ip11,ip23,ip32,isub)+&
&a11a23a32*up(iq)
field_up(ip12,ip23,ip32,isub)=field_up(ip12,ip23,ip32,isub)+&
&a12a23a32*up(iq)
field_up(ip13,ip23,ip32,isub)=field_up(ip13,ip23,ip32,isub)+&
&a13a23a32*up(iq)
field_up(ip10,ip20,ip33,isub)=field_up(ip10,ip20,ip33,isub)+&
&a10a20a33*up(iq)
field_up(ip11,ip20,ip33,isub)=field_up(ip11,ip20,ip33,isub)+&
&a11a20a33*up(iq)
field_up(ip12,ip20,ip33,isub)=field_up(ip12,ip20,ip33,isub)+&
&a12a20a33*up(iq)
field_up(ip13,ip20,ip33,isub)=field_up(ip13,ip20,ip33,isub)+&
&a13a20a33*up(iq)
field_up(ip10,ip21,ip33,isub)=field_up(ip10,ip21,ip33,isub)+&
&a10a21a33*up(iq)
field_up(ip11,ip21,ip33,isub)=field_up(ip11,ip21,ip33,isub)+&
&a11a21a33*up(iq)
field_up(ip12,ip21,ip33,isub)=field_up(ip12,ip21,ip33,isub)+&
&a12a21a33*up(iq)
field_up(ip13,ip21,ip33,isub)=field_up(ip13,ip21,ip33,isub)+&
&a13a21a33*up(iq)
field_up(ip10,ip22,ip33,isub)=field_up(ip10,ip22,ip33,isub)+&
&a10a22a33*up(iq)
field_up(ip11,ip22,ip33,isub)=field_up(ip11,ip22,ip33,isub)+&
&a11a22a33*up(iq)
field_up(ip12,ip22,ip33,isub)=field_up(ip12,ip22,ip33,isub)+&
&a12a22a33*up(iq)
field_up(ip13,ip22,ip33,isub)=field_up(ip13,ip22,ip33,isub)+&
&a13a22a33*up(iq)
field_up(ip10,ip23,ip33,isub)=field_up(ip10,ip23,ip33,isub)+&
&a10a23a33*up(iq)
field_up(ip11,ip23,ip33,isub)=field_up(ip11,ip23,ip33,isub)+&
&a11a23a33*up(iq)
field_up(ip12,ip23,ip33,isub)=field_up(ip12,ip23,ip33,isub)+&
&a12a23a33*up(iq)
field_up(ip13,ip23,ip33,isub)=field_up(ip13,ip23,ip33,isub)+&
&a13a23a33*up(iq)
           ENDDO        ! iq
#endif     
         END DO              ! loop over subs
      CASE(ppm_param_rmsh_kernel_bsp2)
        !-----------------------------------------------------------------------
         !  B-spline 2 (Witch hat)
         !-----------------------------------------------------------------------
        DO isub = 1,topo%nsublist

#if __MODE == __VEC
         !------------------------------------------------------------------------
         !  Unrolled versions for 4-vectors
         !------------------------------------------------------------------------
         IF(lda.EQ.4) THEN
             DO ip = 1,store_info(isub)
               
                isubl = topo%isublist(isub)
                iq    = list_sub(isub,ip)

                !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
                !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
                !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

                x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
                !ip11 = INT(x0(1))+2-istart(1,isubl)
                !ip21 = INT(x0(2))+2-istart(2,isubl)
                !ip31 = INT(x0(3))+2-istart(3,isubl)

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

                field_up(1,ip10,ip20,ip30,isub) = field_up(1,ip10,ip20,ip30,isub) + &
                    & a10*a20*a30*up(1,iq)
                field_up(1,ip10,ip20,ip31,isub) = field_up(1,ip10,ip20,ip31,isub) + &
                    & a10*a20*a31*up(1,iq)
                field_up(1,ip10,ip21,ip30,isub) = field_up(1,ip10,ip21,ip30,isub) + &
                    & a10*a21*a30*up(1,iq)
                field_up(1,ip10,ip21,ip31,isub) = field_up(1,ip10,ip21,ip31,isub) + &
                    & a10*a21*a31*up(1,iq)
                field_up(1,ip11,ip20,ip30,isub) = field_up(1,ip11,ip20,ip30,isub) + &
                    & a11*a20*a30*up(1,iq)
                field_up(1,ip11,ip20,ip31,isub) = field_up(1,ip11,ip20,ip31,isub) + &
                    & a11*a20*a31*up(1,iq)
                field_up(1,ip11,ip21,ip30,isub) = field_up(1,ip11,ip21,ip30,isub) + &
                    & a11*a21*a30*up(1,iq)
                field_up(1,ip11,ip21,ip31,isub) = field_up(1,ip11,ip21,ip31,isub) + &
                    & a11*a21*a31*up(1,iq)
                
                field_up(2,ip10,ip20,ip30,isub) = field_up(2,ip10,ip20,ip30,isub) + &
                    & a10*a20*a30*up(2,iq)
                field_up(2,ip10,ip20,ip31,isub) = field_up(2,ip10,ip20,ip31,isub) + &
                    & a10*a20*a31*up(2,iq)
                field_up(2,ip10,ip21,ip30,isub) = field_up(2,ip10,ip21,ip30,isub) + &
                    & a10*a21*a30*up(2,iq)
                field_up(2,ip10,ip21,ip31,isub) = field_up(2,ip10,ip21,ip31,isub) + &
                    & a10*a21*a31*up(2,iq)
                field_up(2,ip11,ip20,ip30,isub) = field_up(2,ip11,ip20,ip30,isub) + &
                    & a11*a20*a30*up(2,iq)
                field_up(2,ip11,ip20,ip31,isub) = field_up(2,ip11,ip20,ip31,isub) + &
                    & a11*a20*a31*up(2,iq)
                field_up(2,ip11,ip21,ip30,isub) = field_up(2,ip11,ip21,ip30,isub) + &
                    & a11*a21*a30*up(2,iq)
                field_up(2,ip11,ip21,ip31,isub) = field_up(2,ip11,ip21,ip31,isub) + &
                    & a11*a21*a31*up(2,iq)
                    
                field_up(3,ip10,ip20,ip30,isub) = field_up(3,ip10,ip20,ip30,isub) + &
                    & a10*a20*a30*up(3,iq)
                field_up(3,ip10,ip20,ip31,isub) = field_up(3,ip10,ip20,ip31,isub) + &
                    & a10*a20*a31*up(3,iq)
                field_up(3,ip10,ip21,ip30,isub) = field_up(3,ip10,ip21,ip30,isub) + &
                    & a10*a21*a30*up(3,iq)
                field_up(3,ip10,ip21,ip31,isub) = field_up(3,ip10,ip21,ip31,isub) + &
                    & a10*a21*a31*up(3,iq)
                field_up(3,ip11,ip20,ip30,isub) = field_up(3,ip11,ip20,ip30,isub) + &
                    & a11*a20*a30*up(3,iq)
                field_up(3,ip11,ip20,ip31,isub) = field_up(3,ip11,ip20,ip31,isub) + &
                    & a11*a20*a31*up(3,iq)
                field_up(3,ip11,ip21,ip30,isub) = field_up(3,ip11,ip21,ip30,isub) + &
                    & a11*a21*a30*up(3,iq)
                field_up(3,ip11,ip21,ip31,isub) = field_up(3,ip11,ip21,ip31,isub) + &
                    & a11*a21*a31*up(3,iq)
                
                field_up(4,ip10,ip20,ip30,isub) = field_up(4,ip10,ip20,ip30,isub) + &
                    & a10*a20*a30*up(4,iq)
                field_up(4,ip10,ip20,ip31,isub) = field_up(4,ip10,ip20,ip31,isub) + &
                    & a10*a20*a31*up(4,iq)
                field_up(4,ip10,ip21,ip30,isub) = field_up(4,ip10,ip21,ip30,isub) + &
                    & a10*a21*a30*up(4,iq)
                field_up(4,ip10,ip21,ip31,isub) = field_up(4,ip10,ip21,ip31,isub) + &
                    & a10*a21*a31*up(4,iq)
                field_up(4,ip11,ip20,ip30,isub) = field_up(4,ip11,ip20,ip30,isub) + &
                    & a11*a20*a30*up(4,iq)
                field_up(4,ip11,ip20,ip31,isub) = field_up(4,ip11,ip20,ip31,isub) + &
                    & a11*a20*a31*up(4,iq)
                field_up(4,ip11,ip21,ip30,isub) = field_up(4,ip11,ip21,ip30,isub) + &
                    & a11*a21*a30*up(4,iq)
                field_up(4,ip11,ip21,ip31,isub) = field_up(4,ip11,ip21,ip31,isub) + &
                    & a11*a21*a31*up(4,iq)
                    
            END DO
        !------------------------------------------------------------------------
        !  Unrolled versions for 3-vectors
        !------------------------------------------------------------------------
         ELSEIF (lda .EQ. 3) THEN
            DO ip = 1,store_info(isub)
               
                isubl = topo%isublist(isub)
                iq    = list_sub(isub,ip)

                !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
                !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
                !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

                x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
                !ip11 = INT(x0(1))+2-istart(1,isubl)
                !ip21 = INT(x0(2))+2-istart(2,isubl)
                !ip31 = INT(x0(3))+2-istart(3,isubl)

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
                
#ifdef __NOMICROINSTRUCTIONS
                field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1,iq)
                field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(2,iq)
                field_up(3,ip10,ip20,ip30,isub)=field_up(3,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(3,iq)
                field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1,iq)
                field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(2,iq)
                field_up(3,ip11,ip20,ip30,isub)=field_up(3,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(3,iq)
                field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1,iq)
                field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(2,iq)
                field_up(3,ip10,ip21,ip30,isub)=field_up(3,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(3,iq)
                field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1,iq)
                field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(2,iq)
                field_up(3,ip11,ip21,ip30,isub)=field_up(3,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(3,iq)
                field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1,iq)
                field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(2,iq)
                field_up(3,ip10,ip20,ip31,isub)=field_up(3,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(3,iq)
                field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1,iq)
                field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(2,iq)
                field_up(3,ip11,ip20,ip31,isub)=field_up(3,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(3,iq)
                field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1,iq)
                field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(2,iq)
                field_up(3,ip10,ip21,ip31,isub)=field_up(3,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(3,iq)
                field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1,iq)
                field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(2,iq)
                field_up(3,ip11,ip21,ip31,isub)=field_up(3,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(3,iq)
#else
                field_up(1:3,ip10,ip20,ip30,isub)=field_up(1:3,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1:3,iq)
                field_up(1:3,ip11,ip20,ip30,isub)=field_up(1:3,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1:3,iq)
                field_up(1:3,ip10,ip21,ip30,isub)=field_up(1:3,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1:3,iq)
                field_up(1:3,ip11,ip21,ip30,isub)=field_up(1:3,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1:3,iq)
                field_up(1:3,ip10,ip20,ip31,isub)=field_up(1:3,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1:3,iq)
                field_up(1:3,ip11,ip20,ip31,isub)=field_up(1:3,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1:3,iq)
                field_up(1:3,ip10,ip21,ip31,isub)=field_up(1:3,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1:3,iq)
                field_up(1:3,ip11,ip21,ip31,isub)=field_up(1:3,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1:3,iq)
#endif
            END DO
        !------------------------------------------------------------------------
        !  Unrolled versions for 2-vectors
        !------------------------------------------------------------------------
        ELSEIF (lda .EQ. 2) THEN
            DO ip = 1,store_info(isub)
               
                isubl = topo%isublist(isub)
                iq    = list_sub(isub,ip)


                !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
                !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
                !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

                x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
                !ip11 = INT(x0(1))+2-istart(1,isubl)
                !ip21 = INT(x0(2))+2-istart(2,isubl)
                !ip31 = INT(x0(3))+2-istart(3,isubl)

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
                
#ifdef __NOMICROINSTRUCTIONS
                field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1,iq)
                field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(2,iq)
                field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1,iq)
                field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(2,iq)
                field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1,iq)
                field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(2,iq)
                field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1,iq)
                field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(2,iq)
                field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1,iq)
                field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(2,iq)
                field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1,iq)
                field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(2,iq)
                field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1,iq)
                field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(2,iq)
                field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1,iq)
                field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(2,iq)
#else
                field_up(1:2,ip10,ip20,ip30,isub)=field_up(1:2,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1:2,iq)
                field_up(1:2,ip11,ip20,ip30,isub)=field_up(1:2,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1:2,iq)
                field_up(1:2,ip10,ip21,ip30,isub)=field_up(1:2,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1:2,iq)
                field_up(1:2,ip11,ip21,ip30,isub)=field_up(1:2,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1:2,iq)
                field_up(1:2,ip10,ip20,ip31,isub)=field_up(1:2,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1:2,iq)
                field_up(1:2,ip11,ip20,ip31,isub)=field_up(1:2,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1:2,iq)
                field_up(1:2,ip10,ip21,ip31,isub)=field_up(1:2,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1:2,iq)
                field_up(1:2,ip11,ip21,ip31,isub)=field_up(1:2,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1:2,iq)
#endif
            ENDDO
        !------------------------------------------------------------------------
        !  Unrolled versions for 1-vectors
        !------------------------------------------------------------------------
        ELSEIF (lda .EQ. 1) THEN
            DO ip = 1,store_info(isub)
               
                isubl = topo%isublist(isub)
                iq    = list_sub(isub,ip)

                !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
                !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
                !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

                x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
                !ip11 = INT(x0(1))+2-istart(1,isubl)
                !ip21 = INT(x0(2))+2-istart(2,isubl)
                !ip31 = INT(x0(3))+2-istart(3,isubl)
                
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

                field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1,iq)
                field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1,iq)
                field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1,iq)
                field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1,iq)
                field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1,iq)
                field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1,iq)
                field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1,iq)
                field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1,iq)
            ENDDO
        !------------------------------------------------------------------------
        !  All other lda are NOT UNROLLED. Vectorization will be over lda!
        !------------------------------------------------------------------------
        ELSE
            DO ip = 1,store_info(isub)
               
                isubl = topo%isublist(isub)
                iq    = list_sub(isub,ip)

                !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
                !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
                !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

                x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
                !ip11 = INT(x0(1))+2-istart(1,isubl)
                !ip21 = INT(x0(2))+2-istart(2,isubl)
                !ip31 = INT(x0(3))+2-istart(3,isubl)
                
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

                    field_up(ldn,ip10,ip20,ip30,isub)=field_up(ldn,ip10,ip20,ip30,isub)+&
                        &a10a20a30*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip30,isub)=field_up(ldn,ip11,ip20,ip30,isub)+&
                        &a11a20a30*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip30,isub)=field_up(ldn,ip10,ip21,ip30,isub)+&
                        &a10a21a30*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip30,isub)=field_up(ldn,ip11,ip21,ip30,isub)+&
                        &a11a21a30*up(ldn,iq)
                    field_up(ldn,ip10,ip20,ip31,isub)=field_up(ldn,ip10,ip20,ip31,isub)+&
                        &a10a20a31*up(ldn,iq)
                    field_up(ldn,ip11,ip20,ip31,isub)=field_up(ldn,ip11,ip20,ip31,isub)+&
                        &a11a20a31*up(ldn,iq)
                    field_up(ldn,ip10,ip21,ip31,isub)=field_up(ldn,ip10,ip21,ip31,isub)+&
                        &a10a21a31*up(ldn,iq)
                    field_up(ldn,ip11,ip21,ip31,isub)=field_up(ldn,ip11,ip21,ip31,isub)+&
                        &a11a21a31*up(ldn,iq)
                ENDDO    ! lda
            ENDDO        ! iq
        END IF          ! unrolled lda cases
#elif __MODE == __SCA
            DO ip = 1,store_info(isub)
               
                isubl = topo%isublist(isub)
                iq    = list_sub(isub,ip)

                !x0(1)  = (xp(1,iq)-min_phys(1,topoid))*dxi(1)
                !x0(2)  = (xp(2,iq)-min_phys(2,topoid))*dxi(2)
                !x0(3)  = (xp(3,iq)-min_phys(3,topoid))*dxi(3)

                x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)
               
                !ip11 = INT(x0(1))+2-istart(1,isubl)
                !ip21 = INT(x0(2))+2-istart(2,isubl)
                !ip31 = INT(x0(3))+2-istart(3,isubl)
            
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


                field_up(ip10,ip20,ip30,isub)=field_up(ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(iq)
                field_up(ip11,ip20,ip30,isub)=field_up(ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(iq)
                field_up(ip10,ip21,ip30,isub)=field_up(ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(iq)
                field_up(ip11,ip21,ip30,isub)=field_up(ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(iq)
                field_up(ip10,ip20,ip31,isub)=field_up(ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(iq)
                field_up(ip11,ip20,ip31,isub)=field_up(ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(iq)
                field_up(ip10,ip21,ip31,isub)=field_up(ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(iq)
                field_up(ip11,ip21,ip31,isub)=field_up(ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(iq)
            ENDDO        ! iq
#endif
         END DO              ! loop over subs
      CASE DEFAULT
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',    &
     &       'Only Mp4 and Bsp2 are available. Use ppm_rmsh_remesh for other kernels.', &
     &       __LINE__,info)
      END SELECT         ! kernel type
#else
      !--------------------------------------------------------------------------
      !  --- 2D --- 
      !--------------------------------------------------------------------------
      !  loop over subs
      ndata => p_mesh%nnodes

      DO isub = 1,topo%nsublist
         isubl = topo%isublist(isub)
         
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
         
         !-----------------------------------------------------------------
         ! M Prime Four
         !-----------------------------------------------------------------
         DO isub = 1,topo%nsublist
            
            DO ip = 1,store_info(isub)
               
               isubl = topo%isublist(isub)
               iq    = list_sub(isub,ip)

               x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
               
               ip1 = INT(x0(1))+1
               ip2 = INT(x0(2))+1

               xp1 = x0(1)-AINT(x0(1))
               xp2 = x0(2)-AINT(x0(2))

               DO jj = -1,2
                  
                  x2 = ABS(xp2 - REAL(jj,mk))
                  
                  IF(x2.LT.1.0_mk) THEN
                     wx2 = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                  ELSE
                     wx2 = 2.0_mk + (-4.0_mk + &
                          &(2.5_mk - 0.5_mk * x2)*x2)*x2
                  END IF
                  
                  DO ii    = - 1,2
                     
                     x1 = ABS(xp1 - REAL(ii,mk))
                     
                     IF(x1.LT.1.0_MK) THEN
                        wx1 =  1.0_mk - x1**2*(2.5_mk - &
                             & 1.5_mk*x1)
                     ELSE
                        wx1 =  2.0_mk + (-4.0_mk + &
                             & (2.5_mk - 0.5_mk*x1)*x1)*x1
                     END IF
#if __MODE == __SCA                    
                     field_up(ii+ip1,jj+ip2,isub) &
                          &= field_up(ii+ip1,jj+ip2,isub) &
                          &                 + wx1*wx2*up(iq)
#else
                     DO ldn=1,lda
                        field_up(ldn,ii+ip1,jj+ip2,isub) &
                             &= field_up(ldn,ii+ip1,jj+ip2,isub) &
                             &                 + wx1*wx2*up(ldn,iq)
                     ENDDO
#endif                    
                  END DO
               END DO
            END DO
         END DO          ! loop over subs
      CASE DEFAULT
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',    &
     &       'Only Mp4 is available. Use ppm_rmsh_remesh for other kernels.', &
     &       __LINE__,info)
      END SELECT
#endif

9997 CONTINUE
      !--------------------------------------------------------------------------
      !  Before we map the ghosts of field_reno we reset the corresponding 
      !  ghost cells, otherwise, we dont do one-sided interpolation
      !--------------------------------------------------------------------------
#if  __MODE == __VEC
      DO isub=1,nsubs
         isubl = topo%isublist(isub)
         IF (min_sub(3,isubl).LT.(min_phys(3)+0.5*dx(3))) THEN
             DO k=1-ghostsize(3),0
                 DO j=1,ndata(2,isubl)
                     DO i=1,ndata(1,isubl)
                         field_reno(i,j,k,isub) = 0.0_mk
                     END DO
                 END DO
             END DO
         END IF
         IF (max_sub(3,isubl).GT.(max_phys(3)-0.5*dx(3))) THEN
             DO k=ndata(3,isubl)+1,ndata(3,isubl)+ghostsize(3)
                 DO j=1,ndata(2,isubl)
                     DO i=1,ndata(1,isubl)
                         field_reno(i,j,k,isub) = 0.0_mk
                     END DO
                 END DO
             END DO
         END IF
      END DO
#endif
      !--------------------------------------------------------------------------
      !  Now map the ghosts in order to get consistent values at the border of
      !  the subdomains.
      !--------------------------------------------------------------------------


      CALL ppm_map_field_ghost_put(topoid,meshid,ghostsize,info)
      IF (info .NE. 0) GOTO 9999

#if   __MODE == __SCA     
      CALL ppm_map_field_push(topoid,meshid,field_up,info)
#elif __MODE == __VEC     
      CALL ppm_map_field_push(topoid,meshid,field_up,lda,info)
      CALL ppm_map_field_push(topoid,meshid,field_reno,info)
#endif     
      IF (info .NE. 0) GOTO 9999

      CALL ppm_map_field_send(info)
      IF (info .NE. 0) GOTO 9999

#if   __MODE == __SCA     
      CALL ppm_map_field_pop(topoid,meshid,field_up,ghostsize,info)
#elif __MODE == __VEC
      CALL ppm_map_field_pop(topoid,meshid,field_reno,ghostsize,info)
      CALL ppm_map_field_pop(topoid,meshid,field_up,lda,ghostsize,info)
#endif
      IF (info .NE. 0) GOTO 9999
      !-----------------------------------------------------
      !  apply the renormalization at the boundary
      !-----------------------------------------------------
#if  __MODE == __VEC
      DO isub=1,nsubs
         isubl = topo%isublist(isub)
         IF (min_sub(3,isubl).LT.(min_phys(3)+0.5*dx(3))) THEN
             DO k=1,ghostsize(3)
                 DO j=1,ndata(2,isubl)
                     DO i=1,ndata(1,isubl)
                     !-----------------------------------------------------
                     !  TODO: have to finish this
                     !-----------------------------------------------------
                     ! TODO: What if the unity samples to something finite
                     ! but negative, this could happen in some cases...
                     !IF(field_reno(i,j,k,isub).GT.0.0_mk) THEN
                     IF(ABS(field_reno(i,j,k,isub)).GT.0.0_mk) THEN
                         field_up(1,i,j,k,isub) = field_up(1,i,j,k,isub) / &
                          & field_reno(i,j,k,isub)
                         field_up(2,i,j,k,isub) = field_up(2,i,j,k,isub) / &
                          & field_reno(i,j,k,isub)
                         field_up(3,i,j,k,isub) = field_up(3,i,j,k,isub) / &
                          & field_reno(i,j,k,isub)
                     END IF
                     END DO
                 END DO
             END DO
         END IF
      END DO
      DO isub=1,nsubs
         isubl = topo%isublist(isub)
         IF (max_sub(3,isubl).GT.(max_phys(3)-0.5*dx(3))) THEN
             DO k=ndata(3,isubl)-ghostsize(3)+1,ndata(3,isubl)
                 DO j=1,ndata(2,isubl)
                     DO i=1,ndata(1,isubl)
                     !-----------------------------------------------------
                     !  TODO: have to finish this
                     !-----------------------------------------------------
                     ! TODO: What if the unity samples to something finite
                     ! but negative, this could happen in some cases...
                     !IF(field_reno(i,j,k,isub).GT.0.0_mk) THEN
                     IF(ABS(field_reno(i,j,k,isub)).GT.0.0_mk) THEN
                         field_up(1,i,j,k,isub) = field_up(1,i,j,k,isub) / &
                             & field_reno(i,j,k,isub)
                         field_up(2,i,j,k,isub) = field_up(2,i,j,k,isub) / &
                             & field_reno(i,j,k,isub)
                         field_up(3,i,j,k,isub) = field_up(3,i,j,k,isub) / &
                             & field_reno(i,j,k,isub)
                     END IF
                     END DO
                 END DO
             END DO
         END IF
      END DO
#endif
      DEALLOCATE(field_reno)
      NULLIFY(field_reno)
      IF (info .NE. 0) GOTO 9999
      IF(np.EQ.0) GOTO 9999
      !------------------------------------------------------------------------
      ! Deallocation of the arrays....
      !------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      ldu(1) = 0
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
              & 'ppm_interp_p2m_renorm', &
              & 'pb in ilist1 deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
              & 'ppm_interp_p2m_renorm',  &
              & 'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(store_info,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
              & 'ppm_interp_p2m_renorm',  &
              & 'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(list_sub,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
              & 'ppm_interp_p2m_renorm',  &
              & 'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      END IF

      !--------------------------------------------------------------------------
      !  Return 
      !--------------------------------------------------------------------------
9999 CONTINUE
      CALL substop('ppm_interp_p2m_renorm',t0,info)
      RETURN
      
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_ppm_noinit,'ppm_interp_p2m_renorm',  &
     &               'Please call ppm_init first!',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (Np .GT. 0) THEN
           IF (SIZE(xp,2) .LT. Np) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',  &
     &                     'not enough particles contained in xp',__LINE__,info)
            GOTO 8888
           ENDIF
           IF (SIZE(xp,1) .LT.dim) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',  &
     &                     'leading dimension of xp insufficient',__LINE__,info)
            GOTO 8888
           ENDIF
        ENDIF
        IF (Np .LE. 0) THEN
           IF (Np .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',  &
     &                     'particles must be specified',__LINE__,info)
            GOTO 8888
           END IF
           GOTO 8888
        END IF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',  &
     &                     'wrong kernel definition',__LINE__,info)
           GOTO 8888
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
     &               .OR.(kernel_support.EQ.6))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',  &
     &                     'wrong kernel support',__LINE__,info)
           GOTO 8888
        END IF
        CALL ppm_check_topoid(topoid,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',  &
     &                 'topo_id is invalid!',__LINE__,info)
           GOTO 8888
        ENDIF
        CALL ppm_check_meshid(topoid,meshid,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m_renorm',  &
     &                 'mesh_id is invalid!',__LINE__,info)
           GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_renorm_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_renorm_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_renorm_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_renorm_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_renorm_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_renorm_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_renorm_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_renorm_dv_3d
#endif
#endif
#endif       

