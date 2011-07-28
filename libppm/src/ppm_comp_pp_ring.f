      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_comp_pp_ring
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Subroutine which computes kernel interactions by
      !                 direct global (N**2) particle-particle
      !                 interactions using the ring topology.
      !
      !  Input        : xp(:,:)     (F) particle positions
      !                 Np          (I) number of particles on local proc.
      !                 pdata(:,:)  (O) particle data (e.g. strengths).
      !                                 Overloaded types:
      !                                 single,double,single complex,double
      !                                 complex.
      !                 lda         (I) leading dimension of pdata.
      !                 lsymm       (L) Whether to use symmetry or not:
      !                                    .F.  PP interaction w/o symmetry
      !                                    .T.  PP interaction w/  symmetry
      !                 kernel      (O) kernel to be used for PP
      !                                 interactions. To use ppm-internal
      !                                 kernels, specify one of:
      !                                   ppm_param_kerel_laplace2d_2p
      !                                     (2nd order Laplacian,
      !                                     polynomial in 2D)
      !                                   ppm_param_kerel_laplace3d_2p
      !                                     (2nd order Laplacian,
      !                                     polynomial in 3D)
      !                                 To use your own kernel function,
      !                                 pass the function pointer here. Your
      !                                 function should take one argument
      !                                 and return one value. The third
      !                                 possibility is to pass a lookup
      !                                 table with tabulated kernel values.
      !                                 Such a table can be created using
      !                                 ppm_comp_pp_mk_table.
      !                 kpar(:)     (O) Kernel parameters. See documentation
      !                                 or ppm_comp_pp_kernels.inc for
      !                                 description. Type can be single,
      !                                 double, single complex or double
      !                                 complex. When using a lookup table,
      !                                 pass dxtableinv (the inverse of the
      !                                 table spacing) as a scalar here.
      !
      !  Input/output : 
      !
      !  Output       : dpd(:,:)    (O) Change of particle data (pdata) due
      !                                 to interaction. Overloaded types:
      !                                 single,double,single complex,double
      !                                 complex.
      !                 info        (I) return status
      !
      !  Remarks      : dpd needs to be allocated to proper size before
      !                 calling this routine. Also, this routine is not
      !                 resetting dpd to zero before doing the PP
      !                 interactions. This allows contributions from
      !                 different kernels to be accumulated. If needed,
      !                 set it to zero before calling this routine the 
      !                 first time.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_comp_pp_ring.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/07/29 12:29:00  oingo
      !  Adjusted the parameter list of ppm_map_part_ring_shift
      !
      !  Revision 1.3  2004/07/26 13:37:39  ivos
      !  bugfix: function pointer argument does not need to be declared
      !  EXTERNAL if there is an explicit interface.
      !
      !  Revision 1.2  2004/07/26 07:45:24  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.1  2004/07/23 12:57:06  ivos
      !  Preliminary checkin due to module clean-up. Not tested yet, but
      !  it compiles.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_ring_si(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_ring_di(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_ring_sci(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_ring_dci(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_ring_su(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_ring_du(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_ring_scu(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_ring_dcu(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_ring_st(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_ring_dt(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_ring_sct(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_ring_dct(xp,Np,pdata,lda,lsymm,kernel,kpar,   &
     &    dpd,info)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_map_part
      USE ppm_module_map_part_ring_shift
      USE ppm_module_comp_pp_doring
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: xp
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: pdata
      REAL(MK)   , DIMENSION(:,:), INTENT(  OUT) :: dpd
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ) :: pdata
      COMPLEX(MK), DIMENSION(:,:), INTENT(  OUT) :: dpd
#endif
      LOGICAL                    , INTENT(IN   ) :: lsymm
      INTEGER                    , INTENT(IN   ) :: Np,lda
#if   __KERNEL == __INTERNAL
      INTEGER                    , INTENT(IN   ) :: kernel
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
#endif
#elif __KERNEL == __LOOKUP_TABLE
      REAL(MK)                   , INTENT(IN   ) :: kpar
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:  ), INTENT(IN   ) :: kernel
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  ), INTENT(IN   ) :: kernel
#endif
#endif
      INTEGER                    , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                              :: i,j
      INTEGER                              :: hops,isource,itarget
      INTEGER                              :: ll,lu,rl,ru
      INTEGER                              :: Lpart,iopt
      INTEGER    , DIMENSION(2)            :: ldu
      REAL(MK)   , DIMENSION(:,:), POINTER :: xp2
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:), POINTER :: pdata2,dpd2
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), POINTER :: pdata2,dpd2
#endif
      REAL(MK)                             :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
#if   __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              REAL(MK), DIMENSION(:), INTENT(IN) :: kpar
              REAL(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:)  , INTENT(IN   ) :: kpar
      INTERFACE
          FUNCTION kernel(x,kpar)
              IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0E0)
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
              INTEGER, PARAMETER :: MK = KIND(1.0D0)
#endif
              REAL(MK), INTENT(IN) :: x
              COMPLEX(MK), DIMENSION(:), INTENT(IN) :: kpar
              COMPLEX(MK) :: kernel
          END FUNCTION kernel
      END INTERFACE
#endif
#endif

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_comp_pp_ring',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_comp_pp_ring',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_ring',  &
     &            'Np must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_ring',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __KERNEL == __LOOKUP_TABLE
          IF (kpar .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_ring',  &
     &            'kpar (dxtableinv) must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the copy of the local particles
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_dim 
      ldu(2) = Np
      CALL ppm_alloc(xp2,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_comp_pp_ring',     &
     &        'local copy of xp XP2',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = lda
      ldu(2) = Np
      CALL ppm_alloc(pdata2,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_comp_pp_ring',     &
     &        'local copy of pdata PDATA2',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(dpd2,ldu,iopt,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_comp_pp_ring',     &
     &        'local copy of dpd DPD2',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Make a copy of all local particles including all additional info
      !  that they carry
      !-------------------------------------------------------------------------
      Lpart = Np
      xp2(1:ppm_dim,1:Np) = xp(1:ppm_dim,1:Np)
      pdata2(1:lda,1:Np) = pdata(1:lda,1:Np)
      dpd2(1:lda,1:Np) = dpd(1:lda,1:Np)

      !-------------------------------------------------------------------------
      !  Calculate the interactions of the local particles with themselves.
      !-------------------------------------------------------------------------
      CALL ppm_comp_pp_doring(xp,pdata,dpd,Np,xp2,pdata2,dpd2,Lpart,lda, &
     &    lsymm,kernel,kpar,1,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_comp_pp_ring',     &
     &        'Interaction calculation failed.',__LINE__,info)
          GOTO 9999
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Compute whom to send the copy and who to receive the copy from
      !-------------------------------------------------------------------------
      itarget = MOD(ppm_nproc+ppm_rank-1,ppm_nproc)
      isource = MOD(ppm_rank+1,ppm_nproc)
      
      !-------------------------------------------------------------------------
      !  Compute how often we have to shift a copy of the particles to the
      !  neighbour processor
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
          IF (MOD(ppm_nproc,2) .EQ. 0) THEN
              hops = (ppm_nproc/2)-1
          ELSE
              hops = ((ppm_nproc+1)/2)-1
          ENDIF
      ELSE
          hops = ppm_nproc-1
      ENDIF

      DO i = 1,hops
         CALL ppm_map_part_ring_shift(xp2,ppm_dim,Lpart,itarget,isource,info)
         CALL ppm_map_part(pdata2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(dpd2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(xp2,ppm_dim,Lpart,Lpart,-1,ppm_param_map_send,info)
         CALL ppm_map_part(dpd2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(pdata2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(xp2,ppm_dim,Lpart,Lpart,-1,ppm_param_map_pop,info)

         CALL ppm_comp_pp_doring(xp,pdata,dpd,Np,xp2,pdata2,dpd2,Lpart,  &
     &       lda,lsymm,kernel,kpar,0,info)
         IF (info .NE. ppm_param_success) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_comp_pp_ring',     &
     &          'Interaction calculation failed.',__LINE__,info)
            GOTO 9999
         ENDIF
      ENDDO
      
      !-------------------------------------------------------------------------
      !  If we have symmetry we only have send the particles half the round,
      !  so we have to check the case of an even number of processors
      !-------------------------------------------------------------------------
      IF (lsymm .AND. (MOD(ppm_nproc,2) .EQ. 0)) THEN
          CALL ppm_map_part_ring_shift(xp2,ppm_dim,Lpart,itarget,isource,info)
          CALL ppm_map_part(pdata2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
          CALL ppm_map_part(dpd2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
          CALL ppm_map_part(xp2,ppm_dim,Lpart,Lpart,-1,ppm_param_map_send,info)
          CALL ppm_map_part(dpd2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
          CALL ppm_map_part(pdata2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
          CALL ppm_map_part(xp2,ppm_dim,Lpart,Lpart,-1,ppm_param_map_pop,info)

          !---------------------------------------------------------------------
          !  The processor with the higher ppm_rank computes the upper half of
          !  the particles and the opposite processor the lower half
          !---------------------------------------------------------------------
          IF (ppm_rank .GT. hops) THEN
              ll = (Np / 2) + 1
              lu = Np
              rl = 1
              ru = Lpart
          ELSE
              ll = 1
              lu = Np
              rl = 1
              ru = Lpart / 2
          ENDIF

          CALL ppm_comp_pp_doring(xp(:,ll:lu),pdata(:,ll:lu),dpd(:,ll:lu),&
     &                          lu-ll+1,xp2(:,rl:ru),pdata2(:,rl:ru),     &
     &                          dpd2(:,rl:ru),ru-rl+1,lda,lsymm,kernel,kpar, &
     &                          0,info)
          IF (info .NE. ppm_param_success) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_comp_pp_ring',     &
     &           'Interaction calculation failed.',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Send the particles where they belong to
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
          itarget = MOD(ppm_rank+(ppm_nproc/2),ppm_nproc)
          isource = MOD(ppm_nproc+ppm_rank-(ppm_nproc/2),ppm_nproc)
      ENDIF

      IF (itarget .NE. ppm_rank) THEN
          CALL ppm_map_part_ring_shift(xp2,ppm_dim,Lpart,itarget,isource,info)
          CALL ppm_map_part(pdata2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
          CALL ppm_map_part(dpd2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
          CALL ppm_map_part(xp2,ppm_dim,Lpart,Lpart,-1,ppm_param_map_send,info)
          CALL ppm_map_part(dpd2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
          CALL ppm_map_part(pdata2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
          CALL ppm_map_part(xp2,ppm_dim,Lpart,Lpart,-1,ppm_param_map_pop,info)
      
          IF (Lpart .NE. Np) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_part_lost,'ppm_comp_pp_ring',     &
     &           'Not all particles came back!',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Add the particle changes of the local and the long traveled copy
      !-------------------------------------------------------------------------
      IF (lda .EQ. 1) THEN
         DO i = 1,Np
            dpd(1,i) = dpd(1,i) + dpd2(1,i)
         ENDDO
      ELSEIF (lda .EQ. 2) THEN
         DO i = 1,Np
            dpd(1,i) = dpd(1,i) + dpd2(1,i)
            dpd(2,i) = dpd(2,i) + dpd2(2,i)
         ENDDO
      ELSE
         DO i = 1,Np
            DO j = 1,lda
               dpd(j,i) = dpd(j,i) + dpd2(j,i)
            ENDDO
         ENDDO
      ENDIF

 9999 CONTINUE

      !-------------------------------------------------------------------------
      !  Deallocate memory of copies
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(dpd2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_comp_pp_ring',     &
     &        'local copy of dpd DPD2',__LINE__,info)
      ENDIF
      CALL ppm_alloc(pdata2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_comp_pp_ring',     &
     &        'local copy of pdata PDATA2',__LINE__,info)
      ENDIF
      CALL ppm_alloc(xp2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_comp_pp_ring',     &
     &        'local copy of xp XP2',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      CALL substop('ppm_comp_pp_ring',t0,info)
      RETURN
#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_ring_si
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_ring_di
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_ring_sci
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_ring_dci
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_ring_su
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_ring_du
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_ring_scu
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_ring_dcu
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_ring_st
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_ring_dt
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_ring_sct
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_ring_dct
#endif
#endif
