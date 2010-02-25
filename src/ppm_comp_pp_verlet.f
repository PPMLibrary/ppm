      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_comp_pp_verlet
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Subroutine which computes kernel interactions by
      !                 direct particle-particle interactions using Verlet
      !                 lists. 
      !
      !  Input        : xp(:,:)    (F) particle co-ordinates
      !                 Np         (I) numer of particles on local proc.
      !                 pdata(:,:) (O) particle data (strengths).
      !                                Overloaded types:
      !                                single,double,single complex,double
      !                                complex
      !                 lda        (I) leading dimension of pdata.
      !                 lsymm      (L) using symmetry or not?
      !                 kernel     (O) kernel to be used for PP
      !                                interactions. To use ppm-internal
      !                                kernels, specify one of:
      !                                   ppm_param_kerel_laplace2d_2p
      !                                     (2nd order Laplacian,
      !                                     polynomial in 2D)
      !                                   ppm_param_kerel_laplace3d_2p
      !                                     (2nd order Laplacian,
      !                                     polynomial in 3D)
      !                                To use your own kernel function,
      !                                pass the function pointer here. Your
      !                                function should take one argument
      !                                and return one value. The third
      !                                possibility is to pass a lookup
      !                                table with tabulated kernel values.
      !                                Such a table can be created using
      !                                ppm_comp_pp_mk_table.
      !                 kpar(:)    (O) Kernel parameters. See documentation
      !                                or ppm_comp_pp_kernels.inc for
      !                                description. Type can be single,
      !                                double, single complex or double
      !                                complex. When using a lookup table,
      !                                pass dxtableinv (the inverse of the
      !                                table spacing) as a scalar here.
      !                 nvlist(:)  (I) length of the Verlet lists of all
      !                                particles
      !                 vlist(:,:) (I) Verlet lists of all particles. 1st
      !                                index: particles to interact with
      !                                (1..nvlist), 2nd: particle (1..Np).
      !
      !  Input/output : 
      !
      !  Output       : dpd(:,:)   (O) Change of particle data (strengths).
      !                                Overloaded types:
      !                                single,double,single complex,double
      !                                complex.
      !                 info       (I) return status. =0 if no error.
      !
      !  Remarks      : The loops for lda.EQ.1 and lda.EQ.2 are explicitly
      !                 unfolded to allow vectorization in those cases. For
      !                 lda.GT.2 this routine will not vectorize.
      !
      !                 Both nvlist and vlist can be created using
      !                 ppm_neighlist_vlist. This should be done by the
      !                 user program since there could be several sets of
      !                 Verlet lists. This routine can then be called with
      !                 the appropriate one. Do not forget to DEALLOCATE
      !                 nvlist and vlist when they are not needed any more.
      !
      !                 dpd needs to be allocated to proper size before
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
      !  $Log: ppm_comp_pp_verlet.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/10/13 16:42:11  davidch
      !  added support for 3d sph kernels
      !
      !  Revision 1.3  2004/07/26 13:37:39  ivos
      !  bugfix: function pointer argument does not need to be declared
      !  EXTERNAL if there is an explicit interface.
      !
      !  Revision 1.2  2004/07/26 07:45:24  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.1  2004/07/23 12:57:05  ivos
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
      SUBROUTINE ppm_comp_pp_verlet_si(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_di(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_sci(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_dci(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_su(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_du(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_scu(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_dcu(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_st(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_verlet_dt(xp,Np,pdata,lda,lsymm,kernel,kpar,  &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_sct(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_verlet_dct(xp,Np,pdata,lda,lsymm,kernel,kpar, &
     &    nvlist,vlist,dpd,info)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
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
      INTEGER    , DIMENSION(:,:), INTENT(IN   ) :: vlist
      INTEGER    , DIMENSION(  :), INTENT(IN   ) :: nvlist
      INTEGER                    , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                                    :: jpart,ip,jp,ispec,idx
      REAL(MK)                                   :: dx,dy,dz
      REAL(MK)                                   :: factor,factor2
      REAL(MK)                                   :: dij,dij2,dij4,dij5
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                                   :: summ,summ2
      REAL(MK)                                   :: eta,dm
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)                                :: summ,summ2
      COMPLEX(MK)                                :: eta,dm
#endif
      ! timer 
      REAL(MK)                                   :: t0
      CHARACTER(LEN=ppm_char)                    :: mesg
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
      CALL substart('ppm_comp_pp_verlet',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_comp_pp_verlet',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_verlet',  &
     &            'Np must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_verlet',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __KERNEL == __LOOKUP_TABLE
          IF (kpar .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_verlet',  &
     &            'kpar (dxtableinv) must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS using symmetry
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
          IF (lda .EQ. 1) THEN
              DO ip=1,Np
                  summ = 0.0_MK
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp and vice versa
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm        = eta*(pdata(1,jp) - pdata(1,ip))
                      summ      = summ      + dm
                      dpd(1,jp) = dpd(1,jp) - dm
                  ENDDO
                  dpd(1,ip) = dpd(1,ip) + summ
              ENDDO
          ELSEIF (lda .EQ. 2) THEN
              DO ip=1,Np
                  summ = 0.0_MK
                  summ2 = 0.0_MK
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp and vice versa
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm        = eta*(pdata(1,jp) - pdata(1,ip))
                      summ      = summ      + dm
                      dpd(1,jp) = dpd(1,jp) - dm
                      dm        = eta*(pdata(2,jp) - pdata(2,ip))
                      summ2     = summ2     + dm
                      dpd(2,jp) = dpd(2,jp) - dm
                  ENDDO
                  dpd(1,ip) = dpd(1,ip) + summ
                  dpd(2,ip) = dpd(2,ip) + summ2
              ENDDO
          ELSE
              DO ip=1,Np
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp and vice versa
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      DO ispec=1,lda
                          dm   = eta*(pdata(ispec,jp) - pdata(ispec,ip))
                          dpd(ispec,ip) = dpd(ispec,ip) + dm
                          dpd(ispec,jp) = dpd(ispec,jp) - dm
                      ENDDO
                  ENDDO
              ENDDO
          ENDIF

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS not using symmetry
      !-------------------------------------------------------------------------
      ELSE
          IF (lda .EQ. 1) THEN
              DO ip=1,Np
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm = eta*(pdata(1,jp) - pdata(1,ip))
                      dpd(1,ip) = dpd(1,ip) + dm
                  ENDDO
              ENDDO
          ELSEIF (lda .EQ. 2) THEN
              DO ip=1,Np
#ifdef __SXF90
!CDIR NODEP
#endif
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      dm = eta*(pdata(1,jp) - pdata(1,ip))
                      dpd(1,ip) = dpd(1,ip) + dm
                      dm = eta*(pdata(2,jp) - pdata(2,ip))
                      dpd(2,ip) = dpd(2,ip) + dm
                  ENDDO
              ENDDO
          ELSE
              DO ip=1,Np
                  DO jpart=1,nvlist(ip)
                      jp = vlist(jpart,ip)
                      !---------------------------------------------------------
                      !  Calculate the square of the distance between the two 
                      !  particles. It will always be .LE. (cutoff+skin)**2 by 
                      !  construction of the Verlet list.
                      !---------------------------------------------------------
                      dx  = xp(1,ip) - xp(1,jp)
                      dy  = xp(2,ip) - xp(2,jp)
                      IF (ppm_dim .GT. 2) THEN
                          dz  = xp(3,ip) - xp(3,jp)
                          dij = (dx*dx) + (dy*dy) + (dz*dz)
                      ELSE
                          dz = 0.0_MK
                          dij = (dx*dx) + (dy*dy)
                      ENDIF
                      !---------------------------------------------------------
                      !  Particle ip interacts with jp
                      !---------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                      DO ispec=1,lda
                          dm = eta*(pdata(ispec,jp) - pdata(ispec,ip))
                          dpd(ispec,ip) = dpd(ispec,ip) + dm
                      ENDDO
                  ENDDO
              ENDDO
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_comp_pp_verlet',t0,info)
      RETURN
#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_si
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_di
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_sci
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_dci
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_su
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_du
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_scu
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_dcu
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_st
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_verlet_dt
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_sct
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_verlet_dct
#endif
#endif
