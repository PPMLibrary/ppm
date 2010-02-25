      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_comp_pp_doring
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Computes direct particle-particle interactions using the
      !                 ring topology (N**2). This is called by
      !                 ppm_comp_pp_ring and not by the user directly.
      !
      !  Input        : xp(:,:)      (F) particle co-ordinates [Group 1]
      !                 pdata(:,:)   (O) particle data used for interaction
      !                                  (e.g. vorticity, strength, ...)
      !                                  [Group 1]. Overloaded types:
      !                                  single,double,single
      !                                  complex,double complex.
      !                 Mpart        (I) number of particles [Group 1]
      !                 xp2(:,:)     (F) particle co-ordinates [Group 2]
      !                 pdata2(:,:)  (O) particle data used for interaction
      !                                  (e.g. vorticity, strength, ...)
      !                                  [Group 2]. Overloaded types:
      !                                  single,double,single
      !                                  complex,double complex.
      !                 Lpart        (I) number of particles [Group 2]
      !                 lda          (I) leading dimension of pdata,pdata2.
      !                 lsymm        (L) Whether to use symmetry or not:
      !                                   .FALSE. pp interaction w/o symmetry
      !                                   .TRUE.  pp interaction w/  symmetry
      !                 kernel       (O) kernel to be used for PP
      !                                  interactions. To use ppm-internal
      !                                  kernels, specify one of:
      !                                     ppm_param_kerel_laplace2d_2p
      !                                       (2nd order Laplacian,
      !                                       polynomial in 2D)
      !                                     ppm_param_kerel_laplace3d_2p
      !                                       (2nd order Laplacian,
      !                                       polynomial in 3D)
      !                                  To use your own kernel function,
      !                                  pass the function pointer here. Your
      !                                  function should take one argument
      !                                  and return one value. The third
      !                                  possibility is to pass a lookup
      !                                  table with tabulated kernel values.
      !                                  Such a table can be created using
      !                                  ppm_comp_pp_mk_table.
      !                 kpar(:)      (O) Kernel parameters. See documentation
      !                                  or ppm_comp_pp_kernels.inc for
      !                                  description. Type can be single,
      !                                  double, single complex or double
      !                                  complex. When using a lookup table,
      !                                  pass dxtableinv (the inverse of the
      !                                  table spacing) as a scalar here.
      !                 mode         (I) Whether the two groups are the
      !                                  same or not:
      !                                   0 not the same group
      !                                   1 the same group
      !
      !  Input/output : dpd(:,:)     (O) Change of particle data due
      !                                  to interaction [Group 1]
      !                                  Overloaded types:
      !                                  single,double,single
      !                                  complex,double complex.
      !                 dpd2(:,:)    (O) Change of particle data due
      !                                  to interaction [Group 2]
      !                                  Overloaded types:
      !                                  single,double,single
      !                                  complex,double complex.
      !  Output       : info         (I) return status. 0 if no error.
      !
      !  Remarks      : The loops over lda are explicitly unrolled for
      !                 the cases lda.EQ.1 and lda.EQ.2 in order to allow
      !                 vectorization in these cases.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_comp_pp_doring.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/10/14 11:46:43  davidch
      !  added support for 3d sph kernels
      !
      !  Revision 1.3  2004/07/26 13:37:40  ivos
      !  bugfix: function pointer argument does not need to be declared
      !  EXTERNAL if there is an explicit interface.
      !
      !  Revision 1.2  2004/07/26 07:45:24  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.1  2004/07/23 12:57:07  ivos
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
      SUBROUTINE ppm_comp_pp_doring_si(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_doring_di(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_doring_sci(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_doring_dci(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#endif

#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_doring_su(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_doring_du(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_doring_scu(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_doring_dcu(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_comp_pp_doring_st(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_comp_pp_doring_dt(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_doring_sct(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_comp_pp_doring_dct(xp,pdata,dpd,Mpart,xp2,pdata2,dpd2, &
     &    Lpart,lda,lsymm,kernel,kpar,mode,info)
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
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: xp,xp2
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:), INTENT(IN   ) :: pdata,pdata2
      REAL(MK)   , DIMENSION(:,:), INTENT(INOUT) :: dpd,dpd2
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:), INTENT(IN   ) :: pdata,pdata2
      COMPLEX(MK), DIMENSION(:,:), INTENT(INOUT) :: dpd,dpd2
#endif
      INTEGER                    , INTENT(IN   ) :: Mpart,Lpart,mode,lda
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
      LOGICAL                    , INTENT(IN   ) :: lsymm
      INTEGER                    , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      ! counters
      INTEGER                                    :: i,j,ispec,idx
      ! coordinate differences
      REAL(MK)                                   :: dx,dy,dz
      REAL(MK)                                   :: factor,factor2
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                                   :: summ,summ2
      REAL(MK)                                   :: eta,dm
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK)                                :: summ,summ2
      COMPLEX(MK)                                :: eta,dm
#endif
      ! square of inter particle distance
      REAL(MK)                                   :: dij,dij2,dij4,dij5
      REAL(MK)                                   :: t0
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
      CALL substart('ppm_comp_pp_doring',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Mpart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_doring',  &
     &            'Np must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Lpart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_doring',  &
     &            'Np must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_doring',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((mode .NE. 0) .AND. (mode .NE. 1)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_doring',  &
     &            'mode must be either 0 or 1',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __KERNEL == __LOOKUP_TABLE
          IF (kpar .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_comp_pp_doring',  &
     &            'kpar (dxtableinv) must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Check if we are computing the selfinteraction (i.e. Group1 .EQ. Group2)
      !-------------------------------------------------------------------------
      IF (mode .EQ. 1) THEN
          IF (lsymm) THEN
              !-----------------------------------------------------------------
              !  SYMMETRY
              !-----------------------------------------------------------------
              IF (lda .EQ. 1) THEN
                 DO i = 1,Mpart
                    summ = 0.0_MK
                    DO j = i+1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j and vice versa
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       summ = summ + dm
                       dpd2(1,j) = dpd2(1,j) - dm
                    ENDDO
                    dpd(1,i) = dpd(1,i) + summ
                 ENDDO
              ELSEIF (lda .EQ. 2) THEN
                 DO i = 1,Mpart
                    summ = 0.0_MK
                    summ2 = 0.0_MK
                    DO j = i+1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j and vice versa
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       summ = summ + dm
                       dpd2(1,j) = dpd2(1,j) - dm
                       dm   = eta*(pdata2(2,j) - pdata(2,i))
                       summ2 = summ2 + dm
                       dpd2(2,j) = dpd2(2,j) - dm
                    ENDDO
                    dpd(1,i) = dpd(1,i) + summ
                    dpd(2,i) = dpd(2,i) + summ2
                 ENDDO
              ELSE
                 DO i = 1,Mpart
                    DO j = i+1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j and vice versa
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       DO ispec=1,lda
                          dm   = eta*(pdata2(ispec,j) - pdata(ispec,i))
                          dpd(ispec,i) = dpd(ispec,i) + dm
                          dpd2(ispec,j) = dpd2(ispec,j) - dm
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
          ELSE
              !-----------------------------------------------------------------
              !  NO SYMMETRY
              !-----------------------------------------------------------------
              IF (lda .EQ. 1) THEN
                 DO i = 1,Mpart
                    DO j = 1,Lpart
                       IF (i .EQ. j) CYCLE
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j 
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       dpd(1,i) = dpd(1,i) + dm
                    ENDDO
                 ENDDO
              ELSEIF (lda .EQ. 2) THEN
                 DO i = 1,Mpart
                    DO j = 1,Lpart
                       IF (i .EQ. j) CYCLE
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j 
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       dpd(1,i) = dpd(1,i) + dm
                       dm   = eta*(pdata2(2,j) - pdata(2,i))
                       dpd(2,i) = dpd(2,i) + dm
                    ENDDO
                 ENDDO
              ELSE
                 DO i = 1,Mpart
                    DO j = 1,Lpart
                       IF (i .EQ. j) CYCLE
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j 
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       DO ispec=1,lda
                          dm   = eta*(pdata2(ispec,j) - pdata(ispec,i))
                          dpd(ispec,i) = dpd(ispec,i) + dm
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
          ENDIF

      !-------------------------------------------------------------------------
      !  non-self interaction (i.e. Group1 .NE. Group2)
      !-------------------------------------------------------------------------
      ELSE
          IF (lsymm) THEN
              !-----------------------------------------------------------------
              !  SYMMETRY
              !-----------------------------------------------------------------
              IF (lda .EQ. 1) THEN
                 DO i = 1,Mpart
                    summ = 0.0_MK
                    DO j = 1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j and vice versa
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       summ = summ + dm
                       dpd2(1,j) = dpd2(1,j) - dm
                    ENDDO
                    dpd(1,i) = dpd(1,i) + summ
                 ENDDO
              ELSEIF (lda .EQ. 2) THEN
                 DO i = 1,Mpart
                    summ = 0.0_MK
                    summ2 = 0.0_MK
                    DO j = 1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j and vice versa
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       summ = summ + dm
                       dpd2(1,j) = dpd2(1,j) - dm
                       dm   = eta*(pdata2(2,j) - pdata(2,i))
                       summ2 = summ2 + dm
                       dpd2(2,j) = dpd2(2,j) - dm
                    ENDDO
                    dpd(1,i) = dpd(1,i) + summ
                    dpd(2,i) = dpd(2,i) + summ2
                 ENDDO
              ELSE
                 DO i = 1,Mpart
                    DO j = 1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j and vice versa
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       DO ispec=1,lda
                          dm   = eta*(pdata2(ispec,j) - pdata(ispec,i))
                          dpd(ispec,i) = dpd(ispec,i) + dm
                          dpd2(ispec,j) = dpd2(ispec,j) - dm
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
          ELSE
              !-----------------------------------------------------------------
              !  NO SYMMETRY
              !-----------------------------------------------------------------
              IF (lda .EQ. 1) THEN
                 DO i = 1,Mpart
                    DO j = 1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       dpd(1,i) = dpd(1,i) + dm
                    ENDDO
                 ENDDO
              ELSEIF (lda .EQ. 2) THEN
                 DO i = 1,Mpart
                    DO j = 1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j 
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       dm   = eta*(pdata2(1,j) - pdata(1,i))
                       dpd(1,i) = dpd(1,i) + dm
                       dm   = eta*(pdata2(2,j) - pdata(2,i))
                       dpd(2,i) = dpd(2,i) + dm
                    ENDDO
                 ENDDO
              ELSE
                 DO i = 1,Mpart
                    DO j = 1,Lpart
                       dx  = xp(1,i) - xp2(1,j)
                       dy  = xp(2,i) - xp2(2,j)
                       IF (ppm_dim .GT. 2) THEN
                           dz  = xp(3,i) - xp2(3,j)
                           dij = (dx*dx) + (dy*dy) + (dz*dz)
                       ELSE
                           dz = 0.0_MK
                           dij = (dx*dx) + (dy*dy)
                       ENDIF
                       !--------------------------------------------------------
                       !  Particle i interacts with j
                       !--------------------------------------------------------
#include "ppm_comp_pp_kernels.inc"
                       DO ispec=1,lda
                          dm   = eta*(pdata2(ispec,j) - pdata(ispec,i))
                          dpd(ispec,i) = dpd(ispec,i) + dm
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_comp_pp_doring',t0,info)
      RETURN
#if   __KERNEL == __INTERNAL
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_doring_si
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_doring_di
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_doring_sci
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_doring_dci
#endif
      
#elif __KERNEL == __USER_FUNCTION
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_doring_su
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_doring_du
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_doring_scu
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_doring_dcu
#endif

#elif __KERNEL == __LOOKUP_TABLE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_comp_pp_doring_st
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_comp_pp_doring_dt
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_doring_sct
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_comp_pp_doring_dct
#endif
#endif
