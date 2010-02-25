      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_tree_boxcost
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine computes the cost of a set of boxes.
      !
      !  Input        : Nm(:,:)      (I) Number of grid points in each box.
      !                 weights(3)   (F) weights for the three cost
      !                                  contributions: particles, mesh,
      !                                  geometry for defining the total
      !                                  cost of a box.
      !                 min_box(:,:) (F) the minimum coordinate of the
      !                                  boxes. 2nd: box id
      !                 max_box(:,:) (F) the maximum coordinate of the
      !                                  boxes. 2nd: box id.
      !                 nbox         (I) number of boxes to compute cost
      !                                  for.
      !                 lhbx(:)      (I) pointer to first particle in
      !                                  each box to determine cost for.
      !                                  This only needs to be present
      !                                  if there are any particles.
      !                 lpdx(:)      (I) index list of points in
      !                                  each box to determine cost for.
      !                                  This only needs to be present
      !                                  if there are any particles.
      !                 pcost(:)     (F) OPTIONAL argument of length
      !                                  Npart, specifying the
      !                                  computational cost of each
      !                                  particle.
      !
      !  Input/output : boxcost(:)   (F) costs of the boxes 1..nbox.
      !                                  boxcost(i) is the cost of box 
      !                                  boxes(i).
      !
      !  Output       : info         (I) return status
      !
      !  Remarks      : The cost of a particle is counted as pcost (if
      !                 given) or 1. For meshes, the cost is 1 per mesh
      !                 point. For the geometry part, the cost of a
      !                 box is given by its volume. Use weights(3) to
      !                 adjust this if needed.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree_boxcost.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.11  2006/09/04 18:34:57  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.10  2006/04/19 07:21:59  ivos
      !  Cross-check. Thanks Michael for resolving bug 0000021.
      !
      !  Revision 1.9  2006/04/19 06:41:30  michaebe
      !  resolved bug 0000021
      !
      !  Revision 1.8  2005/08/31 13:33:50  ivos
      !  bugfix: removed doubly-declared variables and unused arguments.
      !
      !  Revision 1.7  2005/08/31 12:43:44  ivos
      !  Shark optimizations.
      !
      !  Revision 1.6  2005/08/31 11:24:30  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.5  2005/08/30 13:17:26  ivos
      !  Sharked the routines and unrolled all loops over ppm_dim.
      !
      !  Revision 1.4  2004/12/03 17:14:59  ivos
      !  Switched to the use of particle lists lhbx and lpdx.
      !
      !  Revision 1.3  2004/12/02 16:31:54  ivos
      !  Optimized some loops.
      !
      !  Revision 1.2  2004/09/22 17:25:57  ivos
      !  bugfix: added ppm_module_write.
      !
      !  Revision 1.1  2004/09/22 10:32:02  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_boxcost_s(Nm,weights,min_box,max_box,   &
     &    nbox,lhbx,lpdx,boxcost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_boxcost_d(Nm,weights,min_box,max_box,   &
     &    nbox,lhbx,lpdx,boxcost,info,pcost)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_tree
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_box,max_box
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: weights
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      INTEGER                 , INTENT(IN   ) :: nbox
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: lhbx,lpdx
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: Nm
      REAL(MK), DIMENSION(:  ), POINTER       :: boxcost
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: len_box
      INTEGER , DIMENSION(2)                  :: ldc
      REAL(MK)                                :: t0,meshtotal,geomtotal,dm
      INTEGER                                 :: i,ip,iopt,j
      REAL(MK), DIMENSION(:), POINTER         :: pcst
#ifdef __MPI
      REAL(MK), DIMENSION(:), POINTER         :: pcsum
      INTEGER                                 :: MPTYPE
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_boxcost',t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         DO ip=1,nbox
            DO i=1,ppm_dim
               IF (min_box(i,ip) .GT. max_box(i,ip)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_tree_boxcost',     &
     &                'min_box must be <= max_box !',__LINE__,info)
                  GOTO 9999
               ENDIF 
            ENDDO
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  Set pointer to work memory
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      pcst => pcst_s
#else
      pcst => pcst_d
#endif
#ifdef   __MPI
#if   __KIND == __SINGLE_PRECISION
      pcsum => pcsum_s
#else
      pcsum => pcsum_d
#endif
#endif

      !-------------------------------------------------------------------------
      !  If we have less than 1 box to compute, exit
      !-------------------------------------------------------------------------
      IF (nbox .LT. 1) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_boxcost',   &
     &            'No boxes to be computed. Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Determine MPI data type
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION 
      MPTYPE = MPI_DOUBLE_PRECISION
#endif
#endif

      !-------------------------------------------------------------------------
      !  Computation of particle contributions to box costs
      !-------------------------------------------------------------------------
      IF (have_particles .AND. weights(1) .NE. 0.0_MK) THEN
          DO i=1,nbox
              pcst(i)  = 0.0_MK
#ifdef __MPI
              pcsum(i) = 0.0_MK
#endif
              DO j=lhbx(i),lhbx(i+1)-1
                  ip = lpdx(j)
                  dm = 1.0_MK
                  IF (PRESENT(pcost)) dm = pcost(ip)
                  pcst(i) = pcst(i) + dm
              ENDDO
          ENDDO

#ifdef __MPI
          !---------------------------------------------------------------------
          !  Allreduce of all particle sums
          !---------------------------------------------------------------------
          CALL MPI_Allreduce(pcst,pcsum,nbox,MPTYPE,MPI_SUM,ppm_comm,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_tree_boxcost',   &
     &            'MPI_Allreduce of particle costs',__LINE__,info)
              GOTO 9999
          ENDIF 
          DO i=1,nbox
              pcst(i) = pcsum(i)
          ENDDO
#endif
      ENDIF   ! have_particles

      !-------------------------------------------------------------------------
      !  Add up all cost contributions to total cost.
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          DO i=1,nbox
              !-----------------------------------------------------------------
              !  Mesh and geometry based costs
              !-----------------------------------------------------------------
              meshtotal = 0.0_MK
              len_box(1) = max_box(1,i)-min_box(1,i)
              len_box(2) = max_box(2,i)-min_box(2,i)
              IF (have_mesh .AND. weights(2) .NE. 0) THEN
                  meshtotal = REAL(Nm(1,i),MK)*REAL(Nm(2,i),MK)
              ENDIF
              geomtotal = len_box(1)*len_box(2)
              !-----------------------------------------------------------------
              !  Compute weighted sum of costs and store it
              !-----------------------------------------------------------------
              IF(have_particles) THEN
                 boxcost(i) = pcst(i)*weights(1) + meshtotal*weights(2) +    &
     &                            geomtotal*weights(3)
              ELSE
                 boxcost(i) = meshtotal*weights(2) + geomtotal*weights(3)
              END IF
          ENDDO
      ELSE
          DO i=1,nbox
              !-----------------------------------------------------------------
              !  Mesh and geometry based costs
              !-----------------------------------------------------------------
              meshtotal = 0.0_MK
              len_box(1) = max_box(1,i)-min_box(1,i)
              len_box(2) = max_box(2,i)-min_box(2,i)
              len_box(3) = max_box(3,i)-min_box(3,i)
              IF (have_mesh .AND. weights(2) .NE. 0) THEN
                  meshtotal = REAL(Nm(1,i),MK)*REAL(Nm(2,i),MK)*REAL(Nm(3,i),MK)
              ENDIF
              geomtotal = len_box(1)*len_box(2)*len_box(3)
              !-----------------------------------------------------------------
              !  Compute weighted sum of costs and store it
              !-----------------------------------------------------------------
              IF(have_particles) THEN
                 boxcost(i) = pcst(i)*weights(1) + meshtotal*weights(2) +    &
     &                            geomtotal*weights(3)
              ELSE
                 boxcost(i) = meshtotal*weights(2) + geomtotal*weights(3)
              END IF
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Free local memory
      !-------------------------------------------------------------------------
 9999 CONTINUE
      NULLIFY(pcst)
#ifdef   __MPI
      NULLIFY(pcsum)
#endif
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
      CALL substop('ppm_tree_boxcost',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_boxcost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_boxcost_d
#endif
