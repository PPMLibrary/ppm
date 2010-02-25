      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_set_decomp_cost
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Sets/updates the internal estimate for the
      !                 computational cost associated with redecomposing
      !                 the problem. The user can choose the update method.
      !
      !  Input        : dcost      (F) elapsed time (as measured by
      !                                ppm_time) for defining the
      !                                topology/ies of interest on the
      !                                local processor.
      !                 method     (I) How to update the internal estimate.
      !                                One of:
      !                                     ppm_param_update_overwrite
      !                                       (overwrite old value with
      !                                       new one)
      !                                     ppm_param_update_average
      !                                       (compute running average)
      !                                     ppm_param_update_expfavg
      !                                       (running average with
      !                                       exponential forgetting)
      !
      !  Input/output : 
      !
      !  Output       : info       (I) 0 on success.
      !
      !  Remarks      : The user should time (using ppm_time) topo_mktopo
      !                 for the topology/topologies considered for dynamic
      !                 remaping (ppm does not know which ones they are).
      !                 The elapsed time is given to this
      !                 routine to update the internal estimate. This
      !                 estimate is currently a scalar, so only one
      !                 topology/set of topologies can be monitored.
      !                 Introducing a topology set ID (and corresponding
      !                 translation lists between internal and external
      !                 numbering of these sets) and making the internal
      !                 estimates a vector, this can later be extended to
      !                 multiple topology sets if needed. The topology set
      !                 ID (in external numbering) for which to update the 
      !                 cost estimate will then be an additional argument
      !                 to this routine.
      !
      !                 This routine does a global MPI operation
      !                 (Allreduce).
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_set_decomp_cost.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2006/07/04 16:01:56  ivos
      !  Corrected typo in comment.
      !
      !  Revision 1.3  2004/07/26 07:45:25  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.2  2004/07/21 09:09:36  ivos
      !  Updated comment headers.
      !
      !  Revision 1.1  2004/07/20 16:48:26  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_set_decomp_cost_s(dcost,method,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_set_decomp_cost_d(dcost,method,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_loadbal
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
      REAL(MK)               , INTENT(IN   ) :: dcost
      INTEGER                , INTENT(IN   ) :: method
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: t0,maxcost
      REAL(ppm_kind_double)  :: alpha,beta
      CHARACTER(LEN=ppm_char):: mesg
#ifdef __MPI
      INTEGER                :: MPTYPE
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_set_decomp_cost',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_set_decomp_cost',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (dcost .LT. 0.0_MK) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_set_decomp_cost', &
     &            'dcost must be >= 0.0',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Reduce contributions from all processors. This is needed to make
      !  sure the redecomposition heuristic in ppm_loadbal_inquire returns
      !  the same result on all processors. Otherwise, some would try to
      !  redo the topology, but others not!
      !-------------------------------------------------------------------------
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Define MPI data type
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION
      MPTYPE = MPI_DOUBLE_PRECISION
#endif
      !-------------------------------------------------------------------------
      !  All have to wait for the slowest one to create a new topology
      !-------------------------------------------------------------------------
      CALL MPI_Allreduce(dcost,maxcost,1,MPTYPE,MPI_MAX,ppm_comm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_mpi_fail,'ppm_set_decomp_cost', &
     &         'MPI_ALLREDUCE for max dcost. Nothing updated.',__LINE__,info)
          GOTO 9999
      ENDIF
#else
      maxcost = dcost
#endif

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,F12.6)') 'Old decomp cost estimate: ',   &
     &        ppm_loadbal_decomp_cost
          CALL ppm_write(ppm_rank,'ppm_set_decomp_cost',mesg,info)
          WRITE(mesg,'(A,F12.6)') 'Adding measured value: ',maxcost
          CALL ppm_write(ppm_rank,'ppm_set_decomp_cost',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Choose update method
      !-------------------------------------------------------------------------
      IF     (method .EQ. ppm_param_update_overwrite) THEN
          !---------------------------------------------------------------------
          !  Overwrite with new value
          !---------------------------------------------------------------------
          ppm_loadbal_dcn = 1
#if   __KIND == __SINGLE_PRECISION
          ppm_loadbal_decomp_cost = REAL(maxcost,ppm_kind_double)
#elif __KIND == __DOUBLE_PRECISION
          ppm_loadbal_decomp_cost = maxcost
#endif
          IF (ppm_debug .GT. 1) THEN
              CALL ppm_write(ppm_rank,'ppm_set_decomp_cost',  &
     &            'Overwriting old value.',info)
          ENDIF
      ELSEIF (method .EQ. ppm_param_update_average) THEN
          !---------------------------------------------------------------------
          !  Running average
          !---------------------------------------------------------------------
          ppm_loadbal_dcn = ppm_loadbal_dcn + 1
#if   __KIND == __SINGLE_PRECISION
          alpha = REAL(maxcost,ppm_kind_double) - ppm_loadbal_decomp_cost
#elif __KIND == __DOUBLE_PRECISION
          alpha = maxcost - ppm_loadbal_decomp_cost
#endif
          alpha = alpha/REAL(ppm_loadbal_dcn,ppm_kind_double)
          ppm_loadbal_decomp_cost = ppm_loadbal_decomp_cost + alpha
          IF (ppm_debug .GT. 1) THEN
              CALL ppm_write(ppm_rank,'ppm_set_decomp_cost',  &
     &            'Updating running average.',info)
          ENDIF
      ELSEIF (method .EQ. ppm_param_update_expfavg) THEN
          !---------------------------------------------------------------------
          !  Running average with exponential forgetting. Alpha = 0.5 is a
          !  good value, but can be changed. Hardcoded for simplicity.
          !---------------------------------------------------------------------
          alpha = 0.5_ppm_kind_double
          beta  = 1.0_ppm_kind_double - alpha
          ppm_loadbal_dcn = ppm_loadbal_dcn + 1
#if   __KIND == __SINGLE_PRECISION
          IF (ppm_loadbal_dcn .EQ. 1) THEN
              ppm_loadbal_decomp_cost = REAL(maxcost,ppm_kind_double)
          ELSE
              ppm_loadbal_decomp_cost = (beta*ppm_loadbal_decomp_cost) +   &
     &            (alpha*REAL(maxcost,ppm_kind_double))
          ENDIF
#elif __KIND == __DOUBLE_PRECISION
          IF (ppm_loadbal_dcn .EQ. 1) THEN
              ppm_loadbal_decomp_cost = maxcost
          ELSE
              ppm_loadbal_decomp_cost = (beta*ppm_loadbal_decomp_cost) +   &
     &            (alpha*maxcost)
          ENDIF
#endif
          IF (ppm_debug .GT. 1) THEN
              CALL ppm_write(ppm_rank,'ppm_set_decomp_cost',  &
     &            'Updating running average with exponential forgetting.',info)
          ENDIF
      ELSE
          !---------------------------------------------------------------------
          !  Unknow update method
          !---------------------------------------------------------------------
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_set_decomp_cost',    &
     &       'Unknown update method specified. Nothing updated.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'Number of samples used: ',ppm_loadbal_dcn
          CALL ppm_write(ppm_rank,'ppm_set_decomp_cost',mesg,info)
          WRITE(mesg,'(A,F12.6)') 'New decomp cost estimate: ',   &
     &        ppm_loadbal_decomp_cost
          CALL ppm_write(ppm_rank,'ppm_set_decomp_cost',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_set_decomp_cost',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_set_decomp_cost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_set_decomp_cost_d
#endif
