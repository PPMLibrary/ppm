      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_set_decomp_cost
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE set_dcost_s(dcost,method,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE set_dcost_d(dcost,method,info)
#endif
      !!! Sets/updates the internal estimate for the computational cost
      !!! associated with redecomposing the problem. The user can choose the
      !!! update method.
      !!!
      !!! [TIP]
      !!! The user should time (using `ppm_time`) topo_mktopo for the 
      !!! topology/topologies considered for dynamic remapping (ppm does not
      !!! know which ones they are). The elapsed time is given to this routine
      !!! to update the internal estimate.
      !!!
      !!! [NOTE]
      !!! This estimate is currently a scalar, so only one topology/set of
      !!! topologies can be monitored. Introducing a topology set ID  and
      !!! making the internal estimates a vector, this can later be extended to
      !!! multiple topology sets if needed. The topology set ID (in external
      !!! numbering) for which to update the cost estimate will then be an
      !!! additional argument to this routine.
      !!!
      !!! [NOTE]
      !!! This routine does a global MPI operation (`Allreduce`).
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
      !!! Elapsed time (as measured by `ppm_time`) 
      !!! for defining the topology/ies of interest on the local processor.
      INTEGER                , INTENT(IN   ) :: method
      !!! How to update the internal estimate.
      !!! One of:
      !!!
      !!! * ppm_param_update_replace                                         +
      !!!        (overwrite old value with new one)
      !!! * ppm_param_update_average                                           +
      !!!        (compute running average)
      !!! * ppm_param_update_expfavg                                           +
      !!!        (running average with  exponential forgetting)
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
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
        CALL check
        IF (info .NE. 0) GOTO 9999
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
      IF     (method .EQ. ppm_param_update_replace) THEN
          !---------------------------------------------------------------------
          !  Replace with new value
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
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_set_decomp_cost',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (dcost .LT. 0.0_MK) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_set_decomp_cost', &
     &            'dcost must be >= 0.0',__LINE__,info)
             GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE set_dcost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE set_dcost_d
#endif
