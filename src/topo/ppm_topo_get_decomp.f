      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_get_decomp
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
      SUBROUTINE ppm_topo_get_decomp_s(topoid,nsub,min_sub,max_sub,sub2proc,&
                                  subs_bc,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_get_decomp_d(topoid,nsub,min_sub,max_sub,sub2proc,&
                                  subs_bc,info)
#endif
      !!! This routine returns the subdomain boundaries and boundary
      !!! conditions

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_check_id
      USE ppm_module_alloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                  , INTENT(IN   ) :: topoid
      !!! Topology ID
      INTEGER                  , INTENT(  OUT) :: nsub
      !!! Number of subdomains
      REAL(MK), DIMENSION(:,:), POINTER        :: min_sub
      !!! Returns the min extent of all subdomains
      !!! index 1: dimension
      !!! index 2: sub ID
      REAL(MK), DIMENSION(:,:), POINTER        :: max_sub
      !!! Returns the max extent of all subdomains
      !!! index 1: dimension
      !!! index 2: sub ID
      INTEGER, DIMENSION(:  ), POINTER         :: sub2proc
      !!! subdomain to processor assignment
      INTEGER, DIMENSION(:,:), POINTER         :: subs_bc
      !!! Returns the boundary condition of each subdomain analogous to
      !!! ppm_t_topo definition
      !!!
      !!! - west  : 1
      !!! - east  : 2
      !!! - south : 3
      !!! - north : 4
      !!! - bottom: 5
      !!! - top   : 6
      !!!
      !!! index 1: the index of the 4 or 6 faces in 2 and 3 D
      !!! index 2: the global sub id
      !!!
      !!! states:
      !!!
      !!! - value: 0 the face is internal
      !!! - value: 1 otherwise
      INTEGER                  , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)     :: t0
      LOGICAL                   :: valid
      TYPE(ppm_t_topo), POINTER :: topo => NULL()
      INTEGER, DIMENSION(2)     :: lda
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_get_decomp',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the topology
      !-------------------------------------------------------------------------

      topo => ppm_topo(topoid)%t

      lda(1) = ppm_dim
      lda(2) = topo%nsubs
      call ppm_alloc(min_sub,lda,ppm_param_alloc_fit,info)
      call ppm_alloc(max_sub,lda,ppm_param_alloc_fit,info)
      if (info .NE. 0) then
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_get_decomp',     &
     &        'failed to allocate min_sub, max_sub',__LINE__,info)
          goto 9999
      endif
      lda(1) = ppm_dim*2
      call ppm_alloc(subs_bc,lda,ppm_param_alloc_fit,info)
      if (info .NE. 0) then
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_get_decomp',     &
     &        'failed to allocate subs_bc',__LINE__,info)
          goto 9999
      endif
      lda(1) = topo%nsubs
      call ppm_alloc(sub2proc,lda,ppm_param_alloc_fit,info)
      if (info .NE. 0) then
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_get_decomp',     &
     &        'failed to allocate sub2proc',__LINE__,info)
          goto 9999
      endif

      nsub = topo%nsubs
      IF (topo%prec.EQ.ppm_kind_single) THEN
        min_sub    = topo%min_subs
        max_sub    = topo%max_subs
      ELSE
        min_sub    = topo%min_subd
        max_sub    = topo%max_subd
      ENDIF
      subs_bc     = topo%subs_bc
      sub2proc    = topo%sub2proc


      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_get_decomp',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_get_decomp',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_get_decomp',  &
     &            'Topology ID is invalid!',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_get_decomp_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_get_decomp_d
#endif
