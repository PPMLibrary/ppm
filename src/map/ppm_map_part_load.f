      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_part_load
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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

      SUBROUTINE ppm_map_part_load(info)
      !!! This routine loads the internally stored particle mapping.
      !!!
      !!! See `ppm_map_part_store` for more information.

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_state
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  This routine is ALL integer, except the timing variable t0
      !  We pick double as the precision - no real reason for that 
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER    :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(OUT)  :: info
      !!! Returns stutus, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: ldu
      INTEGER               :: i,j,k
      INTEGER               :: iopt
      REAL(MK)              :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_load',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_load',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  load the map type (is used in _send() line 229)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_map_type_state 

      !-------------------------------------------------------------------------
      !  load the ppm_nsendlist: the number of processers to send to
      !-------------------------------------------------------------------------
      ppm_nsendlist = ppm_nsendlist_state

      !-------------------------------------------------------------------------
      !  load the ppm_nsendbuffer: the size of the send buffer - here zero
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ppm_nsendbuffer_state 

      !-------------------------------------------------------------------------
      !  load the ppm_buffer_set: the current number of sets stored - here zero
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set_state 

      !-------------------------------------------------------------------------
      !  load the ppm_psendbuffer:
      !  first allocate memory for it
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_load', &
     &        'allocation of psendbuffer failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  then load it: pointer to the first particle in the send buffer
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist+1
         ppm_psendbuffer(k) = ppm_psendbuffer_state(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  load the ppm_buffer2part
      !  first allocate memory for it
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = SIZE(ppm_buffer2part_state)
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_load', &
     &        'allocation of buffer2part failed',__LINE__,info)
          GOTO 9999
      ENDIF
 
      !-------------------------------------------------------------------------
      !  then load it: pointer to the particle id of the j-th entry in the 
      !  send buffer
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist
         DO j=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
            ppm_buffer2part(j) = ppm_buffer2part_state(j)
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Store the ppm_nrecvlist and ppm_sendlist 
      !-------------------------------------------------------------------------
      ppm_nrecvlist = ppm_nrecvlist_state 

      ppm_nsendlist = ppm_nsendlist_state 

      !-------------------------------------------------------------------------
      !  load the receive list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      DO k=1,ppm_nrecvlist
         ppm_irecvlist(k) = ppm_irecvlist_state(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  load the send list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      DO k=1,ppm_nsendlist
         ppm_isendlist(k) = ppm_isendlist_state(k) 
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_load',t0,info)
      RETURN
      END SUBROUTINE ppm_map_part_load
