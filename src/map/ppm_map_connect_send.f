      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_connect_send
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

      SUBROUTINE ppm_map_connect_send(cd,lda,Ncon,id,Npart,info)
      !!! This routine sends and receives the connection data
      !!! according to the information stored in `ppm_buffer2part`
      !!!
      !!! [NOTE]
      !!! This routine directly sends and receives the connection data as
      !!! integers

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mpi
      USE ppm_module_util_invert_list
      USE ppm_module_mapping_typedef, ONLY : ppm_buffer_set,ppm_buffer2part, &
      &   ppm_nsendlist,ppm_psendbuffer,ppm_isendlist,ppm_irecvlist
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), POINTER       :: cd
      !!! Connection data
      INTEGER                , INTENT(IN   ) :: lda
      !!! The leading dimension of cd(:,:)
      INTEGER                , INTENT(INOUT) :: Ncon
      !!! Number of connections
      INTEGER, DIMENSION(:)  , TARGET        :: id
      !!! Local to global particle number mapping
      INTEGER                , INTENT(IN   ) :: Npart
      !!! Number of particles
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,l,msend,mrecv,iopt,ipos,part,qpart,icon
      INTEGER               :: tag,cd_size,ncons_local
      INTEGER               :: hops
      INTEGER, PARAMETER    :: big=HUGE(1)
      INTEGER               :: min_,max_

      CHARACTER(LEN=ppm_char) :: caller='ppm_map_connect_send'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if the send buffer is empty
      !-------------------------------------------------------------------------
      IF (ppm_buffer_set.LT.1) THEN
         fail('No mapping defined!', &
         & ppm_err_buffer_empt,ppm_error=ppm_error_warning)
      ENDIF

#ifndef __MPI
      !-------------------------------------------------------------------------
      !  If we only have one processor there is not much to do
      !-------------------------------------------------------------------------
      GOTO 9999
#else
      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(csend,ldu,iopt,info)
      or_fail_alloc('allocating array CSEND',ppm_error=ppm_error_fatal)

      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(crecv,ldu,iopt,info)
      or_fail_alloc('allocating array CRECV',ppm_error=ppm_error_fatal)

      ldu(1) = ppm_nsendlist
      ldu(2) = Ncon
      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_alloc('allocating array PSEND',ppm_error=ppm_error_fatal)

      ldu(1) = Npart
      CALL ppm_alloc(id_send,ldu,iopt,info)
      or_fail_alloc('allocating list ID_SEND',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Initialize the arrays and counter
      !-------------------------------------------------------------------------
      id_send(:) = 0
      csend(:)   = 0
      crecv(:)   = 0
      psend(:,:) = 0

      msend = 0
      mrecv = 0

      !-------------------------------------------------------------------------
      !  Calculate how many particles will stay on this processor
      !-------------------------------------------------------------------------
      ipos = ppm_psendbuffer(2) - ppm_psendbuffer(1)

      !-------------------------------------------------------------------------
      !  Count how many connections we have to send to the other processors
      !-------------------------------------------------------------------------
      DO k = 2,ppm_nsendlist
         qpart = ppm_psendbuffer(k+1) - ppm_psendbuffer(k)
         !----------------------------------------------------------------------
         !  Loop over all particles that go the k-th processor
         !----------------------------------------------------------------------
         DO i = ipos+1,(ipos+qpart)
            part = id(ppm_buffer2part(i))
            id_send(ppm_buffer2part(i)) = part
            !-------------------------------------------------------------------
            !  Loop over all connections to find the particle part
            !-------------------------------------------------------------------
            DO j = 1,Ncon
               DO l = 1,lda
                  IF (part .EQ. cd(l,j)) THEN
                      !---------------------------------------------------------
                      !  Count and mark this particle if it is not already
                      !  counted and marked (this particles goes to the k-th
                      !  processor)
                      !---------------------------------------------------------
                      IF (psend(k,j) .EQ. 0) THEN
                          csend(k) = csend(k) + 1
                          psend(k,j) = 1
                          !----------------------------------------------------
                          !  Mark all particles that will be send away. This
                          !  is needed to do a faster check for orphan
                          !  connections
                          !----------------------------------------------------
                          psend(1,j) = 1
                      ENDIF
                      EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !  Communicate the number of connections we have to send and to
         !  receive
         !----------------------------------------------------------------------
         tag = 100
         CALL MPI_SendRecv(csend(k),1,MPI_INTEGER,ppm_isendlist(k),tag, &
         &                 crecv(k),1,MPI_INTEGER,ppm_irecvlist(k),tag, &
         &                 ppm_comm,MPI_STATUS_IGNORE,info)
         or_fail_MPI("MPI_SendRecv")

         msend = msend + csend(k)
         mrecv = mrecv + crecv(k)

         ipos = ipos + qpart
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate the connection data array in order to hold the incoming
      !  connections
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = lda
      ldu(2) = Ncon + mrecv
      CALL ppm_alloc(cd_local,ldu,iopt,info)
      or_fail_alloc('allocating connectin array CD_LOCAL',ppm_error=ppm_error_fatal)

      cd_local(:,:) = 0
      cd_local(:,1:Ncon) = cd(:,:)

      !-------------------------------------------------------------------------
      !  Allocate memory for the send buffer and receive buffer
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = MAXVAL(csend,1)*lda
      CALL ppm_alloc(sendbuffer,ldu,iopt,info)
      or_fail_alloc('allocating send buffer SENDBUFFER',ppm_error=ppm_error_fatal)

      ldu(1) = MAXVAL(crecv,1)*lda
      CALL ppm_alloc(recvbuffer,ldu,iopt,info)
      or_fail_alloc('allocating receive buffer RECVBUFFER',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Send and receive the connections
      !-------------------------------------------------------------------------
      icon = 1
      DO k = 2,ppm_nsendlist
         !----------------------------------------------------------------------
         !  Fill the sendbuffer for the k-th processor
         !----------------------------------------------------------------------
         ipos = 1
         DO j = 1,Ncon
            IF (psend(k,j) .EQ. 1) THEN
                DO i = 1,lda
                   sendbuffer(ipos) = cd_local(i,j)
                   ipos = ipos + 1
                ENDDO
            ENDIF
         ENDDO

         tag = 200
         CALL MPI_SendRecv(sendbuffer,csend(k)*lda,MPI_INTEGER,         &
         &    ppm_isendlist(k),tag,recvbuffer,crecv(k)*lda,MPI_INTEGER, &
         &    ppm_irecvlist(k),tag,ppm_comm,MPI_STATUS_IGNORE,info)
         or_fail_MPI("MPI_SendRecv")

         !----------------------------------------------------------------------
         !  Add the received connections to the connection data array
         !----------------------------------------------------------------------
         IF (crecv(k) .NE. 0) THEN
             ipos = 1
             DO j = 1,crecv(k)
                DO i = 1,lda
                   cd_local(i,Ncon+icon) = recvbuffer(ipos)
                   ipos = ipos + 1
                ENDDO
                icon = icon + 1
             ENDDO
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  We now have a new number of connections
      !-------------------------------------------------------------------------
      ncons_local = Ncon + mrecv

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(csend,ldu,iopt,info)
      or_fail_dealloc('deallocating array CSEND',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(crecv,ldu,iopt,info)
      or_fail_dealloc('deallocating array CRECV',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(sendbuffer,ldu,iopt,info)
      or_fail_dealloc('deallocating send buffer SENDBUFFER',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(recvbuffer,ldu,iopt,info)
      or_fail_dealloc('deallocating receive buffer RECVBUFFER',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Invert list of particles that has been sent away
      !-------------------------------------------------------------------------
      CALL ppm_util_invert_list(id_send,id_inv,info)
      or_fail("inverting id list",ppm_error=ppm_error_fatal)

      iopt = ppm_param_dealloc
      CALL ppm_alloc(id_send,ldu,iopt,info)
      or_fail_dealloc('deallocating temporary id list ID_SEND',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Find and mark all orphan connections
      !-------------------------------------------------------------------------
      min_=LBOUND(id_inv,1)
      max_=UBOUND(id_inv,1)

      icon = 0
      DO j = 1,Ncon
         IF (psend(1,j) .EQ. 1) THEN
             l = 0
             DO i = 1,lda
                IF ((cd_local(i,j) .GE. min_) .AND.  &
                &   (cd_local(i,j) .LE. max_)) THEN
                   IF (id_inv(cd_local(i,j)) .GT. -big) THEN
                      l = l + 1
                   ENDIF
                ELSE
                   l = l + 1
                ENDIF
             ENDDO
             IF (l .EQ. lda) THEN
                 cd_local(1,j) = -1
                 icon = icon + 1
             ENDIF
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Mark all connections that we have more than once
      !-------------------------------------------------------------------------
      DO j = 1,ncons_local-1
         IF (cd_local(1,j) .EQ. -1) CYCLE

         DO k = j+1,ncons_local
            l = 0
            DO i = 1,lda
               IF (cd_local(i,j).EQ.cd_local(i,k)) THEN
                  l = l + 1
               ENDIF
            ENDDO

            IF (l .EQ. lda) THEN
                cd_local(1,k) = -1
                icon = icon + 1
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Prepare the final connection array
      !-------------------------------------------------------------------------
      ldu(1) = lda
      ldu(2) = ncons_local - icon
      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(cd,ldu,iopt,info)
      or_fail_alloc('allocating connection array CD',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Go through all connections and save all that are not marked
      !-------------------------------------------------------------------------
      Ncon = 0
      DO j = 1,ncons_local
         IF (cd_local(1,j) .NE. -1) THEN
             Ncon = Ncon + 1
             cd(:,Ncon) = cd_local(:,j)
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_dealloc('deallocating array PSEND',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(cd_local,ldu,iopt,info)
      or_fail_dealloc('deallocating connection array CD_LOCAL',ppm_error=ppm_error_fatal)
#endif
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Ncon.LT.0) THEN
           fail('Number of connections must be >=0', &
           & exit_point=8888,ppm_error=ppm_error_fatal)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_connect_send
