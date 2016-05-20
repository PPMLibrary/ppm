      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                 ppm_dbg_print
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
#if   __CTAG == __SCALAR
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE dbg_print_sca_s(topoid,ghostlayer,step,colortag,info,xp,np,mp,append)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE dbg_print_sca_d(topoid,ghostlayer,step,colortag,info,xp,np,mp,append)
#endif
#elif   __CTAG == __VECTOR
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE dbg_print_vec_s(topoid,ghostlayer,step,colortag,info,xp,np,mp,append)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE dbg_print_vec_d(topoid,ghostlayer,step,colortag,info,xp,np,mp,append)
#endif
#endif
      !!! This routine provides a simple means to visualize particles and
      !!! domain decompositions for debugging and monitoring purposes.
      !!!
      !!! The routine is used in conjunction with the ppmdbg.py script.
      !!! It creates two files with names ppm_dbg_####.sub and ppm_dbg_####.dat
      !!! That contain the domain decomposition and particle data respectively

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      USE ppm_module_data
      USE ppm_module_topo
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      ! arguments
      INTEGER,            INTENT(IN)    :: topoid
      !!! topology ID to which we are currently mapped
      REAL(mk),           INTENT(IN)    :: ghostlayer
      !!! ghostlayer, can be set to 0 if we dont care about it
      INTEGER,            INTENT(IN)    :: step
      !!! parameter can be used to create distinct output dump files for each
      !!! timestep
#if   __CTAG == __SCALAR
      INTEGER,            INTENT(IN)    :: colortag
#elif __CTAG == __VECTOR
      INTEGER, DIMENSION(:), POINTER    :: colortag
#endif
      !!! a tag to be able to print out different groups of particles or
      !!! visualize a property
      INTEGER,            INTENT(OUT)   :: info
      REAL(mk), DIMENSION(:,:), POINTER,OPTIONAL :: xp
      !!! a particle position array, this argument is optional
      INTEGER,            INTENT(IN),OPTIONAL    :: np
      !!! number of particles
      !!! if xp is provided, THEN this one must e provided too
      INTEGER,            INTENT(IN),OPTIONAL    :: mp
      !!! number of particles including ghost particles
      !!! if omitted it is assumed that np=mp
      LOGICAL,            INTENT(IN),OPTIONAL    :: append
      !!! Should the particle positions be appended to an existing file or the
      !!! file overwritten (default).
      ! local vars
      CHARACTER(128)                     :: sfmt,pfmt
      CHARACTER(64)                     :: sfname,pfname
      TYPE(ppm_t_topo), POINTER         :: topo
      INTEGER                           :: i
      INTEGER                           :: iunit
      REAL(mk)                          :: t0
      INTEGER                           :: mpart
#ifdef __MPI
      INTEGER, DIMENSION(:),    POINTER :: allnp
      INTEGER, DIMENSION(:),    POINTER :: allmp
#if   __CTAG == __SCALAR
      INTEGER, DIMENSION(:),    POINTER :: allctag => NULL()
#elif __CTAG == __VECTOR
      INTEGER, DIMENSION(:,:),  POINTER :: allctag => NULL()
#endif
      REAL(mk), DIMENSION(:,:,:),POINTER:: allxp   => NULL()
      INTEGER, DIMENSION(3)             :: lda
      INTEGER                           :: maxmp
      INTEGER                           :: iproc
#endif

      CALL substart('ppm_dbg_print',t0,info)


      !------------------------------------------------------------------------
      ! Prepare format strings and file names for I/O
      !------------------------------------------------------------------------
      iunit = 2342
      WRITE(sfmt,'(A,I1,A,I1,A)') '(',ppm_dim*2,'F12.8,I4,',ppm_dim*2,'I2)'
      WRITE(pfmt,'(A,I1,A)') '(',ppm_dim,'E18.9,I4)'
      WRITE(sfname,'(A,I5.5,A)') 'ppm_dbg_',step,'.sub'
      WRITE(pfname,'(A,I5.5,A)') 'ppm_dbg_',step,'.dat'


      !------------------------------------------------------------------------
      ! Write domain decomp and topology info.
      ! In the MPI case, only rank 0 creates this file
      !------------------------------------------------------------------------

#ifdef __MPI
      IF (ppm_rank.eq.0) THEN
#endif
      NULLIFY(topo)
      CALL ppm_topo_get(topoid,topo,info)
      OPEN(iunit,file=sfname)

      WRITE(iunit,'(I1)') ppm_dim
      WRITE(iunit,'(F12.8)') ghostlayer
      DO i=1,topo%nsubs
         WRITE(iunit,sfmt) topo%min_subd(:,i),topo%max_subd(:,i),topo%sub2proc(i),topo%subs_bc(:,i)
      ENDDO
      CLOSE(iunit)
#ifdef __MPI
      ENDIF
#endif

      IF (PRESENT(xp).AND.PRESENT(np)) THEN
          IF (PRESENT(mp)) THEN
              mpart = mp
          ELSE
              mpart = np
          ENDIF

#ifdef __MPI
          !--------------------------------------------------------------------
          ! Send all data to rank 0
          !--------------------------------------------------------------------
          ! first allocate the size info arrays
          NULLIFY(allnp,allmp)
          lda(1) = ppm_nproc
          CALL ppm_alloc(allnp,lda,ppm_param_alloc_fit,info)
          CALL ppm_alloc(allmp,lda,ppm_param_alloc_fit,info)
#if   __CTAG == __SCALAR
          CALL ppm_alloc(allctag,lda,ppm_param_alloc_fit,info)
#endif
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_dbg_print',     &
     &            'failed to allocate allnp, allmp or allctag',__LINE__,info)
              GOTO 9999
          ENDIF
          ! gather the np and mp at the root
          CALL MPI_Gather(np,1,MPI_INTEGER,allnp,1,MPI_INTEGER,0,ppm_comm,info)
          CALL MPI_Gather(mpart,1,MPI_INTEGER,allmp,1,MPI_INTEGER,0,ppm_comm,info)
#if   __CTAG == __SCALAR
          CALL MPI_Gather(colortag,1,MPI_INTEGER,allctag,1,MPI_INTEGER,0,&
          &               ppm_comm,info)
#endif
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_dbg_print',     &
     &            'failed to gather allnp, allmp or allctag',__LINE__,info)
              GOTO 9999
          ENDIF

          IF (ppm_rank.eq.0) THEN
              ! allocate allxp array
              maxmp = maxval(allmp)
              lda(1) = ppm_dim
              lda(2) = maxmp
              lda(3) = ppm_nproc
              CALL ppm_alloc(allxp,lda,ppm_param_alloc_fit,info)
#if __CTAG == __VECTOR
              lda(1) = maxmp
              lda(2) = ppm_nproc
              CALL ppm_alloc(allctag,lda,ppm_param_alloc_fit,info)
#endif
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_dbg_print',     &
         &            'failed to allocate allxp',__LINE__,info)
                  GOTO 9999
              ENDIF

              DO i=1,allmp(1)
                  allxp(:,i,1) = xp(:,i)
#if __CTAG == __VECTOR
                  allctag(i,1) = colortag(i)
#endif
              ENDDO

              ! now let all procs communicate with rank 0
              DO iproc=1,ppm_nproc-1
                  CALL MPI_Recv(allxp(:,:,iproc+1),allmp(iproc+1)*ppm_dim,&
                  &             ppm_mpi_kind,iproc,  &
                  &             0,ppm_comm,MPI_STATUS_IGNORE,info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_mpi_fail,'ppm_dbg_print',   &
     &                'failed to sendrecv xp',__LINE__,info)
                      GOTO 9999
                  ENDIF
#if __CTAG == __VECTOR
                  CALL MPI_Recv(allctag(:,iproc+1),allmp(iproc+1),&
                  &             MPI_INTEGER,iproc,  &
                  &             0,ppm_comm,MPI_STATUS_IGNORE,info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_mpi_fail,'ppm_dbg_print',   &
     &                'failed to sendrecv ctag',__LINE__,info)
                      GOTO 9999
                  ENDIF
#endif
              ENDDO
              IF (PRESENT(append).AND.append) THEN
                  OPEN(iunit,file=pfname,access='append')
              ELSE
                  OPEN(iunit,file=pfname)
              ENDIF
              DO iproc=1,ppm_nproc
                  DO i=1,allnp(iproc)
#if __CTAG == __SCALAR
                      WRITE(iunit,pfmt) allxp(:,i,iproc),allctag(iproc)
#elif __CTAG == __VECTOR
                      WRITE(iunit,pfmt) allxp(:,i,iproc),allctag(i,iproc)
#endif
                  ENDDO
                  DO i=allnp(iproc)+1,allmp(iproc)
#if __CTAG == __SCALAR
                      WRITE(iunit,pfmt) allxp(:,i,iproc),-1
#elif __CTAG == __VECTOR
                      WRITE(iunit,pfmt) allxp(:,i,iproc),allctag(i,iproc)
#endif
                  ENDDO
              ENDDO
              CLOSE(iunit)
              CALL ppm_alloc(allxp,lda,ppm_param_dealloc,info)
          ELSE
              CALL MPI_Send(xp,mpart*ppm_dim,ppm_mpi_kind,0,0,ppm_comm,info)
#if __CTAG == __VECTOR
              CALL MPI_Send(colortag,mpart,MPI_INTEGER,0,0,ppm_comm,info)
#endif
          ENDIF
#else
          IF (PRESENT(append).AND.append) THEN
              OPEN(iunit,FILE=pfname,ACCESS='APPEND')
          ELSE
              OPEN(iunit,FILE=pfname)
          ENDIF
          DO i=1,np
#if __CTAG == __SCALAR
              WRITE(iunit,pfmt) xp(:,i),colortag
#elif __CTAG == __VECTOR
              WRITE(iunit,pfmt) xp(:,i),colortag(i)
#endif
          ENDDO
          DO i=np+1,mpart
#if __CTAG == __SCALAR
              WRITE(iunit,pfmt) xp(:,i),-1
#elif __CTAG == __VECTOR
              WRITE(iunit,pfmt) xp(:,i),colortag(i)
#endif
          ENDDO
          CLOSE(iunit)
#endif
      ENDIF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL ppm_alloc(allnp,lda,ppm_param_dealloc,info)
      CALL ppm_alloc(allmp,lda,ppm_param_dealloc,info)
#if   __CTAG == __SCALAR
      CALL ppm_alloc(allctag,lda,ppm_param_dealloc,info)
#elif __CTAG == __VECTOR
      CALL ppm_alloc(allctag,lda,ppm_param_dealloc,info)
#endif
      CALL ppm_alloc(allxp,(/0,0,0/),ppm_param_dealloc,info)
#endif

      9999  CONTINUE
      CALL substop('ppm_dbg_print',t0,info)

#if   __CTAG == __SCALAR
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE dbg_print_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE dbg_print_sca_d
#endif
#elif   __CTAG == __VECTOR
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE dbg_print_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE dbg_print_vec_d
#endif
#endif
