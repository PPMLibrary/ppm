      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                 ppm_dbg_print
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
SUBROUTINE ppm_dbg_print_s(topoid,ghostlayer,step,colortag,info,xp,np,mp,append)
#elif __KIND == __DOUBLE_PRECISION
SUBROUTINE ppm_dbg_print_d(topoid,ghostlayer,step,colortag,info,xp,np,mp,append)
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
      
      implicit none
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
      INTEGER,            INTENT(IN)    :: colortag
      !!! a tag to be able to print out different groups of particles
      INTEGER,            INTENT(OUT)   :: info
      REAL(mk), DIMENSION(:,:), POINTER,OPTIONAL :: xp
      !!! a particle position array, this argument is optional
      INTEGER,            INTENT(IN),OPTIONAL    :: np
      !!! number of particles
      !!! if xp is provided, then this one must e provided too
      INTEGER,            INTENT(IN),OPTIONAL    :: mp
      !!! number of particles including ghost particles
      !!! if omitted it is assumed that np=mp
      LOGICAL,            INTENT(IN),OPTIONAL    :: append
      !!! Should the particle positions be appended to an existing file or the
      !!! file overwritten (default).
      ! local vars
      CHARACTER(128)                     :: sfmt,pfmt
      CHARACTER(64)                     :: sfname,pfname
      TYPE(ppm_t_topo), POINTER         :: topo => NULL()
      INTEGER                           :: i
      INTEGER                           :: iunit
      REAL(mk)                          :: t0
      INTEGER                           :: mpart
#ifdef __MPI
      INTEGER, DIMENSION(:),    POINTER :: allnp => NULL()
      INTEGER, DIMENSION(:),    POINTER :: allmp => NULL()
      REAL(mk), DIMENSION(:,:,:),POINTER:: allxp => NULL()
      INTEGER, DIMENSION(:,:),  POINTER :: buf   => NULL()
      INTEGER, DIMENSION(3)             :: lda
      INTEGER                           :: maxmp
      INTEGER                           :: iproc
      INTEGER, DIMENSION(:),  POINTER   :: req => NULL()
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
      CALL ppm_topo_get(topoid,topo,info)
      OPEN(iunit,file=sfname)
   
      WRITE(iunit,'(I1)') ppm_dim
      WRITE(iunit,'(F12.8)') ghostlayer
      DO i=1,topo%nsubs
          WRITE(iunit,sfmt) topo%min_subd(:,i),topo%max_subd(:,i),&
 &                          topo%sub2proc(i),&
 &                          topo%subs_bc(:,i)
      ENDDO
      CLOSE(iunit)
#ifdef __MPI
      ENDIF
#endif
      
      IF (present(xp).AND.present(np)) THEN
          IF (present(mp)) THEN
              mpart = mp
          ELSE
              mpart = np
          ENDIF

#ifdef __MPI
          !--------------------------------------------------------------------
          ! Send all data to rank 0
          !--------------------------------------------------------------------
          ! first allocate the size info arrays
          lda(1) = ppm_nproc
          CALL ppm_alloc(allnp,lda,ppm_param_alloc_fit,info)
          CALL ppm_alloc(allmp,lda,ppm_param_alloc_fit,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_dbg_print',     &
     &            'failed to allocate allnp or allmp',__LINE__,info)
              GOTO 9999
          ENDIF
          ! gather the np and mp at the root
          CALL mpi_gather(np,1,MPI_INTEGER,allnp,1,MPI_INTEGER,0,ppm_comm,info)
          CALL mpi_gather(mpart,1,MPI_INTEGER,allmp,1,MPI_INTEGER,0,ppm_comm,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_dbg_print',     &
     &            'failed to gather allnp or allmp',__LINE__,info)
              GOTO 9999
          ENDIF
          
          IF (ppm_rank.eq.0) THEN
              ! allocate allxp array
              maxmp = maxval(allmp)
              lda(1) = ppm_dim
              lda(2) = maxmp
              lda(3) = ppm_nproc
              CALL ppm_alloc(allxp,lda,ppm_param_alloc_fit,info)
              CALL ppm_alloc(buf,lda,ppm_param_alloc_fit,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_dbg_print',     &
         &            'failed to allocate allxp or buf',__LINE__,info)
                  GOTO 9999
              ENDIF
              
              lda(1) = ppm_nproc
              CALL ppm_alloc(req,lda,ppm_param_alloc_fit,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_dbg_print',     &
         &            'failed to allocate req',__LINE__,info)
                  GOTO 9999
              ENDIF
              
              DO i=1,allmp(1)
                  allxp(:,i,1) = xp(:,i)
              ENDDO

              ! now let all procs communicate with rank 0
              DO iproc=1,ppm_nproc-1
!                  CALL mpi_irecv(buf,allmp(iproc+1)*ppm_dim,ppm_mpi_kind,iproc,  &
! &                              0,ppm_comm,req(iproc+1),info)
                  CALL mpi_recv(allxp(:,:,iproc+1),allmp(iproc+1)*ppm_dim,ppm_mpi_kind,iproc,  &
 &                              0,ppm_comm,MPI_STATUS_IGNORE,info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_mpi_fail,'ppm_dbg_print',   &
     &                'failed to sendrecv xp',__LINE__,info)
                      GOTO 9999
                  ENDIF
              ENDDO
!              do iproc=1,ppm_nproc-1
!                  CALL mpi_wait(req(iproc+1),MPI_STATUS_IGNORE,info)
!                  do i=1,allmp(iproc+1)
!                      allxp(:,i,iproc+1) = buf(:,i)
!                  ENDDO
!                  print *,'xp',iproc,allxp(:,1:allmp(iproc+1),iproc+1)
!              ENDDO
              IF (present(append).AND.append) then
                  OPEN(iunit,file=pfname,access='append')
              ELSE
                  OPEN(iunit,file=pfname)
              ENDIF
              DO iproc=1,ppm_nproc
                  DO i=1,allnp(iproc)
                      WRITE(iunit,pfmt) allxp(:,i,iproc),colortag
                  ENDDO
                  DO i=allnp(iproc)+1,allmp(iproc)
                      WRITE(iunit,pfmt) allxp(:,i,iproc),-1
                  ENDDO
              ENDDO
              CLOSE(iunit)
              CALL ppm_alloc(allxp,lda,ppm_param_dealloc,info)
              CALL ppm_alloc(buf,lda,ppm_param_dealloc,info)
          ELSE
              CALL mpi_send(xp,mpart*ppm_dim,ppm_mpi_kind,0,0,ppm_comm,info)
          ENDIF
#else  
          IF (present(append).AND.append) then
              OPEN(iunit,file=pfname,access='append')
          ELSE
              OPEN(iunit,file=pfname)
          ENDIF
          DO i=1,np
              WRITE(iunit,pfmt) xp(:,i),colortag
          ENDDO
          DO i=np+1,mpart
              WRITE(iunit,pfmt) xp(:,i),-1
          ENDDO
          CLOSE(iunit)
#endif
      ENDIF
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_dbg_print',t0,info)

#if   __KIND == __SINGLE_PRECISION
end SUBROUTINE ppm_dbg_print_s
#elif __KIND == __DOUBLE_PRECISION
end SUBROUTINE ppm_dbg_print_d
#endif
