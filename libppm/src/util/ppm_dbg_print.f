      !-------------------------------------------------------------------------
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
subroutine ppm_dbg_print_s(topoid,ghostlayer,step,colortag,info,xp,np,mp)
#elif __KIND == __DOUBLE_PRECISION
subroutine ppm_dbg_print_d(topoid,ghostlayer,step,colortag,info,xp,np,mp)
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

      use ppm_module_data
      use ppm_module_topo
      use ppm_module_alloc
      use ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      
      implicit none
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif 

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    
      ! arguments
      integer,            intent(in)    :: topoid
      !!! topology ID to which we are currently mapped
      real(mk),           intent(in)    :: ghostlayer
      !!! ghostlayer, can be set to 0 if we don't care about it
      integer,            intent(in)    :: step
      !!! parameter can be used to create distinct output dump files for each
      !!! timestep
      integer,            intent(in)    :: colortag
      !!! a tag to be able to print out different groups of particles
      integer,            intent(out)   :: info
      real(mk), dimension(:,:), pointer,optional :: xp
      !!! a particle position array, this argument is optional
      integer,            intent(in),optional    :: np
      !!! number of particles
      !!! if xp is provided, then this one must e provided too
      integer,            intent(in),optional    :: mp
      !!! number of particles including ghost particles
      !!! if omitted it is assumed that np=mp
      ! local vars
      character(128)                     :: sfmt,pfmt
      character(64)                     :: sfname,pfname
      type(ppm_t_topo), pointer         :: topo => NULL()
      integer                           :: i
      integer                           :: iunit
      real(mk)                          :: t0
      integer                           :: mpart
#ifdef __MPI
      integer, dimension(:),    pointer :: allnp => NULL()
      integer, dimension(:),    pointer :: allmp => NULL()
      real(mk), dimension(:,:,:),pointer:: allxp => NULL()
      integer, dimension(:,:),  pointer :: buf   => NULL()
      integer, dimension(3)             :: lda
      integer                           :: maxmp
      integer                           :: iproc
#endif
       
      CALL substart('ppm_dbg_print',t0,info)


      !------------------------------------------------------------------------
      ! Prepare format strings and file names for I/O
      !------------------------------------------------------------------------
      iunit = 2342
      write(sfmt,'(A,I1,A,I1,A)') '(',ppm_dim*2,'F12.8,I4,',ppm_dim*2,'I2)'
      write(pfmt,'(A,I1,A)') '(',ppm_dim,'E18.9,I4)'
      write(sfname,'(A,I5.5,A)') 'ppm_dbg_',step,'.sub'
      write(pfname,'(A,I5.5,A)') 'ppm_dbg_',step,'.dat'


      !------------------------------------------------------------------------
      ! Write domain decomp and topology info.
      ! In the MPI case, only rank 0 creates this file
      !------------------------------------------------------------------------

#ifdef __MPI
      if (ppm_rank.eq.0) then
#endif
      call ppm_topo_get(topoid,topo,info)
      open(iunit,file=sfname)
   
      write(iunit,'(I)') ppm_dim
      write(iunit,'(F12.8)') ghostlayer
      do i=1,topo%nsubs
          write(iunit,sfmt) topo%min_subd(:,i),topo%max_subd(:,i),&
 &                          topo%sub2proc(i),&
 &                          topo%subs_bc(:,i)
      enddo
      close(iunit)
#ifdef __MPI
      endif
#endif
      
      if (present(xp).and.present(np)) then
          if (present(mp)) then
              mpart = mp
          else
              mpart = np
          endif

#ifdef __MPI
          !--------------------------------------------------------------------
          ! Send all data to rank 0
          !--------------------------------------------------------------------
          ! first allocate the size info arrays
          lda(1) = ppm_nproc
          call ppm_alloc(allnp,lda,ppm_param_alloc_fit,info)
          call ppm_alloc(allmp,lda,ppm_param_alloc_fit,info)
          if (info .NE. 0) then
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_dbg_print',     &
     &            'failed to allocate allnp or allmp',__LINE__,info)
              goto 9999
          endif
          ! gather the np and mp at the root
          call mpi_gather(np,1,MPI_INTEGER,allnp,1,MPI_INTEGER,0,ppm_comm,info)
          call mpi_gather(mpart,1,MPI_INTEGER,allmp,1,MPI_INTEGER,0,ppm_comm,info)
          if (info .NE. 0) then
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_dbg_print',     &
     &            'failed to gather allnp or allmp',__LINE__,info)
              goto 9999
          endif
          
          ! allocate allxp array
          maxmp = maxval(allmp)
          lda(1) = ppm_dim
          lda(2) = maxmp
          lda(3) = ppm_nproc
          call ppm_alloc(allxp,lda,ppm_param_alloc_fit,info)
          call ppm_alloc(buf,lda,ppm_param_alloc_fit,info)
          if (info .NE. 0) then
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_dbg_print',     &
     &            'failed to allocate allxp or buf',__LINE__,info)
              goto 9999
          endif
          
          ! now let all procs communicate with rank 0
          if (ppm_rank.eq.0) then
              do i=1,allmp(1)
                  allxp(:,i,1) = xp(:,i)
              enddo

              do iproc=2,ppm_nproc
                  call mpi_sendrecv(xp,allmp(iproc),ppm_mpi_kind,buf, &
 &                              0,0,allmp(iproc),ppm_mpi_kind,iproc,  &
 &                              0,ppm_comm,info)
                  if (info .NE. 0) then
                      info = ppm_error_fatal
                      call ppm_error(ppm_err_mpi_fail,'ppm_dbg_print',   &
     &                'failed to sendrecv xp',__LINE__,info)
                      goto 9999
                  endif
                  do i=1,allmp(iproc)
                      allxp(:,i,iproc) = buf(:,i)
                  enddo
              enddo
          else
              call mpi_sendrecv(xp,mpart,ppm_mpi_kind,buf,0,0,mpart,&
 &                              ppm_mpi_kind,ppm_rank,0,ppm_comm,info)
          endif
          open(iunit,file=pfname,access='append')
          do iproc=1,ppm_nproc
              do i=1,allnp(iproc)
                  write(iunit,pfmt) allxp(:,i,iproc),colortag
              enddo
              do i=allnp(iproc)+1,allmp(iproc)
                  write(iunit,pfmt) allxp(:,i,iproc),-1
              enddo
          enddo
          close(iunit)
          call ppm_alloc(allxp,lda,ppm_param_dealloc,info)
          call ppm_alloc(buf,lda,ppm_param_dealloc,info)

#else  
          open(iunit,file=pfname,access='append')
          do i=1,np
              write(iunit,pfmt) xp(:,i),colortag
          enddo
          do i=np+1,mpart
              write(iunit,pfmt) xp(:,i),-1
          enddo
          close(iunit)
#endif
      endif
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
9999  continue
      CALL substop('ppm_dbg_print',t0,info)

#if   __KIND == __SINGLE_PRECISION
end subroutine ppm_dbg_print_s
#elif __KIND == __DOUBLE_PRECISION
end subroutine ppm_dbg_print_d
#endif
