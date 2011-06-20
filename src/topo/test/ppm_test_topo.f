      !<<<< haeckic begin >>>>!
      !-------------------------------------------------------------------------
      !     Test Case   :                   ppm_test_topo
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

      program ppm_test_topo
      !-------------------------------------------------------------------------
      !     Test Case   :  topologies, i.e. decomposition and proc assignment
      !-------------------------------------------------------------------------

      use ppm_module_typedef
      USE ppm_module_data
      use ppm_module_mktopo
      use ppm_module_topo_get
      use ppm_module_init
      use ppm_module_finalize
      use ppm_module_user
      USE ppm_module_check_id

      implicit none
      #include "../../ppm_define.h"
      #ifdef __MPI
      INCLUDE 'mpif.h'
      #endif

      integer, parameter              :: debug = 0
      integer, parameter              :: mk = ppm_kind_double
      real(mk),parameter              :: pi = 3.1415926535897931_mk
      real(mk),parameter              :: skin = 0._mk
      integer,parameter               :: ndim=2
      integer,parameter               :: pdim=2
      integer                         :: decomp,assig,tolexp
      real(mk)                        :: tol,min_ghost_req,max_ghost_req
      integer                         :: info,comm,rank
      integer                         :: topoid
      integer                         :: np
      integer                         :: mp
      integer                         :: newnp
      real(mk),dimension(:,:),pointer :: xp
      real(mk),dimension(:,:),pointer :: ghost_req
      real(mk),dimension(:  ),pointer :: min_phys,max_phys,h,p_h
      real(mk),dimension(:  ),pointer :: len_phys
      integer                         :: i,j,k,sum1,sum2, ix,iy
      integer                         :: p_i
      integer, dimension(6)           :: bcdef
      real(mk),dimension(:  ),pointer :: cost
      integer, dimension(:  ),pointer :: nm
      type(ppm_t_topo), pointer       :: topo => NULL()
      integer                         :: seedsize
      integer,  dimension(:),allocatable :: seed
      real(mk), dimension(:),allocatable :: randnb
      integer                          :: isymm = 0
      logical                          :: lsymm = .false.,ok, has_one_way
      real(mk)                         :: t0,t1,t2,t3

      integer                          :: np_sqrt = 16
      integer                          :: np_tot, ipart, idom
      integer                          :: tot_sum
      integer, dimension(:),pointer    :: so_sum

      !----------------
      ! setup
      !----------------
      tol = 10.0_mk*epsilon(1.0_mk)
      tolexp = int(log10(epsilon(1.0_mk)))
      min_ghost_req = 0.25_mk
      max_ghost_req = 1.25_mk
      np = np_sqrt*np_sqrt
      np_tot = 4*np

      allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),nm(ndim),h(ndim),p_h(ndim),stat=info)

      min_phys(1:ndim) = 0.0_mk
      max_phys(1:ndim) = 8.0_mk
      len_phys(1:ndim) = max_phys-min_phys
      bcdef(1:6) = ppm_param_bcdef_periodic

      nullify(xp,ghost_req)

#ifdef __MPI
      comm = mpi_comm_world
      call mpi_init(info)
      call mpi_comm_rank(comm,rank,info)
#else
      rank = 0
#endif
      call ppm_init(ndim,mk,tolexp,0,debug,info,99)

      call random_seed(size=seedsize)
      allocate(seed(seedsize))
      allocate(randnb(ndim*np),stat=info)
      do i=1,seedsize
         seed(i)=10+i*i*(rank+1)
      enddo
      call random_seed(put=seed)
      call random_number(randnb)

      !----------------
      ! create particles
      !----------------

      allocate(xp(ndim,np_tot),ghost_req(ndim,np_tot),stat=info)
      xp = 0.0_mk
      ghost_req = 0.0_mk

      allocate(so_sum(ndim),stat=info)
      tot_sum = np_sqrt*(np_sqrt+1)/2 + 1

      so_sum(1) = 0.0_mk
      do ix=1,np_sqrt
         so_sum(1) = so_sum(1) + ix
         so_sum(2) = 0.0_mk

         do iy=1,np_sqrt
            so_sum(2) = so_sum(2) + iy

            ! set positions of particles, s.t. lower ix,iy are closer together
            ! including a random distortion            
            xp(1,(ix-1)*np_sqrt + iy) = min_phys(1)+REAL(so_sum(1),MK)* &
     &                                  len_phys(1)/REAL(tot_sum,MK) +  &
     &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
     &                                  - (ix*len_phys(1)/(tot_sum))/2


            do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
               xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
            enddo
            do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
               xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
            enddo

            xp(2,(ix-1)*np_sqrt + iy) = min_phys(2)+REAL(so_sum(2),MK)* &
     &                                  len_phys(2)/REAL(tot_sum,MK) +  &
     &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np) & 
     &                                  - (iy*len_phys(2)/(tot_sum))/2

            do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
               xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
            enddo
            do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
               xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
            enddo
            ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
            ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
            ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

!                 print *,(ix-1)*np_sqrt + iy,' ',xp(1,(ix-1)*np_sqrt + iy),' ',xp(2,(ix-1)*np_sqrt + iy), &
!         &                  ' ', ghost_req(1,(ix-1)*np_sqrt + iy), ' ', ghost_req(2,(ix-1)*np_sqrt + iy)


         enddo

      enddo

      ! take mirror
      do ix=1,np_sqrt
         do iy=1,np_sqrt
            xp(1,(ix-1)*np_sqrt + iy + np) = - xp(1,(ix-1)*np_sqrt + iy)
            xp(2,(ix-1)*np_sqrt + iy + np) = xp(2,(ix-1)*np_sqrt + iy)
            ghost_req(1,(ix-1)*np_sqrt + iy + np) = ghost_req(1,(ix-1)*np_sqrt + iy)
            ghost_req(2,(ix-1)*np_sqrt + iy + np) = ghost_req(2,(ix-1)*np_sqrt + iy)
         enddo
      enddo

      do ix=1,np_sqrt
         do iy=1,np_sqrt
            xp(1,(ix-1)*np_sqrt + iy + np + np) = xp(1,(ix-1)*np_sqrt + iy)
            xp(2,(ix-1)*np_sqrt + iy + np + np) = - xp(2,(ix-1)*np_sqrt + iy)
            ghost_req(1,(ix-1)*np_sqrt + iy + np + np) = ghost_req(1,(ix-1)*np_sqrt + iy)
            ghost_req(2,(ix-1)*np_sqrt + iy + np + np) = ghost_req(2,(ix-1)*np_sqrt + iy)
         enddo
      enddo

      do ix=1,np_sqrt
         do iy=1,np_sqrt
            xp(1,(ix-1)*np_sqrt + iy + np + np + np) = - xp(1,(ix-1)*np_sqrt + iy)
            xp(2,(ix-1)*np_sqrt + iy + np + np + np) = - xp(2,(ix-1)*np_sqrt + iy)
            ghost_req(1,(ix-1)*np_sqrt + iy + np + np + np) = ghost_req(1,(ix-1)*np_sqrt + iy)
            ghost_req(2,(ix-1)*np_sqrt + iy + np + np + np) = ghost_req(2,(ix-1)*np_sqrt + iy)
         enddo
      enddo

      !----------------
      ! make topology
      !----------------

      ! Test all decompositions
!       decomp = ppm_param_decomp_tree
!       decomp = ppm_param_decomp_pruned_cell
       decomp = ppm_param_decomp_bisection
!       decomp = ppm_param_decomp_xpencil
!       decomp = ppm_param_decomp_ypencil
!       decomp = ppm_param_decomp_zpencil
!       decomp = ppm_param_decomp_cuboid
!       decomp = ppm_param_decomp_user_defined
!       decomp = ppm_param_decomp_xy_slab
!       decomp = ppm_param_decomp_xz_slab
!       decomp = ppm_param_decomp_yz_slab
!       decomp = ppm_param_decomp_cartesian

      assig  = ppm_param_assign_internal

      topoid = 0



      min_phys(1:ndim) = -8.0_mk
      max_phys(1:ndim) = 8.0_mk 
!       call ppm_mktopo(topoid,xp,np_tot,decomp,assig,min_phys,max_phys,bcdef,&
!       &               max_ghost_req,cost,info)

      has_one_way = .TRUE.

       call ppm_mktopo(topoid,xp,np_tot,decomp,assig,min_phys,max_phys,bcdef,&
     &               ghost_req,has_one_way,cost,info)

      call ppm_dbg_print(topoid,max_ghost_req,1,1,info,xp,np_tot)

      call ppm_map_part_global(topoid,xp,np,info)
      call ppm_map_part_push(ghost_req,ndim,np,info)
      call ppm_map_part_send(np,newnp,info)
      call ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
      call ppm_map_part_pop(xp,ndim,np,newnp,info)
      np=newnp
      print *,rank,'global map done'


      ! MAKE TESTS:

      ! 1. Are all particle this proc have in its subs
      call ppm_topo_check(topoid,xp,np,ok,info)
      if (.not. ok) write(*,*) '[',rank,'] topo_check failed'

      ! 2. check if all particles to interact with are in neighboring box
      !         for each particle determine interacting particles and
      !         check if they are in own box or neighbors
      call ppm_topo_check_neigh(topoid,xp,ghost_req,np,has_one_way,ok,info)
      if (.not. ok) write(*,*) '[',rank,'] topo_check_neigh failed'

      ! 3. minboxsizes inside, requires all neighbors correct
      !         check for each box if it fulfills ghost req of each particle inside
      !         this function also checks for has_one way case
       call ppm_topo_check_minbox(topoid,xp,ghost_req,np,has_one_way,ok,info)
       if (.not. ok) write(*,*) '[',rank,'] topo_check_minbox failed'




      !----------------
      ! teardown...
      !----------------
 8000 continue

      call ppm_finalize(info)
      #ifdef __MPI
      call MPI_finalize(info)
      #endif

      deallocate(xp,ghost_req,min_phys,max_phys,len_phys,nm)

      if (rank.eq.0) print *, 'done.'


      end program ppm_test_topo
!<<<< haeckic end >>>>!