test_suite ppm_module_loadbal

use ppm_module_particles_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_operator_typedef
use ppm_module_interfaces
use ppm_module_data
use ppm_module_loadbal
use ppm_module_io_vtk

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = ACOS(-1._mk)
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc,topoid
integer                         :: np_global = 3000
real(mk),parameter              :: cutoff = 0.15_mk
real(mk),dimension(:,:),pointer :: xp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys=>NULL(),max_phys=>NULL()
real(mk),dimension(:  ),pointer :: len_phys=>NULL()
integer                         :: i,j,k,ip,wp_id
integer                         :: nstep
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
integer                         :: isymm = 0
real(mk)                        :: t0,t1,t2,t3,mat_time
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()
real(mk)                        :: err,start_time,end_time,elapsed_t
real(mk)                        :: t_comp,t_comm,t_total

integer, dimension(:), pointer                 :: wp_1i => NULL()
integer, dimension(:,:), pointer               :: wp_2i => NULL()
integer(ppm_kind_int64),dimension(:),  pointer :: wp_1li => NULL()
integer(ppm_kind_int64),dimension(:,:),pointer :: wp_2li => NULL()
integer, dimension(:),   pointer              :: wp_1r => NULL()
real(mk), dimension(:,:), pointer              :: wp_2r => NULL()
complex(mk), dimension(:),   pointer           :: wp_1c => NULL()
complex(mk), dimension(:,:), pointer           :: wp_2c => NULL()
logical, dimension(:),   pointer               :: wp_1l => NULL()
REAL(mk), DIMENSION(:,:), POINTER :: minsub
      !!! Mimimum of extension of subs.
      !!! Used when decomp is user defined
      !!! 
      !!! 1st index: x,y,(z)                                                   
      !!! 2nd: subID
REAL(mk), DIMENSION(:,:), POINTER :: maxsub
      !!! maximum of extension of subs
      !!! Used when decomp is user defined
      !!!
      !!! 1st index: x,y,(z)                                                   
      !!! 2nd: subID
INTEGER,  DIMENSION(:  ), POINTER :: sub2proc


integer, dimension(:),allocatable              :: degree,order
real(ppm_kind_double),dimension(:),allocatable :: coeffs
integer                                        :: nterms
real(mk) :: side_len,run_sum,x,y
integer                         :: nredest,heuristic,nsubs,id,nsubs_x
logical                         :: lredecomp
TYPE(ppm_t_topo),POINTER               :: topo => NULL()


    init

        use ppm_module_init
        use ppm_module_mktopo
        use ppm_module_loadbal
        use ppm_module_particles_typedef
                
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
       
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic
        
#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = int(log10(epsilon(1._mk)))+10
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
!        do i=1,seedsize
!            seed(i)=10+i*i*(rank+1)
!        enddo
        call random_seed(put=seed)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
!        decomp = ppm_param_decomp_user_defined
        assig  = ppm_param_assign_internal
!        assig  = ppm_param_assign_user_defined
        
        !----------------
        ! Define the subs
        ! run the test with 4,16 procs
        !----------------
        nsubs = 4*nproc
        allocate(minsub(ndim,nsubs),maxsub(ndim,nsubs),cost(nsubs),sub2proc(nsubs),stat=info)
        
        !allocate(minsub(ndim,nsubs),maxsub(ndim,nsubs),cost(nsubs),stat=info)

        cost = 1._mk
        if (info.ne.0) print*, 'problem with allocate'
        nsubs_x = INT(sqrt(REAL(nsubs)))
        
        side_len = len_phys(1)/REAL(nsubs_x)
        run_sum = 0._mk
        id = 1
        x = min_phys(1)
        y = min_phys(2)
        do i=1,nsubs_x
           do j=1,nsubs_x
              
              minsub(1,id) = x 
              minsub(2,id) = y 

              maxsub(1,id) = x + side_len
              maxsub(2,id) = y + side_len
              
!              print*,rank,minsub(1:2,id),maxsub(1:2,id)
              y = y + side_len
              id = id + 1
           enddo
           y = min_phys(1)
           x = x + side_len
        enddo   
             
        do i=1,nsubs
           sub2proc(i) = MOD(i,nproc)
!           print*,'sub2proc:',sub2proc(i),'i:',i
        enddo
        topoid = 0
        !if (rank.eq.0)
        !print*,'nsubs:',nsubs
!        print*,'fine till here'
        
       call  ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
!        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info, &
!&           minsub, maxsub,nsubs,sub2proc)

        topo    => ppm_topo(topoid)%t
!        print*,rank,'s neighbors ',topo%ineighproc(1:2)

    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys)

    end finalize


    setup
    end setup
        

    teardown
    end teardown
    test balancing_circuit
        use ppm_module_util_dbg

        type(ppm_t_particles_d)               :: Part1
        class(ppm_t_discr_data),POINTER       :: Prop1=>NULL(),Prop2=>NULL()
        type(ppm_t_field)                     :: Field1
        INTEGER,DIMENSION(:),POINTER          :: nvlist => NULL()
        INTEGER,DIMENSION(:,:),POINTER        :: vlist => NULL()
        class(ppm_t_neighlist_d_),POINTER     :: Nlist => NULL()
        integer                               :: np_local
        REAL(mk), DIMENSION(:  ), POINTER     :: cost
        !!! Aggregate cost for each subdomain
        REAL(mk), DIMENSION(:  ), POINTER     :: proc_cost
        integer                               :: i, nsteps,globalID
        real(mk),DIMENSION(2)                 :: dist
    
        t_comp = 0._mk
        t_comm = 0._mk

        start_time =  MPI_Wtime()
        ! ----------------------------
        ! Initialize the particles
        ! ----------------------------
        call Part1%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)
        Assert_True(associated(Part1%xp))

!        print*,'AFTER INIT'
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        
!        call ppm_dbg_print(topoid,Part1%ghostlayer,0,1,info,xp,Part1%Npart,Part1%Mpart)
        
        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)
        
!        print*,'AFTER GET/SET'
        
        end_time = MPI_Wtime()
        t_comp   = end_time-start_time
        start_time   = MPI_Wtime()
        ! ----------------------------
        ! Map the particles
        ! ----------------------------
        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)
!        print*,'GLOBAL MAPPING ---- DONE'
        call Part1%map_ghosts(info)
        Assert_Equal(info,0)
!        print*,'GHOST MAPPING ---- DONE'
        
        end_time = MPI_Wtime()
        t_comm   = end_time-start_time
        start_time   = MPI_Wtime()
        
        ! ----------------------------
        ! Create initial neighbor list
        ! ----------------------------
        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)
!        print*,'DEBUG =0='        
        Nlist => Part1%get_neighlist(Part1)
        Assert_true(associated(Nlist))
        
        Assert_True(Part1%has_neighlist())
        call Part1%get_vlist(nvlist,vlist,info,Nlist)
        Assert_true(associated(vlist))
        Assert_true(associated(nvlist))
        Assert_Equal(info,0)
        
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        
        end_time = MPI_Wtime()
        t_comp   = t_comp + end_time-start_time
        t_total  = t_comp + t_comm
!         print*,'ARRIVED HERE'
        ! ----------------------------
        ! Estimate the initial cost
        ! ----------------------------
        call ppm_loadbal_bc(topoid,-1,Part1,t_total, &
&                           info,vlist,nvlist,t_comp=t_comp)
        Assert_Equal(info,0)
        ! ----------------------------
        ! Set up a velocity field and 
        ! a scalar test function on 
        ! the particles
        ! ----------------------------
        call Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,&
                lda=2,name="velocity")
            Assert_Equal(info,0)
        call Part1%create_prop(info,discr_data=Prop2,dtype=ppm_type_int,&
                lda=1,name="subID")
            Assert_Equal(info,0)

        call Part1%get(Prop2,wp_1r,info)
                        Assert_Equal(info,0)
        wp_1r = -5
        DO k = 1, topo%nsublist
            globalID = topo%isublist(k)
!               !---------------------------------------------------------------
                !  and check which particles are inside the sub
                !---------------------------------------------------------------
            DO j=1,Part1%Npart
                IF (xp(1,j).GE.topo%min_subd(1,globalID).AND. &
 &                  xp(1,j).LE.topo%max_subd(1,globalID).AND. &
 &                  xp(2,j).GE.topo%min_subd(2,globalID).AND. &
 &                  xp(2,j).LE.topo%max_subd(2,globalID)) THEN

                    wp_1r(j) = k !globalID
!                        print*,rank,'minsub',topo%min_subd(1,k)
  
                ENDIF ! if the particle is inside this sub
            ENDDO
!                PRINT*,rank,' SUB-ID',topo%isublist(k)
        ENDDO
        call ppm_vtk_particles('to_the_corner',Part1,info,step=0)
        Assert_Equal(info,0)
        ! ----------------------------
        ! Start the time loop
        ! ---------------------------
        nsteps = 1
        DO i=1,nsteps
            t_comp = 0._mk
            t_comm = 0._mk
            start_time = MPI_Wtime()
            call Part1%set_xp(xp,info,read_only=.true.)
            Assert_Equal(info,0)
             call Part1%get_xp(xp,info)
            Assert_Equal(info,0)
            
            call Part1%get(Prop1,wp_2r,info)
            Assert_Equal(info,0)
            
            DO ip=1,Part1%Npart
                if (xp(1,ip).GT.1._mk .or. xp(2,ip).GT.1_mk) then
                    dist = 0._mk
                else
                    dist(1) = 1 - xp(1,ip)
                    dist(2) = 1 - xp(2,ip)
                endif
                wp_2r(1:ndim,ip) = dist/nsteps!COS((10._MK*xp(1:ndim,ip))**2)
!                wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
            ENDDO
!            wp_2r = cos(wp_2r) * Part1%ghostlayer
            wp_2r  = wp_2r * Part1%ghostlayer
            ! ----------------------------
            ! Move the particles with this 
            ! displacement field
            ! ----------------------------
            call Part1%move(wp_2r,info)
            Assert_Equal(info,0)
            ! ----------------------------
            ! Apply boundary conditions and 
            ! remap the particles
            ! ----------------------------
            call Part1%apply_bc(info)
            Assert_Equal(info,0)
            
!            t_comp = end_time - start_time
            call Part1%map(info)
            Assert_Equal(info,0)
            
            ! ----------------------------
            ! Get the new ghosts
            ! ----------------------------
            call Part1%map_ghosts(info)
            Assert_Equal(info,0)  
            
            call Part1%get(Prop2,wp_1r,info)
                        Assert_Equal(info,0)
            wp_1r = -5
            call Part1%get_xp(xp,info)
            Assert_Equal(info,0)
            DO k = 1, topo%nsublist
                globalID = topo%isublist(k)
                !---------------------------------------------------------------
                !  and check which particles are inside the sub
                !---------------------------------------------------------------
                DO j=1,Part1%Npart
                    IF (xp(1,j).GE.topo%min_subd(1,globalID).AND. &
     &                  xp(1,j).LE.topo%max_subd(1,globalID).AND. &
     &                  xp(2,j).GE.topo%min_subd(2,globalID).AND. &
     &                  xp(2,j).LE.topo%max_subd(2,globalID)) THEN
    
                        wp_1r(j) = k !globalID
    !                        print*,rank,'minsub',topo%min_subd(1,k)
      
                    ENDIF ! if the particle is inside this sub
                ENDDO
!                PRINT*,rank,' SUB-ID',topo%isublist(k)
            ENDDO
            call ppm_vtk_particles('to_the_corner',Part1,info,step=i)
            Assert_Equal(info,0)
            ! ----------------------------
            ! Update (re-create) neighbor list 
            ! ----------------------------
            start_time   = MPI_Wtime()
            call Part1%comp_neighlist(info)
            Assert_Equal(info,0)
            Nlist => Part1%get_neighlist(Part1)
            Assert_true(associated(Nlist))
            
            Assert_True(Part1%has_neighlist())
            call Part1%get_vlist(nvlist,vlist,info,Nlist)
            Assert_true(associated(vlist))
            Assert_true(associated(nvlist))
            Assert_Equal(info,0)
    
            ! ----------------------------
            ! Estimate the updated cost
            ! ----------------------------
            call Part1%get_xp(xp,info)
            Assert_Equal(info,0)
            end_time   = MPI_Wtime()
            t_comp     = t_comp + end_time - start_time
            
            
            end_time = MPI_Wtime()
            t_total = end_time - start_time
            call ppm_loadbal_bc(topoid,-1,Part1,t_total, &
&                           info,vlist,nvlist)
            Assert_Equal(info,0)
            call Part1%set_xp(xp,info,read_only=.true.)
            Assert_Equal(info,0)
            
            print*,rank,' required this many time:',t_total
            print*,'%%%%%%%%%%%%%%%%%%%%%%%'
!            ! ----------------------------
!            ! Now, move the particles back
!            ! ----------------------------       
!            call Part1%get(Prop1,wp_2r,info)
!            wp_2r = -wp_2r
!            call Part1%move(wp_2r,info)
!            Assert_Equal(info,0)
!    
!            ! ----------------------------
!            ! Re-apply boundary conditions 
!            ! & remap the particles
!            ! ----------------------------
!            call Part1%apply_bc(info)
!            Assert_Equal(info,0)
!            call Part1%map(info)
!            Assert_Equal(info,0)
!    
!            call Part1%map_ghosts(info)
!            Assert_Equal(info,0)
!            
!            ! ----------------------------
!            ! Update (re-create) neighbor list 
!            ! ----------------------------            
!            call Part1%comp_neighlist(info)
!            Assert_Equal(info,0)
!            Assert_True(Part1%has_neighlist())
!
!            ! ----------------------------
!            ! Compare values and check that
!            ! they are still the same 
!            ! ----------------------------
!            call Part1%get_xp(xp,info)
!            Assert_Equal(info,0)
!            call Part1%get(Prop2,wp_1r,info)
!            Assert_Equal(info,0)
!            err = 0._mk
!            DO ip=1,Part1%Npart
!                err = max(err,abs(wp_1r(ip) - f0_test(xp(1:ndim,ip),ndim)))
!            ENDDO
!            Assert_Equal_Within(err,0,tol)
!            
!            ! ----------------------------
!            ! Update (re-create) neighbor list 
!            ! ----------------------------
!            call Part1%comp_neighlist(info)
!            Assert_Equal(info,0)
!            Nlist => Part1%get_neighlist(Part1)
!            Assert_true(associated(Nlist))
!            
!            Assert_True(Part1%has_neighlist())
!            call Part1%get_vlist(nvlist,vlist,info,Nlist)
!            Assert_true(associated(vlist))
!            Assert_true(associated(nvlist))
!            Assert_Equal(info,0)
!    
!            ! ----------------------------
!            ! Estimate the updated cost
!            ! ----------------------------
!            call Part1%get_xp(xp,info)
!            Assert_Equal(info,0)
!            call ppm_loadbal_bc(topoid,-1,xp,Part1%Npart,t_total, &
!&                           info,vlist,nvlist,t_comp=t_comp,t_comm=t_comm)
!            Assert_Equal(info,0)
!            call Part1%set_xp(xp,info,read_only=.true.)
!            Assert_Equal(info,0)
!            
       enddo

    
    end test
!    test random_motion
!        use ppm_module_util_dbg
!
!        ! test initialization of particles on a grid
!        type(ppm_t_particles_d)               :: Part1
!        class(ppm_t_discr_data),POINTER       :: Prop1=>NULL(),Prop2=>NULL()
!        character(len=32)                     :: fname
!        integer                               :: nsteps,step,npart_this,i
!        real(mk)                              :: random_num,t_comm_old,t_comp_old
!        integer,dimension(:),pointer          :: colortag
!        
!        
!
!        
!        
!        call ppm_dbg_print(topoid,0.0_mk,0,1,info)
!        print*,'ppm_dbg_print'
!        call Part1%initialize(np_global,info,topoid=topoid)
!        Assert_Equal(info,0)
!        print*,'Part1 initialize'
!        call Part1%get_xp(xp,info)
!        Assert_Equal(info,0)
!        print*,'Get Xp'
!        ! Global mapping
!        call Part1%map(info,global=.true.,topoid=topoid)
!        Assert_Equal(info,0)
!        print*,rank,'Global mapping'
!        ! Ghost mapping
!        call Part1%map_ghosts(info)
!        Assert_Equal(info,0)
!        print*,rank,'Ghost mapping'
!        
!        mat_time  = 0._MK
!        lredecomp = .FALSE.
!        nredest   = -1
!        heuristic = ppm_param_loadbal_dlb
!        nsteps    = 4
!        
!        ! initialize the current iteration's timings
!        t_comp = 0._MK
!        t_comm = 0._MK
!
!        ! A sample time loop
!        DO step=1,nsteps
!
!            start_time = 0._MK
!            end_time = 0._MK
!            ! store old time step's timings
!            t_comm_old = t_comm
!            t_comp_old = t_comp
!          
!            start_time = MPI_Wtime()
!            !Set up a velocity field and a scalar test function on the particles
!            call Part1%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,&
!                lda=2,name="velocity")
!            Assert_Equal(info,0)
!            call Part1%get(Prop1,wp_2r,info)
!            Assert_Equal(info,0)
!          
!          !call Part1%get_xp(xp,info)
!          !Assert_Equal(info,0)
! 
!          ! RANDOM MOTION
!          DO ip=1,Part1%Npart
!             call random_number(random_num) 
!             wp_2r(1:ndim,ip) = random_num* Part1%ghostlayer
!          ENDDO
!
!          ! CONSTANT MOTION
          !DO ip=1,Part1%Npart

          !        wp_2r(1:2,ip) = 0.0

          !ENDDO
          
!          !wp_2r = cos(wp_2r) * Part1%ghostlayer
!          
!          !Move the particles with this displacement field
!          call Part1%move(wp_2r,info)
!          Assert_Equal(info,0)
!
!          
!          
!          ! Computation time is being calculated
!          t_comp = MPI_Wtime()
!          t_comp = t_comp - start_time 
!          
!          
!          !write(*,'(2(A,F12.4),I3)')'last time:',elapsed_t,' mat_time:',mat_time,rank
!          ! call ppm_dbg_print(topoid,0.0_mk,step,1,info,xp,Part1%Npart)
!               
!          if (step.gt.2) then
!            print*,mat_time,Part1%Npart
!            call ppm_loadbal_inquire(t_comp_old,t_comm_old,step,Part1%Npart,.FALSE.,&
!&              lredecomp,nredest,info,heuristic,mat_time,topoid)
!            Assert_Equal(info,0)
!            print*,rank,'loadbal inquire is fine!!!'
!            if (lredecomp)   call ppm_loadbal_do_dlb(topoid,Part1,info)
!            Assert_Equal(info,0)
!
!            if (lredecomp) mat_time = 0._MK
!          endif
!          call Part1%set_xp(xp,info)
!          Assert_Equal(info,0)
!          start_time = MPI_Wtime()   
!          !Apply boundary conditions and remap the particles
!          call Part1%apply_bc(info)
!          Assert_Equal(info,0)       
!          call Part1%map(info)
!          Assert_Equal(info,0)
!          
!          !Get the new ghosts
!          call Part1%map_ghosts(info)
!          Assert_Equal(info,0)
!          
!          end_time = MPI_Wtime()
!          t_comm = end_time-start_time
!
!          mat_time = mat_time*(step-1)
!          mat_time = (mat_time + t_comp_old + t_comm_old) / step
!          npart_this = Part1%Npart
!         ! call Part1%get_xp(xp,info)
!         ! Assert_Equal(info,0)
!
!         
!          print*,'***********************'
!          write(*,*) 'step:',step,' do dlb:',lredecomp
!          print*,'***********************'
!
!          !write(*,*) 'fine till here'
!
!        ENDDO
!        !call Part1%print_info(info)
!        !Assert_Equal(info,0)
!        call Part1%destroy(info)
!        Assert_Equal(info,0)
!    end test

!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function f0_test(pos,ndim)

    real(mk)                              :: f0_test
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    f0_test =  sin(2._mk*pi*pos(1)) * cos(2._mk*pi*pos(2)) * &
        & sin(2._mk*pi*pos(ndim))

    return

end function f0_test

end test_suite
