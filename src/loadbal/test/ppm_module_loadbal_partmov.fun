test_suite ppm_module_loadbal_partmov

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
integer                         :: decomp,assig,tolexp,temp
integer                         :: info,comm,rank,nproc,topoid
integer                         :: np_global = 140000
real(mk),parameter              :: cutoff = 0.1_mk
real(mk),dimension(:,:),pointer :: xp=>NULL(),xp_temp=>NULL()
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
real(mk)                        :: t_comp,t_comm,t_total,random_num
real(mk)                        :: one_real,two_reals

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
real(mk) :: side_len,run_sum,x,y,t_total_noDLB
integer                         :: nredest,heuristic,nsubs,id,nsubs_x
logical                         :: lredecomp
TYPE(ppm_t_topo),POINTER               :: topo => NULL()
character(len=19)             :: output_buf
character(len=17)             :: output_dlb


    init

        use ppm_module_init
        use ppm_module_mktopo
        use ppm_module_loadbal
        use ppm_module_particles_typedef
                
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
       
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        if (nproc.eq.7) max_phys(1:ndim) = 1.4_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6)       = ppm_param_bcdef_periodic
        
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
!        decomp = ppm_param_decomp_cuboid
        decomp = ppm_param_decomp_user_defined
!        assig  = ppm_param_assign_internal
        assig  = ppm_param_assign_user_defined
        
        !----------------
        ! Define the subs
        ! run the test with 4,16 procs
        !----------------
        IF (nproc.eq.2) then
            nsubs = 16*nproc
        elseif (nproc.eq.4) then
            nsubs = 16*nproc
        elseif (nproc.eq.5) then        
            nsubs = 20*nproc
        elseif (nproc.eq.7) then
            nsubs = 196
        endif
        allocate(minsub(ndim,nsubs),maxsub(ndim,nsubs),cost(nsubs),sub2proc(nsubs),stat=info)
        cost = 1._mk
        if (info.ne.0) print*, 'problem with allocate'
        cost = 1._mk
        if (info.ne.0) print*, 'problem with allocate'
        nsubs_x = floor(sqrt(REAL(nsubs)))
        if (nproc.eq.7) nsubs_x=14
        side_len = len_phys(1)/REAL(nsubs_x)
!        if (nproc.eq.7) side_len=0.1_mk
        print*,'nsubs_x',nsubs_x,'side_len',side_len,'nsubs',nsubs
        
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
        IF (nproc.eq.2) then
            sub2proc(1:16) = 0
            sub2proc(17:32)= 1
        elseif (nproc.eq.4) then
            do i=1,nsubs
             if (i.LE.16) then
                sub2proc(i) = 0
             elseif (i.LE.32) then
                sub2proc(i) = 1
             elseif (i.LE.48) then
                sub2proc(i) = 2
             else
                sub2proc(i) = 3
             endif
!            print*,'sub2proc:',sub2proc(i),'i:',i
            enddo
        elseif (nproc.eq.5) then        
            do i=1,nsubs
             if (i.LE.20) then
                sub2proc(i) = 0
             elseif (i.LE.40) then
                sub2proc(i) = 1
             elseif (i.LE.60) then
                sub2proc(i) = 2
             elseif (i .le.80) then
                sub2proc(i) = 3
             else
                sub2proc(i) = 4
             endif
!            print*,'sub2proc:',sub2proc(i),'i:',i
            enddo
        elseif (nproc.eq.7) then        
            
                sub2proc(1:28) = 0
             
                sub2proc(29:56) = 1
             
                sub2proc(57:84) = 2
             
                sub2proc(85:112) = 3
             
                sub2proc(113:140) = 4
             
                sub2proc(141:168) = 5
             
                sub2proc(169:196) = 6
             
        endif    
        
        topoid = 0
!        print*,rank,sub2proc
        ! ----------------------------
        ! Create topology using particles
        ! ----------------------------
!        call  ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info, &
&           minsub, maxsub,nsubs,sub2proc)
        topo    => ppm_topo(topoid)%t
        
        ! ----------------------------
        !  Particles' positions
        ! ----------------------------
!        allocate(xp_temp(2,np_global),stat=info)
!        do i=1,np_global
!            call random_number(random_num) 
!            xp_temp(1,i) = random_num*0.25
!            xp_temp(2,i) = random_num*0.25
!        enddo
    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)
        print*,rank,'OK till here'
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
        class(ppm_t_neighlist_d_),POINTER     :: Nlist => NULL()
        integer                               :: np_local
        integer                               :: i, nsteps,globalID
        real(mk),DIMENSION(2)                 :: dist
        logical                               :: isDLB
        character(len=5)                      :: vtkfilename
        character(len=1)                      :: procbuf
    
        isDLB = .TRUE.
        write(procbuf,'(I1)')nproc
        IF (isDLB) THEN
            write(vtkfilename,'(A1,A4)')  trim(procbuf),'wDLB'
        ELSE
            write(vtkfilename,'(A1,A4)')  trim(procbuf),'nDLB'
        ENDIF 
        t_comp = 0._mk
        t_comm = 0._mk

        start_time =  MPI_Wtime()
        ! ----------------------------
        ! Initialize the particles
        ! ----------------------------
        call Part1%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)
        Assert_True(associated(Part1%xp))
!        !print*,'after init'
        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)
!        !print*,'after set_xp'
        call Part1%get_xp(xp,info)
        Assert_Equal(info,0)
        Assert_True(associated(xp))
        
        
        end_time = MPI_Wtime()
        t_comp   = end_time-start_time
        start_time   = MPI_Wtime()
               
        ! ----------------------------
        ! Map the particles
        ! ----------------------------
        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)
        call Part1%map_ghosts(info)
        Assert_Equal(info,0)
        
        end_time = MPI_Wtime()
        t_comm   = end_time-start_time
        start_time   = MPI_Wtime()
        print*,rank,'subs',topo%nsubs,'topo%nsublist',topo%nsublist
        !---------------------------
        ! Create initial neighbor list
        ! ----------------------------
        print*,rank,'# particles',Part1%Npart
        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)
                
        Nlist => Part1%get_neighlist(Part1)
        Assert_true(associated(Nlist))
        print*,Nlist%uptodate
        Assert_True(Part1%has_neighlist())
        call Part1%get_vlist(nvlist,vlist,info,Nlist)
        Assert_true(associated(vlist))
        Assert_true(associated(nvlist))
        Assert_Equal(info,0)

!        print*,rank,'isublist',topo%isublist
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
!        !print*,rank,'size of wp_1r',size(wp_1r,1)
        print*,'topo%nsublist',topo%nsublist
!        DO k = 1, topo%nsublist
!            globalID = topo%isublist(k)
                !---------------------------------------------------------------
                !  and check which particles are inside the sub
                !---------------------------------------------------------------
            DO j=1,Part1%Npart
!                IF (xp(1,j).GE.topo%min_subd(1,globalID).AND. &
! &                  xp(1,j).LE.topo%max_subd(1,globalID).AND. &
! &                  xp(2,j).GE.topo%min_subd(2,globalID).AND. &
! &                  xp(2,j).LE.topo%max_subd(2,globalID)) THEN

                    wp_1r(j) = rank 
!                ENDIF ! if the particle is inside this sub
            ENDDO
!        ENDDO
!        !print*, 'arrived here too'
        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)
        call ppm_vtk_particles(vtkfilename,Part1,info,step=0)
        Assert_Equal(info,0)
!       ! print*,'after vtk - before 1st DLB!!!'
        ! ----------------------------
        ! Estimate the initial cost
        ! ----------------------------
        end_time = MPI_Wtime()
        t_comp   = t_comp + end_time-start_time
        t_total  = t_comp + t_comm
        
!        
!        call Part1%map_ghosts(info)
!        Assert_Equal(info,0)
        
!        call ppm_vtk_particles(vtkfilename,Part1,info,step=1)
!        Assert_Equal(info,0)

        ! ----------------------------
        ! Start the time loop
        ! ---------------------------
        nsteps = 40
        DO i=1,nsteps
            print*,'%%%%%%%%%%%%%%%%%%%%% - ',i
            t_comp = 0._mk
            t_comm = 0._mk
            t_total= 0._mk
            start_time = MPI_Wtime()
            
            ! ALL GET CALLS
            call Part1%get_xp(xp,info)
            Assert_Equal(info,0)
!            Assert_true(associated(Nlist))
            call Part1%get(Prop1,wp_2r,info)
            Assert_Equal(info,0)
            
            call Part1%get(Prop2,wp_1r,info)
            Assert_Equal(info,0)
!            print*,rank,'OK',i
            ! Move particles
!            print*,rank,'wp_1r size:',size(wp_1r),'XP size:',size(xp(1,:)),'Npart',Part1%Npart
            
            DO ip=1,Part1%Npart
                if (xp(1,ip).GT.1._mk-Part1%ghostlayer) then
                    dist(1) = 0._mk
                    dist(2) = 0._mk !0.5_mk  - xp(2,ip)
                else
                    dist(1) = max_phys(1) - xp(1,ip)
                    dist(2) = 0._mk !0.5_mk  - xp(2,ip)
                endif
                wp_2r(1:ndim,ip) = dist/1.1_mk!COS((10._MK*xp(1:ndim,ip))**2)
!                !wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
            ENDDO

            wp_2r  = wp_2r * Part1%ghostlayer
            ! ----------------------------
            ! Move the particles with this 
            ! displacement field
            ! ----------------------------
            call Part1%move(wp_2r,info)
            Assert_Equal(info,0)
            ! ----------------------------
            ! Apply boundary conditions and 
            ! remap the particles and get new ghosts
            ! ----------------------------
            call Part1%apply_bc(info)
            Assert_Equal(info,0)
            call Part1%map(info)
            Assert_Equal(info,0)
            call Part1%map_ghosts(info)
            Assert_Equal(info,0)
            ! ----------------------------
            ! Update neighbor lists
            ! ----------------------------
            nvlist => NULL()
            vlist  => NULL()
            Nlist => NULL()
!            print*,rank,'OK222',i 
            call Part1%comp_neighlist(info)
            Assert_Equal(info,0)
!            print*,rank,'OK333',i
            Nlist => Part1%get_neighlist(Part1)
            Assert_true(associated(Nlist))
!            print*,rank,'OK444',i
!            Assert_True(Part1%has_neighlist())
            call Part1%get_vlist(nvlist,vlist,info,Nlist)
            Assert_true(associated(vlist))
            Assert_true(associated(nvlist))
            
    
            call Part1%get_xp(xp,info)
            call Part1%get(Prop1,wp_2r,info)
            Assert_Equal(info,0)
            
            call Part1%get(Prop2,wp_1r,info)
            Assert_Equal(info,0)
            Assert_Equal(info,0)            
            wp_1r = -5
!            print*,rank,'OK222',i
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
    
                        wp_1r(j) = rank !globalID
    !                        print*,rank,'minsub',topo%min_subd(1,k)
      
                    ENDIF ! if the particle is inside this sub
                ENDDO
            ENDDO
!            print*,rank,'OK555',i
            DO k = 1,Part1%Npart,2
                DO j=1,Part1%Npart
                    one_real  = cos(sin(cos(cos (1._mk))))
                    two_reals = cos(sin(cos(sqrt(1._mk))))
                ENDDO
            ENDDO
!            call Part1%set_xp(xp,info,read_only=.true.)
!            Assert_Equal(info,0)
!            print*,rank,'OK666',i
!            call Part1%set_xp(xp,info,read_only=.true.)
!            Assert_Equal(info,0)
            
!            print*,rank,'====wp_1r size:',size(wp_1r),'XP size:',size(xp(1,:)),'Npart',Part1%Npart
!            call ppm_vtk_particles(vtkfilename,Part1,info,step=i+2)
!            Assert_Equal(info,0)
!            print*,rank,'OK777',i
            ! ----------------------------
            ! Update (re-create) neighbor list 
            ! ----------------------------
!            !start_time   = MPI_Wtime()
            
!            Assert_Equal(info,0)
            
            ! ----------------------------
            ! Estimate the updated cost
            ! ----------------------------
!            
            
            end_time = MPI_Wtime()
            
            t_total_noDLB =  end_time - start_time
            
            start_time = MPI_Wtime()
!            print*,rank,'t_total:',t_total
            !!!!! DLB !!!!!
            
            IF (isDLB) THEN
                call ppm_loadbal_bc(topoid,-1,Part1,t_total_noDLB, &
&                           info,vlist,nvlist)
                Assert_Equal(info,0)
!            ! ----------------------------
!            ! Get the new ghosts
!            ! ----------------------------
!                call Part1%map_ghosts(info)
!                Assert_Equal(info,0)
!                print*,rank,'after ghosts after loadbal'
            ENDIF
!            print*,rank,'before get_xp:Npart',Part1%Npart
            end_time = MPI_Wtime()
            t_total = t_total_noDLB + (end_time-start_time)
            IF (.NOT. isDLB) THEN
                write(output_buf,'(A1,A,I1,A)') trim(procbuf),'output_NODLB_',rank,'.txt'
                open (unit=rank, file=output_buf, status='unknown',access='append')
                write(rank,'(I6,F10.6)'), Part1%Npart,t_total
                close(rank)
            ELSE
                  
                write(output_dlb,'(A1,A,I1,A)') trim(procbuf),'output_dlb_',rank,'.txt'
                open (unit=rank, file=output_dlb, status='unknown',access='append')
                write(rank,'(I6,F10.6,F10.6)'), Part1%Npart,t_total,t_total_noDLB
                close(rank)
            ENDIF
            
            call Part1%get_xp(xp,info)
            Assert_Equal(info,0)
!            print*,rank,'before wp_1r'
            call Part1%get(Prop2,wp_1r,info)
            Assert_Equal(info,0)
            wp_1r = -5
!            print*,rank,'topo%nsublist',topo%nsublist
            DO k = 1, topo%nsublist
                globalID = topo%isublist(k)
!                !print*,rank,'globalID',globalID
                !---------------------------------------------------------------
                !  and check which particles are inside the sub
                !---------------------------------------------------------------
                DO j=1,Part1%Npart
                    IF (xp(1,j).GE.topo%min_subd(1,globalID).AND. &
     &                  xp(1,j).LE.topo%max_subd(1,globalID).AND. &
     &                  xp(2,j).GE.topo%min_subd(2,globalID).AND. &
     &                  xp(2,j).LE.topo%max_subd(2,globalID)) THEN
    
                        wp_1r(j) = rank !globalID
    !                        print*,rank,'minsub',topo%min_subd(1,k)
      
                    ENDIF ! if the particle is inside this sub
                ENDDO
!                !PRINT*,rank,' SUB-ID',topo%isublist(k)
            ENDDO
!            !print*,rank,'before set_xp'
            call Part1%set_xp(xp,info,read_only=.true.)
            Assert_Equal(info,0)            
!            !print*,rank,'after set_xp'
            call ppm_vtk_particles(vtkfilename,Part1,info,step=i)
            Assert_Equal(info,0)
!            
            
            end_time = MPI_Wtime()
            t_total = t_total + (end_time-start_time)
            
            print*,rank,' required this much time:',t_total
            ! Number of particles and time spent
!            IF (.NOT. isDLB) THEN
!                write(output_buf,'(A,I1,A)') 'output_NODLB_',rank,'.txt'
!                open (unit=rank, file=output_buf, status='unknown',access='append')
!                write(rank,'(I6,F10.6)'), Part1%Npart,t_total
!                close(rank)
!            ELSE
!                  
!                write(output_dlb,'(A,I1,A)') 'output_dlb_',rank,'.txt'
!                open (unit=rank, file=output_dlb, status='unknown',access='append')
!                write(rank,'(I6,F10.6)'), Part1%Npart,t_total
!                close(rank)
!            ENDIF
!            
            
            print*,rank,'%%%%ENDS%%%%%%%%%%%%%%%%%%'         
       enddo
!
!    
    end test

end test_suite
