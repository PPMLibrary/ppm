test_suite ppm_module_part_modify



#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer,parameter               :: pdim=2
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol,min_rcp,max_rcp
integer                         :: info,comm,rank,nproc
integer                         :: topoid
integer                         :: np = 200
integer                         :: mp
integer                         :: newnp
real(mk),dimension(:,:),pointer :: xp => NULL()
real(mk),dimension(:  ),pointer :: rcp => NULL()
real(mk),dimension(:,:),pointer :: wp => NULL()
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()
real(mk),dimension(:  ),pointer :: h => NULL()
real(mk),dimension(:  ),pointer :: p_h => NULL()
real(mk),dimension(:  ),pointer :: len_phys => NULL()
real(mk),dimension(:  ),pointer :: ghostlayer => NULL()
integer, dimension(:  ),pointer :: ghostsize => NULL()
integer                         :: i,j,k,sum1,sum2
integer                         :: p_i
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
integer, dimension(:  ),pointer :: nm => NULL()
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok
real(mk)                         :: t0,t1,t2,t3

    init

        use ppm_module_typedef
        use ppm_module_init
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostsize(ndim),ghostlayer(2*ndim),&
            &         nm(ndim),h(ndim),p_h(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        ghostlayer(1:2*ndim) = max_rcp
        bcdef(1:6) = ppm_param_bcdef_periodic
        tol = epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))
        
        nullify(xp,rcp,wp)

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys,ghostsize,nm)

    end finalize


    setup

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        allocate(randnb((1+ndim)*np),stat=info)
        do i=1,seedsize
            seed(i)=11+i*i
        enddo
        call random_seed(put=seed)
        call random_number(randnb)
        
        allocate(xp(ndim,np),rcp(np),wp(pdim,np),stat=info)

    end setup
        

    teardown
        
        deallocate(xp,rcp,wp,stat=info)
        deallocate(seed,randnb)

    end teardown

    test ghost_get_add
        ! test ghost get followed by adding particles

        use ppm_module_typedef
        use ppm_module_data
        use ppm_module_data_buffers
        use ppm_module_mktopo
        use ppm_module_topo_check
        real(mk), parameter             :: gl = 0.1_mk
        real(mk),dimension(:,:),pointer :: xpn => NULL()
        real(mk),dimension(:,:),pointer :: wpn => NULL()
        real(mk),dimension(:,:),allocatable :: rand_num
        real(mk),dimension(:),pointer   :: rcpn => NULL()
        integer                         :: np_added=0,mp_added=0
        integer                         :: mp_new=0,np_new=0
        integer                         :: i,j,k,ipart,npart
        integer                         :: isub,iproc
        integer                         :: nb_samp = 200
        integer                         :: dummy = -HUGE(1)
        integer                         :: Nr,Ng
        type(ppm_t_topo),pointer        :: topo => NULL()
        logical                         :: isreal,isghost,ok

        !----------------
        ! create particles
        !----------------

        npart = np
        xp = 0.0_mk
        rcp = 0.0_mk
        min_rcp = 0.01_mk
        max_rcp = 0.2_mk 

        do i=1,npart
            do j=1,ndim
                xp(j,i) = min_phys(j)+&
                len_phys(j)*randnb((ndim+1)*i-(ndim-j))
            enddo
            rcp(i) = -1._mk
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo
!        xp(1,1) = 0.25_mk
!        xp(1,2) = 0.75_mk
!        xp(2,1:2) = 0.25_mk + 0.5_mk*rank


        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_rcp,cost,info)
        topo => ppm_topo(topoid)%t

        call ppm_map_part_global(topoid,xp,npart,info)
        call ppm_map_part_push(rcp,npart,info)
        call ppm_map_part_push(wp,pdim,npart,info)
        call ppm_map_part_send(npart,newnp,info)
        call ppm_map_part_pop(wp,pdim,npart,newnp,info)
        call ppm_map_part_pop(rcp,npart,newnp,info)
        call ppm_map_part_pop(xp,ndim,npart,newnp,info)
        npart=newnp
        call ppm_topo_check(topoid,xp,npart,ok,info)
        assert_true(ok)

        call ppm_map_part_ghost_get(topoid,xp,ndim,npart,0,gl,info)
        call ppm_map_part_push(rcp,npart,info)
        call ppm_map_part_push(wp,pdim,npart,info)
        call ppm_map_part_send(npart,mp,info)
        call ppm_map_part_pop(wp,pdim,npart,mp,info)
        call ppm_map_part_pop(rcp,npart,mp,info)
        call ppm_map_part_pop(xp,ndim,npart,mp,info)
        
        call ppm_topo_check(topoid,xp,npart,ok,info)
        assert_true(ok)

!        do i=1,npart
!            write(10+rank,'(2(E14.7,2X))') xp(1:2,i)
!        enddo
!        do i=npart+1,mp
!            write(20+rank,'(2(E14.7,2X))') xp(1:2,i)
!        enddo

        np_added = nb_samp*topo%nsublist

        allocate(xpn(ndim,np_added),wpn(pdim,np_added),rcpn(np_added),STAT=info)
        assert_true(info.eq.0)
        allocate(rand_num(ndim,nb_samp*topo%nsubs),stat=info)
        assert_true(info.eq.0)

!        if (rank .eq. 0 ) then
!            write(*,*) '[--------------]'
!            do k=1,topo%nsubs
!                write(*,*) ' subdomain ',k
!                write(*,'(F7.4,A,F7.4)') topo%min_subd(1,k), ' | ',topo%max_subd(1,k)
!                write(*,'(F7.4,A,F7.4)') topo%min_subd(2,k), ' | ',topo%max_subd(2,k)
!            enddo
!            write(*,*) '[--------------]'
!        endif

        call random_seed(put=seed)
        call random_number(rand_num)

        do k=1,topo%nsublist !for each subdomain on that proc
            isub = topo%isublist(k)
            do j=1,nb_samp !generate a particle at each corner
                ipart = nb_samp*(k-1)+j
                xpn(1:ndim,ipart)=sample_xp(ndim,topo,gl,isub,j,rand_num,nb_samp)
                wpn(1,ipart)=       real(j,mk)
                wpn(2,ipart)= 7._mk*real(j,mk)
                rcpn(ipart) = 2._mk*real(j,mk)
            enddo
        enddo

!        do i=1,np_added
!            write(100+rank,'(2(E14.7,2X))') xpn(1:2,i)
!        enddo

        call ppm_part_split_compute(topoid,xpn,np_added,0,gl,info)
        assert_true(info.eq.0)
!        do i=1,modify%Nrnew
!            write(160+rank,'(2(E14.7,2X))') xpn(1:2,modify%idx_real_new(i))
!        enddo
!        do i=1,modify%Ngnew
!            write(170+rank,'(2(E14.7,2X))') xpn(1:2,modify%idx_ghost_new(i))
!        enddo
!        write(*,'(A,I0,A,I0,A,I0,1X,I0)') &
!            '[',rank,'] Adding ',np_added,' particles ',modify%Nrnew,modify%Ngnew

        call ppm_part_modify_add(topoid,xpn,ndim,npart,np_added, &
            ppm_param_add_ghost_particles,0,gl,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(wpn,pdim,np_added,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(rcpn,np_added,info)
        assert_true(info.eq.0)
        call ppm_part_modify_send(np_added,mp_added,info)
        assert_true(info.eq.0)

        np_new = npart + np_added
        mp_new = mp + mp_added

        call ppm_part_modify_pop(rcp,rcpn,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(wp,wpn,pdim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(xp,xpn,ndim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)

!        do i=1,modify%Nrnew
!            write(110+rank,'(2(E14.7,2X))') xpn(1:2,modify%idx_real_new(i))
!        enddo
!        write(*,'(A,I0,A,I0,A,3(I0,1X))') &
!            '[',rank,'] Np_a ',np_added,' Mp_a ',mp_added,modify%Nrnew,modify%Ngnew

        call ppm_part_modify_add(topoid,xpn,ndim,npart,np_added, &
            ppm_param_add_real_particles,0,gl,info)
        assert_true(info.eq.0)

        call ppm_part_modify_push(wpn,pdim,np_added,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(rcpn,np_added,info)
        assert_true(info.eq.0)
        call ppm_part_modify_send(np_added,mp_added,info)
        assert_true(info.eq.0)

        np_new = npart + modify%Nrnew
        mp_new = mp    + mp_added - np_added + modify%Nrnew

        call ppm_part_modify_pop(rcp,rcpn,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(wp,wpn,pdim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(xp,xpn,ndim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)

        mp = mp_new

!        do i=1,np_new
!            write(300+rank,'(3(E14.7,2X))') xp(1:2,i),rcp(i)
!        enddo
!        do i=np_new+1,mp
!            write(400+rank,'(3(E14.7,2X))') xp(1:2,i),rcp(i)
!        enddo


        !modify properties
        rcp = 0
        wp = 0
        do i=1,np_new
          rcp(i)  =  3.14_MK
          wp(1,i) = -4._mk
          wp(2,i) = -8._mk
        enddo

        !update ghosts without calling ghost_get
        call ppm_map_part_push(rcp,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_push(wp,pdim,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_send(np_new,mp,info)
        assert_true(info.eq.0)
        call ppm_map_part_pop(wp,pdim,np_new,mp,info)
        assert_true(info.eq.0)
        call ppm_map_part_pop(rcp,np_new,mp,info)
        assert_true(info.eq.0)

!        do i=1,np_new
!            write(310+rank,'(3(E14.7,2X))') xp(1:2,i),rcp(i)
!        enddo
!        do i=np_new+1,mp
!            write(410+rank,'(3(E14.7,2X))') xp(1:2,i),rcp(i)
!        enddo

        assert_true(rcp(np_new+1).eq.3.14_mk)
        assert_true(rcp(mp).eq.3.14_mk)

        deallocate(xpn,wpn,rcpn)

        !Run the check

        ! create the same particles than the ones we added for the test
        ! This time, we do it globally (for all subdomains, not only
        ! those of the current processor).
        np_added = 0
        do iproc=1,topo%ncommseq
        do k=1,topo%nsubs
        if (topo%sub2proc(k).eq.topo%icommseq(iproc) .or. &
            topo%sub2proc(k) .eq. rank) then
            np_added = np_added + nb_samp
        endif
        enddo
        enddo


        allocate(xpn(ndim,np_added),wpn(pdim,np_added),rcpn(np_added))

        ipart = 0
        do iproc=1,topo%ncommseq ! for each neighboring proc
            do k=1,topo%nsubs
            !loop through subdomains that belong to that proc
                if (topo%sub2proc(k).eq.topo%icommseq(iproc) .or. &
                    topo%sub2proc(k) .eq. rank) then
                    do j=1,nb_samp
                        ipart = ipart+1

                        xpn(1:ndim,ipart)=sample_xp(ndim,topo,gl,k,j,rand_num,nb_samp)
                        wpn(1,ipart)= real(j,mk)
                        wpn(2,ipart)= real(k,mk)
                        rcpn(ipart)= 0.33_mk + real(j,mk)


                        !check whether this particle should be on this proc, 
                        !either as ghost or as real particles.
                        isreal = .false.; isghost = .false.
                        do i=1,topo%nsublist !for each subdomain on the current proc
                            isub = topo%isublist(i)
                            !find the real particles
                            if (is_inside_sub(xpn(1:ndim,ipart),isub,ndim,topo)) then
                                isreal = .true.
                            endif
                        enddo
                        do i=1,topo%nsublist !for each subdomain on the current proc
                            isub = topo%isublist(i)
                            !find the ghosts
                            if (is_ghostof_sub(xpn(1:ndim,ipart),&
                                isub,ndim,0,gl,topo)) then
                                if (.not. isreal) isghost = .true.
                            endif
                        enddo
                        if (isreal) then
                            !check that this particle is in xp(1:Np)
                            ok = belongs_to(xpn(1:ndim,ipart),&
                                xp(1:ndim,1:np_new),ndim,np_new)
                            if (.not.ok) then
                                write(*,'(A,I0,A,I0,1X,I0,F7.4,1X,F7.4)') &
                                    '[',rank,'] real particle ',j,ipart,&
                                    xpn(1:ndim,ipart)
                                write(*,'(A,I0,A,F7.4,1X,F7.4)') &
                                    '[',rank,'] min_sudomain ', topo%min_subd(1:2,k)
                                write(*,'(A,I0,A,F7.4,1X,F7.4)') &
                                    '[',rank,'] max_sudomain ', topo%max_subd(1:2,k)
                            endif
                            assert_true(ok)
                        endif

                        if (isghost) then
                            !check that this particle is in xp(Np+1,Mp)
                            ok = belongs_to(xpn(1:ndim,ipart),&
                                xp(1:ndim,np_new+1:mp_new),ndim,mp_new-np_new+1)
                            if (.not.ok) then
                                write(*,'(A,I0,A,I0,1X,I0,F7.4,1X,F7.4)') &
                                    '[',rank,'] ghost particle ',j,ipart,xpn(1:2,ipart)
                                write(*,'(A,I0,A,F7.4,1X,F7.4)') &
                                    '[',rank,'] min_sudomain ', topo%min_subd(1:2,k)
                                write(*,'(A,I0,A,F7.4,1X,F7.4)') &
                                    '[',rank,'] max_sudomain ', topo%max_subd(1:2,k)
                                ok = belongs_to(xpn(1:ndim,ipart),&
                                    xp(1:ndim,1:np_new),ndim,np_new)
                                if (ok) then
                                    write(*,*) &
                                        'but this particle is found as a real particle'
                                endif
                                ok = .false.
                            endif
                            assert_true(ok)
                        endif
                    enddo
                endif
            enddo
        enddo


!        do i=1,np_added
!            write(200+rank,'(2(E14.7,2X))') xpn(1:2,i)
!        enddo
!        call mpi_barrier(comm,info)


        deallocate(xpn,wpn,rcpn)
        deallocate(rand_num)

!        write(*,*) 'stopping now'
!        stop

    end test

    function is_inside_sub(xp,isub,ndim,topo)

        use ppm_module_typedef
        use ppm_module_data

        use ppm_module_mktopo
        real(mk),dimension(ndim)        :: xp
        integer                         :: isub
        integer                         :: ndim
        real(mk),dimension(ndim,2)      :: sub
        type(ppm_t_topo),pointer        :: topo => NULL()
        logical                         :: is_inside_sub

        sub(1:ndim,1)=topo%min_subd(1:ndim,isub)
        sub(1:ndim,2)=topo%max_subd(1:ndim,isub)

        is_inside_sub = .false.
        if (ndim.eq.2) then
            if (xp(1).GE.sub(1,1) .AND. xp(1).LT.sub(1,2) .AND. &
                xp(2).GE.sub(2,1) .AND. xp(2).LT.sub(2,2)) then
                    is_inside_sub = .true.
                    return
            endif
        endif
        if (ndim.eq.3) then
            if (xp(1).GE.sub(1,1) .AND. xp(1).LT.sub(1,2) .AND. &
                xp(2).GE.sub(2,1) .AND. xp(2).LT.sub(2,2) .AND. &
                xp(3).GE.sub(3,1) .AND. xp(3).LT.sub(3,2)) then
                    is_inside_sub = .true.
                    return
            endif
        endif

    end function is_inside_sub

    function is_ghostof_sub(xp,isub,ndim,isymm,gl,topo)

        use ppm_module_typedef
        use ppm_module_data

        use ppm_module_mktopo
        real(mk),dimension(ndim)        :: xp
        integer                         :: isub
        integer                         :: ndim
        integer                         :: isymm
        real(mk)                        :: gl
        real(mk),dimension(ndim,2)      :: sub_gl
        real(mk),dimension(ndim,2)      :: sub
        type(ppm_t_topo),pointer        :: topo => NULL()
        logical                         :: is_ghostof_sub

        is_ghostof_sub = .false.
        if (isymm .eq. 0 ) then
            sub_gl(1:ndim,1)=topo%min_subd(1:ndim,isub)-gl
        else
            sub_gl(1:ndim,1)=topo%min_subd(1:ndim,isub)
        endif
        sub_gl(1:ndim,2)=topo%max_subd(1:ndim,isub)+gl

        sub(1:ndim,1)=topo%min_subd(1:ndim,isub)
        sub(1:ndim,2)=topo%max_subd(1:ndim,isub)

        if (ndim.eq.2) then
            if (xp(1).GE.sub_gl(1,1) .AND. xp(1).LT.sub_gl(1,2) .AND. &
                xp(2).GE.sub_gl(2,1) .AND. xp(2).LT.sub_gl(2,2)) then
                if (xp(1).LT.sub(1,1) .OR. xp(1).GE.sub(1,2) .AND. &
                    xp(2).LT.sub(2,1) .OR. xp(2).GE.sub(2,2)) then
                        is_ghostof_sub = .true.
                        return
                endif
            endif
        endif
        if (ndim.eq.3) then
            if (xp(1).GE.sub_gl(1,1) .AND. xp(1).LT.sub_gl(1,2) .AND. &
                xp(2).GE.sub_gl(2,1) .AND. xp(2).LT.sub_gl(2,2) .AND. &
                xp(3).GE.sub_gl(3,1) .AND. xp(3).LT.sub_gl(3,2)) then
                if (xp(1).LT.sub(1,1) .OR. xp(1).GE.sub(1,2) .AND. &
                    xp(2).LT.sub(2,1) .OR. xp(2).GE.sub(2,2) .AND. &
                    xp(3).LT.sub(3,1) .OR. xp(3).GE.sub(3,2)) then
                        is_ghostof_sub = .true.
                        return
                endif
            endif
        endif

    end function is_ghostof_sub

    function belongs_to(xp,xp_0,ndim,np)

        use ppm_module_typedef
        use ppm_module_data

        real(mk),dimension(ndim)        :: xp
        real(mk),dimension(ndim,np)     :: xp_0
        integer                         :: ndim
        integer                         :: np
        logical                         :: belongs_to
        integer                         :: ip

        belongs_to = .false.
        do ip = 1,np
           if (maxval(abs(xp-xp_0(1:ndim,ip))) .lt. ppm_myepsd) then
               belongs_to = .true.
               return
           endif
        enddo

    end function belongs_to

    function sample_xp(ndim,topo,gl,isub,j,rand_num,nb_samp)
        !generate particles in a given subdomain
        ! either in some specific locations relative to its
        ! corners or borders, or randomly in or near the subdomain.
        ! Useful to generate test cases.
        use ppm_module_typedef
        use ppm_module_data
        use ppm_module_mktopo
        real(mk),dimension(ndim)        :: sample_xp
        integer                         :: isub,j,ndim,nb_samp
        real(mk)                        :: gl
        real(mk),dimension(:,:)         :: rand_num
        type(ppm_t_topo),pointer        :: topo => NULL()
        real(mk)                        :: disp1,disp2,disp3


        disp1 = 0.9_mk * gl
        disp2 = 0.1_mk * gl
        disp3 = -0.1_mk * gl

        do i=1,ndim
            sample_xp(i)= topo%min_subd(i,isub) 
        enddo

        SELECT CASE((j-1)/6+1)
        CASE (1)
        CASE (2)
            sample_xp(1)= topo%max_subd(1,isub)
        CASE (3)
            sample_xp(2)= topo%max_subd(2,isub) 
        CASE (4)
            sample_xp(1)= topo%max_subd(1,isub) 
            sample_xp(2)= topo%max_subd(2,isub) 
        CASE (5)
            sample_xp(ndim)= topo%max_subd(ndim,isub) 
        CASE (6)
            sample_xp(1)= topo%max_subd(1,isub) 
            sample_xp(ndim)= topo%max_subd(ndim,isub) 
        CASE (7)
            sample_xp(2)= topo%max_subd(2,isub) 
            sample_xp(ndim)= topo%max_subd(ndim,isub) 
        CASE (8)
            sample_xp(1)= topo%max_subd(1,isub) 
            sample_xp(2)= topo%max_subd(2,isub) 
            sample_xp(ndim)= topo%max_subd(ndim,isub) 
        CASE DEFAULT !random point inside the subdomain
            sample_xp(1:ndim) = rand_num(1:ndim,(isub-1)*nb_samp + j)
            sample_xp(1:ndim) = &
                topo%min_subd(1:ndim,isub) -gl + sample_xp(1:ndim) * &
                (2._mk*gl + topo%max_subd(1:ndim,isub) - &
                topo%min_subd(1:ndim,isub))
        END SELECT

        !move the particle in one direction (only for nonrandom cases)
        if (j.le.(6*8)) then 
          SELECT CASE(MOD(j,6))
          CASE(0)
          CASE(1)
            sample_xp(1) = sample_xp(1) + disp1
          CASE(2)
            sample_xp(2) = sample_xp(2) + disp2
          CASE(3)
            sample_xp(2) = sample_xp(2) + disp1
          CASE(4)
            sample_xp(2) = sample_xp(2) + disp1
            IF (ndim.eq.3) THEN
                sample_xp(3) = sample_xp(3) + disp3
            ENDIF
          CASE(5)
            sample_xp(1) = sample_xp(1) + disp1
            IF (ndim.eq.3) THEN
                sample_xp(3) = sample_xp(3) + disp1
            ENDIF
          END SELECT
        endif

    end function sample_xp

end test_suite

