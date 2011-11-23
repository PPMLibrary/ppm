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
integer                         :: np = 10
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
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)
        call random_number(randnb)
        
        allocate(xp(ndim,np),rcp(np),wp(pdim,np),stat=info)

    end setup
        

    teardown
        
        deallocate(xp,rcp,wp,stat=info)
        deallocate(seed,randnb)

    end teardown

    test ghost_get_add_generic
        ! test ghost get followed by adding particles

        use ppm_module_typedef
        use ppm_module_data
        use ppm_module_mktopo
        use ppm_module_topo_check
        real(mk), parameter             :: gl = 0.1_mk
        real(mk),dimension(:,:),pointer :: xpn => NULL()
        real(mk),dimension(:,:),pointer :: wpn => NULL()
        real(mk),dimension(:,:),pointer :: rand_num => NULL()
        real(mk),dimension(ndim)        :: disp
        real(mk),dimension(:),pointer   :: rcpn => NULL()
        integer                         :: np_added=0,mp_added=0
        integer                         :: mp_new=0,np_new=0
        integer                         :: i,j,k,ipart,npart
        integer                         :: isub,iproc
        integer                         :: nb_samp = 4
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
            rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo

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

        np_added = nb_samp*topo%nsublist

        allocate(xpn(ndim,np_added),wpn(pdim,np_added),rcpn(np_added))
        allocate(rand_num(ndim,nb_samp*topo%nsubs))

        !set a displacement vector, smaller than the ghost layer size
        disp(1) = gl*0.5_mk
        disp(2) = gl*0.3_mk
        if (ndim .eq. 3) disp(ndim) = gl*0.1_mk
        call random_seed(put=seed)
        call random_number(rand_num)

        do k=1,topo%nsublist !for each subdomain on that proc
            isub = topo%isublist(k)
            do j=1,nb_samp !generate a particle at each corner
                ipart = nb_samp*(k-1)+j
                do i=1,ndim
                    xpn(i,ipart)= topo%min_subd(i,isub) 
                enddo
                SELECT CASE(j)
                CASE (1)
                CASE (2)
                    xpn(1,ipart)= topo%max_subd(1,isub) 
                CASE (3)
                    xpn(2,ipart)= topo%max_subd(2,isub) 
                CASE (4)
                    xpn(ndim,ipart)= topo%max_subd(ndim,isub) 
                CASE (5)
                    xpn(1,ipart)= topo%max_subd(1,isub) 
                    xpn(2,ipart)= topo%max_subd(2,isub) 
                CASE (6)
                    xpn(1,ipart)= topo%max_subd(1,isub) 
                    xpn(ndim,ipart)= topo%max_subd(ndim,isub) 
                CASE (7)
                    xpn(2,ipart)= topo%max_subd(2,isub) 
                    xpn(ndim,ipart)= topo%max_subd(ndim,isub) 
                CASE (8)
                    xpn(1,ipart)= topo%max_subd(1,isub) 
                    xpn(2,ipart)= topo%max_subd(2,isub) 
                    xpn(ndim,ipart)= topo%max_subd(ndim,isub) 
                CASE DEFAULT !random point inside the subdomain
                    xpn(1:ndim,ipart) = rand_num(1:ndim,(isub-1)*nb_samp*+j)
                    do i=1,ndim
                        xpn(i,ipart) = &
                        topo%min_subd(i,isub) + xpn(i,ipart) / &
                        (topo%max_subd(i,isub) - topo%min_subd(i,isub))
                    enddo
                END SELECT
                !move the particle in one direction
                xpn(1,ipart) = xpn(1,ipart) + disp(1)
                xpn(2,ipart) = xpn(2,ipart) + disp(2)
                if (ndim.eq.3) xpn(ndim,ipart) = xpn(ndim,ipart) + gl*disp(ndim)
                wpn(1,ipart)= real(j,mk)
                wpn(2,ipart)= 7._mk*real(j,mk)
                rcpn(ipart)= 0.33_mk + real(j,mk)
            enddo
        enddo

        write(*,*) '[',rank,'] ','Adding ',np_added,' particles'
        do i=1,np_added
            write(100+rank,'(2(E14.7,2X))') xpn(1:2,i)
        enddo


        call ppm_part_split_compute(topoid,xpn,np_added,0,gl,info)
        assert_true(info.eq.0)
        write(*,*) 'Nr = ',modify%Nrnew
        write(*,*) 'Ng = ',modify%Ngnew

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

        write(*,*) 'np_added = ',np_added
        write(*,*) 'mp_added = ',mp_added
        write(*,*) 'np_new = ',np_new
        write(*,*) 'mp_new = ',mp_new

        write(*,*) '*********'
        write(*,*) 'Nr = ',modify%Nrnew
        write(*,*) 'Ng = ',modify%Ngnew

        
        call ppm_part_modify_add(topoid,xpn,ndim,npart,np_added, &
            ppm_param_add_real_particles,0,gl,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(wpn,pdim,np_added,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(rcpn,np_added,info)
        assert_true(info.eq.0)
        call ppm_part_modify_send(np_added,mp_added,info)
        assert_true(info.eq.0)
        np_new = npart + np_added
        mp_new = mp    + mp_added
        write(*,*) 'npart = ',npart
        write(*,*) 'np_added = ',np_added
        write(*,*) 'mp_added = ',mp_added
        write(*,*) 'np_new = ',np_new
        write(*,*) 'mp_new = ',mp_new
        write(*,*) 'AAAA = ',np_new,2*np_new-npart,npart+1,np_new
        call MPI_BARRIER(comm,info)
        write(*,*) 'ok, lets carry on'

        call ppm_part_modify_pop(rcp,rcpn,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(wp,wpn,pdim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(xp,xpn,ndim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)

        mp = mp_new

        !modify properties
        do i=1,np_new
          rcp(i) = 3.14_mk
          wp(1,i) = 7.14_mk
          wp(2,i) = 8.14_mk
        enddo

        !update ghosts without calling ghost_get
        call ppm_map_part_push(rcp,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_push(wp,pdim,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_send(np_new,mp,info)
        assert_true(info.eq.0)
        assert_true(np_new.eq.npart+np_added)
        call ppm_map_part_pop(wp,pdim,np_new,mp,info)
        assert_true(info.eq.0)
        call ppm_map_part_pop(rcp,np_new,mp,info)
        assert_true(info.eq.0)
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
        if (topo%sub2proc(k).eq.topo%icommseq(iproc)) then
            np_added = np_added + nb_samp
        endif
        enddo
        enddo

        allocate(xpn(ndim,np_added),wpn(pdim,np_added),rcpn(np_added))
        allocate(rand_num(ndim,nb_samp*topo%nsubs))

        !set a displacement vector, smaller than the ghost layer size
        disp(1) = gl*0.5_mk
        disp(2) = gl*0.3_mk
        if (ndim .eq. 3) disp(ndim) = gl*0.1_mk
        call random_seed(put=seed)
        call random_number(rand_num)

        ipart = 0
        do iproc=1,topo%ncommseq ! for each neighboring proc
            do k=1,topo%nsubs
               !loop through subdomains that belong to that proc
               if (topo%sub2proc(k).eq.topo%icommseq(iproc)) then
                ipart = ipart+1
                do i=1,ndim
                    xpn(i,ipart)= topo%min_subd(i,k) 
                enddo
                SELECT CASE(j)
                CASE (1)
                CASE (2)
                    xpn(1,ipart)= topo%max_subd(1,k) 
                CASE (3)
                    xpn(2,ipart)= topo%max_subd(2,k) 
                CASE (4)
                    xpn(ndim,ipart)= topo%max_subd(ndim,k) 
                CASE (5)
                    xpn(1,ipart)= topo%max_subd(1,k)
                    xpn(2,ipart)= topo%max_subd(2,k)
                CASE (6)
                    xpn(1,ipart)= topo%max_subd(1,k)
                    xpn(ndim,ipart)= topo%max_subd(ndim,k)
                CASE (7)
                    xpn(2,ipart)= topo%max_subd(2,k)
                    xpn(ndim,ipart)= topo%max_subd(ndim,k)
                CASE (8)
                    xpn(1,ipart)= topo%max_subd(1,k)
                    xpn(2,ipart)= topo%max_subd(2,k)
                    xpn(ndim,ipart)= topo%max_subd(ndim,k)
                CASE DEFAULT !random point inside the subdomain
                    xpn(1:ndim,ipart) = rand_num(1:ndim,(k-1)*nb_samp*+j)
                    do i=1,ndim
                        xpn(i,ipart) = &
                        topo%min_subd(i,k) + xpn(i,ipart) / &
                        (topo%max_subd(i,k) - topo%min_subd(i,k))
                    enddo
                END SELECT
                !move the particle in one direction
                xpn(1,ipart) = xpn(1,ipart) + disp(1)
                xpn(2,ipart) = xpn(2,ipart) + disp(2)
                if (ndim.eq.3) xpn(ndim,ipart) = xpn(ndim,ipart) + gl*disp(ndim)
                wpn(1,ipart)= real(j,mk)
                wpn(2,ipart)= 7._mk*real(j,mk)
                rcpn(ipart)= 0.33_mk + real(j,mk)
              endif
            enddo
          enddo

        do i=1,np_added
            write(200+rank,'(2(E14.7,2X))') xpn(1:2,i)
        enddo
        do i=1,np_new
            write(300+rank,'(2(E14.7,2X))') xp(1:2,i)
        enddo
        !Now loop through all these particles, and check whether they should be on
        !that proc, either as ghost or as real particles.
        !Check if that is the case. 
        do ipart = 1,np_added
            isreal = .false.; isghost = .false.
            do k=1,topo%nsublist !for each subdomain on the current proc
                isub = topo%isublist(k)
                !find the real particles
                if (is_inside_sub(xpn(1:ndim,ipart),isub,ndim,topo)) then
                    isreal = .true.
                endif
            enddo
            do k=1,topo%nsublist !for each subdomain on the current proc
                isub = topo%isublist(k)
                !find the ghosts
                if (is_ghostof_sub(xpn(1:ndim,ipart),isub,ndim,0,gl,topo)) then
                    if (.not. isreal) isghost = .true.
                endif
            enddo
            if (isreal) then
                !check that this particle is in xp(1:Np)
                ok = belongs_to(xpn(1:ndim,ipart),xp(1:ndim,1:np_new),ndim,np_new)
                assert_true(ok)
            endif
            if (isghost) then
                !check that this particle is in xp(Np+1,Mp)
                ok = belongs_to(xpn(1:ndim,ipart),&
                    xp(1:ndim,np_new+1:mp_new),ndim,mp_new-np_new+1)
                assert_true(ok)
            endif
        enddo

        deallocate(xpn,wpn,rcpn)
        deallocate(rand_num)

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
            sub_gl(1:ndim,1)=topo%min_subd(1:ndim,isub)
        else
            sub_gl(1:ndim,1)=topo%min_subd(1:ndim,isub)-gl
        endif
        sub_gl(1:ndim,2)=topo%max_subd(1:ndim,isub)+gl

        sub(1:ndim,1)=topo%min_subd(1:ndim,isub)
        sub(1:ndim,2)=topo%max_subd(1:ndim,isub)

        if (ndim.eq.2) then
            if (xp(1).GE.sub_gl(1,1) .AND. xp(1).LT.sub_gl(1,2) .AND. &
                xp(2).GE.sub_gl(2,1) .AND. xp(2).LT.sub_gl(2,2)) then
                if (xp(1).LT.sub(1,1) .AND. xp(1).GE.sub(1,2) .AND. &
                    xp(2).LT.sub(2,1) .AND. xp(2).GE.sub(2,2)) then
                        is_ghostof_sub = .true.
                        return
                endif
            endif
        endif
        if (ndim.eq.3) then
            if (xp(1).GE.sub_gl(1,1) .AND. xp(1).LT.sub_gl(1,2) .AND. &
                xp(2).GE.sub_gl(2,1) .AND. xp(2).LT.sub_gl(2,2) .AND. &
                xp(3).GE.sub_gl(3,1) .AND. xp(3).LT.sub_gl(3,2)) then
                if (xp(1).LT.sub(1,1) .AND. xp(1).GE.sub(1,2) .AND. &
                    xp(2).LT.sub(2,1) .AND. xp(2).GE.sub(2,2) .AND. &
                    xp(3).LT.sub(3,1) .AND. xp(3).GE.sub(3,2)) then
                        is_ghostof_sub = .true.
                        return
                endif
            endif
        endif

    end function is_ghostof_sub

    function belongs_to(xp,xp_0,ndim,np)

        use ppm_module_typedef
        use ppm_module_data

        use ppm_module_mktopo
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

end test_suite

