      SUBROUTINE DTYPE(vbp_create)(Pc,Npart,info,name)
          !!! create a set of particles
          !!! This allocates the particle positions.
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_vbp))                   :: Pc
          !!! Data structure containing the particles
          INTEGER,                    INTENT(IN   ) :: Npart
          !!! Number of particles
          INTEGER,                    INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: name
          !!! give a name to this Particle set
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: i

          LOGICAL :: lalloc,ldealloc

          start_subroutine("vbp_create")

          !Call the parent function
          CALL Pc%DTYPE(ppm_t_particles)%create(Npart,info,name)
          or_fail("failed to initialize non-vbp particle set")

          !and update the few fields that are specific to VBP
          Pc%adaptive = .FALSE.

          check_false(<#ASSOCIATED(Pc%rcp)#>,&
          & "The rcp property (cutoff radii) is already defined for that particle set. Use destroy() before create()")

          Pc%rcp => NULL()

          end_subroutine()
      END SUBROUTINE DTYPE(vbp_create)

      SUBROUTINE DTYPE(vbp_destroy)(Pc,info)
          !!! destroy a set of particles
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          DEFINE_MK()
          CLASS(DTYPE(ppm_t_vbp))                :: Pc
          !!! Data structure containing the particles
          INTEGER,                INTENT(  OUT)  :: info
          !!! Returns status, 0 upon success.

          start_subroutine("vbp_destroy")

          !re-initialize the few fields that are specific to VBP
          Pc%adaptive = .FALSE.
          !yaser
          Pc%rcp => NULL()

          !Call the parent function
          CALL Pc%DTYPE(ppm_t_particles)%destroy(info)
          or_fail("failed to destroy vbp particle set")

          end_subroutine()
      END SUBROUTINE DTYPE(vbp_destroy)

      SUBROUTINE DTYPE(vbp_set_cutoff)(Pc,cutoff,info,Nlist)
          !!! Set a constant cutoff radius for a Particle set with varying sizes
          !!! and update the ghostlayer sizes.
          !!! The cutoff radius concerns the default neighbor list, unless
          !!! specified otherwise.
          !!! If the cutoff is increased from its previous value, the neighbour
          !!! list is flagged as "not up-to-date" and will have to be recomputed
          !!! before it is used again.
          !!! If the ghostlayer sizes are increased from their previous values,
          !!! the ghosts are flagged as "not up-to-date" and will have to be
          !!! recomputed before they are used again.
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_vbp))                  :: Pc
          REAL(MK),                 INTENT(IN   )  :: cutoff
          !!! cutoff radius (same number of elements as we have particles)
          INTEGER,                  INTENT(   OUT) :: info
          !!! return status. On success, 0
          CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,INTENT(INOUT) :: NList
          !!! Neighbor list for which this cutoff radius
          !!! applies. By default, this is the "standard" Verlet list, with neighbours
          !!! sought within the particle set itself.

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
!           CLASS(DTYPE(ppm_t_neighlist)_),   POINTER :: nl => NULL()
          REAL(MK) :: max_cutoff
!           REAL(MK),DIMENSION(:),POINTER       :: rcp => NULL()
          REAL(MK), DIMENSION(:), ALLOCATABLE :: cutoff_v

          start_subroutine("vbp_set_cutoff")

          IF (ASSOCIATED(Pc%rcp)) THEN
              !Varying-blob particle set already has a varying cutoff radius. Let us
              !overwrite it with the new value.
              ! NOTE: a better thing to do would be to delete it and use the constant
              ! cutoff instead. This is doable but would require every routine to
              ! check which of the two cutoffs are defined (the scalar or the vector
              ! variable) and use the right one.
              ALLOCATE(cutoff_v(1:Pc%Npart),STAT=info)
              or_fail_alloc("cutoff_v")

              cutoff_v = cutoff

              CALL Pc%set_varying_cutoff(cutoff_v,info,Nlist)
              or_fail("Failed to set constant cutoff radius")

              Pc%adaptive = .TRUE.

              DEALLOCATE(cutoff_v,STAT=info)
              or_fail_dealloc("cutoff_v")
          ELSE
              CALL Pc%DTYPE(ppm_t_particles)%set_cutoff(cutoff,info,Nlist)
              or_fail("Failed to set constant cutoff radius")

              Pc%adaptive = .FALSE.
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(vbp_set_cutoff)


      SUBROUTINE DTYPE(vbp_set_varying_cutoff)(Pc,cutoff,info,Nlist)
          !!! Set a cutoff radius for a Particle set with varying sizees
          !!! and update the ghostlayer sizes.
          !!! The cutoff radius concerns the default neighbor list, unless
          !!! specified otherwise.
          !!! If the cutoff is increased from its previous value, the neighbour
          !!! list is flagged as "not up-to-date" and will have to be recomputed
          !!! before it is used again.
          !!! If the ghostlayer sizes are increased from their previous values,
          !!! the ghosts are flagged as "not up-to-date" and will have to be
          !!! recomputed before they are used again.
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          ! Arguments
          !-------------------------------------------------------------------------
          DEFINE_MK()

          CLASS(DTYPE(ppm_t_vbp))                                 :: Pc

          REAL(MK), DIMENSION(:),                   INTENT(IN   ) :: cutoff
          !!! cutoff radius (same number of elements as we have particles)

          INTEGER,                                  INTENT(  OUT) :: info
          !!! return status. On success, 0

          CLASS(DTYPE(ppm_t_neighlist)_), OPTIONAL, INTENT(INOUT) :: NList
          !!! Neighbor list for which this cutoff radius
          !!! applies. By default, this is the "standard" Verlet list, with neighbours
          !!! sought within the particle set itself.
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: Nl

          REAL(MK), DIMENSION(:), POINTER :: rcp
          REAL(MK)                        :: max_cutoff

          INTEGER :: ip,datatype

          start_subroutine("vbp_set_varying_cutoff")

          !-------------------------------------------------------------------------
          !  Set new cutoff
          !-------------------------------------------------------------------------
          IF (.NOT.ASSOCIATED(Pc%rcp)) THEN
             datatype=MERGE(ppm_type_real,ppm_type_real_single,MK.EQ.ppm_kind_double)

             CALL Pc%create_prop(info,part_prop=Pc%rcp,dtype=datatype,name='rcp')
             or_fail("could not create property for varying cutoff radius rcp")
          ENDIF

          NULLIFY(rcp)

          CALL Pc%get(Pc%rcp,rcp,info)
          or_fail("could not access varying cutoff radius rcp")

          DO ip=1,Pc%Npart
             rcp(ip) = cutoff(ip)
          ENDDO

          max_cutoff = MAXVAL(rcp(1:Pc%Npart))

          IF (PRESENT(NList)) THEN
             check_true(<#Pc%neighs%has(NList)#>,&
             "Neighbour list does not concern this particle set")
             IF (max_cutoff .LT. NList%cutoff) NList%uptodate = .FALSE.
             NList%cutoff = max_cutoff
          ELSE
             Nl => Pc%get_neighlist()
             check_associated(Nl,"Compute neighbour lists first")
             IF (max_cutoff .LT. Nl%cutoff) Nl%uptodate = .FALSE.
             Nl%cutoff = max_cutoff
          ENDIF

          ! Compute ghostlayer sizes
          IF (max_cutoff.GT.Pc%ghostlayer) THEN
             !If the new cutoff is larger than the current ghostsize
             ! then the new ghostsize is the new cutoff
             Pc%ghostlayer = max_cutoff
             ! update states
             Pc%flags(ppm_part_ghosts) = .FALSE.
          ELSE IF (max_cutoff .LT. Pc%ghostlayer) THEN
             !Else, we find the new maximum cutoff radius amongst
             !all existing neighbor lists on this Particle set
             Pc%ghostlayer = 0._MK
             Nl => Pc%neighs%begin()
             DO WHILE (ASSOCIATED(Nl))
                IF (Nl%cutoff .GT. Pc%ghostlayer) THEN
                   Pc%ghostlayer = Nl%cutoff
                ENDIF
                Nl => Pc%neighs%next()
             ENDDO
             !no need to update states: ghosts are still ok.
          ENDIF

          Pc%adaptive = .TRUE.

          end_subroutine()
      END SUBROUTINE DTYPE(vbp_set_varying_cutoff)

      SUBROUTINE DTYPE(vbp_neigh_create)(this,Part_src,info,&
      &          name,skin,symmetry,cutoff,Nlist)
          !!! Create a data structure to store a neighbour list
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_vbp))                                 :: this

          CLASS(DTYPE(ppm_t_particles)_),           TARGET        :: Part_src
          !!! Particle set to which the neighbours belong (can be the same as this)

          INTEGER,                                  INTENT(  OUT) :: info

          CHARACTER(LEN=*),               OPTIONAL, INTENT(IN   ) :: name
          !!! name of this neighbour list

          REAL(MK),                       OPTIONAL, INTENT(IN   ) :: skin
          REAL(MK),                       OPTIONAL, INTENT(IN   ) :: cutoff

          LOGICAL,                        OPTIONAL, INTENT(IN   ) :: symmetry

          CLASS(DTYPE(ppm_t_neighlist)_), OPTIONAL, POINTER       :: Nlist
          !!! returns a pointer to the newly created verlet list

          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: Nl

          REAL(MK), DIMENSION(:), POINTER :: rcp

          INTEGER :: vec_size,i,datatype

          start_subroutine("vbp_neigh_create")

          ! Create the neighbour list
          ALLOCATE(DTYPE(ppm_t_neighlist)::Nl,STAT=info)
          or_fail_alloc("Nl")

          IF (PRESENT(name)) THEN
             Nl%name = name
          ELSE
             WRITE(Nl%name,*) 'Nl',TRIM(ADJUSTL(this%name)),'_',&
             & TRIM(ADJUSTL(Part_src%name))
          ENDIF

          check_associated(<#Part_src%xp#>,"Invalid particle set Part_src")

          Nl%Part => Part_src

          ASSOCIATE (ghosts => this%flags(ppm_part_ghosts))
             IF (.NOT.ASSOCIATED(this%rcp)) THEN
                datatype=MERGE(ppm_type_real,ppm_type_real_single,MK.EQ.ppm_kind_double)

                CALL this%create_prop(info,part_prop=this%rcp, &
                &    dtype=datatype,name='rcp',with_ghosts=ghosts)
                or_fail("Creating property for rcp failed")
             ENDIF

             NULLIFY(rcp)
             !yaser I resolved the bug, by if conditional
             !TOCHECK
             IF (PRESENT(cutoff)) THEN
                CALL this%get(this%rcp,rcp,info,with_ghosts=ghosts)
                or_fail("Cannot access this%rcp")

                rcp = cutoff

                CALL this%set(this%rcp,rcp,info,ghosts_ok=ghosts)
                or_fail("Cannot set this%rcp pointer free")
             ELSE
                !yaser: I added read_only=.TRUE. otherwise
                !discr_data%flags(ppm_ppt_ghosts) will be set to false
                !in this case, it should not be, as the vbp cutoff values
                !are available as this%ghostlayer in each domain
                CALL this%get(this%rcp,rcp,info,with_ghosts=ghosts,read_only=.TRUE.)
                or_fail("Cannot access this%rcp")

                rcp = this%ghostlayer

                !yaser: I added read_only=.TRUE. otherwise
                !discr_data%flags(ppm_ppt_ghosts) will be set to false
                CALL this%set(this%rcp,rcp,info,ghosts_ok=ghosts,read_only=.TRUE.)
                or_fail("Cannot set this%rcp pointer free")
             ENDIF

          END ASSOCIATE

          Nl%cutoff = -1._MK
          !this field should not be used with adaptive particles

          Nl%skin = MERGE(skin,0.0_MK,PRESENT(skin))

          IF (PRESENT(symmetry)) THEN
             Nl%isymm = MERGE(1,0,symmetry)
          ELSE
             Nl%isymm = 0
          ENDIF

          Nl%uptodate = .FALSE.
          Nl%nneighmin = 0
          Nl%nneighmax = 0

          !returning a pointer to the neighbour list,
          !before it is pushed into the collection.
          IF (PRESENT(Nlist)) Nlist => Nl

          CALL this%neighs%push(Nl,info)
          or_fail("pushing new neighbour list into collection failed")

          end_subroutine()

      END SUBROUTINE DTYPE(vbp_neigh_create)


      SUBROUTINE DTYPE(vbp_neighlist)(this,info,P_xset,name, &
      &          skin,symmetry,cutoff,lstore,incl_ghosts,knn)
        !!!  Neighbor lists for particles
        !!!  Compute the Verlet lists for the target particles, using neighbours
        !!!  from the particle set P_set (the default is that P_xset is the
        !!!  same set as the target particles)
        !!!-----------------------------------------------------------------
        !!!  Assumptions:
        !!! * Particles positions need to have been mapped onto the topology
        !!! * Ghost positions have been computed
        USE ppm_module_mpi
        USE ppm_module_neighlist
#ifdef __WITH_CNL
        USE ppm_module_cnl
#endif
        USE ppm_module_inl_vlist
        USE ppm_module_inl_xset_vlist
        USE ppm_module_inl_k_vlist
        USE ppm_module_kdtree
        USE ppm_module_util_time
        IMPLICIT NONE

        DEFINE_MK()
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(ppm_t_vbp)),                  TARGET        :: this
        !!! Data structure containing the particles
        INTEGER,                                  INTENT(  OUT) :: info
        !!! Return status, on success 0.
        !-------------------------------------------------------------------------
        !  Optional arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(ppm_t_particles)_), OPTIONAL, TARGET        :: P_xset
        !!! Particle set from which the neighbours are sought
        CHARACTER(LEN=*),               OPTIONAL, INTENT(IN   ) :: name
        !!! name of this neighbour list
        REAL(MK),                       OPTIONAL, INTENT(IN   ) :: skin
        !!! skin
        LOGICAL,                        OPTIONAL, INTENT(IN   ) :: symmetry
        !!! if using symmetry
        REAL(MK),                       OPTIONAL, INTENT(IN   ) :: cutoff
        !!! cutoff radius
        LOGICAL,                        OPTIONAL, INTENT(IN   ) :: lstore
        !!! store verlet lists
        LOGICAL,                        OPTIONAL, INTENT(IN   ) :: incl_ghosts
        !!! if true, then verlet lists are computed for all particles, incl. ghosts.
        !!! Default is false.
        INTEGER,                        OPTIONAL, INTENT(IN   ) :: knn
        !!! if present, neighbour lists are constructed such that each particle
        !!! has at least knn neighbours.
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_t_topo), POINTER :: topo

        TYPE(DTYPE(kdtree)), POINTER :: tree

        TYPE(DTYPE(kdtree_result)), DIMENSION(:), ALLOCATABLE, TARGET :: results

        CLASS(ppm_t_operator_discr_),   POINTER :: op
        CLASS(DTYPE(ppm_t_particles)_), POINTER :: Part_src
        CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: Nlist

        REAL(MK), DIMENSION(:), POINTER :: rcp
        REAL(MK), DIMENSION(2*ppm_dim)  :: ghostlayer
        REAL(KIND(1.D0))                :: t1,t2
        REAL(MK)                        :: tskin

        INTEGER :: topoid
        INTEGER :: nneighmin,nneighmax
        INTEGER :: op_id,np_target,i
        INTEGER :: ip,ineigh
        !!! index variable

        LOGICAL :: lknn,lsymm
        !!! uses a neighbour-finding algorithm that finds enough neighbours
        LOGICAL :: xset_neighlists

        start_subroutine("vbp_comp_neighlist")

        !-----------------------------------------------------------------
        !  Checks
        !-----------------------------------------------------------------
        check_associated(<#this%xp#>,&
        & "Particles structure had not been defined. Call allocate first")

        check_true(<#this%flags(ppm_part_partial)#>,&
        & "Particles not mapped. Do a partial/global mapping")

        xset_neighlists = .FALSE.
        IF (PRESENT(P_xset)) THEN
           Part_src => P_xset

           check_associated(<#Part_src%xp#>,&
           & "Cross-Set particles have not been defined. Call allocate first")

           check_true(<#Part_src%flags(ppm_part_partial)#>,&
           & "Particles not mapped. Do a partial/global mapping")

           IF (.NOT.ASSOCIATED(Part_src,this)) THEN
              xset_neighlists = .TRUE.
           ENDIF
        ELSE
           Part_src => this
        ENDIF

        check_true(<#Part_src%flags(ppm_part_ghosts)#>,&
        & "Ghosts have not been updated. They are needed for neighlists")

        check_associated(<#this%neighs#>)

        !check whether the neighbour list already exists
        IF (this%has_neighlist(Part_src)) THEN
           Nlist => this%get_neighlist(Part_src)
           IF (PRESENT(skin).OR.PRESENT(symmetry).OR.PRESENT(cutoff)) THEN
              stdout("the optional arguments skin,",&
              & "symmetry or cutoff will not be used",&
              & " because the neighbour list already exists. We",&
              & " should perhaps change the API?  ")
              fail("Need to destroy/re-create this neighbour list first")
           ENDIF
           IF (Nlist%cutoff.LT.0._MK) THEN

           ENDIF
        ELSE
           NULLIFY(Nlist)
           CALL this%create_neighlist(Part_src,info,name=name,skin=skin,&
           &    symmetry=symmetry,cutoff=cutoff,Nlist=Nlist)
           or_fail("failed to create neighbour list")
        ENDIF

        check_associated(Nlist)

        !check that we have a cutoff radius
        check_associated(<#this%rcp#>,&
        & "cutoff radii for adaptive particles have not been defined")

        NULLIFY(rcp)
        CALL this%get(this%rcp,rcp,info,with_ghosts=.TRUE.,read_only=.TRUE.)
        or_fail("could not access cutoff radii")

        lsymm =Nlist%isymm.EQ.1
        lknn  =PRESENT(knn)
        tskin =Nlist%skin
        topoid=this%active_topoid

        do_something: IF (Nlist%uptodate .OR. this%Npart.EQ.0) THEN
           !neighbor lists are already up-to-date, or no particles on this proc
           !nothing to do
           IF (Nlist%uptodate) THEN
              fail('neighlists are already up-to-date, NOTHING to do', &
              & 999,exit_point=no,ppm_error=ppm_error_notice)
              info = 0
           ELSE
              Nlist%nneighmin = 0
              Nlist%nneighmax = 0
           ENDIF
        ELSE
           !hack to build (potentially incomplete) neighbour lists even
           !for ghost particles
           np_target = this%Npart
           IF (PRESENT(incl_ghosts)) THEN
              IF (incl_ghosts) THEN
                 np_target = this%Mpart
                 topo => ppm_topo(topoid)%t
                 SELECT CASE (MK)
                 CASE (ppm_kind_single)
                    topo%min_subs(:,:) = topo%min_subs(:,:) - topo%ghostsizes
                    topo%max_subs(:,:) = topo%max_subs(:,:) + topo%ghostsizes

                 CASE (ppm_kind_double)
                    topo%min_subd(:,:) = topo%min_subd(:,:) - topo%ghostsized
                    topo%max_subd(:,:) = topo%max_subd(:,:) + topo%ghostsized

                 END SELECT
              ENDIF
           ENDIF

           IF (lknn) THEN
              this%stats%nb_kdtree = this%stats%nb_kdtree+1

              CALL ppm_util_time(t1)

              ALLOCATE(tree,STAT=info)
              or_fail_alloc("tree")

              CALL tree%create(Part_src%xp(1:ppm_dim,1:Part_src%Mpart), &
              &    info,sort=.TRUE.,rearrange=.TRUE.)
              or_fail("tree%create")

              ALLOCATE(results(knn+1),STAT=info)
              or_fail_alloc("results")

              ldc(1) = knn
              ldc(2) = this%Npart
              CALL ppm_alloc(Nlist%vlist,ldc,ppm_param_alloc_grow,info)
              or_fail_alloc("Nlist%vlist")

              DO ip=1,this%Npart
                 CALL kdtree_n_nearest(tree,this%xp(1:ppm_dim,ip),&
                 &    knn+1,results,info)
                 or_fail("kdtree_n_nearest")

                 ! If the tree is not sorted you need to remove the
                 ! particle ip from the list of neighbors
                 Nlist%vlist(1:knn,ip)=results(2:knn+1)%idx

!                  ineigh=0
!                  DO i=1,knn+1
!                     IF (results(i)%idx.NE.ip) THEN
!                        ineigh=ineigh+1
!                        Nlist%vlist(ineigh,ip)=results(i)%idx
!                     ENDIF
!                  ENDDO
              ENDDO

              ldc(1) = this%Npart
              CALL ppm_alloc(Nlist%nvlist,ldc,ppm_param_alloc_grow,info)
              or_fail_alloc("Nlist%nvlist")

              Nlist%nvlist=knn

              CALL tree%destroy(info)
              or_fail("tree%destroy")

              DEALLOCATE(tree,STAT=info)
              or_fail_dealloc("tree")
              NULLIFY(tree)

              DEALLOCATE(results,STAT=info)
              or_fail_dealloc("results")

              CALL ppm_util_time(t2)
              this%stats%t_kdtree = this%stats%t_kdtree+(t2-t1)
           ELSE
              !FIXME: when adaptive ghost layers are available
              ghostlayer(1:2*ppm_dim)=Part_src%ghostlayer

#ifdef __WITH_CNL
              conventionalinl: IF (this%conventionalinl) THEN
                 this%stats%nb_cinl = this%stats%nb_cinl+1

                 CALL ppm_util_time(t1)

                 !HUGLY HACK to make CNL routines work on a topology with
                 !several subdomains
#if   __KIND == __SINGLE_PRECISION
                 CALL cnl_vlist(this%xp,                                &
                 &    rcp,this%Npart,this%Mpart,                        &
                 &    ppm_topo(topoid)%t%min_subs(:,1)-this%ghostlayer, &
                 &    ppm_topo(topoid)%t%max_subs(:,1)+this%ghostlayer, &
                 &    this%nvlist,this%vlist,ppm_dim,info)
#elif __KIND == __DOUBLE_PRECISION
                 CALL cnl_vlist(this%xp,                                &
                 &    rcp,this%Npart,this%Mpart,                        &
                 &    ppm_topo(topoid)%t%min_subd(:,1)-this%ghostlayer, &
                 &    ppm_topo(topoid)%t%max_subd(:,1)+this%ghostlayer, &
                 &    this%nvlist,this%vlist,ppm_dim,info)
#endif
                 or_fail("ppm_cinl_vlist failed")
                 !end HUGLY HACK

                 CALL ppm_util_time(t2)
                 this%stats%t_cinl = this%stats%t_cinl + (t2 - t1)
              ELSE
#endif
                 !__WITH_CNL
                 IF (xset_neighlists) THEN
                    this%stats%nb_xset_nl = this%stats%nb_xset_nl + 1

                    CALL ppm_util_time(t1)

                    CALL ppm_inl_xset_vlist(topoid,this%xp,                &
                    &    this%Npart,this%Mpart,Part_src%xp,Part_src%Npart, &
                    &    Part_src%Mpart,rcp,                               &
                    &    tskin,ghostlayer,info,Nlist%vlist,                &
                    &    Nlist%nvlist,lstore)
                    or_fail("ppm_inl_xset_vlist failed")

                    CALL ppm_util_time(t2)
                    this%stats%t_xset_nl = this%stats%t_xset_nl + (t2 - t1)
                 ELSE
                    this%stats%nb_inl = this%stats%nb_inl+1

                    CALL ppm_util_time(t1)

                    CALL ppm_inl_vlist(topoid,this%xp,np_target,     &
                    &    this%Mpart,rcp,tskin,lsymm,ghostlayer,info, &
                    &    Nlist%vlist,Nlist%nvlist)
                    or_fail("ppm_inl_vlist failed")

                    CALL ppm_util_time(t2)
                    this%stats%t_inl = this%stats%t_inl + (t2 - t1)
                 ENDIF ! XSET
#ifdef __WITH_CNL
              ENDIF conventionalinl
#endif

           ENDIF

           !restore subdomain sizes (revert hack)
           IF (PRESENT(incl_ghosts)) THEN
              IF (incl_ghosts) THEN
                 SELECT CASE (MK)
                 CASE (ppm_kind_single)
                    topo%min_subs(:,:) = topo%min_subs(:,:) + topo%ghostsizes
                    topo%max_subs(:,:) = topo%max_subs(:,:) - topo%ghostsizes

                 CASE (ppm_kind_double)
                    topo%min_subd(:,:) = topo%min_subd(:,:) + topo%ghostsized
                    topo%max_subd(:,:) = topo%max_subd(:,:) - topo%ghostsized

                 END SELECT
                 topo => NULL()
              ENDIF
           ENDIF

           !-----------------------------------------------------------------------
           !Update state
           !-----------------------------------------------------------------------
           Nlist%uptodate = .TRUE.

           Nlist%nneighmin = MINVAL(Nlist%nvlist(1:this%Npart))
           Nlist%nneighmax = MAXVAL(Nlist%nvlist(1:np_target))

           ! DC operators that do not use a xset neighbour list, if they exist,
           ! are no longer valid (they depend on the neighbour lists)
           IF (ASSOCIATED(this%ops)) THEN
              op => this%ops%begin()
              DO WHILE (ASSOCIATED(op))
                 IF (.NOT.op%flags(ppm_ops_interp)) THEN
                    op%flags(ppm_ops_iscomputed) = .FALSE.
                 ENDIF
                 op => this%ops%next()
              ENDDO
           ENDIF

           ! We want to distinguish between "self" neighbour lists
           ! and cross-set ones.
           IF (ASSOCIATED(Nlist%Part,this)) THEN
              this%flags(ppm_part_neighlists) = .TRUE.
           ENDIF

        ENDIF do_something

        end_subroutine()
        RETURN
      END SUBROUTINE DTYPE(vbp_neighlist)


      SUBROUTINE DTYPE(vbp_updated_cutoff)(this,max_cutoff,info,Nlist)
          !!! Update state variables when cutoff radii (presumably of self-
          !!! organizing particles) have changed externally
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          ! Arguments
          !-------------------------------------------------------------------------
          DEFINE_MK()
          CLASS(DTYPE(ppm_t_vbp))            :: this
          REAL(MK)                                 :: max_cutoff
          !!! New value for the largest cutoff radius
          INTEGER,                  INTENT(   OUT) :: info
          !!! return status. On success, 0
          CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,INTENT(INOUT) :: NList
          !!! Neighbor list for which this cutoff radius
          !!! applies. By default, this is the "standard" Verlet list, with neighbours
          !!! sought within the particle set itself.

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl

          start_subroutine("part_updated_cutoff")

          !-------------------------------------------------------------------------
          !  Set new cutoff
          !-------------------------------------------------------------------------
          IF (PRESENT(NList)) THEN
              check_true(<#this%neighs%has(NList)#>,&
              & "Neighbour list does not concern this particle set")
              NList%uptodate = .FALSE.
              NList%cutoff = max_cutoff
          ELSE
              nl => this%get_neighlist()
              check_associated(nl,"Compute neighbour lists first")
              nl%uptodate = .FALSE.
              nl%cutoff = max_cutoff
          ENDIF

          ! Compute ghostlayer sizes
          IF (max_cutoff.GT.this%ghostlayer) THEN
             !If the new cutoff is larger than the current ghostsize
             ! then the new ghostsize is the new cutoff
             this%ghostlayer = max_cutoff
             ! update states
             this%flags(ppm_part_ghosts) = .FALSE.
          ELSE IF (max_cutoff .LT. this%ghostlayer) THEN
             !Else, we find the new maximum cutoff radius amongst
             !all existing neighbor lists on this Particle set
             this%ghostlayer = 0._MK
             nl => this%neighs%begin()
             DO WHILE (ASSOCIATED(nl))
                IF (nl%cutoff .GT. this%ghostlayer) THEN
                   this%ghostlayer = nl%cutoff
                ENDIF
                nl => this%neighs%next()
             ENDDO
             !no need to update states: ghosts are still ok.
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(vbp_updated_cutoff)
