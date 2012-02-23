!!----------------------------------------------------------------------
!! Particle properties
!!----------------------------------------------------------------------
TYPE DTYPE(ppm_t_part_prop)
    !!! Data structure for particle properties
    INTEGER                                        :: data_type
    !!! type of the property
    !!! One of:
    !!!     ppm_param_...
    !!! 
    CHARACTER(len=ppm_char)                        :: name
    !!! name of the property
    LOGICAL, DIMENSION(ppm_param_length_pptflags)  :: flags
    !!! logical flags
    !!!    ppm_ppt_ghosts
    !!!          true if ghost values are up-to-date
    !!!    ppm_ppt_partial
    !!!          true if there is a one-to-one mapping with the particles
    !!!    ppm_ppt_reqput
    !!!    ppm_ppt_map_parts
    !!!          true if partial mappings are desired for this property (default)
    !!!          (if false, the array for this property is not reallocated when
    !!!           particles move to a different processor or when they are
    !!!           interpolated from one distribution to another)
    !!!    ppm_ppt_map_ghosts
    !!!          true if ghost mappings are desired for this property (default)
    INTEGER                                        :: lda
    !!! leading dimension of the data array
    !!!
    !!! pointers to arrays where the scalar-value properties are stored
    INTEGER, DIMENSION(:), POINTER                 :: data_1d_i =>NULL()
    !!! if the data is 1d integer
    INTEGER, DIMENSION(:,:), POINTER               :: data_2d_i =>NULL()
    !!! if the data is 2d integer
    INTEGER(ppm_kind_int64),DIMENSION(:), POINTER  :: data_1d_li =>NULL()
    !!! if the data is 1d long integer
    INTEGER(ppm_kind_int64),DIMENSION(:,:),POINTER :: data_2d_li =>NULL()
    !!! if the data is 2d long integer
    REAL(MK), DIMENSION(:), POINTER                :: data_1d_r => NULL()
    !!! if the data is 1d real
    REAL(MK), DIMENSION(:,:), POINTER              :: data_2d_r =>NULL()
    !!! if the data is 2d real
    COMPLEX(MK), DIMENSION(:), POINTER             :: data_1d_c => NULL()
    !!! if the data is 1d complex
    COMPLEX(MK), DIMENSION(:,:), POINTER           :: data_2d_c =>NULL()
    !!! if the data is 2d complex
    LOGICAL, DIMENSION(:), POINTER                 :: data_1d_l => NULL()
    !!! if the data is 1d logical
    LOGICAL, DIMENSION(:,:), POINTER               :: data_2d_l =>NULL()
    !!! if the data is 2d logical

    CONTAINS
    PROCEDURE :: create => DTYPE(prop_create)
    PROCEDURE :: destroy => DTYPE(prop_destroy)
    PROCEDURE :: print_info => DTYPE(prop_print_info)

END TYPE DTYPE(ppm_t_part_prop)
!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold neighlists
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_part_prop)
    TYPE(DTYPE(ppm_t_part_prop)), POINTER  :: t => NULL()
END TYPE DTYPE(ppm_ptr_part_prop)
!----------------------------------------------------------------------
! Container for properties
!----------------------------------------------------------------------
TYPE, EXTENDS(DTYPE(ppm_t_container)) ::  DTYPE(ppm_c_props)
    TYPE(DTYPE(ppm_ptr_part_prop)),DIMENSION(:),POINTER :: vec => NULL()
    !!! array with pointers to each property for this particle set
    TYPE(DTYPE(ppm_t_part_prop)),              POINTER  :: iterator => NULL()
    CONTAINS
    PROCEDURE       :: begin => DTYPE(begin_prop)
    PROCEDURE       :: next => DTYPE(next_prop)
    PROCEDURE       :: last => DTYPE(last_prop)
    PROCEDURE       :: prev => DTYPE(prev_prop)
    PROCEDURE       :: destroy => DTYPE(prop_container_destroy)

    PROCEDURE     :: DTYPE(data_1d_i_check)
    PROCEDURE     :: DTYPE(data_2d_i_check)
    PROCEDURE     :: DTYPE(data_1d_li_check)
    PROCEDURE     :: DTYPE(data_2d_li_check)
    PROCEDURE     :: DTYPE(data_1d_r_check)
    PROCEDURE     :: DTYPE(data_2d_r_check)
    PROCEDURE     :: DTYPE(data_1d_c_check)
    PROCEDURE     :: DTYPE(data_2d_c_check)
    PROCEDURE     :: DTYPE(data_1d_l_check)
    PROCEDURE     :: DTYPE(data_2d_l_check)


    !NOTE: for some reason, this does not work
    ! (ifort crashes at compile time if there is a call
    !  to this%checktype(wp,wpid,info) )
    ! One can call the not-overloaded procedures directly,
    ! but thats very annoying...
    CONTAINS
        PROCEDURE :: exists => DTYPE(prop_exists)
    GENERIC       :: checktype => & 
        DTYPE(data_1d_i_check),&
        DTYPE(data_2d_i_check),&
        DTYPE(data_1d_li_check),&
        DTYPE(data_2d_li_check),&
        DTYPE(data_1d_r_check),&
        DTYPE(data_2d_r_check),&
        DTYPE(data_1d_c_check),&
        DTYPE(data_2d_c_check),&
        DTYPE(data_1d_l_check),&
        DTYPE(data_2d_l_check)
END TYPE DTYPE(ppm_c_props)


!!----------------------------------------------------------------------
!! Particle neighbor lists
!!----------------------------------------------------------------------
TYPE DTYPE(ppm_t_neighlist)
    CHARACTER(len=ppm_char)                        :: name
    !!! name of the neighbour list
    INTEGER                                         :: P_id = 0
    !!! Id of the set of particles that this neighbour list refers to
    !!! The default, 0, stands for "self".
    REAL(MK)                                        :: cutoff 
    !!! cutoff radius
    REAL(MK)                                        :: skin
    !!! skin layer around the particles
    INTEGER                                         :: isymm
    !!! using symmetry
    INTEGER              , DIMENSION(:  ), POINTER  :: nvlist=> NULL()
    !!! Number of neighbors of each particles
    INTEGER              , DIMENSION(:,:), POINTER  :: vlist=> NULL()
    !!! Neighbor lists
    LOGICAL                                         :: uptodate = .FALSE.
    !!! true if the neighbor lists have been computed
    INTEGER                                         :: nneighmin = 0
    !!! smallest number of neighbors on this processor
    INTEGER                                         :: nneighmax = 0
    !!! highest number of neighbors on this processor

    CONTAINS
    PROCEDURE :: destroy => DTYPE(neigh_destroy)

END TYPE DTYPE(ppm_t_neighlist)
!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold neighlists
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_t_neighlist)
    TYPE(DTYPE(ppm_t_neighlist)), POINTER  :: t => NULL()
END TYPE DTYPE(ppm_ptr_t_neighlist)
!----------------------------------------------------------------------
! Container for neighbour lists
!----------------------------------------------------------------------
TYPE, EXTENDS(DTYPE(ppm_t_container)) ::  DTYPE(ppm_c_neighlists)
    TYPE(DTYPE(ppm_ptr_t_neighlist)),DIMENSION(:),POINTER :: vec => NULL()
    !!! array with pointers to each neighbor list for this particle set
    TYPE(DTYPE(ppm_t_neighlist)),                 POINTER :: iterator => NULL()
    CONTAINS
    PROCEDURE       :: begin => DTYPE(begin_neigh)
    PROCEDURE       :: next => DTYPE(next_neigh)
    PROCEDURE       :: last => DTYPE(last_neigh)
    PROCEDURE       :: prev => DTYPE(prev_neigh)
    PROCEDURE       :: destroy => DTYPE(neigh_container_destroy)
    PROCEDURE       :: exists => DTYPE(neigh_exists)
END TYPE DTYPE(ppm_c_neighlists)




TYPE DTYPE(particles_stats)
    !!! Data structure containing statistics about a particle set
    !!! 
    INTEGER                                          :: nb_nl = 0
    !!! number of neighbour lists constructions
    INTEGER                                          :: nb_inl = 0
    !!! number of inhomogeneous neighbour lists constructions
    INTEGER                                          :: nb_cinl = 0
    !!! number of conventional (old and depreciated) inl constructions
    INTEGER                                          :: nb_xset_inl = 0
    !!! number of xset inhomogeneous neighbour lists constructions
    INTEGER                                          :: nb_xset_nl = 0
    !!! number of xset neighbour lists constructions
    INTEGER                                          :: nb_dc_comp = 0
    !!! number of DC operators computation (matrix inversions)
    INTEGER                                          :: nb_dc_apply = 0
    !!! number of DC operators evaluation (no matrix inversions)
    INTEGER                                          :: nb_kdtree = 0
    !!! number of kdtree constructions
    INTEGER                                          :: nb_global_map = 0
    !!! number of global mappings
    INTEGER                                          :: nb_part_map = 0
    !!! number of partial mappings
    INTEGER                                          :: nb_ghost_get = 0
    !!! number of partial mappings
    INTEGER                                          :: nb_ghost_push = 0
    !!! number of partial mappings
    REAL(MK)                                       :: t_nl = 0._MK
    !!! time spent for neighbour lists constructions
    REAL(MK)                                       :: t_inl = 0._MK
    !!! time spent for inhomogeneous neighbour lists constructions
    REAL(MK)                                       :: t_cinl = 0._MK
    !!! time spent for conventional (old and depreciated) inl constructions
    REAL(MK)                                       :: t_xset_inl = 0._MK
    !!! time spent for xset inhomogeneous neighbour lists constructions
    REAL(MK)                                       :: t_xset_nl = 0._MK
    !!! time spent for xset neighbour lists constructions
    REAL(MK)                                       :: t_dc_comp = 0._MK
    !!! time spent for DC operators computation (matrix inversions)
    REAL(MK)                                       :: t_dc_apply = 0._MK
    !!! time spent for DC operators evaluation (no matrix inversions)
    REAL(MK)                                       :: t_kdtree = 0._MK
    !!! time spent for kdtree constructions
    REAL(MK)                                       :: t_global_map = 0._MK
    !!! time spent for global mappings
    REAL(MK)                                       :: t_part_map = 0._MK
    !!! time spent for partial mappings
    REAL(MK)                                       :: t_ghost_get = 0._MK
    !!! time spent for partial mappings
    REAL(MK)                                       :: t_ghost_push = 0._MK
    !!! time spent for partial mappings

    INTEGER                                          :: nb_ls = 0
    REAL(MK)                                       :: t_ls = 0._MK
    REAL(MK)                                       :: t_add = 0._MK
    REAL(MK)                                       :: t_del = 0._MK
    REAL(MK)                                       :: t_compD = 0._MK

END TYPE DTYPE(particles_stats)


TYPE DTYPE(ppm_t_particles)
    !!! Data structure for a particle set

    CHARACTER(LEN=ppm_char)                         :: name
    !!! name for this particle set
    REAL(MK), DIMENSION(:,:), POINTER               :: xp => NULL()
    !!! positions of the particles
    INTEGER                                         :: Npart
    !!! Number of real particles on this processor
    INTEGER                                         :: Mpart
    !!! Number of particles (including ghosts) on this processor
    LOGICAL, DIMENSION(ppm_param_length_partflags)  :: flags
    !!! logical flags
    !!!    ppm_part_ghosts
    !!!          true if ghost values are up-to-date
    !!!    ppm_part_partial
    !!!          true if the particles have been mapped on the active topology
    !!!    ppm_part_reqput
    !!!    ppm_part_areinside
    !!!          true if all the particles are inside the comp. domain
    !!!    ppm_part_cartesian
    !!!          true if the particles form a Cartesian grid
    !!!    ppm_part_neighlists
    !!!          true if neighbour lists within this particle cloud are
    !!!          up-to-date
    !!!    ppm_part_global_index
    !!!          true if particles have an uptodate global index


    ! Topology (links between Particles and topologies)
    INTEGER                                         :: active_topoid
    !!! Topology on which particles are currently mapped

    ! Special data IDs
    INTEGER                                         :: gi_ID
    !!! Global index


    !ghost layers
    REAL(MK)                                        :: ghostlayer
    !!! ghost layer size used for ghost mapping
    !!! (note: could in principle be different from the cutoff used to 
    !!!  compute the neighbour lists)
    INTEGER                                         :: isymm

    ! Container for Particles properties
    TYPE(DTYPE(ppm_c_props))                        :: props


    ! Container for Neighbor lists
    TYPE(DTYPE(ppm_c_neighlists))                   :: neighs


    ! Container for differential operators
    TYPE(DTYPE(ppm_c_operators))                    :: ops


    ! Container for particle mappings
    TYPE(DTYPE(ppm_c_pmaps))                        :: maps

    ! Container for lists of (pointers to) other particle clouds
    TYPE(DTYPE(ppm_c_Pc))                           :: set_Pc


    ! Load balancing
    REAL(MK), DIMENSION(:), POINTER                 :: pcost=> NULL()
    !!! cost associated to each particle 


    ! stats
    TYPE(DTYPE(particles_stats))                    :: stats
    !!! runtime statistics (e.g. timings, memory)


    ! clocks
    REAL(MK)                                        :: time
    !!! clock for this set of particles
    INTEGER                                         :: itime
    !!! iteration number for this set of particles


    !------------miscellaneous----------

    ! spacings between particles
    REAL(MK)                                      :: h_avg
    !!! (global) average distance between particles
    REAL(MK)                                      :: h_min
    !!! (global) minimum distance between particles

    CONTAINS
    PROCEDURE     :: create => DTYPE(part_create)
    PROCEDURE     :: destroy => DTYPE(part_destroy)
    PROCEDURE     :: initialize => DTYPE(part_initialize)
    PROCEDURE     :: del_parts => DTYPE(part_del_parts)

    PROCEDURE     :: create_prop => DTYPE(part_prop_create)
    PROCEDURE     :: destroy_prop => DTYPE(part_prop_destroy)
    PROCEDURE     :: realloc_prop => DTYPE(part_prop_realloc)

    PROCEDURE     :: create_neighlist => DTYPE(part_neigh_create)
    PROCEDURE     :: set_cutoff => DTYPE(part_set_cutoff)
    PROCEDURE     :: destroy_neighlist => DTYPE(part_neigh_destroy)
    PROCEDURE     :: comp_neighlist => DTYPE(part_neighlist)
    PROCEDURE     :: get_nvlist => DTYPE(get_nvlist)
    PROCEDURE     :: get_vlist => DTYPE(get_vlist)

    PROCEDURE     :: create_op => DTYPE(part_op_create)
    PROCEDURE     :: destroy_op => DTYPE(part_op_destroy)
    PROCEDURE     :: DTYPE(ppm_dcop_compute2d)
    PROCEDURE     :: DTYPE(ppm_dcop_compute3d)
    PROCEDURE     :: comp_op => DTYPE(part_op_compute)
    PROCEDURE     :: apply_op => DTYPE(part_op_apply)
    !       PROCEDURE     :: updated_positions => DTYPE(part_updated_positions)

    PROCEDURE     :: map_part_push_legacy => DTYPE(part_prop_push)
    PROCEDURE     :: map_part_pop_legacy => DTYPE(part_prop_pop)
    PROCEDURE     :: map => DTYPE(part_mapping)
    PROCEDURE     :: map_ghosts => DTYPE(part_mapping_ghosts)

    PROCEDURE     :: move => DTYPE(part_move)
    PROCEDURE     :: apply_bc => DTYPE(part_apply_bc)


    PROCEDURE     :: print_info => DTYPE(part_print_info)

    PROCEDURE     :: comp_global_index => DTYPE(part_comp_global_index)

    PROCEDURE     :: get_xp => DTYPE(get_xp)
    PROCEDURE     :: set_xp => DTYPE(set_xp)
    PROCEDURE     :: get_dcop => DTYPE(get_dcop)
    PROCEDURE     :: set_dcop => DTYPE(set_dcop)
    PROCEDURE     :: DTYPE(data_1d_i_get)
    PROCEDURE     :: DTYPE(data_2d_i_get)
    PROCEDURE     :: DTYPE(data_1d_li_get)
    PROCEDURE     :: DTYPE(data_2d_li_get)
    PROCEDURE     :: DTYPE(data_1d_r_get)
    PROCEDURE     :: DTYPE(data_2d_r_get)
    PROCEDURE     :: DTYPE(data_1d_c_get)
    PROCEDURE     :: DTYPE(data_2d_c_get)
    PROCEDURE     :: DTYPE(data_1d_l_get)
    PROCEDURE     :: DTYPE(data_2d_l_get)

    PROCEDURE     :: DTYPE(data_1d_i_set)
    PROCEDURE     :: DTYPE(data_2d_i_set)
    PROCEDURE     :: DTYPE(data_1d_li_set)
    PROCEDURE     :: DTYPE(data_2d_li_set)
    PROCEDURE     :: DTYPE(data_1d_r_set)
    PROCEDURE     :: DTYPE(data_2d_r_set)
    PROCEDURE     :: DTYPE(data_1d_c_set)
    PROCEDURE     :: DTYPE(data_2d_c_set)
    PROCEDURE     :: DTYPE(data_1d_l_set)
    PROCEDURE     :: DTYPE(data_2d_l_set)
    GENERIC       :: get =>  &
        DTYPE(data_1d_i_get),&
        DTYPE(data_2d_i_get),&
        DTYPE(data_1d_li_get),&
        DTYPE(data_2d_li_get),&
        DTYPE(data_1d_r_get),&
        DTYPE(data_2d_r_get),&
        DTYPE(data_1d_c_get),&
        DTYPE(data_2d_c_get),&
        DTYPE(data_1d_l_get),&
        DTYPE(data_2d_l_get)
    GENERIC       :: set =>  &
        DTYPE(data_1d_i_set),&
        DTYPE(data_2d_i_set),&
        DTYPE(data_1d_li_set),&
        DTYPE(data_2d_li_set),&
        DTYPE(data_1d_r_set),&
        DTYPE(data_2d_r_set),&
        DTYPE(data_1d_c_set),&
        DTYPE(data_2d_c_set),&
        DTYPE(data_1d_l_set),&
        DTYPE(data_2d_l_set)

    PROCEDURE       :: create_map => DTYPE(part_map_create)
    PROCEDURE       :: destroy_map => DTYPE(part_map_destroy)
    PROCEDURE       :: DTYPE(map_part_send)

    PROCEDURE       :: DTYPE(map_part_pop_1d)
    PROCEDURE       :: DTYPE(map_part_pop_1dc)
    PROCEDURE       :: DTYPE(map_part_pop_1di)
    PROCEDURE       :: DTYPE(map_part_pop_1dl)
    PROCEDURE       :: DTYPE(map_part_pop_2d)
    PROCEDURE       :: DTYPE(map_part_pop_2dc)
    PROCEDURE       :: DTYPE(map_part_pop_2di)
    PROCEDURE       :: DTYPE(map_part_pop_2dl)
    PROCEDURE       :: DTYPE(map_part_push_1d)
    PROCEDURE       :: DTYPE(map_part_push_1dc)
    PROCEDURE       :: DTYPE(map_part_push_1di)
    PROCEDURE       :: DTYPE(map_part_push_1dl)
    PROCEDURE       :: DTYPE(map_part_push_2d)
    PROCEDURE       :: DTYPE(map_part_push_2dc)
    PROCEDURE       :: DTYPE(map_part_push_2di)
    PROCEDURE       :: DTYPE(map_part_push_2dl)

    GENERIC         ::  map_part_pop => &
        DTYPE(map_part_pop_1d),&
        DTYPE(map_part_pop_1dc),&
        DTYPE(map_part_pop_1di),&
        DTYPE(map_part_pop_1dl),&
        DTYPE(map_part_pop_2d),&
        DTYPE(map_part_pop_2dc),&
        DTYPE(map_part_pop_2di),&
        DTYPE(map_part_pop_2dl)
    GENERIC         ::  map_part_push => &
        DTYPE(map_part_push_1d),&
        DTYPE(map_part_push_1dc),&
        DTYPE(map_part_push_1di),&
        DTYPE(map_part_push_1dl),&
        DTYPE(map_part_push_2d),&
        DTYPE(map_part_push_2dc),&
        DTYPE(map_part_push_2di),&
        DTYPE(map_part_push_2dl)

END TYPE DTYPE(ppm_t_particles)
!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold ppm_t_particles
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_t_particles)
    TYPE(DTYPE(ppm_t_particles)), POINTER  :: t => NULL()
END TYPE DTYPE(ppm_ptr_t_particles)
!----------------------------------------------------------------------
! Container for lists of (pointers to) Particle clouds
!----------------------------------------------------------------------
TYPE, EXTENDS(DTYPE(ppm_t_container)) ::  DTYPE(ppm_c_Pc)
    TYPE(DTYPE(ppm_ptr_t_particles)),DIMENSION(:), POINTER :: vec => NULL()
    !!! array with pointers to each neighbor list for this particle set
    CLASS(DTYPE(ppm_t_particles)),                 POINTER :: iterator => NULL()
    !    CONTAINS
    !        PROCEDURE       :: begin => DTYPE(begin_pc)
    !        PROCEDURE       :: next => DTYPE(next_pc)
    !        PROCEDURE       :: last => DTYPE(last_pc)
    !        PROCEDURE       :: prev => DTYPE(prev_pc)
    !        PROCEDURE       :: destroy => DTYPE(pc_container_destroy)
END TYPE DTYPE(ppm_c_Pc)


TYPE, EXTENDS(DTYPE(ppm_t_particles)) :: DTYPE(ppm_t_sop)
    !!! an extension of the Particle cloud data structure
    !!! for Self-Organizing Particles

    INTEGER                                         :: nn_sq_id
    !!! index of the wps array where nearest-neighbour distances are stored

    ! Adaptive particles
    LOGICAL                                         :: adaptive
    !!! true if the particles have their own cutoff radii
    !!! in this case, the cutoff will be stored in wps(rcp_id)%vec
    INTEGER                                         :: rcp_id
    !!! index of the wps array where the cutoff radius is stored
    INTEGER                                         :: D_id
    !!! index of the wps array where D is stored
    INTEGER                                         :: Dtilde_id
    !!! index of the wps array where D_tilde is stored
    INTEGER                                         :: adapt_wpid
    !!! index of the wps array where is stored the property on 
    !!! which adaptation is based 
    !!! default is first 1d property that is not rcp_id (if any)
    !!! otherwise, it is rcp_id
    !    INTEGER                                         :: adapt_wpgradid
    !    !!! index of the wpv array where is stored the gradient of the property 
    !    !!! on which adaptation is based (if needed)
    LOGICAL                                         :: level_set
    !!! true if particles carry a level-set function
    INTEGER                                         :: level_id
    !!! index of the wps array where the level-set is stored
    !    INTEGER                                         :: level_old_id
    !!! index of the wps array where the level-set is backed up before adapt

    INTEGER                                         :: level_grad_id
    !!! index of the wps array where the gradient of the level-set is stored
    !    INTEGER                                         :: level_grad_old_id
    !!! index of the wps array where the gradient of the level-set 
    !!! is backed up before adapt


    ! Container for lists of (pointers to) other adaptive particle clouds
    TYPE(DTYPE(ppm_c_aPc))                          :: set_aPc


    ! Anisotropic particles
    LOGICAL                                         :: anisotropic
    !!! true if the particles have their own cutoff radii
    !!! in this case, the G tensor will be stored in wpv(G_id)%vec
    INTEGER                                         :: G_id
    !!! index where G is stored

    !    CONTAINS
    !        PRIVATE
    !        PROCEDURE     :: create => DTYPE(sop_part_create)


END TYPE DTYPE(ppm_t_sop)
!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold ppm_t_sop
!----------------------------------------------------------------------
TYPE DTYPE(ppm_ptr_t_sop)
    TYPE(DTYPE(ppm_t_sop)), POINTER  :: t => NULL()
END TYPE DTYPE(ppm_ptr_t_sop)
!----------------------------------------------------------------------
! Container for lists of (pointers to) adaptive Particle clouds
!----------------------------------------------------------------------
TYPE, EXTENDS(DTYPE(ppm_t_container)) ::  DTYPE(ppm_c_aPc)
    TYPE(DTYPE(ppm_ptr_t_sop)),DIMENSION(:),POINTER  :: vec => NULL()
    !!! array with pointers to each neighbor list for this particle set
    CLASS(DTYPE(ppm_t_sop)),                 POINTER :: iterator => NULL()
    !    CONTAINS
    !        PROCEDURE       :: begin => DTYPE(begin_pc)
    !        PROCEDURE       :: next => DTYPE(next_pc)
    !        PROCEDURE       :: last => DTYPE(last_pc)
    !        PROCEDURE       :: prev => DTYPE(prev_pc)
    !        PROCEDURE       :: destroy => DTYPE(pc_container_destroy)
END TYPE DTYPE(ppm_c_aPc)


#undef   MK
#undef   _MK
#undef   DTYPE
