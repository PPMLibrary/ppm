!!----------------------------------------------------------------------
!! Particle properties
!!----------------------------------------------------------------------
TYPE,ABSTRACT,EXTENDS(ppm_t_discr_data):: DTYPE(ppm_t_part_prop)_
    !!! Data structure for particle properties
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
    PROCEDURE(DTYPE(prop_create)_),    DEFERRED :: create
    PROCEDURE(DTYPE(prop_destroy)_),   DEFERRED :: destroy 
    PROCEDURE(DTYPE(prop_print_info)_),DEFERRED :: print_info

    PROCEDURE(DTYPE(data_1d_i_check)_),DEFERRED :: DTYPE(data_1d_i_check)
    PROCEDURE(DTYPE(data_2d_i_check)_),DEFERRED :: DTYPE(data_2d_i_check)
    PROCEDURE(DTYPE(data_1d_li_check)_),DEFERRED :: DTYPE(data_1d_li_check)
    PROCEDURE(DTYPE(data_2d_li_check)_),DEFERRED :: DTYPE(data_2d_li_check)
    PROCEDURE(DTYPE(data_1d_r_check)_),DEFERRED :: DTYPE(data_1d_r_check)
    PROCEDURE(DTYPE(data_2d_r_check)_),DEFERRED :: DTYPE(data_2d_r_check)
    PROCEDURE(DTYPE(data_1d_c_check)_),DEFERRED :: DTYPE(data_1d_c_check)
    PROCEDURE(DTYPE(data_2d_c_check)_),DEFERRED :: DTYPE(data_2d_c_check)
    PROCEDURE(DTYPE(data_1d_l_check)_),DEFERRED :: DTYPE(data_1d_l_check)
    PROCEDURE(DTYPE(data_2d_l_check)_),DEFERRED :: DTYPE(data_2d_l_check)
    !NOTE: for some reason, this does not work
    ! (ifort crashes at compile time if there is a call
    !  to this%checktype(wp,wpid,info) )
    ! One can call the not-overloaded procedures directly,
    ! but thats very annoying...
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

END TYPE DTYPE(ppm_t_part_prop)_
minclude ppm_create_collection(DTYPE(part_prop)_,DTYPE(part_prop)_,generate="abstract")

!!----------------------------------------------------------------------
!! Particle neighbor lists
!!----------------------------------------------------------------------
TYPE,ABSTRACT :: DTYPE(ppm_t_neighlist)_
    CHARACTER(len=ppm_char)                         :: name
    !!! name of the neighbour list
    CLASS(ppm_t_discr_kind),POINTER                 :: Part => NULL()
    !!! Pointer to the set of particles to which the neighbours 
    !!! (in the neighbor lists) belongs to.
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
    PROCEDURE(DTYPE(neigh_destroy)_),DEFERRED :: destroy 

END TYPE DTYPE(ppm_t_neighlist)_
minclude ppm_create_collection(DTYPE(neighlist)_,DTYPE(neighlist)_,generate="abstract")


TYPE,ABSTRACT :: DTYPE(particles_stats)_
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

END TYPE DTYPE(particles_stats)_


TYPE,ABSTRACT,EXTENDS(ppm_t_discr_kind) :: DTYPE(ppm_t_particles)_
    !!! Data structure for a particle set

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
    !!!          true if neighbour lists within this Particle set are
    !!!          up-to-date
    !!!    ppm_part_global_index
    !!!          true if particles have an uptodate global index


    ! Topology (links between Particles and topologies)
    INTEGER                                         :: active_topoid
    !!! Topology on which particles are currently mapped

    CLASS(DTYPE(ppm_t_part_prop)_),POINTER          :: gi => NULL()
    !!! Global index


    !ghost layers
    REAL(MK)                                        :: ghostlayer
    !!! ghost layer size used for ghost mapping
    !!! (note: could in principle be different from the cutoff used to 
    !!!  compute the neighbour lists)
    INTEGER                                         :: isymm
    !!! Symmetric interactions
    ! Container for Particles properties
    CLASS(DTYPE(ppm_c_part_prop)_),POINTER          :: props => NULL()


    ! Container for Neighbor lists
    CLASS(DTYPE(ppm_c_neighlist)_),POINTER          :: neighs => NULL()



    ! Container for particle mappings
    CLASS(DTYPE(ppm_c_part_mapping)_),POINTER       :: maps => NULL()

    !! List of IDs of other Particle sets
    !TYPE(idList)                                    :: set_Pc


    ! Load balancing
    REAL(MK), DIMENSION(:), POINTER                 :: pcost=> NULL()
    !!! cost associated to each particle 


    CLASS(ppm_v_main_abstr),POINTER                 :: field_ptr => NULL()
    !!! Pointers to the fields that are currently discretized on this 
    !!! Particle set

    CLASS(ppm_c_operator_discr_),POINTER            :: ops => NULL()
    !!! Pointers to the operators that are currently discretized for this 
    !!! Particle set


    ! stats
    CLASS(DTYPE(particles_stats)_),ALLOCATABLE      :: stats
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
    PROCEDURE(DTYPE(part_create)_),      DEFERRED :: create 
    PROCEDURE(DTYPE(part_destroy)_),     DEFERRED :: destroy 
    PROCEDURE(DTYPE(part_initialize)_),  DEFERRED :: initialize 
    PROCEDURE(DTYPE(part_del_parts)_),   DEFERRED :: del_parts 

    PROCEDURE(DTYPE(part_prop_create)_), DEFERRED :: create_prop 
    PROCEDURE(DTYPE(part_prop_destroy)_),DEFERRED :: destroy_prop 
    PROCEDURE(DTYPE(part_prop_realloc)_),DEFERRED :: realloc_prop 
    PROCEDURE(DTYPE(part_get_discr)_),   DEFERRED :: get_discr
    PROCEDURE(DTYPE(part_prop_zero)_),   DEFERRED :: zero

    PROCEDURE(DTYPE(part_neigh_create)_),DEFERRED :: create_neighlist 
    PROCEDURE(DTYPE(part_set_cutoff)_),  DEFERRED :: set_cutoff 
    PROCEDURE(DTYPE(part_neigh_destroy)_),DEFERRED:: destroy_neighlist 
    PROCEDURE(DTYPE(part_neighlist)_),   DEFERRED :: comp_neighlist 
    PROCEDURE(DTYPE(get_nvlist)_),       DEFERRED :: get_nvlist 
    PROCEDURE(DTYPE(get_vlist)_),        DEFERRED :: get_vlist 
    PROCEDURE(DTYPE(get_neighlist)_),    DEFERRED :: get_neighlist
    PROCEDURE(DTYPE(has_neighlist)_),    DEFERRED :: has_neighlist
    PROCEDURE(DTYPE(has_ghosts)_),       DEFERRED :: has_ghosts

!    PROCEDURE(DTYPE(part_dcop_create)_),DEFERRED :: create_op 
!    PROCEDURE(DTYPE(part_dcop_destroy)_),DEFERRED :: destroy_op 
!    PROCEDURE(DTYPE(ppm_dcop_compute2d)_),DEFERRED :: DTYPE(ppm_dcop_compute2d)
!    PROCEDURE(DTYPE(ppm_dcop_compute3d)_),DEFERRED :: DTYPE(ppm_dcop_compute3d)
!    PROCEDURE(DTYPE(part_op_compute)_),DEFERRED :: comp_op 
!    PROCEDURE(DTYPE(part_op_apply)_),DEFERRED :: apply_op 
    !       PROCEDURE(DTYPE(part_updated_positions)_),DEFERRED :: updated_positions 

    PROCEDURE(DTYPE(part_prop_push)_),         DEFERRED :: map_part_push_legacy 
    PROCEDURE(DTYPE(part_prop_pop)_),          DEFERRED  :: map_part_pop_legacy 
    PROCEDURE(DTYPE(part_map)_),               DEFERRED  :: map 
    PROCEDURE(DTYPE(part_map_ghost_get)_),     DEFERRED :: map_ghost_get
    PROCEDURE(DTYPE(part_map_ghost_push)_),    DEFERRED :: map_ghost_push
    PROCEDURE(DTYPE(part_map_ghost_send)_),    DEFERRED :: map_ghost_send
    PROCEDURE(DTYPE(part_map_ghost_pop)_),     DEFERRED :: map_ghost_pop
    PROCEDURE(DTYPE(part_map_ghost_pop_pos)_), DEFERRED :: map_ghost_pop_pos
    PROCEDURE(DTYPE(part_map_ghosts)_),        DEFERRED :: map_ghosts 

    PROCEDURE(DTYPE(part_move)_),DEFERRED :: move 
    PROCEDURE(DTYPE(part_apply_bc)_),DEFERRED :: apply_bc 

    PROCEDURE(DTYPE(part_p2m)_), DEFERRED :: interp_to_mesh

    PROCEDURE(DTYPE(part_print_info)_),DEFERRED :: print_info 

    PROCEDURE(DTYPE(part_comp_global_index)_),DEFERRED :: comp_global_index 

    PROCEDURE(DTYPE(get_xp)_),DEFERRED :: get_xp 
    PROCEDURE(DTYPE(set_xp)_),DEFERRED :: set_xp 

    PROCEDURE(DTYPE(data_1d_i_get_prop)_),DEFERRED :: DTYPE(data_1d_i_get_prop)
    PROCEDURE(DTYPE(data_2d_i_get_prop)_),DEFERRED :: DTYPE(data_2d_i_get_prop)
    PROCEDURE(DTYPE(data_1d_li_get_prop)_),DEFERRED :: &
        DTYPE(data_1d_li_get_prop)
    PROCEDURE(DTYPE(data_2d_li_get_prop)_),DEFERRED :: &
        DTYPE(data_2d_li_get_prop)
    PROCEDURE(DTYPE(data_1d_r_get_prop)_),DEFERRED :: DTYPE(data_1d_r_get_prop)
    PROCEDURE(DTYPE(data_2d_r_get_prop)_),DEFERRED :: DTYPE(data_2d_r_get_prop)
    PROCEDURE(DTYPE(data_1d_c_get_prop)_),DEFERRED :: DTYPE(data_1d_c_get_prop)
    PROCEDURE(DTYPE(data_2d_c_get_prop)_),DEFERRED :: DTYPE(data_2d_c_get_prop)
    PROCEDURE(DTYPE(data_1d_l_get_prop)_),DEFERRED :: DTYPE(data_1d_l_get_prop)
    PROCEDURE(DTYPE(data_2d_l_get_prop)_),DEFERRED :: DTYPE(data_2d_l_get_prop)
    !GENERIC       :: get_prop =>  &
        !DTYPE(data_1d_i_get_prop),&
        !DTYPE(data_2d_i_get_prop),&
        !DTYPE(data_1d_li_get_prop),&
        !DTYPE(data_2d_li_get_prop),&
        !DTYPE(data_1d_r_get_prop),&
        !DTYPE(data_2d_r_get_prop),&
        !DTYPE(data_1d_c_get_prop),&
        !DTYPE(data_2d_c_get_prop),&
        !DTYPE(data_1d_l_get_prop),&
        !DTYPE(data_2d_l_get_prop)

    PROCEDURE(DTYPE(data_1d_i_set_prop)_),DEFERRED :: DTYPE(data_1d_i_set_prop)
    PROCEDURE(DTYPE(data_2d_i_set_prop)_),DEFERRED :: DTYPE(data_2d_i_set_prop)
    PROCEDURE(DTYPE(data_1d_li_set_prop)_),DEFERRED :: &
        DTYPE(data_1d_li_set_prop)
    PROCEDURE(DTYPE(data_2d_li_set_prop)_),DEFERRED :: &
        DTYPE(data_2d_li_set_prop)
    PROCEDURE(DTYPE(data_1d_r_set_prop)_),DEFERRED :: DTYPE(data_1d_r_set_prop)
    PROCEDURE(DTYPE(data_2d_r_set_prop)_),DEFERRED :: DTYPE(data_2d_r_set_prop)
    PROCEDURE(DTYPE(data_1d_c_set_prop)_),DEFERRED :: DTYPE(data_1d_c_set_prop)
    PROCEDURE(DTYPE(data_2d_c_set_prop)_),DEFERRED :: DTYPE(data_2d_c_set_prop)
    PROCEDURE(DTYPE(data_1d_l_set_prop)_),DEFERRED :: DTYPE(data_1d_l_set_prop)
    PROCEDURE(DTYPE(data_2d_l_set_prop)_),DEFERRED :: DTYPE(data_2d_l_set_prop)
    !GENERIC       :: set_prop =>  &
        !DTYPE(data_1d_i_set_prop),&
        !DTYPE(data_2d_i_set_prop),&
        !DTYPE(data_1d_li_set_prop),&
        !DTYPE(data_2d_li_set_prop),&
        !DTYPE(data_1d_r_set_prop),&
        !DTYPE(data_2d_r_set_prop),&
        !DTYPE(data_1d_c_set_prop),&
        !DTYPE(data_2d_c_set_prop),&
        !DTYPE(data_1d_l_set_prop),&
        !DTYPE(data_2d_l_set_prop)

    PROCEDURE(DTYPE(data_1d_i_get_field)_),DEFERRED :: DTYPE(data_1d_i_get_field)
    PROCEDURE(DTYPE(data_2d_i_get_field)_),DEFERRED :: DTYPE(data_2d_i_get_field)
    PROCEDURE(DTYPE(data_1d_li_get_field)_),DEFERRED :: &
        DTYPE(data_1d_li_get_field)
    PROCEDURE(DTYPE(data_2d_li_get_field)_),DEFERRED :: &
        DTYPE(data_2d_li_get_field)
    PROCEDURE(DTYPE(data_1d_r_get_field)_),DEFERRED :: DTYPE(data_1d_r_get_field)
    PROCEDURE(DTYPE(data_2d_r_get_field)_),DEFERRED :: DTYPE(data_2d_r_get_field)
    PROCEDURE(DTYPE(data_1d_c_get_field)_),DEFERRED :: DTYPE(data_1d_c_get_field)
    PROCEDURE(DTYPE(data_2d_c_get_field)_),DEFERRED :: DTYPE(data_2d_c_get_field)
    PROCEDURE(DTYPE(data_1d_l_get_field)_),DEFERRED :: DTYPE(data_1d_l_get_field)
    PROCEDURE(DTYPE(data_2d_l_get_field)_),DEFERRED :: DTYPE(data_2d_l_get_field)
    GENERIC       :: get_field =>  &
        DTYPE(data_1d_i_get_field),&
        DTYPE(data_2d_i_get_field),&
        DTYPE(data_1d_li_get_field),&
        DTYPE(data_2d_li_get_field),&
        DTYPE(data_1d_r_get_field),&
        DTYPE(data_2d_r_get_field),&
        DTYPE(data_1d_c_get_field),&
        DTYPE(data_2d_c_get_field),&
        DTYPE(data_1d_l_get_field),&
        DTYPE(data_2d_l_get_field)

    PROCEDURE(DTYPE(data_1d_i_set_field)_),DEFERRED :: DTYPE(data_1d_i_set_field)
    PROCEDURE(DTYPE(data_2d_i_set_field)_),DEFERRED :: DTYPE(data_2d_i_set_field)
    PROCEDURE(DTYPE(data_1d_li_set_field)_),DEFERRED :: &
        DTYPE(data_1d_li_set_field)
    PROCEDURE(DTYPE(data_2d_li_set_field)_),DEFERRED :: &
        DTYPE(data_2d_li_set_field)
    PROCEDURE(DTYPE(data_1d_r_set_field)_),DEFERRED :: DTYPE(data_1d_r_set_field)
    PROCEDURE(DTYPE(data_2d_r_set_field)_),DEFERRED :: DTYPE(data_2d_r_set_field)
    PROCEDURE(DTYPE(data_1d_c_set_field)_),DEFERRED :: DTYPE(data_1d_c_set_field)
    PROCEDURE(DTYPE(data_2d_c_set_field)_),DEFERRED :: DTYPE(data_2d_c_set_field)
    PROCEDURE(DTYPE(data_1d_l_set_field)_),DEFERRED :: DTYPE(data_1d_l_set_field)
    PROCEDURE(DTYPE(data_2d_l_set_field)_),DEFERRED :: DTYPE(data_2d_l_set_field)
    GENERIC       :: set_field =>  &
        DTYPE(data_1d_i_set_field),&
        DTYPE(data_2d_i_set_field),&
        DTYPE(data_1d_li_set_field),&
        DTYPE(data_2d_li_set_field),&
        DTYPE(data_1d_r_set_field),&
        DTYPE(data_2d_r_set_field),&
        DTYPE(data_1d_c_set_field),&
        DTYPE(data_2d_c_set_field),&
        DTYPE(data_1d_l_set_field),&
        DTYPE(data_2d_l_set_field)

    PROCEDURE(DTYPE(data_1d_i_get)_),DEFERRED :: DTYPE(data_1d_i_get)
    PROCEDURE(DTYPE(data_2d_i_get)_),DEFERRED :: DTYPE(data_2d_i_get)
    PROCEDURE(DTYPE(data_1d_li_get)_),DEFERRED :: DTYPE(data_1d_li_get)
    PROCEDURE(DTYPE(data_2d_li_get)_),DEFERRED :: DTYPE(data_2d_li_get)
    PROCEDURE(DTYPE(data_1d_r_get)_),DEFERRED :: DTYPE(data_1d_r_get)
    PROCEDURE(DTYPE(data_2d_r_get)_),DEFERRED :: DTYPE(data_2d_r_get)
    PROCEDURE(DTYPE(data_1d_c_get)_),DEFERRED :: DTYPE(data_1d_c_get)
    PROCEDURE(DTYPE(data_2d_c_get)_),DEFERRED :: DTYPE(data_2d_c_get)
    PROCEDURE(DTYPE(data_1d_l_get)_),DEFERRED :: DTYPE(data_1d_l_get)
    PROCEDURE(DTYPE(data_2d_l_get)_),DEFERRED :: DTYPE(data_2d_l_get)
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
        DTYPE(data_2d_l_get),&
        DTYPE(data_1d_i_get_prop),&
        DTYPE(data_2d_i_get_prop),&
        DTYPE(data_1d_li_get_prop),&
        DTYPE(data_2d_li_get_prop),&
        DTYPE(data_1d_r_get_prop),&
        DTYPE(data_2d_r_get_prop),&
        DTYPE(data_1d_c_get_prop),&
        DTYPE(data_2d_c_get_prop),&
        DTYPE(data_1d_l_get_prop),&
        DTYPE(data_2d_l_get_prop),&
        DTYPE(data_1d_i_get_field),&
        DTYPE(data_2d_i_get_field),&
        DTYPE(data_1d_li_get_field),&
        DTYPE(data_2d_li_get_field),&
        DTYPE(data_1d_r_get_field),&
        DTYPE(data_2d_r_get_field),&
        DTYPE(data_1d_c_get_field),&
        DTYPE(data_2d_c_get_field),&
        DTYPE(data_1d_l_get_field),&
        DTYPE(data_2d_l_get_field)


    PROCEDURE(DTYPE(data_1d_i_set)_),DEFERRED :: DTYPE(data_1d_i_set)
    PROCEDURE(DTYPE(data_2d_i_set)_),DEFERRED :: DTYPE(data_2d_i_set)
    PROCEDURE(DTYPE(data_1d_li_set)_),DEFERRED :: DTYPE(data_1d_li_set)
    PROCEDURE(DTYPE(data_2d_li_set)_),DEFERRED :: DTYPE(data_2d_li_set)
    PROCEDURE(DTYPE(data_1d_r_set)_),DEFERRED :: DTYPE(data_1d_r_set)
    PROCEDURE(DTYPE(data_2d_r_set)_),DEFERRED :: DTYPE(data_2d_r_set)
    PROCEDURE(DTYPE(data_1d_c_set)_),DEFERRED :: DTYPE(data_1d_c_set)
    PROCEDURE(DTYPE(data_2d_c_set)_),DEFERRED :: DTYPE(data_2d_c_set)
    PROCEDURE(DTYPE(data_1d_l_set)_),DEFERRED :: DTYPE(data_1d_l_set)
    PROCEDURE(DTYPE(data_2d_l_set)_),DEFERRED :: DTYPE(data_2d_l_set)
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
        DTYPE(data_2d_l_set),&
        DTYPE(data_1d_i_set_prop),&
        DTYPE(data_2d_i_set_prop),&
        DTYPE(data_1d_li_set_prop),&
        DTYPE(data_2d_li_set_prop),&
        DTYPE(data_1d_r_set_prop),&
        DTYPE(data_2d_r_set_prop),&
        DTYPE(data_1d_c_set_prop),&
        DTYPE(data_2d_c_set_prop),&
        DTYPE(data_1d_l_set_prop),&
        DTYPE(data_2d_l_set_prop)

    PROCEDURE(DTYPE(part_map_create)_),DEFERRED :: create_map 
    PROCEDURE(DTYPE(part_map_destroy)_),DEFERRED :: destroy_map 

    !PROCEDURE(DTYPE(map_part_send)_),DEFERRED :: DTYPE(map_part_send)

    !PROCEDURE(DTYPE(map_part_pop_1d)_),DEFERRED :: DTYPE(map_part_pop_1d)
    !PROCEDURE(DTYPE(map_part_pop_1dc)_),DEFERRED :: DTYPE(map_part_pop_1dc)
    !PROCEDURE(DTYPE(map_part_pop_1di)_),DEFERRED :: DTYPE(map_part_pop_1di)
    !PROCEDURE(DTYPE(map_part_pop_1dl)_),DEFERRED :: DTYPE(map_part_pop_1dl)
    !PROCEDURE(DTYPE(map_part_pop_2d)_),DEFERRED :: DTYPE(map_part_pop_2d)
    !PROCEDURE(DTYPE(map_part_pop_2dc)_),DEFERRED :: DTYPE(map_part_pop_2dc)
    !PROCEDURE(DTYPE(map_part_pop_2di)_),DEFERRED :: DTYPE(map_part_pop_2di)
    !PROCEDURE(DTYPE(map_part_pop_2dl)_),DEFERRED :: DTYPE(map_part_pop_2dl)
    !PROCEDURE(DTYPE(map_part_push_1d)_),DEFERRED :: DTYPE(map_part_push_1d)
    !PROCEDURE(DTYPE(map_part_push_1dc)_),DEFERRED :: DTYPE(map_part_push_1dc)
    !PROCEDURE(DTYPE(map_part_push_1di)_),DEFERRED :: DTYPE(map_part_push_1di)
    !PROCEDURE(DTYPE(map_part_push_1dl)_),DEFERRED :: DTYPE(map_part_push_1dl)
    !PROCEDURE(DTYPE(map_part_push_2d)_),DEFERRED :: DTYPE(map_part_push_2d)
    !PROCEDURE(DTYPE(map_part_push_2dc)_),DEFERRED :: DTYPE(map_part_push_2dc)
    !PROCEDURE(DTYPE(map_part_push_2di)_),DEFERRED :: DTYPE(map_part_push_2di)
    !PROCEDURE(DTYPE(map_part_push_2dl)_),DEFERRED :: DTYPE(map_part_push_2dl)

    !GENERIC         ::  map_part_pop => &
        !DTYPE(map_part_pop_1d),&
        !DTYPE(map_part_pop_1dc),&
        !DTYPE(map_part_pop_1di),&
        !DTYPE(map_part_pop_1dl),&
        !DTYPE(map_part_pop_2d),&
        !DTYPE(map_part_pop_2dc),&
        !DTYPE(map_part_pop_2di),&
        !DTYPE(map_part_pop_2dl)
    !GENERIC         ::  map_part_push => &
        !DTYPE(map_part_push_1d),&
        !DTYPE(map_part_push_1dc),&
        !DTYPE(map_part_push_1di),&
        !DTYPE(map_part_push_1dl),&
        !DTYPE(map_part_push_2d),&
        !DTYPE(map_part_push_2dc),&
        !DTYPE(map_part_push_2di),&
        !DTYPE(map_part_push_2dl)


END TYPE DTYPE(ppm_t_particles)_
minclude ppm_create_collection(DTYPE(particles)_,DTYPE(particles)_,generate="abstract")
!minclude define_abstract_collection_type(DTYPE(ppm_t_sop)_)

!TYPE,ABSTRACT,EXTENDS(DTYPE(ppm_t_particles)_) :: DTYPE(ppm_t_sop)_
    !!!! an extension of the Particle set data structure
    !!!! for Self-Organizing Particles

    !INTEGER                                         :: nn_sq_id
    !!!! index of the wps array where nearest-neighbour distances are stored

    !! Adaptive particles
    !LOGICAL                                         :: adaptive
    !!!! true if the particles have their own cutoff radii
    !!!! in this case, the cutoff will be stored in wps(rcp_id)%vec
    !INTEGER                                         :: rcp_id
    !!!! index of the wps array where the cutoff radius is stored
    !INTEGER                                         :: D_id
    !!!! index of the wps array where D is stored
    !INTEGER                                         :: Dtilde_id
    !!!! index of the wps array where D_tilde is stored
    !INTEGER                                         :: adapt_wpid
    !!!! index of the wps array where is stored the property on 
    !!!! which adaptation is based 
    !!!! default is first 1d property that is not rcp_id (if any)
    !!!! otherwise, it is rcp_id
    !!    INTEGER                                         :: adapt_wpgradid
    !!    !!! index of the wpv array where is stored the gradient of the property 
    !!    !!! on which adaptation is based (if needed)
    !LOGICAL                                         :: level_set
    !!!! true if particles carry a level-set function
    !INTEGER                                         :: level_id
    !!!! index of the wps array where the level-set is stored
    !!    INTEGER                                         :: level_old_id
    !!!! index of the wps array where the level-set is backed up before adapt

    !INTEGER                                         :: level_grad_id
    !!!! index of the wps array where the gradient of the level-set is stored
    !!    INTEGER                                         :: level_grad_old_id
    !!!! index of the wps array where the gradient of the level-set 
    !!!! is backed up before adapt


    !! List of IDs of other adaptive Particle sets
    !TYPE(idList)                                    :: set_aPc


    !! Anisotropic particles
    !LOGICAL                                         :: anisotropic
    !!!! true if the particles have their own cutoff radii
    !!!! in this case, the G tensor will be stored in wpv(G_id)%vec
    !INTEGER                                         :: G_id
    !!!! index where G is stored

    !!    CONTAINS
    !!        PRIVATE
    !!        PROCEDURE     :: create => DTYPE(sop_part_create)


!END TYPE DTYPE(ppm_t_sop)_


