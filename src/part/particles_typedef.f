!!----------------------------------------------------------------------
!! Particle properties
!!----------------------------------------------------------------------
TYPE,EXTENDS(DTYPE(ppm_t_part_prop)_) :: DTYPE(ppm_t_part_prop)
    !!! Data structure for particle properties
    CONTAINS
    PROCEDURE :: create => DTYPE(prop_create)
    PROCEDURE :: destroy => DTYPE(prop_destroy)
    PROCEDURE :: print_info => DTYPE(prop_print_info)

    PROCEDURE :: DTYPE(data_1d_i_check)
    PROCEDURE :: DTYPE(data_2d_i_check)
    PROCEDURE :: DTYPE(data_1d_li_check)
    PROCEDURE :: DTYPE(data_2d_li_check)
    PROCEDURE :: DTYPE(data_1d_r_check)
    PROCEDURE :: DTYPE(data_2d_r_check)
    PROCEDURE :: DTYPE(data_1d_c_check)
    PROCEDURE :: DTYPE(data_2d_c_check)
    PROCEDURE :: DTYPE(data_1d_l_check)
    PROCEDURE :: DTYPE(data_2d_l_check)

END TYPE
minclude ppm_create_collection(DTYPE(part_prop),DTYPE(part_prop),generate="extend")

!!----------------------------------------------------------------------
!! Particle neighbor lists
!!----------------------------------------------------------------------
TYPE,EXTENDS(DTYPE(ppm_t_neighlist)_) :: DTYPE(ppm_t_neighlist)
    CONTAINS
    PROCEDURE :: destroy => DTYPE(neigh_destroy)
END TYPE
minclude ppm_create_collection(DTYPE(neighlist),DTYPE(neighlist),generate="extend")

TYPE,EXTENDS( DTYPE(particles_stats)_) :: DTYPE(particles_stats)
END TYPE


TYPE,EXTENDS(DTYPE(ppm_t_particles)_) :: DTYPE(ppm_t_particles)
    !!! Data structure for a particle set
    CONTAINS
    PROCEDURE     :: create => DTYPE(part_create)
    PROCEDURE     :: destroy => DTYPE(part_destroy)
    PROCEDURE     :: initialize => DTYPE(part_initialize)
    PROCEDURE     :: del_parts => DTYPE(part_del_parts)

    PROCEDURE     :: create_prop  => DTYPE(part_prop_create)
    PROCEDURE     :: destroy_prop => DTYPE(part_prop_destroy)
    PROCEDURE     :: realloc_prop => DTYPE(part_prop_realloc)
    PROCEDURE     :: get_discr    => DTYPE(part_get_discr)
    PROCEDURE     :: zero         => DTYPE(part_prop_zero)

    PROCEDURE     :: create_neighlist => DTYPE(part_neigh_create)
    PROCEDURE     :: set_cutoff => DTYPE(part_set_cutoff)
    PROCEDURE     :: destroy_neighlist => DTYPE(part_neigh_destroy)
    PROCEDURE     :: comp_neighlist => DTYPE(part_neighlist)
    PROCEDURE     :: get_nvlist => DTYPE(get_nvlist)
    PROCEDURE     :: get_vlist => DTYPE(get_vlist)
    PROCEDURE     :: get_neighlist => DTYPE(get_neighlist)
    PROCEDURE     :: has_neighlist => DTYPE(has_neighlist)
    PROCEDURE     :: has_ghosts => DTYPE(has_ghosts)

!    PROCEDURE     :: create_op  => DTYPE(part_dcop_create)
!    PROCEDURE     :: destroy_op => DTYPE(part_dcop_destroy)
!    PROCEDURE     :: DTYPE(ppm_dcop_compute2d)
!    PROCEDURE     :: DTYPE(ppm_dcop_compute3d)
!    PROCEDURE     :: comp_op => DTYPE(part_op_compute)
!    PROCEDURE     :: apply_op => DTYPE(part_op_apply)
    !       PROCEDURE     :: updated_positions => DTYPE(part_updated_positions)

    PROCEDURE     :: map_part_push_legacy => DTYPE(part_prop_push)
    PROCEDURE     :: map_part_pop_legacy  => DTYPE(part_prop_pop)

    PROCEDURE     :: map_ghost_get     => DTYPE(part_map_ghost_get)
    PROCEDURE     :: map_ghost_push    => DTYPE(part_map_ghost_push)
    PROCEDURE     :: map_ghost_send    => DTYPE(part_map_ghost_send)
    PROCEDURE     :: map_ghost_pop     => DTYPE(part_map_ghost_pop)
    PROCEDURE     :: map_ghost_pop_positions => DTYPE(part_map_ghost_pop_pos)
    PROCEDURE     :: map_ghost_push_positions => DTYPE(part_map_ghost_push_pos)
    PROCEDURE     :: map_ghosts        => DTYPE(part_map_ghosts)
    PROCEDURE     :: map               => DTYPE(part_map)
    PROCEDURE     :: map_positions     => DTYPE(part_map_positions)
    PROCEDURE     :: map_push          => DTYPE(part_map_push)
    PROCEDURE     :: map_send          => DTYPE(part_map_send)
    PROCEDURE     :: map_pop           => DTYPE(part_map_pop)
    PROCEDURE     :: map_pop_positions => DTYPE(part_map_pop_positions)

    PROCEDURE     :: move              => DTYPE(part_move)
    PROCEDURE     :: apply_bc          => DTYPE(part_apply_bc)

    PROCEDURE     :: interp_to_mesh    => DTYPE(part_p2m)
    PROCEDURE     :: remesh            => DTYPE(part_remesh)


    PROCEDURE     :: print_info        => DTYPE(part_print_info)

    PROCEDURE     :: comp_global_index => DTYPE(part_comp_global_index)

    PROCEDURE     :: get_xp => DTYPE(get_xp)
    PROCEDURE     :: set_xp => DTYPE(set_xp)

    PROCEDURE     :: DTYPE(data_1d_i_get_prop)
    PROCEDURE     :: DTYPE(data_2d_i_get_prop)
    PROCEDURE     :: DTYPE(data_1d_li_get_prop)
    PROCEDURE     :: DTYPE(data_2d_li_get_prop)
    PROCEDURE     :: DTYPE(data_1d_r_get_prop)
    PROCEDURE     :: DTYPE(data_2d_r_get_prop)
    PROCEDURE     :: DTYPE(data_1d_c_get_prop)
    PROCEDURE     :: DTYPE(data_2d_c_get_prop)
    PROCEDURE     :: DTYPE(data_1d_l_get_prop)
    PROCEDURE     :: DTYPE(data_2d_l_get_prop)

    PROCEDURE     :: DTYPE(data_1d_i_set_prop)
    PROCEDURE     :: DTYPE(data_2d_i_set_prop)
    PROCEDURE     :: DTYPE(data_1d_li_set_prop)
    PROCEDURE     :: DTYPE(data_2d_li_set_prop)
    PROCEDURE     :: DTYPE(data_1d_r_set_prop)
    PROCEDURE     :: DTYPE(data_2d_r_set_prop)
    PROCEDURE     :: DTYPE(data_1d_c_set_prop)
    PROCEDURE     :: DTYPE(data_2d_c_set_prop)
    PROCEDURE     :: DTYPE(data_1d_l_set_prop)
    PROCEDURE     :: DTYPE(data_2d_l_set_prop)

    PROCEDURE     :: DTYPE(data_1d_i_get_field)
    PROCEDURE     :: DTYPE(data_2d_i_get_field)
    PROCEDURE     :: DTYPE(data_1d_li_get_field)
    PROCEDURE     :: DTYPE(data_2d_li_get_field)
    PROCEDURE     :: DTYPE(data_1d_r_get_field)
    PROCEDURE     :: DTYPE(data_2d_r_get_field)
    PROCEDURE     :: DTYPE(data_1d_c_get_field)
    PROCEDURE     :: DTYPE(data_2d_c_get_field)
    PROCEDURE     :: DTYPE(data_1d_l_get_field)
    PROCEDURE     :: DTYPE(data_2d_l_get_field)

    PROCEDURE     :: DTYPE(data_1d_i_set_field)
    PROCEDURE     :: DTYPE(data_2d_i_set_field)
    PROCEDURE     :: DTYPE(data_1d_li_set_field)
    PROCEDURE     :: DTYPE(data_2d_li_set_field)
    PROCEDURE     :: DTYPE(data_1d_r_set_field)
    PROCEDURE     :: DTYPE(data_2d_r_set_field)
    PROCEDURE     :: DTYPE(data_1d_c_set_field)
    PROCEDURE     :: DTYPE(data_2d_c_set_field)
    PROCEDURE     :: DTYPE(data_1d_l_set_field)
    PROCEDURE     :: DTYPE(data_2d_l_set_field)

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

    PROCEDURE       :: create_map => DTYPE(part_map_create)
    PROCEDURE       :: destroy_map => DTYPE(part_map_destroy)
!    PROCEDURE       :: DTYPE(map_part_send)

!    PROCEDURE       :: DTYPE(map_part_pop_1d)
!    PROCEDURE       :: DTYPE(map_part_pop_1dc)
!    PROCEDURE       :: DTYPE(map_part_pop_1di)
!    PROCEDURE       :: DTYPE(map_part_pop_1dl)
!    PROCEDURE       :: DTYPE(map_part_pop_2d)
!    PROCEDURE       :: DTYPE(map_part_pop_2dc)
!    PROCEDURE       :: DTYPE(map_part_pop_2di)
!    PROCEDURE       :: DTYPE(map_part_pop_2dl)
!    PROCEDURE       :: DTYPE(map_part_push_1d)
!    PROCEDURE       :: DTYPE(map_part_push_1dc)
!    PROCEDURE       :: DTYPE(map_part_push_1di)
!    PROCEDURE       :: DTYPE(map_part_push_1dl)
!    PROCEDURE       :: DTYPE(map_part_push_2d)
!    PROCEDURE       :: DTYPE(map_part_push_2dc)
!    PROCEDURE       :: DTYPE(map_part_push_2di)
!    PROCEDURE       :: DTYPE(map_part_push_2dl)
!
!    GENERIC         ::  map_part_pop => &
!        DTYPE(map_part_pop_1d),&
!        DTYPE(map_part_pop_1dc),&
!        DTYPE(map_part_pop_1di),&
!        DTYPE(map_part_pop_1dl),&
!        DTYPE(map_part_pop_2d),&
!        DTYPE(map_part_pop_2dc),&
!        DTYPE(map_part_pop_2di),&
!        DTYPE(map_part_pop_2dl)
!    GENERIC         ::  map_part_push => &
!        DTYPE(map_part_push_1d),&
!        DTYPE(map_part_push_1dc),&
!        DTYPE(map_part_push_1di),&
!        DTYPE(map_part_push_1dl),&
!        DTYPE(map_part_push_2d),&
!        DTYPE(map_part_push_2dc),&
!        DTYPE(map_part_push_2di),&
!        DTYPE(map_part_push_2dl)

END TYPE DTYPE(ppm_t_particles)
minclude ppm_create_collection(DTYPE(particles),DTYPE(particles),generate="extend")


#undef   MK
#undef   _MK
#undef   DTYPE
#undef   CTYPE
