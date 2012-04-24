FUNCTION DTYPE(has_neighlist)(this,Part) RESULT(res)
    !!! Check whether there exists a neighbour list between 
    !!! one particle set and another 
    !!! (default is that the two sets are the same)
    CLASS(DTYPE(ppm_t_particles)),TARGET           :: this
    !!! particle set
    CLASS(DTYPE(ppm_t_particles)_),OPTIONAL,TARGET :: Part
    !!! particle set within which the neighbours are sought
    LOGICAL                                        :: res
    !!! Neighbour list

    !Local variables
    INTEGER                                        :: info
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER         :: NList => NULL()

    res = .TRUE.
    IF (PRESENT(Part)) THEN
        NList => this%neighs%begin()
        DO WHILE(ASSOCIATED(NList))
            !TODO FIXME
                !SELECT TYPE(Part_src => NList%Part)
                !CLASS IS (DTYPE(ppm_t_particles)_)
                    !IF (ASSOCIATED(Part_src,Part)) RETURN
                !END SELECT
            NList => this%neighs%next()
        ENDDO
    ELSE
        NList => this%neighs%begin()
        DO WHILE(ASSOCIATED(NList))
            !TODO FIXME
            !SELECT TYPE(Part_src => NList%Part)
            !CLASS IS (DTYPE(ppm_t_particles)_)
                !IF (ASSOCIATED(Part_src,this)) RETURN
            !END SELECT
            NList => this%neighs%next()
        ENDDO
    ENDIF
    res = .FALSE.
    RETURN

END FUNCTION DTYPE(has_neighlist)


FUNCTION DTYPE(get_neighlist)(this,Part) RESULT(NList)
    !!! Returns the neighbour list between one particle set 
    !!! and another (default is of course that the two sets
    !!! are the same)
    CLASS(DTYPE(ppm_t_particles)),TARGET           :: this
    !!! particle set
    CLASS(DTYPE(ppm_t_particles)_),OPTIONAL,TARGET :: Part
    !!! particle set within which the neighbours are sought
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER         :: NList
    !!! Neighbour list

    !Local variables
    INTEGER                                        :: info

    IF (PRESENT(Part)) THEN
        NList => this%neighs%begin()
        DO WHILE(ASSOCIATED(NList))
            !TODO FIXME
            !SELECT TYPE(Part_src => NList%Part)
            !CLASS IS (DTYPE(ppm_t_particles)_)
                !IF (ASSOCIATED(Part_src,Part)) RETURN
            !END SELECT
            NList => this%neighs%next()
        ENDDO
    ELSE
        NList => this%neighs%begin()
        DO WHILE(ASSOCIATED(NList))
            !TODO FIXME
            !SELECT TYPE(Part_src => NList%Part)
            !CLASS IS (DTYPE(ppm_t_particles)_)
                !IF (ASSOCIATED(Part_src,this)) RETURN
            !END SELECT
            NList => this%neighs%next()
        ENDDO
    ENDIF

    CALL ppm_error(ppm_err_argument,&
        "Could not find neighbour list, returning null pointer",&
        "get_neighlist",__LINE__,info)

    NList => NULL()
END FUNCTION DTYPE(get_neighlist)

SUBROUTINE DTYPE(get_vlist)(this,nvlist,vlist,info,NList)
    !!! returns pointers to the arrays nvlist and vlist
    !!! that contain the Verlet lists for the neighbour list
    !!! of ID nlid.
    CLASS(DTYPE(ppm_t_particles))                :: this
    INTEGER,DIMENSION(:),POINTER,  INTENT( OUT)  :: nvlist
    !!! number of neighbours for each particle
    INTEGER,DIMENSION(:,:),POINTER,INTENT( OUT)  :: vlist
    !!! verlet list
    INTEGER,            INTENT(INOUT)     :: info
    !!! return status. On success, 0
    CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,TARGET :: NList
    !!! Neighbour list (if not the default one)

    CLASS(DTYPE(ppm_t_neighlist)_),POINTER  :: nl => NULL()
    start_subroutine("get_vlist")

    
    IF (PRESENT(NList)) THEN
        nl => NList
    ELSE
        nl => this%get_neighlist()
    ENDIF

    check_associated(nl,"Could not find neighbour list. Make sure they are already computed")
    check_true(nl%uptodate,"Neighbour lists need to be updated")

    nvlist => nl%nvlist
    vlist  => nl%vlist

    end_subroutine()
END SUBROUTINE DTYPE(get_vlist)

SUBROUTINE DTYPE(get_nvlist)(this,nvlist,info,NList)
    !!! returns pointers to the arrays nvlist
    !!! that contain the Verlet lists for the neighbour list
    !!! of ID nlid.
    CLASS(DTYPE(ppm_t_particles))                :: this
    INTEGER,DIMENSION(:),POINTER,  INTENT( OUT)  :: nvlist
    !!! number of neighbours for each particle
    INTEGER,                       INTENT(INOUT) :: info
    !!! return status. On success, 0
    CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,TARGET :: NList
    !!! Neighbour list (if not the default one)

    CLASS(DTYPE(ppm_t_neighlist)_),POINTER  :: nl => NULL()
    start_subroutine("get_vlist")

    
    IF (PRESENT(NList)) THEN
        nl => NList
    ELSE
        nl => this%get_neighlist()
    ENDIF

    check_associated(nl,"Could not find neighbour list. Make sure they are already computed")
    check_true(nl%uptodate,"Neighbour lists need to be updated") 

    nvlist => nl%nvlist

    end_subroutine()
END SUBROUTINE DTYPE(get_nvlist)
