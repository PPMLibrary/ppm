SUBROUTINE DTYPE(get_vlist)(Pc,nvlist,vlist,nlid)
    !!! returns pointers to the arrays nvlist and vlist
    !!! that contain the Verlet lists for the neighbour list
    !!! of ID nlid.
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,DIMENSION(:),        POINTER  :: nvlist
    !!! number of neighbours for each particle
    INTEGER,DIMENSION(:,:),      POINTER  :: vlist
    !!! verlet list
    INTEGER                               :: nlid
    !!! id of the neighbour list 
    INTEGER                               :: info


    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    CHARACTER(LEN = ppm_char)             :: caller = 'get_vlist'

    
    nvlist => NULL()
    vlist => NULL()

    IF (.NOT.Pc%neighs%exists(nlid)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'neighbour list is invalid or not allocated',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Pc%neighs%vec(nlid)%t%uptodate) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'neighbour lists have not been computed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    nvlist => Pc%neighs%vec(nlid)%t%nvlist
    vlist => Pc%neighs%vec(nlid)%t%vlist


    9999 CONTINUE

END SUBROUTINE DTYPE(get_vlist)

SUBROUTINE DTYPE(get_nvlist)(Pc,nvlist,nlid)
    !!! returns a pointer to the array of nb of neighbors
    !!! for the neighbour list of id nlID
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,DIMENSION(:),        POINTER  :: nvlist
    !!! number of neighbours for each particle
    INTEGER                               :: nlid
    !!! id of the neighbour list 
    INTEGER                               :: info

    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    CHARACTER(LEN = ppm_char)             :: caller = 'get_nvlist'
    
    nvlist => NULL()

    IF (.NOT.Pc%neighs%exists(nlid)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'neighbour list is invalid or not allocated',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Pc%neighs%vec(nlid)%t%uptodate) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'neighbour lists have not been computed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    nvlist => Pc%neighs%vec(nlid)%t%nvlist

    9999 CONTINUE

END SUBROUTINE DTYPE(get_nvlist)
