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

    start_subroutine("get_vlist")

    
    nvlist => NULL()
    vlist => NULL()

    IF (.NOT.Pc%neighs%exists(nlid)) THEN
        fail("neighbour list is invalid or not allocated")
    ENDIF

    IF (.NOT.Pc%neighs%vec(nlid)%t%uptodate) THEN
        fail("neighbour lists have not been computed")
    ENDIF

    nvlist => Pc%neighs%vec(nlid)%t%nvlist
    vlist => Pc%neighs%vec(nlid)%t%vlist

    end_subroutine()
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

    start_subroutine("get_nvlist")
    
    nvlist => NULL()

    IF (.NOT.Pc%neighs%exists(nlid)) THEN
        fail("neighbour list is invalid or not allocated")
    ENDIF

    IF (.NOT.Pc%neighs%vec(nlid)%t%uptodate) THEN
        fail("neighbour lists have not been computed")
    ENDIF

    nvlist => Pc%neighs%vec(nlid)%t%nvlist

    end_subroutine()
END SUBROUTINE DTYPE(get_nvlist)
