      FUNCTION DTYPE(has_neighlist)(this,Part) RESULT(res)
        !!! Check whether there exists a neighbour list between
        !!! one particle set and another
        !!! (default is that the two sets are the same)
        IMPLICIT NONE

        CLASS(DTYPE(ppm_t_particles)),     TARGET :: this
        !!! particle set
        CLASS(ppm_t_discr_kind), OPTIONAL, TARGET :: Part
        !!! particle set within which the neighbours are sought
        LOGICAL                                   :: res
        !!! Neighbour list

        !Local variables
        CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: NList

        IF (PRESENT(Part)) THEN
           NList => this%neighs%begin()
           DO WHILE(ASSOCIATED(NList))
              IF (ASSOCIATED(NList%Part,Part)) THEN
                 res = .TRUE.
                 RETURN
              ENDIF
              NList => this%neighs%next()
           ENDDO
        ELSE
           NList => this%neighs%begin()
           DO WHILE(ASSOCIATED(NList))
              IF (ASSOCIATED(NList%Part,this)) THEN
                 res = .TRUE.
                 RETURN
              ENDIF
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
        IMPLICIT NONE

        CLASS(DTYPE(ppm_t_particles)),     TARGET :: this
        !!! particle set
        CLASS(ppm_t_discr_kind), OPTIONAL, TARGET :: Part
        !!! particle set within which the neighbours are sought
        !!! default is the same particle set as "this"
        CLASS(DTYPE(ppm_t_neighlist)_),   POINTER :: NList
        !!! Neighbour list

        !Local variables
        INTEGER :: info

        IF (PRESENT(Part)) THEN
           NList => this%neighs%begin()
           DO WHILE(ASSOCIATED(NList))
              IF (ASSOCIATED(NList%Part,Part)) RETURN
              NList => this%neighs%next()
           ENDDO
        ELSE
           NList => this%neighs%begin()
           DO WHILE(ASSOCIATED(NList))
              IF (ASSOCIATED(NList%Part,this)) RETURN
              NList => this%neighs%next()
           ENDDO
        ENDIF

        info=ppm_error_notice
        CALL ppm_error(ppm_err_argument,&
        &    "Could not find neighbour list, returning null pointer",&
        &    "get_neighlist",__LINE__,info)

        NList => NULL()
      END FUNCTION DTYPE(get_neighlist)

      SUBROUTINE DTYPE(get_vlist)(this,nvlist,vlist,info,NList)
        !!! returns pointers to the arrays nvlist and vlist
        !!! that contain the Verlet lists for the neighbour list
        !!! NList
        IMPLICIT NONE

        CLASS(DTYPE(ppm_t_particles))                           :: this

        INTEGER, DIMENSION(:),                    POINTER       :: nvlist
        !!! number of neighbours for each particle
        INTEGER, DIMENSION(:,:),                  POINTER       :: vlist
        !!! verlet list
        INTEGER,                                  INTENT(INOUT) :: info
        !!! return status. On success, 0
        CLASS(DTYPE(ppm_t_neighlist)_), OPTIONAL, TARGET        :: NList
        !!! Neighbour list (if not the default one)

        !Local variables
        CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: NL

        start_subroutine("get_vlist")

        IF (PRESENT(NList)) THEN
           NL => NList
        ELSE
           NL => this%get_neighlist()
        ENDIF

        check_associated(NL,"Could not find neighbour list. Make sure they are already computed")

        check_true(NL%uptodate,"Neighbour lists need to be updated")

        nvlist => NL%nvlist

        vlist  => NL%vlist

        end_subroutine()
      END SUBROUTINE DTYPE(get_vlist)

      SUBROUTINE DTYPE(get_nvlist)(this,nvlist,info,NList)
        !!! returns pointers to the arrays nvlist
        !!! that contain the Verlet lists for the neighbour list
        !!! of ID nlid.
        IMPLICIT NONE

        CLASS(DTYPE(ppm_t_particles))                    :: this

        INTEGER, DIMENSION(:),             POINTER       :: nvlist
        !!! number of neighbours for each particle
        INTEGER,                           INTENT(INOUT) :: info
        !!! return status. On success, 0
        CLASS(DTYPE(ppm_t_neighlist)_), OPTIONAL, TARGET :: NList
        !!! Neighbour list (if not the default one)

        !Local variables
        CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: NL

        start_subroutine("get_vlist")

        IF (PRESENT(NList)) THEN
           NL => NList
        ELSE
           NL => this%get_neighlist()
        ENDIF

        check_associated(NL,"Could not find neighbour list. Make sure they are already computed")

        check_true(<#NL%uptodate#>,"Neighbour lists need to be updated")

        nvlist => NL%nvlist

        end_subroutine()

      END SUBROUTINE DTYPE(get_nvlist)
