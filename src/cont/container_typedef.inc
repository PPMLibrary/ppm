!----------------------------------------------------------------------
! Container class
!----------------------------------------------------------------------
TYPE, ABSTRACT :: ppm_t_container
    INTEGER                                         :: nb = 0
    !!! number of different instances
    INTEGER                                         :: iter_id = 0
    !!! index at which the iterator is currently pointing
    INTEGER                                         :: min_id = HUGE(1)
    !!! smallest non-null index in the array
    INTEGER                                         :: max_id = 0
    !!! highest non-null index in the array
    INTEGER                                         :: vec_size = 0
    !!! current allocated size of the array
END TYPE ppm_t_container

!ABSTRACT INTERFACE
!    FUNCTION DTYPE(begin)(this) RESULT(iterator)
!        IMPORT DTYPE(ppm_t_container)
!        CLASS(DTYPE(ppm_t_container)), INTENT(IN) :: this
!        CLASS(*), POINTER    :: iterator
!    END FUNCTION DTYPE(begin)
!END INTERFACE
!ABSTRACT INTERFACE
!    FUNCTION DTYPE(next)(this) RESULT(iterator)
!        IMPORT DTYPE(ppm_t_container)
!        CLASS(DTYPE(ppm_t_container)), INTENT(IN) :: this
!        CLASS(*), POINTER  :: iterator
!    END FUNCTION DTYPE(next)
!END INTERFACE




#undef   DTYPE
