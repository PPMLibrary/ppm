#define WRAP(a) a
#define __FUNCNAME DTYPE(WRAP(DATANAME)_get)

SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,with_ghosts)
    CLASS(DTYPE(ppm_t_particles))   :: Pc
    INTEGER                         :: ppt_id
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER     :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER   :: wp
#endif
    INTEGER                         :: info
    LOGICAL,OPTIONAL                :: with_ghosts

    IF (ppt_id .LE. 0) THEN
        write(cbuf,*) 'ERROR: failed to get DATANAME for property ',& 
            'ppt_id = ',ppt_id
        CALL ppm_write(ppm_rank,'DATANAME',cbuf,info)
        wp => NULL()
        RETURN
    ENDIF

    IF (ppt_id .LE. Pc%max_wpid) THEN
        IF (Pc%props(ppt_id)%t%flags(ppm_ppt_partial)) THEN
            IF (PRESENT(with_ghosts)) THEN
                IF (with_ghosts) THEN
                    IF (Pc%props(ppt_id)%t%flags(ppm_ppt_ghosts)) THEN
                        wp => &
#if   __DIM == 1
                     Pc%props(ppt_id)%t%WRAP(DATANAME)(1:Pc%Mpart)
#elif __DIM == 2
                     Pc%props(ppt_id)%t%WRAP(DATANAME)(:,1:Pc%Mpart)
#endif
                    ELSE
                        write(*,*) line_of_stars
                        write(*,*) 'ERROR: tried to get DATANAME (name = ',&
                            & TRIM(ADJUSTL(Pc%props(ppt_id)%t%name)),&
                            & ') with ghosts when ghosts are not up-to-date. ',&
                            & 'Returning NULL pointer'
                        write(*,*) 'Run with traceback option to debug'
                        write(*,*) line_of_stars
                        wp => NULL()
#ifdef __crash_on_null_pointers
                        !segfault the program. Compile with appropriate compiler
                        !options to check for array bounds and provide traceback
                        write(*,*) WRAP(DATANAME)(1)
#endif
                    ENDIF
                    RETURN
                ENDIF
            ENDIF
            wp => &
#if   __DIM == 1
                Pc%props(ppt_id)%t%WRAP(DATANAME)(1:Pc%Npart)
#elif __DIM == 2
                Pc%props(ppt_id)%t%WRAP(DATANAME)(:,1:Pc%Npart)
#endif
            RETURN
        ENDIF
    ENDIF
    write(*,*) line_of_stars
    write(*,*) 'ERROR: tried to get DATANAME (name = ',&
        & TRIM(ADJUSTL(Pc%props(ppt_id)%t%name)),&
        & ') when mapping is not up-to-date. ',&
        & 'Returning NULL pointer'
    write(*,*) 'Run with traceback option to debug'
    write(*,*) line_of_stars
    wp => NULL()
#ifdef __crash_on_null_pointers
    !segfault the program. Compile with appropriate compiler
    !options to check for array bounds and provide traceback
    write(*,*) WRAP(DATANAME)(1)
#endif

END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set)

SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,read_only,ghosts_ok)
    CLASS(DTYPE(ppm_t_particles))    :: Pc
    INTEGER                          :: ppt_id
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER      :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER    :: wp
#endif

    !If read_only was not explicitely set to true, then assume
    !that ghosts are no longer up to date, unless ghosts_ok was
    ! explicitely set to true
    IF (PRESENT(ghosts_ok)) THEN
        IF (ghosts_ok) THEN
            wp => NULL()
            RETURN
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (read_only) THEN
            wp => NULL()
            RETURN
        ENDIF
    ENDIF

    !Assume that the ghost values are now incorrect
    Pc%props(ppt_id)%t%flags(ppm_ppt_ghosts) = .FALSE.
    wp => NULL()

END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#undef DATANAME
#undef __TYPE
