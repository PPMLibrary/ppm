#define WRAP(a) a


#define __FUNCNAME DTYPE(WRAP(DATANAME)_check)
SUBROUTINE __FUNCNAME(Pc,wp,id,info) 
    !!!------------------------------------------------------------------------!
    !!! Check whether a Data Structure exists and can be accessed at this id
    !!!------------------------------------------------------------------------!

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles))                       :: Pc
    !!! Data structure containing the particles
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER,        INTENT(IN   )   :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER,      INTENT(IN   )   :: wp
#endif
    INTEGER,                            INTENT(IN   )   :: id
    !!! id where the data is stored
    INTEGER,                            INTENT(   OUT)  :: info
    !!! Return status, on success 0.
    CHARACTER(LEN=ppm_char)                             :: caller ='prop_check'
    INTEGER, DIMENSION(:),POINTER :: nullv=>NULL()
    
    info = 0
    !-------------------------------------------------------------------------
    ! Check arguments
    !-------------------------------------------------------------------------
    IF (.NOT.Pc%props%exists(id)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'This property id does not exist.',&
            __LINE__,info)
        RETURN
    ENDIF

#if   __DIM == 1
    IF (Pc%props%vec(id)%t%lda.NE.1) THEN
#elif __DIM == 2
    IF (Pc%props%vec(id)%t%lda.LT.2) THEN
#endif
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'Argument has wrong dimension for this property id',&
            __LINE__,info)
        info = nullv(1)
        RETURN
    ENDIF


    IF (Pc%props%vec(id)%t%data_type.NE. &
#if   __MYTYPE == __INTEGER
        ppm_type_int&
#elif __MYTYPE == __LONGINT
        ppm_type_longint& 
#elif __MYTYPE == __REAL
        ppm_type_real& 
#elif __MYTYPE == __COMPLEX
        ppm_type_comp& 
#elif __MYTYPE == __LOGICAL
        ppm_type_logical& 
#else
        ppm_type_none&
#endif
       &  ) THEN

        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'Argument has the wrong type for this property id',&
            __LINE__,info)
        info = nullv(1)
        RETURN
    ENDIF

END SUBROUTINE __FUNCNAME
#undef __FUNCNAME
#undef __MYTYPE

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get)
#define __CHECKTYPE DTYPE(WRAP(DATANAME)_check)
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
    CHARACTER (LEN = ppm_char)      :: caller = '__FUNCNAME'

    IF (ppt_id .LE. 0) THEN
        write(cbuf,*) 'ERROR: failed to get DATANAME for property ',& 
            'ppt_id = ',ppt_id
        CALL ppm_write(ppm_rank,'DATANAME',cbuf,info)
        wp => NULL()
        RETURN
    ENDIF

    !CALL Pc%props%checktype(wp,ppt_id,info)
    CALL Pc%__CHECKTYPE(wp,ppt_id,info)
    !IF (info.NE.0) THEN
        !info = ppm_error_error
        !WRITE(cbuf,'(A,I0,A)') & 
            !'Type conflict between the requested property id:',ppt_id,&
            !& ' and the target array (1st argument)'  
        !CALL ppm_error(ppm_err_alloc,caller,cbuf,__LINE__,info)
        !wp => NULL()
        !RETURN
    !ENDIF

    IF (ppt_id .LE. Pc%props%max_id) THEN
        ASSOCIATE (prop => Pc%props%vec(ppt_id)%t)
        IF (prop%flags(ppm_ppt_partial)) THEN
            IF (PRESENT(with_ghosts)) THEN
                IF (with_ghosts) THEN
                    IF (prop%flags(ppm_ppt_ghosts)) THEN
                        wp => &
#if   __DIM == 1
                     prop%WRAP(DATANAME)(1:Pc%Mpart)
#elif __DIM == 2
                     prop%WRAP(DATANAME)(:,1:Pc%Mpart)
#endif
                    ELSE
                        write(*,*) line_of_stars
                        write(*,*) 'ERROR: tried to get DATANAME (name = ',&
                            & TRIM(ADJUSTL(prop%name)),&
                            & ') with ghosts when ghosts are not up-to-date. ',&
                            & 'Returning NULL pointer'
                        write(*,*) 'Run with traceback option to debug'
                        write(*,*) line_of_stars
                        wp => NULL()
#ifdef __crash_on_null_pointers
                        !segfault the program. Compile with appropriate compiler
                        !options to check for array bounds and provide traceback
#if   __DIM == 1
                        write(*,*) prop%WRAP(DATANAME)(1)
#elif   __DIM == 2
                        write(*,*) prop%WRAP(DATANAME)(1,1)
#endif
#endif
                    ENDIF
                    RETURN
                ENDIF
            ENDIF
            wp => &
#if   __DIM == 1
                prop%WRAP(DATANAME)(1:Pc%Npart)
#elif __DIM == 2
                prop%WRAP(DATANAME)(:,1:Pc%Npart)
#endif
            RETURN
        ENDIF
        END ASSOCIATE
    ENDIF

    write(cbuf,*) 'ERROR: tried to get DATANAME (name = ',&
        & TRIM(ADJUSTL(Pc%props%vec(ppt_id)%t%name)),&
        & ') when mapping is not up-to-date. ',&
        & 'Returning NULL pointer', &
        'Run with traceback option to debug'
    CALL ppm_write(ppm_rank,cbuf,caller,info)

    wp => NULL()
#ifdef __crash_on_null_pointers
    !segfault the program. Compile with appropriate compiler
    !options to check for array bounds and provide traceback
#if   __DIM == 1
    write(*,*) Pc%props%vec(ppt_id)%t%WRAP(DATANAME)(1)
#elif   __DIM == 2
    write(*,*) Pc%props%vec(ppt_id)%t%WRAP(DATANAME)(1,1)
#endif
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
    Pc%props%vec(ppt_id)%t%flags(ppm_ppt_ghosts) = .FALSE.
    wp => NULL()

END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#undef DATANAME
#undef __TYPE
#undef __MYTYPE
