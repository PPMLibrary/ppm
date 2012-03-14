#define CONTAINER DTYPE(ppm_c_mdat)
#define __CONTAINER(a) DTYPE(ppm_c_mdat)_/**/a
#define VEC_TYPE DTYPE(ppm_t_ptr_mesh_data)
#define ITERATOR_TYPE DTYPE(ppm_t_mesh_data)
#include "cont/extended_container_typeproc.f"

!CREATE ENTRY
SUBROUTINE DTYPE(mdat_create)(mdat,meshid,datatype,lda,name,flags,info,&
        zero,patchid)
    !!! Constructor for mesh data data structure
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_mesh_data))      :: mdat
    INTEGER,                INTENT(IN) :: meshid
    INTEGER,                INTENT(IN) :: datatype
    INTEGER,                INTENT(IN) :: lda
    !!! number of data components per mesh node
    CHARACTER(LEN=*)                   :: name
    !!! name to this mesh data
    LOGICAL, DIMENSION(ppm_param_length_mdatflags),INTENT(IN) :: flags
    INTEGER,               INTENT(OUT) :: info
    LOGICAL, OPTIONAL,     INTENT( IN) :: zero
    !!! if true, then initialize the data to zero
    INTEGER, OPTIONAL,      INTENT(IN) :: patchid
    !!! patch data structure to be used. Default is 0, i.e. one 
    !!! patch per subdomain, covering the entire subdomain and its ghost
    !!! layer.

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt, ndim,topoid,i
    CHARACTER(LEN=ppm_char)            :: caller = 'mdat_create'
    TYPE(ppm_t_topo),POINTER    :: topo => NULL()


    CALL substart(caller,t0,info)

    topoid = ppm_mesh%vec(meshid)%t%topoid
    topo => ppm_topo(topoid)%t

    mdat%lda       = lda
    mdat%data_type = datatype
    mdat%meshid = meshid
    mdat%name      = name
    mdat%flags     = flags
    mdat%topoid = topoid

    ALLOCATE(mdat%subs(topo%nsublist),STAT=info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'allocating subs array failed',__LINE__,info)
        GOTO 9999
    ENDIF

    DO i=1,topo%nsublist
        CALL mdat%subs(i)%create(topoid,meshid,i,datatype,lda,flags,info,&
            zero,patchid)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'creating subs_data failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDDO

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE DTYPE(mdat_create)

SUBROUTINE DTYPE(sub_create)(sub,topoid,meshid,subid,datatype,lda,&
        flags,info,zero,patchid)
    !!! Constructor for subdomain data data structure
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_sub_data))      :: sub
    INTEGER,                INTENT(IN) :: topoid
    INTEGER,                INTENT(IN) :: meshid
    INTEGER,                INTENT(IN) :: subid
    INTEGER,                INTENT(IN) :: datatype
    INTEGER,                INTENT(IN) :: lda
    !!! number of data components per mesh node
    LOGICAL, DIMENSION(ppm_param_length_mdatflags),INTENT(IN) :: flags
    INTEGER,               INTENT(OUT) :: info
    LOGICAL, OPTIONAL,     INTENT( IN) :: zero
    !!! if true, then initialize the data to zero
    INTEGER, OPTIONAL,      INTENT(IN) :: patchid
    !!! patch data structure to be used. Default is 0, i.e. one 
    !!! patch per subdomain, covering the entire subdomain and its ghost
    !!! layer.

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt,ndim,lpatchid,i,npatch
    CHARACTER(LEN=ppm_char)            :: caller = 'sub_create'
    LOGICAL                            :: zero_data
    TYPE(ppm_t_topo),POINTER    :: topo => NULL()



    CALL substart(caller,t0,info)

    IF (PRESENT(patchid)) THEN
        lpatchid = patchid
    ELSE
        lpatchid = 1
    ENDIF

    IF (lpatchid.EQ.1) THEN
        npatch = 1
    ELSE
        npatch = ppm_topo(topoid)%t%patches(lpatchid)%npatch(subid)
    ENDIF

    sub%topoid = topoid
    sub%subid  = subid
    sub%npatch = npatch

    IF (ASSOCIATED(sub%p)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'ppm_t_sub_data structure not clean. Use destroy() before create()',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    ALLOCATE(sub%p(npatch),STAT=info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'allocating patch array failed',__LINE__,info)
        GOTO 9999
    ENDIF

    DO i=1,npatch
        CALL sub%p(i)%create(lpatchid,topoid,subid,datatype,lda,info,zero)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'creating patch data failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDDO

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE DTYPE(sub_create)

SUBROUTINE DTYPE(sub_destroy)(sub,info)
    !!! Destructor for subdomain data data structure
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_sub_data))      :: sub
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt,i
    CHARACTER(LEN=ppm_char)            :: caller = 'sub_create'
    LOGICAL                            :: zero_data
    TYPE(ppm_t_topo),POINTER    :: topo => NULL()



    CALL substart(caller,t0,info)


    DO i=1,sub%npatch
        CALL sub%p(i)%destroy(info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'destroying patch data failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDDO
    DEALLOCATE(sub%p,STAT=info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'deallocating patch array failed',__LINE__,info)
        GOTO 9999
    ENDIF
    NULLIFY(sub%p)
    sub%topoid = -1
    sub%subid  = -1
    sub%npatch = -1

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE DTYPE(sub_destroy)


!DESTROY ENTRY
SUBROUTINE DTYPE(mdat_destroy)(mdat,info)
    CLASS(DTYPE(ppm_t_mesh_data))      :: mdat
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'mdat_destroy'
    INTEGER                            :: i,topoid

    CALL substart(caller,t0,info)

    IF (ASSOCIATED(mdat%subs)) THEN
        DO i=1,ppm_topo(mdat%topoid)%t%nsublist
            CALL mdat%subs(i)%destroy(info)
        ENDDO
    ENDIF
    NULLIFY(mdat%subs)

    mdat%data_type = ppm_type_none
    mdat%lda = 0
    mdat%flags = .FALSE.
    mdat%name = ''

    CALL substop(caller,t0,info)

END SUBROUTINE DTYPE(mdat_destroy)

SUBROUTINE DTYPE(mdat_print_info)(mdat,info,level,fileunit,mdatid)
    !-----------------------------------------------------------------------
    ! Print out summary information about this mdaterty
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_mesh_data))                          :: mdat
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
    !!! indentation level at which to printout the info. Default = 0
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
    !!! Already open file unit for printout. Default = stdout
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: mdatid
    !!! id of this mdaterty in the parent struct
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                              :: lev,fileu,id
    REAL(KIND(1.D0))                     :: t0
    CHARACTER(LEN = ppm_char)            :: caller = 'mdat_print_info'
    CHARACTER(LEN = ppm_char)            :: myformat
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)


    IF (PRESENT(fileunit)) THEN
        fileu = fileunit
    ELSE
        fileu = 6
    ENDIF
    IF (PRESENT(level)) THEN
        lev = MAX(level,1)
    ELSE
        lev = 1
    ENDIF
    IF (PRESENT(mdatid)) THEN
        id = mdatid
    ELSE
        id = 1
    ENDIF

    WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,I0,A,A,A,I0)'

    WRITE(fileu,myformat) 'mdaterty ',id,': ',TRIM(color_print(mdat%name,33)),&
        ' Type: ',mdat%data_type

    lev = lev + 1

    WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,I0)'
    WRITE(fileu,myformat) 'lda: ',mdat%lda

    WRITE(myformat,'(A,I0,A,I0,A)') '(',4*lev,'X,A,',&
        ppm_param_length_mdatflags,'L)'
    WRITE(fileu,myformat) 'flags: ',mdat%flags


    CALL substop(caller,t0,info)

END SUBROUTINE DTYPE(mdat_print_info)





#undef DTYPE
#undef DEFINE_MK
