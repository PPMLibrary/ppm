SUBROUTINE DTYPE(dcop_create)(this,Part_src,Part_to,info,nterms,with_ghosts,&
        vector,interp,order)
    !!! Create a differential operator
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_dcop))              :: this
    CLASS(ppm_t_discr_kind),INTENT(IN   ),TARGET :: Part_src
    !!! Particle set that this operator takes data from.
    CLASS(ppm_t_discr_kind),INTENT(IN   ),TARGET :: Part_to
    !!! Particle set that this operator returns data on.
    INTEGER,                INTENT(OUT)   :: info
    !!! Returns status, 0 upon success.
    INTEGER,                INTENT(IN):: nterms
    !!! Number of terms in the differential operator
    LOGICAL,OPTIONAL,       INTENT(IN   ) :: with_ghosts
    !!! True if the operator should be computed for ghost particles too. 
    !!! Note that the resulting values will be wrong for the ghost particles
    !!! that have some neighbours outside the ghost layers. Default is false.
    LOGICAL,OPTIONAL,       INTENT(IN   ) :: vector
    !!! True if the operator is a vector field. Default is false.
    LOGICAL,OPTIONAL,       INTENT(IN   ) :: interp
    !!! True if the operator interpolates data from one set of particles to
    !!! another. Default is false.
    INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN):: order
    !!! Order of approximation for each term of the differential operator

    start_subroutine("op_create")


    this%flags = .FALSE.
    IF (PRESENT(with_ghosts)) this%flags(ppm_ops_inc_ghosts) = with_ghosts
    IF (PRESENT(interp))      this%flags(ppm_ops_interp) = interp
    IF (PRESENT(vector))      this%flags(ppm_ops_vector) = vector
    ldc(1)=nterms
    CALL ppm_alloc(this%order,ldc,ppm_param_alloc_fit,info)
        or_fail_alloc("this%order")
    IF (PRESENT(order)) THEN
        IF (MINVAL(order).LT.0) THEN
            fail("invalid approx order: must be positive")
        ENDIF
        IF (SIZE(order).NE.nterms) THEN
            fail("Wrong size for Order parameter. Should be nterms")
        ENDIF
        this%order = order
    ELSE
        this%order = 2
    ENDIF
    this%discr_src => Part_src
    this%discr_to => Part_to

    end_subroutine()
END SUBROUTINE DTYPE(dcop_create)

!DESTROY ENTRY
SUBROUTINE DTYPE(dcop_destroy)(this,info)
    !!! Destroy the description for a differential operator
    CLASS(DTYPE(ppm_t_dcop))                  :: this
    INTEGER,                   INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    start_subroutine("op_destroy")

    CALL ppm_alloc(this%ker,ldc,ppm_param_dealloc,info)
        or_fail_dealloc("op%ker")

    this%flags = .FALSE.
    this%discr_src => NULL()
    this%discr_to => NULL()

    dealloc_pointer(this%order)

    end_subroutine()
END SUBROUTINE DTYPE(dcop_destroy)



!COMPUTE DC OPERATOR ON A PARTICLE SET
SUBROUTINE DTYPE(dcop_compute)(this,Field_src,Field_to,info)
#ifdef __MPI
    INCLUDE "mpif.h"
#endif
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_dcop))                     :: this
    CLASS(ppm_t_field_),TARGET,INTENT(IN)    :: Field_src
    CLASS(ppm_t_field_),TARGET,INTENT(INOUT) :: Field_to
    INTEGER,                       INTENT(OUT)   :: info

    INTEGER                                    :: ip,iq,j,ineigh,lda,np_target
    REAL(KIND(1.D0))                           :: t1,t2
    REAL(MK),DIMENSION(:,:),POINTER            :: eta => NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: wps1 => NULL(),wps2=>NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: wpv1 => NULL(),wpv2=>NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: dwps => NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: dwpv => NULL()
    INTEGER, DIMENSION(:),  POINTER            :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER            :: vlist => NULL()
    REAL(MK)                                   :: sig
    LOGICAL                                    :: vector_output
    LOGICAL                                    :: vector_input
    LOGICAL                                    :: with_ghosts

    CLASS(ppm_t_discr_data),POINTER            :: data_src => NULL()
    CLASS(ppm_t_discr_data),POINTER            :: data_to  => NULL()
    start_subroutine("dcop_compute")

    !-------------------------------------------------------------------------
    ! Check arguments
    !-------------------------------------------------------------------------
    check_true("this%flags(ppm_ops_iscomputed)",&
        "Operator is not correctly discretized. Maybe the particles have moved?")
    check_associated("this%discr_src")
    check_associated("this%discr_to")

    CALL Field_src%get_discr(this%discr_src,data_src,info)
        or_fail("Field_src%get_discr failed")

    SELECT TYPE (Part_src => this%discr_src)
    CLASS IS (DTYPE(ppm_t_particles)_)
    SELECT TYPE (Part_to  => this%discr_to)
    CLASS IS (DTYPE(ppm_t_particles)_)


    vector_output = this%flags(ppm_ops_vector)
    vector_input = (data_src%lda .NE. 1)

    check_true("Part_src%has_neighlist(Part_to)",&
        "Please compute xset neighbor lists first")

    check_true("Part_src%has_ghosts(Field_src)",&
        "Ghost for source field need to be updated")

    !If with_ghosts has been set to true (only used in some special cases)
    !the operator is computed for all particles, including ghosts. The
    !normal usage is to loop from 1 to Npart only.
    with_ghosts = this%flags(ppm_ops_inc_ghosts)
    IF (with_ghosts) THEN
        np_target = Part_src%Mpart
    ELSE
        np_target = Part_src%Npart
    ENDIF


    !allocate output field if needed
    !otherwise simply check that the output array had been allocated
    !to the right size
    IF (.NOT. Field_to%is_discretized_on(Part_to)) THEN
        CALL Field_to%discretize_on(Part_to,info,with_ghosts=with_ghosts)
            or_fail("Failed to initialize destination field discretization")
    ENDIF

    CALL Field_to%get_discr(this%discr_to,data_to,info)
        or_fail("Field_to%get_to failed")

    IF (.NOT. data_to%flags(ppm_ppt_partial) .OR. &
            &       with_ghosts .AND. .NOT. data_to%flags(ppm_ppt_ghosts)) THEN
        CALL Field_to%discretize_on(Part_to,info,with_ghosts=with_ghosts)
            or_fail("Failed to reinitialize destination field discretization")
        CALL Field_to%get_discr(this%discr_to,data_to,info)
            or_fail("Field_to%get_to failed")
    ENDIF

    IF (vector_output) THEN
        CALL Part_to%get_field(field_to,dwpv,info,with_ghosts=with_ghosts)
        DO ip = 1,np_target
            DO j=1,this%op_ptr%nterms
                dwpv(1:lda,ip) = 0._MK
            ENDDO
        ENDDO
    ELSE
        CALL Part_to%get_field(field_to,dwps,info,with_ghosts=with_ghosts)
        DO ip = 1,np_target
            dwps(ip) = 0._MK
        ENDDO
    ENDIF

    eta => this%ker(:,1:np_target)


    IF (this%flags(ppm_ops_interp)) THEN
        CALL Part_to%get_vlist(nvlist,vlist,info)
            or_fail("could not access neighbour lists")
        IF (vector_output) THEN
            IF(vector_input) THEN
                CALL Part_src%get_field(field_src,wpv2,info,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            wpv2(1:lda,iq) * eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Part_src%set_field(field_src,wpv2,info,read_only=.TRUE.)
            ELSE
                CALL Part_src%get_field(field_src,wps2,info,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            wps2(iq) * eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Part_src%set_field(field_src,wps2,info,read_only=.TRUE.)
            ENDIF
        ELSE
            CALL Part_src%get_field(field_src,wps2,info,with_ghosts=.TRUE.)
            DO ip = 1,np_target
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwps(ip) = dwps(ip) + wps2(iq) * eta(ineigh,ip)
                ENDDO
            ENDDO
            CALL Part_src%set_field(field_src,wps2,info,read_only=.TRUE.)
        ENDIF
    ELSE
        CALL Part_src%get_vlist(nvlist,vlist,info)
            or_fail("could not access neighbour lists")
        sig = -1._mk !TODO FIXME
        IF (vector_output) THEN
            IF(vector_input) THEN
                CALL Part_src%get_field(field_src,wpv1,info,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            (wpv1(1:lda,iq) + sig*(wpv1(1:lda,ip)))* &
                            eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Part_src%set_field(field_src,wpv1,info,read_only=.TRUE.)
            ELSE
                CALL Part_src%get_field(field_src,wps1,info,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            (wps1(iq) + sig*(wps1(ip)))* &
                            eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Part_src%set_field(field_src,wps1,info,read_only=.TRUE.)
            ENDIF
        ELSE
            CALL Part_src%get_field(field_src,wps1,info,with_ghosts=.TRUE.)
            DO ip = 1,np_target
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwps(ip) = dwps(ip) + &
                        (wps1(iq)+sig*(wps1(ip))) * eta(ineigh,ip)
                ENDDO
            ENDDO
            CALL Part_src%set_field(field_src,wps1,info,read_only=.TRUE.)
        ENDIF
    ENDIF

    eta => NULL()

    IF (vector_output) THEN
        IF (with_ghosts) THEN
            !we assume that the ghosts are up-to-date even though
            !they clearly are not. we assume you know what you are
            !doing when using this option.
            CALL Part_to%set_field(field_to,dwpv,info,ghosts_ok=.TRUE.)
        ELSE
            CALL Part_to%set_field(field_to,dwpv,info)
        ENDIF
    ELSE
        IF (with_ghosts) THEN
            CALL Part_to%set_field(field_to,dwps,info,ghosts_ok=.TRUE.)
        ELSE
            CALL Part_to%set_field(field_to,dwps,info)
        ENDIF
    ENDIF
    nvlist => NULL()
    vlist => NULL()

    Part_src%stats%nb_dc_apply = Part_src%stats%nb_dc_apply + 1
#ifdef __MPI
    t2 = MPI_WTIME(info)
    Part_src%stats%t_dc_apply = Part_src%stats%t_dc_apply+(t2-t1)
#endif


    CLASS DEFAULT
        fail("Wrong type. Operator should be discretized on a particle set")
    END SELECT !discr_to
    CLASS DEFAULT
        fail("Wrong type. Operator should be discretized on a particle set")
    END SELECT !discr_src

    end_subroutine()
END SUBROUTINE

#undef DTYPE
#undef DEFINE_MK
