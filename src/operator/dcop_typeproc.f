      SUBROUTINE DTYPE(dcop_create)(this,Op,Part_src,Part_to,info, &
      &          with_ghosts,vector,interp,order,order_v,prop)
      !!! Create a discretized differential operator (For DC-PSE)
          DEFINE_MK()
          CLASS(DTYPE(ppm_t_dcop))                                       :: this
          CLASS(ppm_t_operator_),                          TARGET        :: Op
          !!! Generic differential operator that this is a discretization of.
          CLASS(ppm_t_discr_kind),                         TARGET        :: Part_src
          !!! Particle set that this operator takes data from.
          CLASS(ppm_t_discr_kind),                         TARGET        :: Part_to
          !!! Particle set that this operator returns data on.
          INTEGER,                                         INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          LOGICAL,                 OPTIONAL,               INTENT(IN   ) :: with_ghosts
          !!! True if the operator should be computed for ghost particles too.
          !!! Note that the resulting values will be wrong for the ghost particles
          !!! that have some neighbours outside the ghost layers. Default is false.
          LOGICAL,                 OPTIONAL,               INTENT(IN   ) :: vector
          !!! True if the operator is a vector field. Default is false.
          LOGICAL,                 OPTIONAL,               INTENT(IN   ) :: interp
          !!! True if the operator interpolates data from one set of particles to
          !!! another. Default is false.
          !!! WARNING: on output of this routine, one has to update the ghost
          !!! values for the Part_src particle set (because it will have a
          !!! new property that stores nearest-neighbour distances, the ghost
          !!! values of which will be needed)
          INTEGER,                 OPTIONAL,               INTENT(IN   ) :: order
          !!! Order of approximation for all terms of the differential operator
          INTEGER,                 OPTIONAL, DIMENSION(:), POINTER       :: order_v
          !!! Order of approximation for each term of the differential operator
          CLASS(ppm_t_discr_data), OPTIONAL,               TARGET        :: prop

          REAL(MK), DIMENSION(:),   POINTER :: nn2
          REAL(MK), DIMENSION(:,:), POINTER :: xp => NULL()
          REAL(MK)                          :: dist2
          REAL(MK), PARAMETER               :: big=HUGE(1.0_MK)

          INTEGER, DIMENSION(:),   POINTER :: nvlist
          INTEGER, DIMENSION(:,:), POINTER :: vlist
          INTEGER                          :: ip,iq,ineigh
          INTEGER                          :: datatype

          !nearest-neighbor distance, for interpolating operators. Either they
          !were computed already, or this has to be done now.

          start_subroutine("op_create")

          this%flags(ppm_ops_iscomputed)=.FALSE.
          this%flags(ppm_ops_inc_ghosts)=MERGE(with_ghosts,.FALSE.,PRESENT(with_ghosts))
          this%flags(ppm_ops_interp)    =MERGE(interp,     .FALSE.,PRESENT(interp))
          this%flags(ppm_ops_vector)    =MERGE(vector,     .FALSE.,PRESENT(vector))

          ldc(1)=Op%nterms
          CALL ppm_alloc(this%order,ldc,ppm_param_alloc_fit,info)
          or_fail_alloc("this%order")

          IF (PRESENT(order_v)) THEN
             IF (ASSOCIATED(order_v)) THEN
                IF (MINVAL(order_v).LT.0) THEN
                   fail("invalid approx order: must be positive")
                ENDIF
                IF (SIZE(order_v).NE.Op%nterms) THEN
                   fail("Wrong size for Order parameter. Should be Op%nterms")
                ENDIF
                this%order = order_v
             ELSE
                this%order = MERGE(order,2,PRESENT(order))
             ENDIF
          ELSE
             this%order = MERGE(order,2,PRESENT(order))
          ENDIF

          this%discr_src => Part_src
          this%discr_to  => Part_to
          this%op_ptr    => Op

          !!---------------------------------------------------------------------!
          !!There is some work to do for interpolating operators
          !!The square of the distance to nearest neighbours has to be stored
          !!for each particle, including ghost ones.
          !!---------------------------------------------------------------------!
          IF (this%flags(ppm_ops_interp)) THEN
             IF (.NOT.PRESENT(prop)) THEN
                SELECT TYPE(Part_src)
                CLASS IS (DTYPE(ppm_t_particles)_)
                   datatype=MERGE(ppm_type_real,ppm_type_real_single,ppm_kind.EQ.ppm_kind_double)

                   !Create a new property and make this%nn_sq point to it
                   CALL Part_src%create_prop(info,part_prop=this%nn_sq,&
                   &    dtype=datatype,name='nearest_neighb_squared')
                   or_fail("could not create property for nn_sq")

                   NULLIFY(nvlist,vlist,nn2)

                   CALL Part_src%get(this%nn_sq,nn2,info)
                   or_fail("could not get pointer to nn_sq")

                   CALL Part_src%get_xp(xp,info,with_ghosts=.TRUE.)
                   or_fail("could not get pointer to xp")

                   CALL Part_src%get_vlist(nvlist,vlist,info)
                   or_fail("could not get pointer to Verlet lists")

                   DO ip=1,Part_src%Npart
                      nn2(ip) = big
                      DO ineigh=1,nvlist(ip)
                         iq=vlist(ineigh,ip)
                         dist2 = SUM((xp(1:ppm_dim,iq)-xp(1:ppm_dim,ip))**2)
                         nn2(ip) = MIN(nn2(ip),dist2)
                      ENDDO
                   ENDDO

                   CALL Part_src%set_xp(xp,info,read_only=.TRUE.)

                CLASS DEFAULT
                   fail("does not work for this type of discretization")

                END SELECT
             ELSE
                SELECT TYPE(prop)
                CLASS IS (DTYPE(ppm_t_part_prop)_)
                   this%nn_sq => prop

                CLASS DEFAULT
                   fail("wrong type. Prop argument should be of class ppm_t_part_prop")

                END SELECT
             ENDIF
          ENDIF !this%flags(ppm_ops_interp)

          end_subroutine()
      END SUBROUTINE DTYPE(dcop_create)

      !DESTROY ENTRY
      SUBROUTINE DTYPE(dcop_destroy)(this,info)
          !!! Destroy the discretized differential operator
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_dcop)) :: this

          INTEGER,   INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.

          start_subroutine("op_destroy")

          CALL ppm_alloc(this%ker,ldc,ppm_param_dealloc,info)
          or_fail_dealloc("op%ker")

          IF (this%flags(ppm_ops_interp)) THEN
             SELECT TYPE(Part_src => this%discr_src)
             CLASS IS (DTYPE(ppm_t_particles)_)
                CALL Part_src%props%remove(info,this%nn_sq)
                or_fail("Failed to remove property")

                DEALLOCATE(this%nn_sq,STAT=info)
                or_fail_dealloc("Failed to destroy property")

                this%nn_sq => NULL()

             END SELECT
          ENDIF

          this%flags     = .FALSE.

          this%discr_src => NULL()
          this%discr_to  => NULL()
          this%op_ptr    => NULL()

          dealloc_pointer(this%order)

          end_subroutine()
      END SUBROUTINE DTYPE(dcop_destroy)

      !DESTROY ENTRY
      SUBROUTINE DTYPE(dcop_comp_weights)(this,info,c,min_sv)
          !!! Compute the weights for a DC-PSE operator
          !!! (this is an expensive step and has to
          !!! be re-done everytime the particles move)
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_dcop)) :: this

          INTEGER,   INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !---------------------------------------------------------
          ! Optional arguments
          !---------------------------------------------------------
          REAL(MK), OPTIONAL, INTENT(IN   ) :: c
          !!! ratio h/epsilon
          REAL(MK), OPTIONAL, INTENT(  OUT) :: min_sv
          !!! if present, compute the singular value decomposition of the
          !!! vandermonde matrix for each operator and return the smallest one
          start_subroutine("op_comp_weights")

          IF (ppm_dim.EQ.2) THEN
             CALL this%comp_weights_2d(info,c,min_sv)
          ELSE
             CALL this%comp_weights_3d(info,c,min_sv)
          ENDIF
          or_fail("Failed to compute weights for DC operator")

          SELECT TYPE(P => this%discr_src)
          CLASS IS (DTYPE(ppm_t_particles))
             P%stats%nb_dc_comp = P%stats%nb_dc_comp + 1

          END SELECT

          this%flags(ppm_ops_iscomputed) = .TRUE.

          end_subroutine()
      END SUBROUTINE DTYPE(dcop_comp_weights)

      !COMPUTE DC OPERATOR ON A PARTICLE SET
      SUBROUTINE DTYPE(dcop_compute)(this,Field_src,Field_to,info)
        !!! Evaluate a discretized operator (DC-PSE) on a field
        !!! (Field_src) and put the result in another field (Field_to)
        !!! The weights of the operator must have been previously computed.
        !!! The destination field, Field_to, must have been defined (using
        !!! Field_to%create()) but it does not have to be be discretized or
        !!! initialized. This will be done if necessary.
        !!! The dimension of Field_to needs to conform the output of the
        !!! operator (vector or scalar, type, etc...).
        USE ppm_module_mpi
        USE ppm_module_util_time
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Include
        !-------------------------------------------------------------------------

        DEFINE_MK()
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(ppm_t_dcop))               :: this
        !!! Discretized DC-PSE operator, with pre-computed weights
        CLASS(ppm_t_main_abstr), TARGET        :: Field_src
        !!! Input field (or particle property)
        CLASS(ppm_t_main_abstr), TARGET        :: Field_to
        !!! Output field (or particle property)
        INTEGER,                 INTENT(  OUT) :: info
        !!! Return status, on success 0.

        !-------------------------------------------------------------------------
        ! Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_discr_data), POINTER :: data_src
        CLASS(ppm_t_discr_data), POINTER :: data_to
!         CLASS(ppm_t_discr_kind), POINTER :: null_discr  => NULL()

        REAL(MK), DIMENSION(:,:), POINTER :: eta
        REAL(MK), DIMENSION(:),   POINTER :: wps1
        REAL(MK), DIMENSION(:),   POINTER :: wps2
        REAL(MK), DIMENSION(:,:), POINTER :: wpv1
        REAL(MK), DIMENSION(:,:), POINTER :: wpv2
        REAL(MK), DIMENSION(:),   POINTER :: dwps
        REAL(MK), DIMENSION(:,:), POINTER :: dwpv
        REAL(MK)                          :: sig
        REAL(KIND(1.D0))                  :: t1,t2

        INTEGER, DIMENSION(:),   POINTER :: nvlist
        INTEGER, DIMENSION(:,:), POINTER :: vlist
        INTEGER                          :: i,j
        INTEGER                          :: ip,iq,jdim
        INTEGER                          :: ineigh,lda,np_target

        LOGICAL :: vector_operator
        !true if operator returns a vector field
        LOGICAL :: vector_output
        LOGICAL :: vector_input
        LOGICAL :: with_ghosts

        CHARACTER(LEN=ppm_char) :: fname

        start_subroutine("dcop_compute")

        CALL ppm_util_time(t1)

        !-------------------------------------------------------------------------
        ! Check arguments
        !-------------------------------------------------------------------------
        check_true(<#this%flags(ppm_ops_iscomputed)#>,&
        & "Operator is not correctly discretized. Maybe the particles have moved?")
        check_associated(<#this%discr_src#>)
        check_associated(<#this%discr_to#>)

        SELECT TYPE (Part_src => this%discr_src)
        CLASS IS (DTYPE(ppm_t_particles)_)
           SELECT TYPE (Part_to => this%discr_to)
           CLASS IS (DTYPE(ppm_t_particles)_)
              !Get a pointer to the discretized data (2 cases, depending on whether
              ! the input is a Field or a Particle Property object).
              SELECT TYPE(Field_src)
              CLASS IS (ppm_t_field_)
                 NULLIFY(data_src)
                 CALL Field_src%get_discr(Part_src,data_src,info)
                 or_fail("Failed to access the discretization of the source field")

              CLASS IS (DTYPE(ppm_t_part_prop)_)
                 data_src => Field_src

              CLASS DEFAULT
                 fail("Wrong type for Field_src")

              END SELECT

              vector_operator=this%flags(ppm_ops_vector)
              vector_input   =data_src%lda.NE.1

              ! Check that the destination field has been defined
              ! if not, create new one with default dimensions.
              SELECT TYPE(Field_to)
              CLASS IS (ppm_t_field_)
                 IF (Field_to%lda.LE.0) THEN
                    check_false(<#ASSOCIATED(Field_to%discr_info)#>,&
                    & "Destination field seems to be corrupted - Try Field%destroy()?")

                    WRITE(fname,'(A,A)') "Output_from_",TRIM(ADJUSTL(this%op_ptr%name))

                    IF (vector_operator) THEN
                       CALL Field_to%create(this%op_ptr%nterms,info,name=fname)
                    ELSE
                       CALL Field_to%create(data_src%lda,info,name=fname)
                    ENDIF
                    or_fail("could not create new destination field")
                 ENDIF

              CLASS IS (DTYPE(ppm_t_part_prop)_)

              CLASS DEFAULT
                fail("Wrong type for Field_to")

              END SELECT

              !set lda, the leading dimension of the output data
              IF (vector_operator) THEN
                 !the operator returns vector values
                 lda = this%op_ptr%nterms
                 ! check for compatibility with the input
                 IF (vector_input) THEN
                    check_true(<#lda.EQ.data_src%lda#>,&
                    & "Leading dimension of input field does not match nterms of operator")
                 ENDIF
              ELSE
                 ! Each component of the input field will be fed
                 ! to the operator.
                 lda = data_src%lda
              ENDIF

              check_true(<#Part_src%has_neighlist(Part_to)#>,&
              & "Please compute (maybe xset?) neighbor lists first")

              check_true(<#data_src%has_ghosts()#>,&
              & "Ghost for source field need to be updated")

              !If with_ghosts has been set to true (only used in some special cases)
              !the operator is computed for all particles, including ghosts. The
              !normal usage is to loop from 1 to Npart only.
              with_ghosts = this%flags(ppm_ops_inc_ghosts)
              np_target = MERGE(Part_src%Mpart,Part_src%Npart,with_ghosts)

              !Get a pointer to the discretized data (2 cases, depending on whether
              ! the input is a Field or a Particle Property object).
              SELECT TYPE(Field_to)
              CLASS IS (ppm_t_field_)
                 IF (.NOT. ASSOCIATED(Field_to%discr_info)) THEN
                    CALL Field_to%create(lda,info,name=this%op_ptr%name,&
                    &    dtype=data_src%data_type)
                    or_fail("Failed to create field for the output of the operator")
                 ENDIF
                 !allocate output field if needed
                 !otherwise simply check that the output array had been allocated
                 !to the right size
                 IF (.NOT. Field_to%is_discretized_on(Part_to)) THEN
                    CALL Field_to%discretize_on(Part_to,info,with_ghosts=with_ghosts)
                    or_fail("Failed to initialize destination field discretization")
                 ENDIF

                 NULLIFY(data_to)
                 CALL Field_to%get_discr(this%discr_to,data_to,info)
                 or_fail("Field_to%get_to failed")

                 IF (.NOT.data_to%flags(ppm_ppt_partial).OR. &
                 &  with_ghosts.AND..NOT.data_to%flags(ppm_ppt_ghosts)) THEN
                    CALL Field_to%discretize_on(Part_to,info,with_ghosts=with_ghosts)
                    or_fail("Failed to reinitialize destination field discretization")

                    CALL Field_to%get_discr(this%discr_to,data_to,info)
                    or_fail("Field_to%get_to failed")
                 ENDIF

              CLASS IS (DTYPE(ppm_t_part_prop)_)
                 data_to => Field_to

              CLASS DEFAULT
                 fail("Wrong type for Field_to")

              END SELECT


              ! check for compatibility with the output
              IF (.NOT.vector_operator) THEN
                 check_true(<#data_to%lda.EQ.data_src%lda#>,&
                 & "With this operator, output and input fields should have the same dimensions.")
              ENDIF

              vector_output = (data_to%lda.GT.1)

              IF (vector_output) THEN
                 NULLIFY(dwpv)
                 CALL Part_to%get(data_to,dwpv,info,with_ghosts=with_ghosts)
                 or_fail("Cannot access data_to on this particle set")

                 FORALL (i=1:lda,j=1:np_target) dwpv(i,j)=0.0_MK
              ELSE
                 NULLIFY(dwps)
                 CALL Part_to%get(data_to,dwps,info,with_ghosts=with_ghosts)
                 or_fail("Cannot access data_to on this particle set")

                 FORALL (i=1:np_target) dwps(i)=0.0_MK
              ENDIF

              eta => this%ker(:,1:np_target)

              NULLIFY(nvlist,vlist)
              IF (this%flags(ppm_ops_interp)) THEN
                 CALL Part_to%get_vlist(nvlist,vlist,info)
                 or_fail("could not access neighbour lists")

                 IF (vector_operator) THEN
                    IF (vector_input) THEN
                       NULLIFY(wpv2)
                       CALL Part_src%get(data_src,wpv2,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wpv2")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             DO jdim=1,lda
                                dwpv(jdim,ip) = dwpv(jdim,ip)+ &
                                & wpv2(jdim,iq) * eta(jdim+(ineigh-1)*lda,ip)
                             ENDDO
                          ENDDO
                       ENDDO
                    ELSE
                       NULLIFY(wps2)
                       CALL Part_src%get(data_src,wps2,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wps2")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             DO jdim=1,lda
                                dwpv(jdim,ip) = dwpv(jdim,ip)+ &
                                & wps2(iq) * eta(jdim+(ineigh-1)*lda,ip)
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDIF
                 ELSE
                    IF (vector_output) THEN
                       NULLIFY(wpv2)
                       CALL Part_src%get(data_src,wpv2,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wpv2")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             DO jdim=1,lda
                                dwpv(jdim,ip) = dwpv(jdim,ip)+ &
                                & wpv2(jdim,iq) * eta(ineigh,ip)
                             ENDDO
                          ENDDO
                       ENDDO
                    ELSE
                       NULLIFY(wps2)
                       CALL Part_src%get(data_src,wps2,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wps2")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             dwps(ip) = dwps(ip) + wps2(iq) * eta(ineigh,ip)
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDIF
              ELSE
                 CALL Part_src%get_vlist(nvlist,vlist,info)
                 or_fail("could not access neighbour lists")
                 sig = -1._mk
                 IF (vector_operator) THEN
                    IF (vector_input) THEN
                       NULLIFY(wpv1)
                       CALL Part_src%get(data_src,wpv1,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wpv1")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             DO jdim=1,lda
                                dwpv(jdim,ip) = dwpv(jdim,ip)+(wpv1(jdim,iq)+ &
                                & sig*(wpv1(jdim,ip)))*eta(jdim+(ineigh-1)*lda,ip)
                             ENDDO
                          ENDDO
                       ENDDO
                    ELSE
                       NULLIFY(wps1)
                       CALL Part_src%get(data_src,wps1,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wps1")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             DO jdim=1,lda
                                dwpv(jdim,ip)=dwpv(jdim,ip)+(wps1(iq)+sig*(wps1(ip)))* &
                                & eta(jdim+(ineigh-1)*lda,ip)
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDIF
                 ELSE
                    IF (vector_output) THEN
                       NULLIFY(wpv1)
                       CALL Part_src%get(data_src,wpv1,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wpv1")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             DO jdim=1,lda
                                dwpv(jdim,ip) = dwpv(jdim,ip)+ &
                                & (wpv1(jdim,iq)+sig*(wpv1(jdim,ip)))*eta(ineigh,ip)
                             ENDDO
                          ENDDO
                       ENDDO
                    ELSE
                       NULLIFY(wps1)
                       CALL Part_src%get(data_src,wps1,info,&
                       &    with_ghosts=.TRUE.,read_only=.TRUE.)
                       or_fail("could not access wps1")

                       DO ip = 1,np_target
                          DO ineigh = 1,nvlist(ip)
                             iq = vlist(ineigh,ip)
                             dwps(ip) = dwps(ip)+(wps1(iq)+sig*(wps1(ip)))*eta(ineigh,ip)
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDIF
              ENDIF

              NULLIFY(eta)

              IF (lda.GT.1) THEN
                 IF (with_ghosts) THEN
                    !we assume that the ghosts are up-to-date even though
                    !they clearly are not. we assume you know what you are
                    !doing when using this option.
                    CALL Part_to%set(data_to,dwpv,info,ghosts_ok=.TRUE.)
                 ELSE
                    CALL Part_to%set(data_to,dwpv,info)
                 ENDIF
              ELSE
                 IF (with_ghosts) THEN
                    CALL Part_to%set(data_to,dwps,info,ghosts_ok=.TRUE.)
                 ELSE
                    CALL Part_to%set(data_to,dwps,info)
                 ENDIF
              ENDIF

              Part_src%stats%nb_dc_apply = Part_src%stats%nb_dc_apply + 1

              CALL ppm_util_time(t2)
              Part_src%stats%t_dc_apply = Part_src%stats%t_dc_apply+(t2-t1)

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
