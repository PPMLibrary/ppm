#define WRAP(a) a

#define __FUNCNAME DTYPE(WRAP(DATANAME)_check)
#define __FUNCTYPE WRAP(DATANAME)
      SUBROUTINE __FUNCNAME(this,wp,info)
          !!!------------------------------------------------------------------------!
          !!! Check whether a Data Structure exists and can be accessed
          !!!------------------------------------------------------------------------!
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop))                                 :: this
          !!! Data structure containing the particles
#if   __DIM == 1
          __TYPE, DIMENSION(:  ), POINTER, INTENT(IN   ) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(IN   ) :: wp
#endif

          INTEGER,                                        INTENT(  OUT) :: info
          !!! Return status, on success 0.

          start_subroutine(__FUNCNAME)
          !-------------------------------------------------------------------------
          ! Check arguments
          !-------------------------------------------------------------------------
#if   __DIM == 1
          IF (this%lda.NE.1) THEN
#elif __DIM == 2
          IF (this%lda.LT.2) THEN
#endif
             fail ("Argument has wrong dimension for this property")
          ENDIF

          info=ppm_error_notice

          SELECT CASE (this%data_type)
#if   __MYTYPE == __INTEGER
          CASE (ppm_type_int)
#elif __MYTYPE == __LONGINT
          CASE (ppm_type_longint)
#elif __MYTYPE == __REAL
          CASE (ppm_type_real,ppm_type_real_single)
#elif __MYTYPE == __COMPLEX
          CASE (ppm_type_comp,ppm_type_comp_single)
#elif __MYTYPE == __LOGICAL
          CASE (ppm_type_logical)
#endif
             IF (ASSOCIATED(this%__FUNCTYPE)) info=0
          CASE DEFAULT
             fail ("Argument has wrong datatype for this property")
          END SELECT

          end_subroutine()

      END SUBROUTINE __FUNCNAME
#undef __FUNCNAME
#undef __FUNCTYPE

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get_prop)
      SUBROUTINE __FUNCNAME(this,discr_data,wp,info,with_ghosts,read_only,skip_checks)

          IMPLICIT NONE

          CLASS(DTYPE(ppm_t_particles))                                 :: this
          CLASS(ppm_t_discr_data),                        INTENT(INOUT) :: discr_data

#if   __DIM == 1
          __TYPE, DIMENSION(:  ), POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          !!! data array

          INTEGER,                                        INTENT(  OUT) :: info
          !!! Return status, on success 0.

          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: with_ghosts
          !!! returns array between 1:Mpart (default is 1:Npart)
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          !!! swear on your favourite book that you will not modify the data
          !!! that you access. Default is false (which implies that the state
          !!! variables will be modify with the assumption that the data
          !!! accessed has been changed, e.g. so that a subsequent ghost update
          !!! will also update this property)
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: skip_checks
          !!! Only for users with huge cojones
          !!! [NOTE]
          !!! wp points to the data array

          INTEGER :: np

          LOGICAL :: skip

          start_subroutine(__FUNCNAME)

          NULLIFY(wp)

          SELECT TYPE(discr_data)
          CLASS IS (DTYPE(ppm_t_part_prop)_)
             IF (ppm_debug.GE.1) THEN
                CALL discr_data%checktype(wp,info)
                or_fail("Argument is not accessible for this property")
             ENDIF

             skip=MERGE(skip_checks,.FALSE.,PRESENT(skip_checks))

             IF (skip) THEN
#if   __DIM == 1
                wp => discr_data%WRAP(DATANAME)
#elif __DIM == 2
                wp => discr_data%WRAP(DATANAME)
#endif
             ELSE
                np=this%Npart
                IF (PRESENT(with_ghosts)) THEN
                   IF (with_ghosts) THEN
                      IF (discr_data%flags(ppm_ppt_ghosts)) THEN
                         np=this%Mpart
                      ELSE
                         stdout("ERROR: tried to get DATANAME (name = ", &
                         & 'TRIM(ADJUSTL(discr_data%name))',             &
                         & ") with ghosts when ghosts are not up-to-date. Returning NULL pointer")

                         fail("ghosts are not up-to-date")
                      ENDIF
                   ENDIF
                ENDIF

                IF (discr_data%flags(ppm_ppt_partial)) THEN
#if   __DIM == 1
                   wp => discr_data%WRAP(DATANAME)(1:np)
#elif __DIM == 2
                   wp => discr_data%WRAP(DATANAME)(:,1:np)
#endif
                ELSE
                   stdout("ERROR: tried to get DATANAME (name = ", &
                   & 'TRIM(ADJUSTL(discr_data%name))',             &
                   & ") when mapping is not up-to-date.",          &
                   & "Returning NULL pointer Run with traceback option to debug")

                   fail("unmapped particles")
                ENDIF
             ENDIF !skip

             !Assume that the ghost values are now incorrect, unless explicitely
             !told otherwise
             IF (PRESENT(read_only)) THEN
                IF (.NOT.read_only) THEN
                   discr_data%flags(ppm_ppt_ghosts)=.FALSE.
                ENDIF
             ELSE
                discr_data%flags(ppm_ppt_ghosts)=.FALSE.
             ENDIF
          CLASS DEFAULT
             fail("wrong type. Discretized data should be a ppm_t_part_prop")
          END SELECT

          check_associated(wp,"Get_Prop returned a NULL pointer",ppm_err_sub_failed)

          end_subroutine()
      END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set_prop)
      SUBROUTINE __FUNCNAME(this,discr_data,wp,info,read_only,ghosts_ok)

          IMPLICIT NONE

          CLASS(DTYPE(ppm_t_particles))                                 :: this
          CLASS(ppm_t_discr_data),                        INTENT(INOUT) :: discr_data

#if   __DIM == 1
          __TYPE, DIMENSION(:  ), POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          !!! data array

          INTEGER,                                        INTENT(  OUT) :: info
          !!! Return status, on success 0.

          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: ghosts_ok

          start_subroutine(__FUNCNAME)

          !If read_only was not explicitely set to true, then assume
          !that ghosts are no longer up to date, unless ghosts_ok was
          ! explicitely set to true
          IF (PRESENT(ghosts_ok)) THEN
             IF (.NOT.ghosts_ok) THEN
                !Assume that the ghost values are now incorrect
                discr_data%flags(ppm_ppt_ghosts)=.FALSE.
             ENDIF
          ENDIF

          IF (PRESENT(read_only)) THEN
             IF (.NOT.read_only) THEN
                !Assume that the ghost values are now incorrect
                discr_data%flags(ppm_ppt_ghosts)=.FALSE.
             ENDIF
          ENDIF

          wp => NULL()

          end_subroutine()
      END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get_field)
      SUBROUTINE __FUNCNAME(this,Field,wp,info,with_ghosts,read_only,skip_checks)
          !!! Returns a pointer to the data array where that contains
          !!! the discretized elements of Field on this particle set.
          CLASS(DTYPE(ppm_t_particles))                                 :: this
          CLASS(ppm_t_field_),                            INTENT(IN   ) :: Field

#if   __DIM == 1
          __TYPE, DIMENSION(:  ), POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          !!! data array

          INTEGER,                                        INTENT(  OUT) :: info
          !!! Return status, on success 0.

          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: with_ghosts
          !!! returns array between 1:Mpart (default is 1:Npart)
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          !!! swear on your favourite book that you will not modify the data
          !!! that you access. Default is false (which implies that the state
          !!! variables will be modify with the assumption that the data
          !!! accessed has been changed, e.g. so that a subsequent ghost update
          !!! will also update this property)
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: skip_checks
          !!! Only for users with huge cojones
          !!! [NOTE]
          !!! wp points to the data array

          CLASS(ppm_t_discr_data), POINTER :: discr_data

          INTEGER :: np

          LOGICAL :: skip

          start_subroutine(__FUNCNAME)

          NULLIFY(wp,discr_data)
          CALL Field%get_discr(this,discr_data,info)
          or_fail("could not get discr data for this field on that particle set")

          SELECT TYPE(prop => discr_data)
          CLASS IS (DTYPE(ppm_t_part_prop))
              IF (ppm_debug.GE.1) THEN
                 CALL prop%checktype(wp,info)
                 or_fail("Argument is not accessible for this property")
              ENDIF

              skip=MERGE(skip_checks,.FALSE.,PRESENT(skip_checks))

              IF (skip) THEN
#if   __DIM == 1
                 wp => prop%WRAP(DATANAME)
#elif __DIM == 2
                 wp => prop%WRAP(DATANAME)
#endif
              ELSE
                 np=this%Npart
                 IF (PRESENT(with_ghosts)) THEN
                    IF (with_ghosts) THEN
                       IF (prop%flags(ppm_ppt_ghosts)) THEN
                          np=this%Mpart
                       ELSE
                          stdout("ERROR: tried to get DATANAME (name = ",     &
                          & 'TRIM(ADJUSTL(prop%name))',                       &
                          & ") with ghosts when ghosts are not up-to-date. ", &
                          & "Returning NULL pointer")

                          fail("ghosts are not up-to-date")
                       ENDIF
                    ENDIF
                 ENDIF

                 IF (prop%flags(ppm_ppt_partial)) THEN
#if   __DIM == 1
                    wp => prop%WRAP(DATANAME)(1:np)
#elif __DIM == 2
                    wp => prop%WRAP(DATANAME)(:,1:np)
#endif
                 ELSE
                    stdout("ERROR: tried to get DATANAME (name = ", &
                    & 'TRIM(ADJUSTL(prop%name))',                   &
                    & ") when mapping is not up-to-date. ",         &
                    & "Returning NULL pointer",                     &
                    & "Run with traceback option to debug")

                    fail("unmapped particles")
                 ENDIF
              ENDIF !skip

              !Assume that the ghost values are now incorrect, unless explicitely
              !told otherwise
              IF (PRESENT(read_only)) THEN
                 IF (.NOT.read_only) THEN
                    prop%flags(ppm_ppt_ghosts)=.FALSE.
                 ENDIF
              ELSE
                 prop%flags(ppm_ppt_ghosts)=.FALSE.
              ENDIF
          CLASS DEFAULT
             fail("wrong type. Discretized data should be a ppm_t_part_prop")
          END SELECT

          check_associated(wp,"Get_Field returned a NULL pointer",ppm_err_sub_failed)

          end_subroutine()
      END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set_field)
      SUBROUTINE __FUNCNAME(this,Field,wp,info,read_only,ghosts_ok)
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles))                                 :: this
          CLASS(ppm_t_field_),                            INTENT(IN   ) :: Field

#if   __DIM == 1
          __TYPE, DIMENSION(:  ), POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          !!! data array

          INTEGER,                                        INTENT(  OUT) :: info
          !!! Return status, on success 0.

          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: ghosts_ok

          CLASS(ppm_t_discr_data), POINTER :: discr_data => NULL()

          start_subroutine(__FUNCNAME)

          CALL Field%get_discr(this,discr_data,info)
          or_fail("could not get discr data for this field on that particle set")

          SELECT TYPE(prop => discr_data)
          CLASS IS (DTYPE(ppm_t_part_prop))
             !If read_only was not explicitely set to true, then assume
             !that ghosts are no longer up to date, unless ghosts_ok was
             ! explicitely set to true
             IF (PRESENT(ghosts_ok)) THEN
                IF (.NOT.ghosts_ok) THEN
                   !Assume that the ghost values are now incorrect
                   prop%flags(ppm_ppt_ghosts)=.FALSE.
                ENDIF
             ENDIF

             IF (PRESENT(read_only)) THEN
                IF (.NOT.read_only) THEN
                   !Assume that the ghost values are now incorrect
                   prop%flags(ppm_ppt_ghosts)=.FALSE.
                ENDIF
             ENDIF

          END SELECT

          wp => NULL()

          end_subroutine()
      END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get)
#define __CHECKTYPE DTYPE(WRAP(DATANAME)_check)
      SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,info,with_ghosts,read_only,skip_checks)
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_particles))                                 :: Pc

#if   __DIM == 1
          __TYPE, DIMENSION(:  ), POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          !!! data array

          INTEGER,                                        INTENT(IN   ) :: ppt_id
          INTEGER,                                        INTENT(  OUT) :: info
          !!! Return status, on success 0.

          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: with_ghosts
          !!! returns array between 1:Mpart (default is 1:Npart)
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          !!! swear on your favourite book that you will not modify the data
          !!! that you access. Default is false (which implies that the state
          !!! variables will be modify with the assumption that the data
          !!! accessed has been changed, e.g. so that a subsequent ghost update
          !!! will also update this property)
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: skip_checks
          !!! Only for users with huge cojones

          LOGICAL :: lghosts
          LOGICAL :: skip

          start_subroutine(__FUNCNAME)

          NULLIFY(wp)

          IF (ppt_id.LE.0) THEN
             stdout("ERROR: failed to get DATANAME for property with ppt_id = ",ppt_id)
             fail("Cannot get property. Returning Null pointer")
          ENDIF

          IF (ppt_id.LE.Pc%props%max_id) THEN
             ASSOCIATE (prop => Pc%props%vec(ppt_id)%t)
                skip=MERGE(skip_checks,.FALSE.,PRESENT(skip_checks))
                IF (skip) THEN
#if   __DIM == 1
                   wp => prop%WRAP(DATANAME)
#elif __DIM == 2
                   wp => prop%WRAP(DATANAME)
#endif
                ELSE IF (prop%flags(ppm_ppt_partial)) THEN
                   lghosts =MERGE(with_ghosts,.FALSE.,PRESENT(with_ghosts))
                   IF (lghosts) THEN
                      IF (prop%flags(ppm_ppt_ghosts)) THEN
#if   __DIM == 1
                         wp => prop%WRAP(DATANAME)(1:Pc%Mpart)
#elif __DIM == 2
                         wp => prop%WRAP(DATANAME)(:,1:Pc%Mpart)
#endif
                      ELSE
                         stdout("ERROR: tried to get DATANAME (name = ",     &
                         & 'TRIM(ADJUSTL(prop%name))',                       &
                         & ") with ghosts when ghosts are not up-to-date. ", &
                         & "Returning NULL pointer")

                         fail("Ghosts not up-to-date. Call map_ghosts()?")
                      ENDIF
                   ELSE
#if   __DIM == 1
                      wp => prop%WRAP(DATANAME)(1:Pc%Npart)
#elif __DIM == 2
                      wp => prop%WRAP(DATANAME)(:,1:Pc%Npart)
#endif
                   ENDIF
                ELSE
                   stdout("ERROR: tried to get DATANAME (name = ", &
                   & 'TRIM(ADJUSTL(Pc%props%vec(ppt_id)%t%name))', &
                   & ") when mapping is not up-to-date. ",         &
                   & "Returning NULL pointer")

                   fail("unmapped particles")
                ENDIF

                IF (PRESENT(read_only)) THEN
                   IF (.NOT.read_only) prop%flags(ppm_ppt_ghosts)=.FALSE.
                ELSE
                   prop%flags(ppm_ppt_ghosts)=.FALSE.
                ENDIF
             END ASSOCIATE
          ELSE
             fail("Invalid id for particle property, returning NULL pointer")
          ENDIF

          end_subroutine()
      END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set)
      SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,info,read_only,ghosts_ok)
          CLASS(DTYPE(ppm_t_particles))                                 :: Pc
#if   __DIM == 1
          __TYPE, DIMENSION(:  ), POINTER, INTENT(INOUT) :: wp
#elif __DIM == 2
          __TYPE, DIMENSION(:,:), POINTER, INTENT(INOUT) :: wp
#endif
          !!! data array

          INTEGER,                                        INTENT(IN   ) :: ppt_id
          INTEGER,                                        INTENT(  OUT) :: info

          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                              OPTIONAL, INTENT(IN   ) :: ghosts_ok

          start_subroutine(__FUNCNAME)

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
          Pc%props%vec(ppt_id)%t%flags(ppm_ppt_ghosts)=.FALSE.

          wp => NULL()

          end_subroutine()
      END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#undef DATANAME
#undef __TYPE
#undef __MYTYPE
