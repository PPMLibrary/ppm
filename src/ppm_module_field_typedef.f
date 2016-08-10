! minclude ppm_header(ppm_module_field_typedef)

      MODULE ppm_module_field_typedef
      !!! Declares field data type

      !----------------------------------------------------------------------
      !  Modules
      !----------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_topo_typedef
      USE ppm_module_interfaces
      USE ppm_module_util_functions
      IMPLICIT NONE

      PRIVATE

      !----------------------------------------------------------------------
      ! Internal parameters
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------------------
      INTEGER               :: ppm_nb_fields = 0

      !----------------------------------------------------------------------
      ! Module variables
      !----------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldc

      !----------------------------------------------------------------------
      ! Type declaration
      !----------------------------------------------------------------------
      TYPE, EXTENDS(ppm_t_discr_info_):: ppm_t_discr_info
      CONTAINS
          PROCEDURE :: create => discr_info_create
          PROCEDURE :: destroy => discr_info_destroy
      END TYPE
minclude ppm_create_collection(discr_info,discr_info,generate="extend")

      TYPE,EXTENDS(ppm_t_field_) :: ppm_t_field
      CONTAINS
          PROCEDURE :: create            => field_create
          PROCEDURE :: destroy           => field_destroy
          PROCEDURE :: set_rel_discr     => field_set_rel_discr

          PROCEDURE :: map_ghost_push    => field_map_ghost_push
          PROCEDURE :: map_ghost_pop     => field_map_ghost_pop
          PROCEDURE :: map_push          => field_map_push
          PROCEDURE :: map_pop           => field_map_pop

          PROCEDURE :: is_discretized_on => field_is_discretized_on
          PROCEDURE :: discretize_on     => field_discretize_on
          PROCEDURE :: get_pid           => field_get_pid
          PROCEDURE :: get_discr         => field_get_discr
      END TYPE ppm_t_field
minclude ppm_create_collection(field,field,generate="extend")
minclude ppm_create_collection(field,field,generate="extend",vec=true)

      !----------------------------------------------------------------------
      ! PUBLIC from ppm_module_interfaces
      !----------------------------------------------------------------------
      PUBLIC :: ppm_t_discr_info_
      PUBLIC :: ppm_t_field_

      !----------------------------------------------------------------------
      ! PUBLIC
      !----------------------------------------------------------------------
      PUBLIC :: ppm_t_discr_info
      PUBLIC :: ppm_c_discr_info
      PUBLIC :: ppm_t_field
      PUBLIC :: ppm_c_field
      PUBLIC :: ppm_v_field


      !----------------------------------------------------------------------
      !  Type-bound procedures
      !----------------------------------------------------------------------
      CONTAINS

minclude ppm_create_collection_procedures(discr_info,discr_info_)
minclude ppm_create_collection_procedures(field,field_)
minclude ppm_create_collection_procedures(field,field_,vec=true)

      !CREATE
      SUBROUTINE discr_info_create(this,discr,discr_data,lda,flags,info,p_idx)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_info)                                        :: this
          CLASS(ppm_t_discr_kind),                         TARGET        :: discr
          CLASS(ppm_t_discr_data),                         TARGET        :: discr_data
          INTEGER,                                         INTENT(IN   ) :: lda
          !!! number of components
          LOGICAL, DIMENSION(ppm_param_length_mdataflags), INTENT(IN   ) :: flags
          INTEGER,                                         INTENT(  OUT) :: info
          INTEGER, OPTIONAL,                               INTENT(IN   ) :: p_idx
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          start_subroutine("discr_info_create")

          this%discrID    = discr%ID
          this%discr_ptr  => discr
          this%discr_data => discr_data
          this%lda        = lda
          this%flags      = flags
          IF (PRESENT(p_idx)) this%p_idx = p_idx

          end_subroutine()
      END SUBROUTINE discr_info_create
      !DESTROY
      SUBROUTINE discr_info_destroy(this,info)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_info) :: this

          INTEGER,  INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_),      POINTER :: sbp
          CLASS(ppm_t_subpatch_data_), POINTER :: sbpdat

          start_subroutine("discr_info_destroy")

          IF (this%discrID.GT.0) THEN
             ! yaser: I think this was a major bug!
             ! Without this all the allocated data still would be there
             ! even though you destroy them.
             SELECT TYPE (discr => this%discr_ptr)
             CLASS IS (ppm_t_equi_mesh_)
                IF (ASSOCIATED(this%discr_data)) THEN
                   SELECT TYPE (ddata => this%discr_data)
                   CLASS IS (ppm_t_mesh_discr_data_)
                      IF (ASSOCIATED(ddata%field_ptr)) THEN
                         IF (ASSOCIATED(discr%field_ptr)) THEN
                            CALL discr%field_ptr%remove(info,ddata%field_ptr)
                            or_fail("Failed to remove the element")
                         ENDIF !ASSOCIATED(discr%field_ptr)
                      ENDIF !ASSOCIATED(ddata%field_ptr)

                      IF (ASSOCIATED(discr%mdata)) THEN
                         CALL discr%mdata%remove(info,ddata)
                         or_fail("Failed to remove the element")
                      ENDIF !ASSOCIATED(discr%field_ptr)

                      CALL ddata%destroy(info)
                      or_fail("Failed to destroy ddata")

                      DEALLOCATE(this%discr_data,STAT=info)
                      or_fail_dealloc("Failed to deallocate ddata")
                   END SELECT

                   IF (ASSOCIATED(discr%subpatch)) THEN
                      IF (ASSOCIATED(discr%subpatch%vec)) THEN
                         sbp => discr%subpatch%begin()
                         sbp_loop: DO WHILE (ASSOCIATED(sbp))
                            IF (ASSOCIATED(sbp%subpatch_data)) THEN
                               IF (ASSOCIATED(sbp%subpatch_data%vec)) THEN
                                  sbpdat => sbp%subpatch_data%begin()
                                  DO WHILE (ASSOCIATED(sbpdat))
                                     IF (ASSOCIATED(sbpdat%discr_data)) THEN
                                        IF (ASSOCIATED(sbpdat%discr_data,this%discr_data)) THEN
                                           CALL sbp%subpatch_data%remove(info)
                                           or_fail("Failed to remove and destroy element")

                                           DEALLOCATE(sbpdat,STAT=info)
                                           or_fail_dealloc("Failed to deallocate sbpdat")

                                           sbp => discr%subpatch%next()
                                           CYCLE sbp_loop
                                        ENDIF !ASSOCIATED(sbpdat%discr_data,this%discr_data)
                                     ENDIF !ASSOCIATED(sbpdat%discr_data)
                                     sbpdat => sbp%subpatch_data%next()
                                  ENDDO
                               ENDIF !ASSOCIATED(sbp%subpatch_data%vec)
                            ENDIF !ASSOCIATED(sbp%subpatch_data)
                            sbp => discr%subpatch%next()
                         ENDDO sbp_loop
                      ENDIF !ASSOCIATED(discr%subpatch%vec)
                   ENDIF !ASSOCIATED(discr%subpatch)
                ENDIF !ASSOCIATED(this%discr_data)

             CLASS IS (ppm_t_particles_s_)
                IF (ASSOCIATED(discr%props)) THEN
                   SELECT TYPE (ddata => this%discr_data)
                   CLASS IS (ppm_t_discr_data)
                      CALL discr%field_ptr%remove(info,ddata%field_ptr)
                      or_fail("Failed to remove the bookkeeping element")

                      CALL discr%destroy_prop(ddata,info)
                      or_fail("particle discr data destroy failed.")
                   END SELECT
                ELSE
                   SELECT TYPE (ddata => this%discr_data)
                   CLASS IS (ppm_t_part_prop_s_)
                      CALL ddata%destroy(info)
                      or_fail("discr data destroy failed.")

                   END SELECT
                ENDIF

                DEALLOCATE(this%discr_data,STAT=info)
                or_fail_dealloc("this%discr_data")

             CLASS IS (ppm_t_particles_d_)
                IF (ASSOCIATED(discr%props)) THEN
                   SELECT TYPE (ddata => this%discr_data)
                   CLASS IS (ppm_t_discr_data)
                      CALL discr%field_ptr%remove(info,ddata%field_ptr)
                      or_fail("Failed to remove the bookkeeping element")

                      CALL discr%destroy_prop(ddata,info)
                      or_fail("particle discr data destroy failed.")

                   END SELECT
                ELSE
                   SELECT TYPE (ddata => this%discr_data)
                   CLASS IS (ppm_t_part_prop_d_)
                      CALL ddata%destroy(info)
                      or_fail("discr data destroy failed.")

                   END SELECT
                ENDIF
                DEALLOCATE(this%discr_data,STAT=info)
                or_fail_dealloc("this%discr_data")

             END SELECT
          ENDIF !(this%discrID.GT.0)

          this%discrID    = 0
          this%discr_ptr  => NULL()
          this%discr_data => NULL()
          this%lda        = 0
          this%flags      = .FALSE.

          end_subroutine()
      END SUBROUTINE discr_info_destroy

      !CREATE
      SUBROUTINE field_create(this,lda,info,dtype,name,init_func)
          !!! Constructor for fields
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                             :: this

          INTEGER,                         INTENT(IN   ) :: lda
          !!! number of components
          INTEGER,                         INTENT(  OUT) :: info
          !!!
          INTEGER,               OPTIONAL, INTENT(IN   ) :: dtype
          !!! data type (ppm_type_int, ppm_type_real, ppm_type_comp, ppm_type_logical)
          CHARACTER(LEN=*),      OPTIONAL, INTENT(IN   ) :: name
          !!!
          REAL(ppm_kind_double), OPTIONAL, POINTER, EXTERNAL :: init_func
          !!! support for initialization function not finished (need to think
          !!! about data types...)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          start_subroutine("field_create")

          !Use a global ID (this may not be needed, after all...)
          ppm_nb_fields = ppm_nb_fields + 1
          this%ID = ppm_nb_fields

          check_true(<#lda.GT.0#>,&
          & "Trying to create field with zero components. LDA must be > 0")

          this%lda = lda
          IF (PRESENT(name)) THEN
             this%name = TRIM(ADJUSTL(name))
          ELSE
             WRITE(this%name,'(A,I0)') "dft_field_",ppm_nb_fields
          ENDIF
          IF (PRESENT(dtype)) THEN
             check_true(<#dtype.GT.0#>,"dtype must be > 0")
             this%data_type = dtype
          ELSE
             this%data_type = MERGE(ppm_type_real,ppm_type_real_single,ppm_kind.EQ.ppm_kind_double)
          ENDIF

          check_false(<#ASSOCIATED(this%discr_info)#>,&
          & "Seems like this field was already allocated - Call destroy() first?")

          ALLOCATE(ppm_c_discr_info::this%discr_info,STAT=info)
          or_fail_alloc("this%discr_info")

          end_subroutine()
      END SUBROUTINE field_create
      !DESTROY
      SUBROUTINE field_destroy(this,info)
          !!! Destructor for subdomain data data structure
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field), TARGET        :: this

          INTEGER,            INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_info_),      POINTER :: tdi
          CLASS(ppm_t_subpatch_),        POINTER :: sbp
          CLASS(ppm_t_subpatch_data_),   POINTER :: sbpdat
          CLASS(ppm_t_mesh_discr_data_), POINTER :: mddata
          CLASS(ppm_t_part_prop_d_),     POINTER :: pddata_d
          CLASS(ppm_t_part_prop_s_),     POINTER :: pddata_s

          start_subroutine("field_destroy")

          this%ID = 0
          this%name = ''
          this%lda = 0
          this%data_type = 0

          !Destroy the bookkeeping entries in the fields that are
          !discretized on this mesh or particle
          !TODO !! yaser I think this is all
          IF (ASSOCIATED(this%discr_info)) THEN
             tdi => this%discr_info%begin()
             DO WHILE (ASSOCIATED(tdi))
                SELECT TYPE (dp => tdi%discr_ptr)
                CLASS IS (ppm_t_equi_mesh_)
                   IF (this%is_discretized_on(dp)) THEN
                      CALL dp%field_ptr%remove(info,this)
                      or_fail("Failed to remove the bookkeeping field entries in Mesh")

                      IF (ASSOCIATED(dp%subpatch)) THEN
                         sbp => dp%subpatch%begin()
                         sbp_loop: DO WHILE (ASSOCIATED(sbp))
                            IF (ASSOCIATED(sbp%subpatch_data)) THEN
                               sbpdat => sbp%subpatch_data%begin()
                               sbpdat_loop: DO WHILE (ASSOCIATED(sbpdat))
                                  IF (ASSOCIATED(sbpdat%discr_data)) THEN
                                     IF (ASSOCIATED(sbpdat%discr_data%field_ptr)) THEN
                                        IF (ASSOCIATED(sbpdat%discr_data%field_ptr,this)) THEN
                                           CALL sbp%subpatch_data%remove(info,sbpdat)
                                           or_fail("Failed to destroy and remove the discretized field entries in Mesh")

                                           DEALLOCATE(sbpdat,STAT=info)
                                           or_fail_dealloc("Failed to deallocate sbpdat")
                                        ENDIF !ASSOCIATED(sbpdat%discr_data%field_ptr,this)
                                     ENDIF !ASSOCIATED(sbpdat%discr_data%field_ptr)
                                  ENDIF !ASSOCIATED(sbpdat%discr_data)
                                  sbpdat => sbp%subpatch_data%next()
                               ENDDO sbpdat_loop
                            ENDIF !(ASSOCIATED(sbp%subpatch_data))
                            sbp => dp%subpatch%next()
                         ENDDO sbp_loop
                      ENDIF !(ASSOCIATED(dp%subpatch))

                      IF (ASSOCIATED(dp%mdata)) THEN
                         mddata => dp%mdata%begin()
                         mddata_loop: DO WHILE (ASSOCIATED(mddata))
                            IF (ASSOCIATED(mddata%field_ptr,this)) THEN
                               CALL dp%mdata%remove(info,mddata)
                               or_fail("Failed to remove the bookkeeping field entries in Mesh mdata")

                               CALL mddata%destroy(info)
                               or_fail("Failed to destroy mddata")

                               DEALLOCATE(mddata,STAT=info)
                               or_fail_dealloc("Failed to deallocate mddata")
                            ENDIF
                            mddata => dp%mdata%next()
                         ENDDO mddata_loop
                      ENDIF !(ASSOCIATED(dp%mdata))
                   ENDIF !(this%is_discretized_on(dp))

                CLASS IS (ppm_t_particles_d_)
                   IF (this%is_discretized_on(dp)) THEN
                      CALL dp%field_ptr%remove(info,this)
                      or_fail("Failed to remove the bookkeeping field entries in Particle")

                      IF (ASSOCIATED(dp%props)) THEN
                         pddata_d => dp%props%begin()
                         pddata_d_loop: DO WHILE (ASSOCIATED(pddata_d))
                            IF (ASSOCIATED(pddata_d%field_ptr,this)) THEN
                                CALL dp%props%remove(info,pddata_d)
                                or_fail("Failed to remove and destroy the discretized field entries in Particle")

                                DEALLOCATE(pddata_d,STAT=info)
                                or_fail_dealloc("Failed to deallocate pddata_d")
                            ENDIF !(ASSOCIATED(sbp%subpatch_data))
                            pddata_d => dp%props%next()
                         ENDDO pddata_d_loop
                      ENDIF !(ASSOCIATED(dp%props))
                   ENDIF !(this%is_discretized_on(dp))

                CLASS IS (ppm_t_particles_s_)
                   IF (this%is_discretized_on(dp)) THEN
                      CALL dp%field_ptr%remove(info,this)
                      or_fail("Failed to remove the bookkeeping field entries in Particle")

                      IF (ASSOCIATED(dp%props)) THEN
                         pddata_s => dp%props%begin()
                         pddata_s_loop: DO WHILE (ASSOCIATED(pddata_s))
                            IF (ASSOCIATED(pddata_s%field_ptr,this)) THEN
                                CALL dp%props%remove(info,pddata_s)
                                or_fail("Failed to remove the discretized field entries in Particle")

                                DEALLOCATE(pddata_s,STAT=info)
                                or_fail_dealloc("Failed to deallocate pddata_s")
                            ENDIF !(ASSOCIATED(sbp%subpatch_data))
                            pddata_s => dp%props%next()
                         ENDDO pddata_s_loop
                      ENDIF !(ASSOCIATED(dp%props))
                   ENDIF !(this%is_discretized_on(dp))

                END SELECT !TYPE (dp => tdi%discr_ptr)

                tdi%discrID=-1

                tdi => this%discr_info%next()
             ENDDO
          ENDIF

          destroy_collection_ptr(this%discr_info)

          end_subroutine()
      END SUBROUTINE field_destroy

      SUBROUTINE field_discretize_on(this,discr,info,datatype,with_ghosts,discr_info)
          !!! Allocate field on a mesh
          !!! If the field has a procedure for initialization (e.g. an
          !!! initial condition), then the field is also initialized.
          !!! If the mesh has patches, the field is allocated only on these
          !!! patches. If no patches have been defined, it is assumed
          !!! that the user expects the field to be allocated on the whole domain.
          !!! A single patch is then defined, covering all the subdomains.
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field),                 TARGET        :: this
          CLASS(ppm_t_discr_kind),            TARGET        :: discr
          !!! mesh or Particle set onto which this field is to be discretized
          INTEGER,                            INTENT(  OUT) :: info
          INTEGER,                  OPTIONAL, INTENT(IN   ) :: datatype
          !!! By default, the type is assumed to be real, double-precision.
          LOGICAL,                  OPTIONAL, INTENT(IN   ) :: with_ghosts
          !!! By default, the data arrays are allocated with Mpart iif the ghost
          !!! particles are up-to-date. Otherwise, they are allocated with size Npart.
          !!! Setting with_ghosts to true or false forces the allocation to be done
          !!! with size Mpart or Npart, respectively.
          CLASS(ppm_t_discr_info_), OPTIONAL, POINTER       :: discr_info
          !!! Pointer to the new data and bookkeeping information for the mesh
          !!! or particle set on which this field has been discretized.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
!           CLASS(ppm_t_subpatch_),     POINTER :: p => NULL()
!           CLASS(ppm_t_subpatch_data_),POINTER :: subpdat => NULL()
          INTEGER                             :: dtype,p_idx
          LOGICAL                             :: lghosts
          CLASS(ppm_t_mesh_discr_data_),POINTER :: mddata
          !CLASS(ppm_t_part_prop_d_),    POINTER :: pddata => NULL()
          CLASS(ppm_t_discr_data),    POINTER :: pddata
          CLASS(ppm_t_main_abstr),      POINTER :: el

          start_subroutine("field_discretize_on")

          IF (PRESENT(datatype)) THEN
             dtype = datatype
          ELSE
             check_true(<#this%data_type.NE.0#>,"field data type has not been initialized. Fix constructor.")
             dtype = this%data_type
          ENDIF

          !Check whether this field has already been initialized
          check_true(<#this%ID.GT.0.AND.this%lda.GT.0#>,&
          & "Field needs to be initialized before calling discretized. Call ThisField%create() first")

          check_false(<#this%is_discretized_on(discr)#>,&
          & "Method to re-discretize a field on an existing discretization (overwriting, reallocation of data) is not yet implemented. TODO!")

          SELECT TYPE(discr)
          CLASS IS (ppm_t_equi_mesh_)
             !Check that the mesh contains patches onto which the data
             ! should be allocated. If not, create a single patch that
             ! covers the whole domain.
             IF (.NOT.ASSOCIATED(discr%subpatch)) THEN
                fail("Mesh not allocated. Use mesh%create() first.")
             ELSE
                IF (discr%npatch.LE.0) THEN
                   CALL discr%def_uniform(info)
                   or_fail("failed to create a uniform patch data structure")
                ENDIF
             ENDIF

             NULLIFY(mddata)
             CALL discr%create_prop(this,mddata,info,p_idx)
             or_fail("discr%create_prop()")

             !Update the bookkeeping table to store the relationship between
             !the mesh and the field.
             CALL this%set_rel_discr(discr,mddata,info,p_idx)
             or_fail("failed to log the relationship between this field and that mesh")

             !CALL discr%set_rel(this,info)
             IF (PRESENT(discr_info)) THEN
                discr_info => this%discr_info%last()
             ENDIF

             el => this
             CALL discr%field_ptr%push(el,info)
             or_fail("failed to log the relationship between this mesh and that field")

          CLASS IS (ppm_t_particles_d_)
             lghosts =MERGE(with_ghosts,discr%flags(ppm_part_ghosts),PRESENT(with_ghosts))

             NULLIFY(pddata)
             !Create a new property data structure in the particle set to store this field
             CALL discr%create_prop(info,this,discr_data=pddata,with_ghosts=lghosts)
             or_fail("discr%create_prop")

             !Update the bookkeeping table to store the relationship between
             ! the particle and the field.
             CALL this%set_rel_discr(discr,pddata,info)
             or_fail("failed to log the relationship between this field and that particle set")
             !CALL this%set_rel(discr,p_idx,info)

             IF (PRESENT(discr_info)) THEN
                discr_info => this%discr_info%last()
             ENDIF

             el => this
             CALL discr%field_ptr%push(el,info)
             or_fail("failed to log the relationship between this particle set and that field")

          CLASS IS (ppm_t_particles_s_)
             lghosts =MERGE(with_ghosts,discr%flags(ppm_part_ghosts),PRESENT(with_ghosts))

             NULLIFY(pddata)
             !Create a new property data structure in the particle set to store this field
             CALL discr%create_prop(info,this,discr_data=pddata,with_ghosts=lghosts)
             or_fail("discr%create_prop")

             !Update the bookkeeping table to store the relationship between
             ! the particle and the field.
             CALL this%set_rel_discr(discr,pddata,info)
             or_fail("failed to log the relationship between this field and that particle set")
             !CALL this%set_rel(discr,p_idx,info)

             IF (PRESENT(discr_info)) THEN
                discr_info => this%discr_info%last()
             ENDIF

             el => this
             CALL discr%field_ptr%push(el,info)
             or_fail("failed to log the relationship between this particle set and that field")

          CLASS DEFAULT
             fail("support for this discretization type is missing")

          END SELECT

          end_subroutine()
      END SUBROUTINE field_discretize_on

      SUBROUTINE field_set_rel_discr(this,discr,discr_data,info,p_idx)
          !!! Create bookkeeping data structure to log the relationship between
          !!! the field and a particle set
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                             :: this
          CLASS(ppm_t_discr_kind),         TARGET        :: discr
          !!! mesh, or particle set, that this field is discretized on
          CLASS(ppm_t_discr_data),         TARGET        :: discr_data
          !!! data (or mesh or on particle set) for this discretization

          INTEGER,                         INTENT(  OUT) :: info
          INTEGER,OPTIONAL,                INTENT(IN   ) :: p_idx
          !!! storage index in the collection of discretization
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_info_),            POINTER    :: dinfo => NULL()

          LOGICAL, DIMENSION(ppm_param_length_mdataflags) :: flags

          start_subroutine("field_set_rel_discr")

          check_false(<#this%is_discretized_on(discr)#>,&
          "This field has already been discretized here. Cannot do that twice")

          flags = .FALSE.

          ALLOCATE(ppm_t_discr_info::dinfo,STAT=info)
          or_fail_alloc("could not allocate new discr_info pointer")

          IF (PRESENT(p_idx)) THEN
             CALL dinfo%create(discr,discr_data,this%lda,flags,info,p_idx)
          ELSE
             CALL dinfo%create(discr,discr_data,this%lda,flags,info)
          ENDIF
          or_fail("could not create new discr_info object")

          IF (.NOT.ASSOCIATED(this%discr_info)) THEN
             ALLOCATE(ppm_c_discr_info::this%discr_info,STAT=info)
             or_fail_alloc("could not allocate discr_info collection")
          ENDIF

          CALL this%discr_info%push(dinfo,info)
          or_fail("could not add dinfo to the discr_info collection")

          end_subroutine()
      END SUBROUTINE field_set_rel_discr

      SUBROUTINE field_map_ghost_push(this,mesh,info)
          !!! Push field data into the buffers of a mesh for ghost mappings
          !!! The field must of course be stored on this mesh
          !!! (for now) we assume that mesh%map_ghost_get() has already been called
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                     :: this
          CLASS(ppm_t_equi_mesh_), TARGET        :: mesh
          !!! mesh that this field is discretized on
          INTEGER,                 INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          start_subroutine("field_map_ghost_push")

          CALL mesh%map_ghost_push(this,info)
          or_fail("mesh%map_ghost_push")

          end_subroutine()
      END SUBROUTINE field_map_ghost_push

      SUBROUTINE field_map_ghost_pop(this,mesh,info,poptype)
          !!! Pop field data from the buffers of a mesh for ghost mappings
          !!! The field must of course be stored on this mesh
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                     :: this
          CLASS(ppm_t_equi_mesh_), TARGET        :: mesh
          !!! mesh that this field is discretized on

          INTEGER,                 INTENT(  OUT) :: info
          INTEGER,      OPTIONAL,  INTENT(IN   ) :: poptype
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          start_subroutine("field_map_ghost_pop")

          SELECT CASE (PRESENT(poptype))
          CASE (.FALSE.)
             CALL mesh%map_ghost_pop(this,info)
          CASE (.TRUE.)
             CALL mesh%map_ghost_pop(this,info,poptype)
          END SELECT
          or_fail("mesh%map_ghost_pop")

          end_subroutine()
      END SUBROUTINE field_map_ghost_pop

      SUBROUTINE field_map_push(this,mesh,info)
          !!! Push field data into the buffers of a mesh for global mappings
          !!! The field must of course be stored on this mesh
          !!! (for now) we assume that mesh%map() has already been called
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                     :: this
          CLASS(ppm_t_equi_mesh_), TARGET        :: mesh
          !!! mesh that this field is discretized on

          INTEGER,                 INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          start_subroutine("field_map_push")

          CALL mesh%map_push(this,info)
          or_fail("mesh%map_push")

          end_subroutine()
      END SUBROUTINE field_map_push

      SUBROUTINE field_map_pop(this,mesh,info)
          !!! Pop field data from the buffers of a mesh for global mappings
          !!! The field must of course be stored on this mesh
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                     :: this
          CLASS(ppm_t_equi_mesh_), TARGET        :: mesh
          !!! mesh that this field is discretized on

          INTEGER,                 INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          start_subroutine("field_map_pop")

          CALL mesh%map_pop(this,info)
          or_fail("mesh%map_pop")

          end_subroutine()
      END SUBROUTINE field_map_pop

      FUNCTION field_get_pid(this,discr_kind,tstep) RESULT(p_idx)
          !!! Returns a pointer to the discretization (mesh or particles)
          !!! on which the field is discretized.
          !!! TODO: Optionally, can retrieve discretizations at different time points
          !!! (like n-1, n-2, etc...).
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                     :: this
          CLASS(ppm_t_discr_kind), TARGET        :: discr_kind
          !!! discretization

          INTEGER,       OPTIONAL, INTENT(IN   ) :: tstep
          !!! TODO
          !!! If the current time step is n, discretizations at previous times
          !!! can be accessed using the tstep argument
          !!!     (tstep =  0, default) => step n
          !!!     (tstep = -1         ) => step n-1
          !!!     (tstep = -2         ) => step n-2
          !!!     etc...
          !!! The number of steps that are stored is set when a field is first
          !!! discretized on a mesh or on a particle set.
          !!! The data management and book-keeping between the data of different
          !!! time steps are done by the time integrator.
          INTEGER                                      :: p_idx
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_info_), POINTER :: dinfo

          start_function("field_get_pid")

          dinfo => this%discr_info%begin()
          loop: DO WHILE(ASSOCIATED(dinfo))
              check_associated(<#dinfo%discr_ptr#>)
              IF (ASSOCIATED(dinfo%discr_ptr,discr_kind)) THEN
                  p_idx = dinfo%p_idx
                  RETURN
              ENDIF
              dinfo => this%discr_info%next()
          ENDDO loop

          p_idx = -10

          end_function()
      END FUNCTION field_get_pid

      SUBROUTINE field_get_discr(this,discr_kind,discr_data,info,tstep)
          !!! Returns a pointer to the discretization (mesh or particles)
          !!! on which the field is discretized.
          !!! TODO: Optionally, can retrieve discretizations at different time points
          !!! (like n-1, n-2, etc...).
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                     :: this
          CLASS(ppm_t_discr_kind), TARGET        :: discr_kind
          CLASS(ppm_t_discr_data), POINTER       :: discr_data
          !!! discretization

          INTEGER,                 INTENT(  OUT) :: info
          INTEGER,       OPTIONAL, INTENT(IN   ) :: tstep
          !!! If the current time step is n, discretizations at previous times
          !!! can be accessed using the tstep argument
          !!!     (tstep =  0, default) => step n
          !!!     (tstep = -1         ) => step n-1
          !!!     (tstep = -2         ) => step n-2
          !!!     etc...
          !!! The number of steps that are stored is set when a field is first
          !!! discretized on a mesh or on a particle set.
          !!! The data management and book-keeping between the data of different
          !!! time steps are done by the time integrator.
          !!! TODO
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_info_), POINTER :: dinfo

          start_subroutine("field_get_discr")

          dinfo => this%discr_info%begin()
          loop: DO WHILE(ASSOCIATED(dinfo))
              IF (ASSOCIATED(dinfo%discr_ptr,discr_kind)) THEN
                 discr_data => dinfo%discr_data
                 EXIT loop
              ENDIF
              dinfo => this%discr_info%next()
          ENDDO loop

          check_associated(discr_data, &
          & "Field seems to not be distretized on this particle set")

          end_subroutine()
      END SUBROUTINE field_get_discr

      FUNCTION field_is_discretized_on(this,discr_kind,discr_info,tstep) RESULT(res)
          !!! Check if this field has been discretized on this
          !!! discretization kind (mesh or particles)
          !!! TODO: Optionally, check for discretizations at different time points
          !!! (like n-1, n-2, etc...).
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field)                                :: this

          CLASS(ppm_t_discr_kind),            TARGET        :: discr_kind
          !!! discretization
          CLASS(ppm_t_discr_info_), OPTIONAL, POINTER       :: discr_info

          INTEGER,                  OPTIONAL, INTENT(IN   ) :: tstep
          !!! If the current time step is n, discretizations at previous times
          !!! can be accessed using the tstep argument
          !!!     (tstep =  0, default) => step n
          !!!     (tstep = -1         ) => step n-1
          !!!     (tstep = -2         ) => step n-2
          !!!     etc...
          !!! The number of steps that are stored is set when a field is first
          !!! discretized on a mesh or on a particle set.
          !!! The data management and book-keeping between the data of different
          !!! time steps are done by the time integrator.
          !!! TODO
          LOGICAL                                           :: res
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_info_), POINTER :: dinfo

          start_function("field_is_discretized_on")

          IF (ASSOCIATED(this%discr_info)) THEN
             dinfo => this%discr_info%begin()
             DO WHILE(ASSOCIATED(dinfo))
                IF (ASSOCIATED(dinfo%discr_ptr)) THEN
                   IF (ASSOCIATED(dinfo%discr_ptr,discr_kind)) THEN
                      res = .TRUE.
                      IF (PRESENT(discr_info)) THEN
                         discr_info => dinfo

                         check_associated(discr_info, &
                         & "Field seems to not be distretized on this set")
                      ENDIF
                      RETURN
                   ENDIF
                ENDIF
                dinfo => this%discr_info%next()
             ENDDO
          ENDIF

          res = .FALSE.
          RETURN

          end_function()
      END FUNCTION field_is_discretized_on

      END MODULE ppm_module_field_typedef
