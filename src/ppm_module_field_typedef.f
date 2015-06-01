! minclude ppm_header(ppm_module_field_typedef)

      MODULE ppm_module_field_typedef
      !!! Declares field data type

      !----------------------------------------------------------------------
      !  Modules
      !----------------------------------------------------------------------
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_topo_typedef
      USE ppm_module_interfaces
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_util_functions
      IMPLICIT NONE

      !----------------------------------------------------------------------
      ! Internal parameters
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------------------
      INTEGER                        :: ppm_nb_fields = 0

      !----------------------------------------------------------------------
      ! Module variables
      !----------------------------------------------------------------------
      INTEGER, PRIVATE, DIMENSION(3) :: ldc

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
      !  Type-bound procedures
      !----------------------------------------------------------------------
      CONTAINS

minclude ppm_create_collection_procedures(discr_info,discr_info_)
minclude ppm_create_collection_procedures(field,field_)
minclude ppm_create_collection_procedures(field,field_,vec=true)

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
          !!! support for initialisation function not finished (need to think
          !!! about data types...)

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
             this%data_type = ppm_type_real
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
          CLASS(ppm_t_field)     :: this

          INTEGER, INTENT(  OUT) :: info

          CLASS(ppm_t_discr_info_), POINTER :: tdi

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
                      or_fail("Failed to destroy the bookkeeping field entries in Mesh")
                   ENDIF

                CLASS IS (ppm_t_particles_d_)
                   IF (this%is_discretized_on(dp)) THEN
                      CALL dp%field_ptr%remove(info,this)
                      or_fail("Failed to destroy the bookkeeping field entries in Particle")
                   ENDIF

                CLASS IS (ppm_t_particles_s_)
                   IF (this%is_discretized_on(dp)) THEN
                      CALL dp%field_ptr%remove(info,this)
                      or_fail("Failed to destroy the bookkeeping field entries in Particle")
                   ENDIF

                END SELECT
                tdi => this%discr_info%next()
             ENDDO
          ENDIF

          destroy_collection_ptr(this%discr_info)

          end_subroutine()
      END SUBROUTINE field_destroy
      !CREATE
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

          CLASS(ppm_t_subpatch_data_), POINTER :: subpdat

          start_subroutine("discr_info_destroy")

          ! yaser: I think this was a major bug!
          ! Without this all the allocated data still would be there
          ! even though you destroy them.
          SELECT TYPE (discr => this%discr_ptr)
          CLASS IS (ppm_t_equi_mesh_)

             !Ugly measure to make sure that we are destroying available data
             !not the one that has been destroyed previously
             !TODO find st nicer
             IF (ASSOCIATED(discr%subpatch)) THEN
                IF (ASSOCIATED(discr%subpatch%vec)) THEN

                   SELECT TYPE (ddata => this%discr_data)
                   CLASS IS (ppm_t_mesh_discr_data_)
                      subpdat => ddata%subpatch%begin()
                      DO WHILE (ASSOCIATED(subpdat))
                         CALL subpdat%destroy(info)
                         or_fail("subpdat%destroy")

                         subpdat => ddata%subpatch%next()
                      ENDDO

                      CALL ddata%destroy(info)
                      or_fail("mesh discr data destroy failed.")

                   END SELECT

                ENDIF
             ENDIF

          CLASS IS (ppm_t_particles_s_)
             IF (ASSOCIATED(discr%props)) THEN
                SELECT TYPE (ddata => this%discr_data)
                CLASS IS (ppm_t_discr_data)
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

          CLASS IS (ppm_t_particles_d_)
             IF (ASSOCIATED(discr%props)) THEN
                SELECT TYPE (ddata => this%discr_data)
                CLASS IS (ppm_t_discr_data)
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

          END SELECT

          this%discrID    = 0
          this%discr_ptr  => NULL()
          this%discr_data => NULL()
          this%lda        = 0
          this%flags      = .FALSE.

          end_subroutine()
      END SUBROUTINE discr_info_destroy

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
          END IF

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

          CLASS(ppm_t_discr_info_),            POINTER    :: dinfo => NULL()
          LOGICAL, DIMENSION(ppm_param_length_mdataflags) :: flags

          start_subroutine("field_set_rel_discr")

          check_false(<#this%is_discretized_on(discr)#>,&
          "This field has already been discretized here. Cannot do that twice")

          flags = .FALSE.

          ALLOCATE(ppm_t_discr_info::dinfo,STAT=info)
          or_fail_alloc("could not allocate new discr_info pointer")

          CALL dinfo%create(discr,discr_data,this%lda,flags,info,p_idx)
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
